# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('sommer', 'dplyr')

set.seed(20260808)

# 03. What this script does -----------------------------------------------------
# Dosage-aware GBLUP genomic prediction, evaluated by five-fold
# cross-validation, for all five simulated traits, so that the comparison
# against the Stage-1 MGIDI index, itself a multi-trait composite, matches in
# predictive scope. Table 5 allows two ways to make the Stage-1-vs-Stage-3
# comparison like-for-like: correlate both methods' scores against "a specific
# target trait, or a pre-specified composite performance score computed
# identically for every genotype." This script takes the composite route. Every
# trait gets its own out-of-fold GBLUP prediction from one shared fold
# assignment across all five traits, so a genotype is held out of every trait
# model at the same time and the individual-level composite score is out-of-fold
# for all five traits. Those five per-trait GEBVs are combined into a single
# composite predicted-merit score with the same formula (per-trait z-score, then
# mean) that 05 and 06 apply to the observed phenotypes, so Stage 3's composite
# score and Stage 1's MGIDI score are evaluated against the same shared outcome
# Table 5 specifies.
#
# Simplification: each trait gets its own univariate GBLUP model (5 separate
# sommer fits per fold) rather than a single joint multi-trait mixed model with
# an estimated cross-trait genetic covariance structure. A joint multi-trait
# GBLUP would let traits borrow strength from their genetic correlations and is
# the next refinement; the independent-per-trait fit keeps the two methods'
# predictive scope matched, since both draw on all five traits and are evaluated
# against one shared composite outcome.
#
# API notes (sommer 4.4.6): vsr(GEN, Gu = Gmat) is the custom-
# relationship-matrix random effect; fit$U[["u:GEN"]] holds the GEBV for
# every genotype level in Gu, including ones whose phenotype was masked to
# NA for the held-out fold; fit$convergence is a documented boolean field,
# checked for all 25 fits below.
#
# GBLUP is fitted on genotype-level trait means (already averaged over the 2
# field replicates) rather than on raw plot-level records with an explicit
# replicate term. This two-step approach (collapse repeated measures to a
# genotype summary, then fit the genomic model on that summary) loses precision
# in the variance-component estimates relative to a one-step model, and matches
# the summary step 05's MGIDI baseline takes by construction (metan::gamem() is
# itself a genotype-summary step before mgidi() runs).

IN  <- "../data_simulated"
OUT <- "../outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

Gmat  <- readRDS(file.path(IN, "G_matrix_hexaploid.rds"))
pheno <- read.csv(file.path(IN, "plot_level_phenotypes.csv"))

TRAITS  <- c("Yield", "DryMatter", "Carotenoid", "RootShapeScore", "StorageScore")
N_FOLDS <- 5

# 04. Genotype-level trait means, aligned to the G-matrix's genotype order ------
geno_mean <- pheno %>%
  group_by(GEN) %>%
  summarise(across(all_of(TRAITS), mean), .groups = "drop") %>%
  filter(GEN %in% rownames(Gmat)) %>%
  mutate(GEN = factor(GEN, levels = rownames(Gmat))) %>%
  arrange(GEN)

stopifnot(nrow(geno_mean) == nrow(Gmat))  # every genotype in G has a phenotype here (simulated, complete)

n <- nrow(geno_mean)

# 05. One shared fold assignment for all five traits ---------------------------
# Genotype i is held out of every trait's training set in the same round, so its
# five GEBVs (combined into one composite score below) are out-of-fold at the
# same time, rather than mixing in-sample predictions for some traits with
# out-of-sample predictions for others in the same individual.
fold_id <- sample(rep(seq_len(N_FOLDS), length.out = n))

gebv_matrix <- matrix(NA_real_, nrow = n, ncol = length(TRAITS),
                      dimnames = list(as.character(geno_mean$GEN), TRAITS))
trait_accuracy <- matrix(NA_real_, nrow = N_FOLDS, ncol = length(TRAITS),
                         dimnames = list(NULL, TRAITS))
convergence_log <- matrix(NA, nrow = N_FOLDS, ncol = length(TRAITS),
                          dimnames = list(NULL, TRAITS))

# 06. Per-trait GBLUP, cross-validated on the shared fold assignment -------------
for (tr in TRAITS) {
  for (k in seq_len(N_FOLDS)) {
    test_idx <- which(fold_id == k)
    y_masked <- geno_mean[[tr]]
    y_masked[test_idx] <- NA

    df_fit <- data.frame(GEN = geno_mean$GEN, y_masked = y_masked)
    fit <- mmer(y_masked ~ 1, random = ~ vsr(GEN, Gu = Gmat), rcov = ~ units,
                data = df_fit, verbose = FALSE, dateWarning = FALSE)

    convergence_log[k, tr] <- isTRUE(fit$convergence)

    gebv <- fit$U[["u:GEN"]][["y_masked"]]
    gebv <- gebv[as.character(geno_mean$GEN[test_idx])]

    gebv_matrix[as.character(geno_mean$GEN[test_idx]), tr] <- gebv
    trait_accuracy[k, tr] <- cor(gebv, geno_mean[[tr]][test_idx])
  }
  cat(sprintf("Trait %-15s cross-validated accuracy: mean = %.3f, sd = %.3f\n",
              tr, mean(trait_accuracy[, tr]), sd(trait_accuracy[, tr])))
}

n_failed <- sum(!convergence_log)
if (n_failed > 0) {
  cat(sprintf("\nWARNING: %d / %d GBLUP fits (trait x fold) did NOT report convergence:\n",
              n_failed, length(convergence_log)))
  print(which(!convergence_log, arr.ind = TRUE))
} else {
  cat(sprintf("\nConvergence check: all %d GBLUP fits (5 traits x %d folds) converged.\n",
              length(convergence_log), N_FOLDS))
}

# 07. Combine into ONE composite predicted-merit score, out-of-fold ---------------
# Same formula 05/06 apply to the observed phenotypes: z-score each trait
# across genotypes (all five traits are "higher is better," matching the
# same all-"h" ideotype MGIDI uses in 05), then average. This makes the
# composite GEBV and the composite observed outcome directly comparable on
# the same standardised scale.
observed_matrix <- as.matrix(geno_mean[, TRAITS])
rownames(observed_matrix) <- as.character(geno_mean$GEN)

compute_composite <- function(mat) rowMeans(scale(mat))

composite_gebv     <- compute_composite(gebv_matrix)
composite_observed <- compute_composite(observed_matrix)

composite_accuracy <- cor(composite_gebv, composite_observed)
cat(sprintf("\nComposite (5-trait, standardise-then-average) GBLUP accuracy: %.3f\n",
            composite_accuracy))

# 08. Write output -----------------------------------------------------------------
saveRDS(list(gebv_matrix = gebv_matrix, observed_matrix = observed_matrix,
             composite_gebv = composite_gebv, composite_observed = composite_observed,
             composite_accuracy = composite_accuracy,
             trait_accuracy = trait_accuracy, convergence_log = convergence_log,
             fold_id = setNames(fold_id, as.character(geno_mean$GEN))),
        file.path(OUT, "gblup_cv_results.rds"))
write.csv(as.data.frame(trait_accuracy) %>% mutate(fold = seq_len(N_FOLDS), .before = 1),
          file.path(OUT, "gblup_cv_accuracy_by_fold.csv"), row.names = FALSE)
cat("Wrote:", file.path(OUT, "gblup_cv_results.rds"), "\n")
