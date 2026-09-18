# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('dplyr')

set.seed(20260808)

# 03. What this script does -----------------------------------------------------
# Runs two validations of the Stage-3 pipeline against Stage-1 MGIDI, for two
# complementary purposes:
#
# (A) Oracle validation (simulation only). Correlates each method's score
# against the true, noise-free simulated breeding value (01's
# true_breeding_value_matrix.rds), a quantity no real breeding programme has
# access to and which is exactly known here because the population is simulated.
# This is how simulation studies validate genomic-prediction pipelines, and it
# avoids a construction artifact in option (B): MGIDI and the observed composite
# are both near-linear combinations of the same five observed traits, so a high
# correlation between them is close to a mathematical certainty rather than
# evidence that MGIDI predicts breeding merit better than GBLUP. The observed
# composite's own correlation with the true breeding value, the noise ceiling
# printed below, bounds what any method can reach against that target. Against
# the true breeding value the two methods separate: GBLUP's accuracy exceeds
# MGIDI's and the bootstrap CI on the difference excludes zero (values printed
# below).
#
# (B) Table 5's literal operational gate. The bootstrap-CI test as Table 5
# specifies it for real deployment, where no true breeding value is available:
# each method's score correlated against the shared observed composite outcome
# (04/05's composite_observed), gated by whether the accuracy difference's
# bootstrap CI excludes zero. Its limitation: MGIDI's score here is fitted using
# each genotype's own observed data, since metan::mgidi() has no mechanism to
# score a genotype blind the way GBLUP's held-out fold prediction can, so this
# comparison structurally favours MGIDI beyond predictive skill. (A) is the
# artifact-free comparison; (B) is what Table 5 asks a real programme, which has
# no true breeding value to fall back on, to compute.

IN  <- "../outputs"
IN_SIM <- "../data_simulated"
OUT <- "../outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_BOOT <- 10000
CI_LEVEL <- 0.95
TRAITS <- c("Yield", "DryMatter", "Carotenoid", "RootShapeScore", "StorageScore")

gblup  <- readRDS(file.path(IN, "gblup_cv_results.rds"))
stage1 <- readRDS(file.path(IN, "stage1_mgidi_results.rds"))
true_bv <- readRDS(file.path(IN_SIM, "true_breeding_value_matrix.rds"))  # genotype x trait, pre-noise

# 04. Assemble the paired per-genotype comparison table ---------------------------
gblup_df <- data.frame(
  GEN = names(gblup$composite_gebv),
  Composite_GEBV = as.numeric(gblup$composite_gebv),
  Composite_Observed = as.numeric(gblup$composite_observed[names(gblup$composite_gebv)])
)

composite_true_bv <- setNames(rowMeans(scale(true_bv[, TRAITS])), rownames(true_bv))

cmp <- gblup_df %>%
  inner_join(stage1$mgidi_df %>% select(GEN, MGIDI), by = "GEN") %>%
  mutate(Composite_TrueBV = composite_true_bv[GEN])

stopifnot(nrow(cmp) == length(gblup$composite_gebv))  # every genotype held out exactly once (04), matched here
cat("Paired comparison table:", nrow(cmp), "genotypes\n\n")

# 05. (A) Oracle validation: accuracy against the TRUE simulated breeding value ---
acc_gblup_oracle  <- cor(cmp$Composite_GEBV, cmp$Composite_TrueBV)
acc_stage1_oracle <- cor(-cmp$MGIDI, cmp$Composite_TrueBV)
noise_ceiling     <- cor(cmp$Composite_Observed, cmp$Composite_TrueBV)

cat("=== (A) Oracle validation (simulation only, vs. TRUE breeding value) ===\n")
cat(sprintf("GBLUP (composite, out-of-fold) accuracy:  %.3f\n", acc_gblup_oracle))
cat(sprintf("Stage-1 (MGIDI) accuracy:                 %.3f\n", acc_stage1_oracle))
cat(sprintf("Noise ceiling (observed composite itself): %.3f  (upper bound any noisy-phenotype-based method could reach)\n",
            noise_ceiling))
cat(sprintf("Difference (GBLUP - Stage-1): %.3f\n\n", acc_gblup_oracle - acc_stage1_oracle))

boot_diff_oracle <- numeric(N_BOOT)
n <- nrow(cmp)
for (b in seq_len(N_BOOT)) {
  idx <- sample.int(n, n, replace = TRUE)
  boot_diff_oracle[b] <- cor(cmp$Composite_GEBV[idx], cmp$Composite_TrueBV[idx]) -
    cor(-cmp$MGIDI[idx], cmp$Composite_TrueBV[idx])
}
alpha <- 1 - CI_LEVEL
ci_oracle <- quantile(boot_diff_oracle, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
cat(sprintf("Bootstrap %.0f%% CI on the oracle accuracy difference: [%.3f, %.3f]\n\n",
            CI_LEVEL * 100, ci_oracle[1], ci_oracle[2]))

# 06. (B) Table 5's literal operational gate: accuracy against the OBSERVED composite ---
acc_gblup_point  <- cor(cmp$Composite_GEBV, cmp$Composite_Observed)
acc_stage1_point <- cor(-cmp$MGIDI, cmp$Composite_Observed)

cat("=== (B) Table 5's literal gate (operational: vs. OBSERVED composite; MGIDI not held out, see header note) ===\n")
cat(sprintf("GBLUP (composite, out-of-fold) accuracy: %.3f\n", acc_gblup_point))
cat(sprintf("Stage-1 (MGIDI, in-sample) accuracy:      %.3f\n", acc_stage1_point))
cat(sprintf("Difference (GBLUP - Stage-1): %.3f\n\n", acc_gblup_point - acc_stage1_point))

boot_diff <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
  idx <- sample.int(n, n, replace = TRUE)
  boot_diff[b] <- cor(cmp$Composite_GEBV[idx], cmp$Composite_Observed[idx]) -
    cor(-cmp$MGIDI[idx], cmp$Composite_Observed[idx])
}
ci <- quantile(boot_diff, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
cat(sprintf("Bootstrap %.0f%% CI on (GBLUP - Stage-1) accuracy: [%.3f, %.3f]\n", CI_LEVEL * 100, ci[1], ci[2]))

verdict <- if (ci[1] > 0) {
  "GO"
} else if (ci[2] < 0) {
  "NO-GO (Stage-1 significantly better; see the MGIDI in-sample caveat above)"
} else {
  "NO-GO (improvement not distinguishable from zero; CI spans zero)"
}
cat(sprintf("\nTable 5 Stage-3 literal gate verdict: %s\n", verdict))
cat("((A) above, the oracle comparison, is the artifact-free check that the pipeline\n")
cat(" itself is recovering real signal; this gate is what Table 5 asks a real\n")
cat(" programme, which never has a true breeding value, to compute operationally.)\n")

# 07. Write output -----------------------------------------------------------------
result <- list(
  comparison_table    = cmp,
  # (A) oracle
  acc_gblup_oracle    = acc_gblup_oracle,
  acc_stage1_oracle   = acc_stage1_oracle,
  noise_ceiling       = noise_ceiling,
  boot_diff_oracle    = boot_diff_oracle,
  ci_oracle           = ci_oracle,
  # (B) Table 5 literal
  acc_gblup_point     = acc_gblup_point,
  acc_stage1_point    = acc_stage1_point,
  boot_diff           = boot_diff,
  ci                  = ci,
  ci_level            = CI_LEVEL,
  n_boot              = N_BOOT,
  verdict             = verdict
)
saveRDS(result, file.path(OUT, "stage_gate_decision.rds"))
write.csv(cmp, file.path(OUT, "stage_gate_comparison_table.csv"), row.names = FALSE)
cat("\nWrote:", file.path(OUT, "stage_gate_decision.rds"), "\n")
