# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('dplyr', 'tidyr', 'purrr')

set.seed(20260808)

# 03. What this script does -----------------------------------------------------
# Simulates a biparental hexaploid (2n = 6x = 90) sweetpotato population to
# validate the Stage-3 genomic-selection pipeline end to end, since no real
# GBS/SNP genotyping data exists for the authors' sweetpotato panel (only
# phenotypic MGIDI data and SSR fingerprinting do).
#
# Simplification (see README.md): inheritance is modeled as random hexasomic
# (polysomic) segregation with no preferential pairing, and each of the 15
# linkage groups (matching sweetpotato's basic chromosome number, 15 x 6 = 90)
# is treated as a single non-recombining block per homolog. Within-linkage-group
# recombination is not simulated. The article records that sweetpotato's meiotic
# behaviour "does not fit auto- or allopolyploid models cleanly" (Table 2; Gao
# et al. 2020, PLoS ONE 15(3) e0229624). This simulator validates pipeline
# mechanics and the Table-5 gate logic; it does not model real sweetpotato
# transmission genetics.
#
# Population size (315 progeny) matches the precedent the article cites:
# Gemenet et al. (2020), the one real sweetpotato genomic-selection
# feasibility study, used a 315-progeny hexaploid biparental population.

N_PROGENY   <- 315
N_LG        <- 15     # linkage groups = sweetpotato's basic chromosome number
MARKERS_PER_LG <- 100  # -> 1500 markers total, GBS-scale
N_COPIES    <- 6       # hexaploid: 6 homologous copies per locus
N_REPS      <- 2       # field replicates per genotype, for BLUP estimation

OUT <- "../data_simulated"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 04. Simulate the two parents' marker genotypes ---------------------------------
# Each marker: population allele frequency drawn from Beta(2,2) (moderate MAF
# spread), then each parent's 6 copies at that marker drawn independently as
# Bernoulli(p) -- approximates two unrelated founders with realistic
# heterozygosity, not a specific real accession.
n_markers <- N_LG * MARKERS_PER_LG
marker_id <- sprintf("SNP_%04d", seq_len(n_markers))
lg_id     <- rep(seq_len(N_LG), each = MARKERS_PER_LG)
allele_freq <- rbeta(n_markers, 2, 2)

simulate_parent_copies <- function(af) {
  matrix(rbinom(length(af) * N_COPIES, 1, af), nrow = length(af), ncol = N_COPIES)
}
P1_copies <- simulate_parent_copies(allele_freq)  # n_markers x 6
P2_copies <- simulate_parent_copies(allele_freq)

# 05. Generate 315 progeny via random hexasomic gamete formation -----------------
# For each parent, for each linkage group, choose 3-of-6 homolog indices
# (without replacement) ONCE per linkage group -- all markers in that group
# inherit together from the chosen copies (single-block-per-LG simplification).
make_gamete <- function(parent_copies) {
  gamete <- integer(nrow(parent_copies))
  for (lg in seq_len(N_LG)) {
    idx <- which(lg_id == lg)
    chosen <- sample.int(N_COPIES, 3)
    # gamete contribution per marker = sum of the 3 chosen homolog alleles
    gamete[idx] <- rowSums(parent_copies[idx, chosen, drop = FALSE])
  }
  gamete
}

progeny_dosage <- matrix(NA_integer_, nrow = n_markers, ncol = N_PROGENY,
                          dimnames = list(marker_id, sprintf("SP%03d", seq_len(N_PROGENY))))
for (j in seq_len(N_PROGENY)) {
  progeny_dosage[, j] <- make_gamete(P1_copies) + make_gamete(P2_copies)
}
# dosage now 0-6 copies of the alternate allele per marker per progeny, as
# expected from a hexaploid GBS dosage call

# 06. Simulate a realistic multi-trait phenotype architecture --------------------
# Traits mirror the jointly-constrained breeding targets the article names for
# sweetpotato: yield, dry matter, carotenoid content, root shape, storage
# behaviour. All five traits enter both the Stage-3 composite (04) and the
# Stage-1 MGIDI index (05), matching the comparison mode Section 7 specifies:
# "a specific target trait, or a pre-specified composite performance score
# computed identically for every genotype".
TRAITS <- c("Yield", "DryMatter", "Carotenoid", "RootShapeScore", "StorageScore")

# assign each trait a set of causal QTL, with effect sizes scaled to a
# target heritability
n_qtl_per_trait <- 40
qtl_effects <- map(TRAITS, function(tr) {
  qtl_idx <- sample.int(n_markers, n_qtl_per_trait)
  effects <- rnorm(n_qtl_per_trait, mean = 0, sd = 1)
  list(idx = qtl_idx, effects = effects)
})
names(qtl_effects) <- TRAITS

# Shared QTL are injected to induce trait correlations. Incidental overlap
# between each trait's independently-drawn QTL set leaves the realised
# correlations near zero. QTL are shared between chosen trait pairs with a
# controlled effect-sign relationship: Yield/Carotenoid
# share QTL with OPPOSITE-sign effects (approximating the real
# beta-carotene/starch trade-off Gemenet et al. document, cited in
# Table 1 of the article); DryMatter/StorageScore share QTL with
# SAME-sign effects (both traits relate to real storage-root starch/
# dry-matter accumulation, so a positive association is the more
# realistic simplification).
N_SHARED_QTL <- 20
inject_shared_qtl <- function(qtl_effects, trait_a, trait_b, n_shared, sign) {
  shared_idx <- sample.int(n_markers, n_shared)
  eff_a <- rnorm(n_shared, mean = 0, sd = 1)
  eff_b <- sign * eff_a
  qtl_effects[[trait_a]]$idx     <- c(qtl_effects[[trait_a]]$idx, shared_idx)
  qtl_effects[[trait_a]]$effects <- c(qtl_effects[[trait_a]]$effects, eff_a)
  qtl_effects[[trait_b]]$idx     <- c(qtl_effects[[trait_b]]$idx, shared_idx)
  qtl_effects[[trait_b]]$effects <- c(qtl_effects[[trait_b]]$effects, eff_b)
  qtl_effects
}
qtl_effects <- inject_shared_qtl(qtl_effects, "Yield", "Carotenoid", N_SHARED_QTL, sign = -1)
qtl_effects <- inject_shared_qtl(qtl_effects, "DryMatter", "StorageScore", N_SHARED_QTL, sign = +1)

H2_TARGET <- 0.40  # broad-sense heritability on a genotype-MEAN basis (see resid_sd formula below)

breeding_value <- sapply(TRAITS, function(tr) {
  q <- qtl_effects[[tr]]
  bv <- t(progeny_dosage[q$idx, , drop = FALSE]) %*% q$effects
  setNames(as.numeric(bv), colnames(progeny_dosage))
})
rownames(breeding_value) <- colnames(progeny_dosage)
colnames(breeding_value) <- TRAITS

# scale each trait's genetic variance, then add plot-level residual noise so
# that gamem()-style REML recovers ~H2_TARGET on GENOTYPE MEANS (not on a
# single plot observation). With N_REPS field replicates per genotype, the
# genotype-mean heritability is h2_mean = Vg / (Vg + Ve_plot/N_REPS), not
# Vg / (Vg + Ve_plot): with N_REPS = 2, Ve_plot/N_REPS is 2x smaller than
# Ve_plot itself, so the plot-level formula would inflate the realised
# genotype-mean h2. Solving h2_mean = Vg / (Vg + Ve_plot/N_REPS) = H2_TARGET
# for the plot-level residual SD gives the sqrt(N_REPS) correction factor
# below; gamem() reports the realised genotype-mean h2 for each trait.
geno_sd_target <- c(Yield = 8, DryMatter = 3, Carotenoid = 5,
                    RootShapeScore = 1.2, StorageScore = 1.0)
plot_raw <- expand.grid(GEN = rownames(breeding_value), REP = seq_len(N_REPS)) %>%
  arrange(GEN, REP)

for (tr in TRAITS) {
  bv <- breeding_value[, tr]
  bv <- (bv - mean(bv)) / sd(bv) * geno_sd_target[[tr]]
  resid_sd <- geno_sd_target[[tr]] * sqrt((1 - H2_TARGET) / H2_TARGET) * sqrt(N_REPS)
  trait_mean <- switch(tr, Yield = 28, DryMatter = 30, Carotenoid = 12,
                        RootShapeScore = 5, StorageScore = 5)
  bv_named <- setNames(bv, rownames(breeding_value))
  plot_raw[[tr]] <- trait_mean + bv_named[as.character(plot_raw$GEN)] +
    rnorm(nrow(plot_raw), 0, resid_sd)
}

# 07. Write outputs ---------------------------------------------------------------
saveRDS(progeny_dosage, file.path(OUT, "progeny_dosage_matrix.rds"))
saveRDS(list(marker_id = marker_id, lg_id = lg_id, allele_freq = allele_freq),
        file.path(OUT, "marker_map.rds"))
write.csv(plot_raw, file.path(OUT, "plot_level_phenotypes.csv"), row.names = FALSE)
saveRDS(qtl_effects, file.path(OUT, "true_qtl_effects.rds"))
# Noise-free per-genotype breeding value, before the residual noise added in
# section 06. A real breeding programme does not have this quantity; 06 uses it
# for the oracle validation check of whether the pipeline recovers genetic
# signal, which is possible in simulation and not on real data. Monotonic
# scaling does not affect correlation-based accuracy metrics, so the raw
# (pre-scaling) values here are equivalent to the scaled version used to build
# the phenotypes.
saveRDS(breeding_value, file.path(OUT, "true_breeding_value_matrix.rds"))

cat(sprintf(
  "Simulated %d progeny x %d markers (%d linkage groups), %d traits, %d reps/genotype.\n",
  N_PROGENY, n_markers, N_LG, length(TRAITS), N_REPS))
cat("Dosage range check:", range(progeny_dosage), "(expect 0-6)\n")

# 08. Diagnostics: realised genetic correlations -------------------------------
cat("\nTrue genetic correlation matrix (breeding values; verifies the injected\n")
cat("Yield/Carotenoid trade-off and DryMatter/StorageScore association):\n")
print(round(cor(breeding_value[, TRAITS]), 2))
cat("\n(expect Yield-Carotenoid negative and DryMatter-StorageScore positive,\n")
cat(" matching the sign of the injected shared-QTL effects; other pairs near\n")
cat(" zero, with no injected relationship between them.)\n")

cat("Wrote:", OUT, "\n")
