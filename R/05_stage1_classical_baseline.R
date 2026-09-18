# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('metan', 'dplyr')

# 03. What this script does -----------------------------------------------------
# Stage-1 classical multi-trait selection baseline (Table 5), computed the
# same way the authors' companion germplasm workflow does it
# (https://github.com/ntducphd/vnua_sweetpotato_improvement: gamem() then a
# multivariate index), applied here to the SAME simulated phenotype data used
# for the Stage-3 GBLUP pipeline (04), so the Stage-1-vs-Stage-3 gate
# comparison in 06 is evaluated on one shared population, not two different
# ones.
#
# Asymmetry between the two methods (see README): MGIDI is a composite
# ideotype-distance score computed from a genotype's own observed multi-trait
# phenotypes, and has no mechanism to extrapolate to a genotype it was given no
# phenotypes for, whereas GBLUP extrapolates through the G-matrix to genotypes
# without their own trial data. Matching how MGIDI is used in breeding practice
# (ranking genotypes that have already been phenotyped), its score here is
# computed once on the full observed dataset and is not cross-validated the way
# GBLUP is in 04. The quantity compared in 06 is each method's correlation with
# the same shared observed outcome (05's composite_observed, computed the same
# way as 04's and evaluated for both methods at the same held-out individuals
# GBLUP's folds use). MGIDI's score for those individuals was informed by their
# own observed data; GBLUP's prediction for its held-out fold was not. Table 5
# sets out this comparison logic, and Section 7 of the article discusses the
# asymmetry.
#
# Scope matching: MGIDI is a 5-trait composite index, so 04 runs GBLUP across
# all five traits and combines them into a composite score with the same
# standardise-then-average formula used here. Both methods are therefore
# evaluated against the same multi-trait composite outcome, the "pre-specified
# composite performance score computed identically for every genotype"
# comparison mode Table 5 allows.
#
# Representation matching: metan::mgidi() defaults to use_data = "blup", so
# MGIDI's ideotype-distance is computed from gamem()'s BLUP-shrunk genotype
# estimates rather than the raw phenotypic means, while composite_observed below
# and 04's GBLUP training target are both raw genotype means. This script sets
# use_data = "pheno" so MGIDI is computed from the same raw phenotypic means it
# is compared against (metan's mgidi.Rd documents use_data = "pheno" as the
# option for "phenotypic means instead [of] BLUPs").

IN  <- "../data_simulated"
OUT <- "../outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("Yield", "DryMatter", "Carotenoid", "RootShapeScore", "StorageScore")

pheno <- read.csv(file.path(IN, "plot_level_phenotypes.csv"))
pheno$GEN <- as.character(pheno$GEN)

# 04. Genotype-mean BLUPs across all 5 traits -------------------------------------
model <- gamem(pheno, gen = GEN, rep = REP, resp = all_of(TRAITS))

# 05. MGIDI multi-trait selection index -------------------------------------------
# All 5 traits treated as "higher is better" (ideotype = all "h").
# RootShapeScore/StorageScore direction in a
# real programme would depend on the actual scoring convention used, which
# this simulator does not encode.
mgidi_index <- mgidi(model, use_data = "pheno", ideotype = rep("h", length(TRAITS)), SI = 15, verbose = FALSE)

mgidi_scores <- mgidi_index$MGIDI  # tibble with columns Genotype, MGIDI (not rownames)
cat("MGIDI computed for", nrow(mgidi_scores), "genotypes.\n")
cat("Selected (top SI%):", length(mgidi_index$sel_gen), "genotypes\n")

# 06. Same shared composite observed outcome as 04 (recomputed independently ------
# here from the same raw phenotype file, so this script runs standalone. The
# computation is deterministic, so it matches 04's value.)
geno_traits <- pheno %>% group_by(GEN) %>% summarise(across(all_of(TRAITS), mean), .groups = "drop")
observed_matrix <- as.matrix(geno_traits[, TRAITS])
rownames(observed_matrix) <- geno_traits$GEN
composite_observed <- setNames(rowMeans(scale(observed_matrix)), rownames(observed_matrix))

mgidi_df <- mgidi_scores %>%
  rename(GEN = Genotype) %>%
  mutate(composite_observed = composite_observed[GEN])

# MGIDI is a distance (lower = closer to ideotype = better), so its negative is
# correlated with the composite outcome and "higher predicted merit" points the
# same direction for both methods in 06's comparison.
stage1_full_cor <- cor(-mgidi_df$MGIDI, mgidi_df$composite_observed, use = "complete.obs")
cat(sprintf("Stage-1 (MGIDI) whole-population correlation with the shared composite outcome: %.3f\n", stage1_full_cor))

# 07. Write output -----------------------------------------------------------------
saveRDS(list(mgidi_df = mgidi_df, stage1_full_cor = stage1_full_cor),
        file.path(OUT, "stage1_mgidi_results.rds"))
write.csv(mgidi_df, file.path(OUT, "stage1_mgidi_scores.csv"), row.names = FALSE)
cat("Wrote:", file.path(OUT, "stage1_mgidi_results.rds"), "\n")
