# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('polyRAD', 'dplyr')

set.seed(20260808)

# 03. What this script does -----------------------------------------------------
# Takes the true simulated dosage from 01, which a real breeding programme does
# not observe directly, then (a) simulates noisy GBS read counts from it and
# (b) calls dosage back out of those read counts with polyRAD, the
# dosage-uncertain genotype-calling tool the article cites (Clark et al. 2019)
# as the precedent for this step. Downstream scripts (03 onwards) use polyRAD's
# estimated dosage, not the true simulated dosage, matching what a real
# programme has available.
#
# Uses polyRAD's mapping-population pipeline (PipelineMapping2Parents), not the
# population-level IterateHWE, since this is a biparental F1 from two known,
# unrelated hexaploid parents rather than a population in Hardy-Weinberg
# equilibrium. The pipeline expects the parent1/parent2 plus progeny naming
# convention used below; SetDonorParent() and SetRecurrentParent() take a taxon
# name and return a modified object, so they are called as ordinary functions
# and not as `<-` replacement functions.

IN  <- "../data_simulated"
OUT <- "../data_simulated"

progeny_dosage <- readRDS(file.path(IN, "progeny_dosage_matrix.rds"))  # markers x progeny, TRUE dosage 0-6
mm <- readRDS(file.path(IN, "marker_map.rds"))
marker_id <- mm$marker_id
lg_id     <- mm$lg_id
n_markers <- length(marker_id)
n_progeny <- ncol(progeny_dosage)
PLOIDY    <- 6

# 04. Simulate realistic GBS read counts from the true dosage --------------------
# Depth: polyploid dosage calling needs deeper coverage than diploid SNP calling
# to resolve the 7 dosage classes (0-6), so depth is Poisson with mean 40x per
# marker x sample (typical GBS depth variability).
# Read sampling: alt-allele reads ~ Binomial(depth, p), where p is the true
# allele ratio (dosage/6) adjusted for a small sequencing/PCR error rate. The
# read counts therefore carry sampling and error noise rather than re-encoding
# the true dosage, which would make polyRAD's estimation step circular.
MEAN_DEPTH  <- 40
ERROR_RATE  <- 0.01

depth <- matrix(rpois(n_markers * n_progeny, MEAN_DEPTH), nrow = n_markers, ncol = n_progeny)
true_ratio <- progeny_dosage / PLOIDY
obs_ratio  <- true_ratio * (1 - ERROR_RATE) + (1 - true_ratio) * ERROR_RATE
alt_reads  <- matrix(rbinom(n_markers * n_progeny, depth, obs_ratio), nrow = n_markers, ncol = n_progeny)
ref_reads  <- depth - alt_reads
dimnames(alt_reads) <- dimnames(ref_reads) <- dimnames(progeny_dosage)

# 05. Parental read counts (PipelineMapping2Parents genotypes the parents to
# anchor allele frequencies) ---------------------------------------------------
# 01 does not save the parental copy matrices, so parent read counts are
# reconstructed from the population-level allele frequency: each parent's
# expected dosage under the allele frequency used to found the cross.
parent_dose <- matrix(rbinom(n_markers * 2, PLOIDY, mm$allele_freq), nrow = n_markers, ncol = 2,
                       dimnames = list(marker_id, c("parent1", "parent2")))
parent_ratio <- parent_dose / PLOIDY
parent_depth <- matrix(rpois(n_markers * 2, MEAN_DEPTH), nrow = n_markers, ncol = 2)
parent_obs_ratio <- parent_ratio * (1 - ERROR_RATE) + (1 - parent_ratio) * ERROR_RATE
parent_alt <- matrix(rbinom(n_markers * 2, parent_depth, parent_obs_ratio), nrow = n_markers, ncol = 2,
                      dimnames = list(marker_id, c("parent1", "parent2")))
parent_ref <- parent_depth - parent_alt

alt_reads_all <- cbind(parent_alt, alt_reads)
ref_reads_all <- cbind(parent_ref, ref_reads)
all_taxa <- colnames(alt_reads_all)

# 06. Build the RADdata object ----------------------------------------------------
# alleleDepth: taxa in rows, alleles in columns; 2 allele columns per SNP
# (ref, alt), interleaved per locus, matching polyRAD's expected layout.
alleleDepth <- matrix(0L, nrow = length(all_taxa), ncol = 2 * n_markers,
                       dimnames = list(all_taxa, NULL))
alleleDepth[, seq(1, 2 * n_markers, by = 2)] <- t(ref_reads_all)
alleleDepth[, seq(2, 2 * n_markers, by = 2)] <- t(alt_reads_all)
colnames(alleleDepth) <- paste0(rep(marker_id, each = 2), c("_ref", "_alt"))

alleles2loc <- rep(seq_len(n_markers), each = 2)
locTable <- data.frame(Chr = lg_id, Pos = seq_len(n_markers), row.names = marker_id)

rad <- RADdata(
  alleleDepth       = alleleDepth,
  alleles2loc       = alleles2loc,
  locTable          = locTable,
  possiblePloidies  = list(PLOIDY),          # pure autohexaploid, no ambiguity
  contamRate        = 0.001,
  alleleNucleotides = rep(c("A", "T"), n_markers),  # placeholder bases (no real sequence available)
  # polyRAD internally computes effective ploidy as
  # sum(possiblePloidies[[i]]) * taxaPloidy / 2: with
  # possiblePloidies = list(6), taxaPloidy must be 2 to get effective
  # ploidy 6*2/2 = 6, matching the real exampleRAD_mapping's own diploid
  # case (possiblePloidies = list(2), taxaPloidy = 2 -> 2*2/2 = 2). Setting
  # taxaPloidy = PLOIDY (6) here would give effective ploidy 18.
  taxaPloidy        = 2L
)

rad <- SetDonorParent(rad, "parent1")
rad <- SetRecurrentParent(rad, "parent2")

# 07. Run polyRAD's mapping-population pipeline -----------------------------------
# Simple F1 (no backcrossing, intermating or selfing generations); useLinkage =
# FALSE because locTable$Pos here is an arbitrary marker index rather than a
# physical or genetic map position, so linkage-based prior updating is not
# meaningful with these placeholder coordinates.
# freqAllowedDeviation is set to 0.01 instead of the 0.05 default.
# PipelineMapping2Parents builds its expected-frequency grid as
# seq(0, 1, length.out = (pld.donor + pld.recurrent) * max(ploidy)/2 + 1); for
# two hexaploid (ploidy 6) parents that is (6+6)*6/2 + 1 = 37 points, i.e. steps
# of 1/36 (~0.0278), against 5 points and steps of 0.25 in the diploid case the
# 0.05 default is tuned for. Half the minimum gap is ~0.0139, so
# freqAllowedDeviation must sit below that.
rad <- PipelineMapping2Parents(
  rad,
  n.gen.backcrossing  = 0,
  n.gen.intermating   = 0,
  n.gen.selfing       = 0,
  useLinkage          = FALSE,
  freqAllowedDeviation = 0.01
)

# 08. Extract estimated dosage (posterior mean genotypes, rescaled to 0-6) -------
# omit1allelePerLocus/omitCommonAllele default to TRUE and keep whichever
# allele is RARER per locus, which mixes ref- and alt-allele dosage across
# markers. Both allele columns are requested here and "_alt" selected by
# name, so the comparison to progeny_dosage (alt-allele dosage in 01) uses
# one allele throughout.
wmg <- GetWeightedMeanGenotypes(rad, minval = 0, maxval = PLOIDY,
                                 omit1allelePerLocus = FALSE)
progeny_rows <- setdiff(rownames(wmg), c("parent1", "parent2"))
alt_cols <- grep("_alt$", colnames(wmg), value = TRUE)
est_dosage <- t(wmg[progeny_rows, alt_cols, drop = FALSE])
rownames(est_dosage) <- sub("_alt$", "", rownames(est_dosage))
est_dosage <- est_dosage[marker_id, , drop = FALSE]

# 09. Sanity-check estimated dosage against the ground truth (diagnostic only) ---
common_progeny <- intersect(colnames(est_dosage), colnames(progeny_dosage))
r_check <- cor(as.vector(est_dosage[, common_progeny]),
               as.vector(progeny_dosage[, common_progeny]), use = "complete.obs")
cat(sprintf("polyRAD estimated vs. true simulated dosage, marker-wise correlation: %.3f\n", r_check))
cat("(expected well below 1.0; this is the dosage uncertainty the article\n")
cat(" discusses. A value close to 0 would indicate a calling problem rather\n")
cat(" than noise.)\n")

# 10. Write output -----------------------------------------------------------------
saveRDS(est_dosage, file.path(OUT, "estimated_dosage_matrix.rds"))
cat("Estimated dosage matrix:", nrow(est_dosage), "markers x", ncol(est_dosage), "progeny\n")
cat("Dosage value range:", round(range(est_dosage, na.rm = TRUE), 2), "\n")
cat("Wrote:", file.path(OUT, "estimated_dosage_matrix.rds"), "\n")
