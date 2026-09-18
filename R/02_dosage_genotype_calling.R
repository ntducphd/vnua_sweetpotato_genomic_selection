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
# freqAllowedDeviation is set to 0.01 instead of the 0.05 default, because the
# default is too large for this grid and PipelineMapping2Parents stops with
# "allowedDeviation is too large given intervals within expectedFreqs".
# The grid is built as
#   allelesin <- (pld.don + pld.rec) * pld.max / 2
#   possfreq  <- seq(0, 1, length.out = (n.gen.backcrossing + 1) * allelesin + 1)
# where pld.don and pld.rec are the two parents' taxaPloidy and pld.max is
# max(sum(possiblePloidies)). Here that is (2 + 2) * 6 / 2 = 12, so 13 grid
# points with a spacing of 1/12 (about 0.0833) and an upper limit of half that
# spacing, 1/24 (about 0.0417). 0.01 therefore sits well inside the limit; it is
# about four times tighter than necessary, and that choice has a cost, recorded
# in section 09 below and in README's Known limitations: the tighter the
# tolerance, the more loci whose observed allele frequency matches no expected
# mapping frequency and are left without any dosage call.
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

# 09. Call rate, then sanity-check the called dosage against the ground truth ------
# The call rate is reported first and explicitly, because the correlation below
# is computed only over cells polyRAD actually called and would otherwise read
# as if every locus had been genotyped. Loci left uncalled here are removed
# later by 03's missing-rate filter, so this is where that loss originates.
uncalled_alt  <- sum(rowMeans(is.na(est_dosage)) == 1)
wmg_progeny   <- wmg[progeny_rows, , drop = FALSE]
uncalled_both <- sum(colMeans(is.na(wmg_progeny[, paste0(marker_id, "_ref"), drop = FALSE])) == 1 &
                     colMeans(is.na(wmg_progeny[, paste0(marker_id, "_alt"), drop = FALSE])) == 1)
cat(sprintf("Call rate: %d / %d loci uncalled for BOTH alleles; %d more called only on the\n",
            uncalled_both, nrow(est_dosage), uncalled_alt - uncalled_both))
cat(sprintf(" reference allele, so %d / %d markers (%.1f%% of cells) leave this script with no\n",
            uncalled_alt, nrow(est_dosage), 100 * mean(is.na(est_dosage))))
cat(" alt-allele dosage at all. See the freqAllowedDeviation note in section 07.\n")

common_progeny <- intersect(colnames(est_dosage), colnames(progeny_dosage))
r_check <- cor(as.vector(est_dosage[, common_progeny]),
               as.vector(progeny_dosage[, common_progeny]), use = "complete.obs")
cat(sprintf("polyRAD estimated vs. true simulated dosage, marker-wise correlation over the\n"))
cat(sprintf(" CALLED cells only: %.3f\n", r_check))
cat("(expected well below 1.0; this is the dosage uncertainty the article\n")
cat(" discusses. A value close to 0 would indicate a calling problem rather\n")
cat(" than noise.)\n")

# 10. Write output -----------------------------------------------------------------
saveRDS(est_dosage, file.path(OUT, "estimated_dosage_matrix.rds"))
cat("Estimated dosage matrix:", nrow(est_dosage), "markers x", ncol(est_dosage), "progeny\n")
cat("Dosage value range:", round(range(est_dosage, na.rm = TRUE), 2), "\n")
cat("Wrote:", file.path(OUT, "estimated_dosage_matrix.rds"), "\n")
