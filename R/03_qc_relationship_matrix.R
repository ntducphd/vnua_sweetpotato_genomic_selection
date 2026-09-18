# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('AGHmatrix', 'dplyr')

# 03. What this script does -----------------------------------------------------
# QC-filters polyRAD's estimated dosage matrix (02) and builds a hexaploid
# additive genomic relationship (G) matrix with AGHmatrix, a CRAN package with
# autopolyploid support (VanRaden method extended to arbitrary ploidy). Also
# runs a duplicate and near-identical-individual check, the genomic-marker
# analogue of the SSR-fingerprinting congruence check in the authors' companion
# germplasm workflow (https://github.com/ntducphd/vnua_sweetpotato_improvement):
# the same purpose, catching mislabelled or duplicate genotypes before they
# enter downstream selection decisions, with a different marker system.

IN  <- "../data_simulated"
OUT <- "../data_simulated"
PLOIDY <- 6

est_dosage <- readRDS(file.path(IN, "estimated_dosage_matrix.rds"))  # markers x progeny
mm <- readRDS(file.path(IN, "marker_map.rds"))

# 04. QC filtering ------------------------------------------------------------------
# Missing-rate and MAF filtering: drop markers with >20% missing genotype
# calls, and markers with minor allele frequency < 0.05 among the progeny
# (uninformative for within-family prediction).
# On the committed simulated population the missing-rate filter is not a
# formality. 02 leaves 675 of the 1,500 loci with no dosage call for either
# allele and 75 more called only on the reference allele, so 750 markers reach
# this script entirely NA and are removed here; a further 67 go on MAF, leaving
# 683. That is where the 1,500 -> 683 reduction comes from, and it is a property
# of how polyRAD's mapping pipeline is configured in 02 (see its
# freqAllowedDeviation note), not of MAF filtering a complete matrix.
miss_rate <- rowMeans(is.na(est_dosage))
maf <- pmin(rowMeans(est_dosage, na.rm = TRUE) / PLOIDY,
            1 - rowMeans(est_dosage, na.rm = TRUE) / PLOIDY)

keep <- miss_rate <= 0.20 & maf >= 0.05
cat(sprintf("QC: %d / %d markers retained (missing<=20%%, MAF>=0.05)\n",
            sum(keep), length(keep)))

dosage_qc <- est_dosage[keep, , drop = FALSE]
# Defensive clip to the valid [0, ploidy] range. It is inert on output from
# GetWeightedMeanGenotypes(minval = 0, maxval = PLOIDY), which rescales into
# that interval, and is kept for a real-data loader that might supply dosage
# from another caller.
dosage_qc[dosage_qc < 0] <- 0
dosage_qc[dosage_qc > PLOIDY] <- PLOIDY

# 05. Duplicate / near-identical individual check -----------------------------------
# Genomic analogue of the project's SSR-based identity check: flag progeny pairs
# whose marker dosage correlation exceeds 0.98 (near-clonal identity, which in a
# real dataset would indicate accidental replication or a labelling error). None
# is expected here: the 315 progeny are full sibs of one biparental cross, drawn
# independently of each other, with no duplication introduced.
# The correlation is computed on marker-CENTRED dosage. On raw 0..6 dosage all
# individuals share the same marker-frequency profile, which compresses every
# pair into r = 0.88 to 0.94 on this population and leaves only a narrow margin
# below the threshold, so anything short of a near-perfect duplicate passes.
# Centring each marker removes that shared profile: the off-diagonal then runs
# from -0.30 to 0.29 here, while an injected exact clone still scores 1.00.
dosage_centred <- dosage_qc - rowMeans(dosage_qc)
cor_mat <- cor(dosage_centred, use = "pairwise.complete.obs")
diag(cor_mat) <- NA
dup_pairs <- which(cor_mat > 0.98, arr.ind = TRUE)
dup_pairs <- dup_pairs[dup_pairs[, 1] < dup_pairs[, 2], , drop = FALSE]
if (nrow(dup_pairs) > 0) {
  cat(sprintf("WARNING: %d near-identical genotype pair(s) found (r > 0.98):\n", nrow(dup_pairs)))
  for (k in seq_len(nrow(dup_pairs))) {
    cat(" ", rownames(cor_mat)[dup_pairs[k, 1]], "<->", colnames(cor_mat)[dup_pairs[k, 2]], "\n")
  }
} else {
  cat("Duplicate/identity check: no near-identical genotype pairs found.\n")
}

# 06. Build the hexaploid additive genomic relationship (G) matrix ------------------
# AGHmatrix expects individuals in rows, markers in columns, dosage coded
# 0..ploidy, the transpose of the marker x progeny convention used here.
SNPmatrix <- t(dosage_qc)

Gmat <- Gmatrix(SNPmatrix, method = "VanRaden", ploidy = PLOIDY,
                maf = 0, missingValue = NA, integer = FALSE)

cat("G-matrix dimensions:", dim(Gmat), "\n")
cat("G-matrix diagonal range (expect centered around 1):", round(range(diag(Gmat)), 3), "\n")
cat("G-matrix off-diagonal range:", round(range(Gmat[upper.tri(Gmat)]), 3), "\n")
stopifnot(nrow(Gmat) == ncol(dosage_qc), all(rownames(Gmat) == colnames(dosage_qc)))

# 07. Write output --------------------------------------------------------------------
saveRDS(Gmat, file.path(OUT, "G_matrix_hexaploid.rds"))
saveRDS(dosage_qc, file.path(OUT, "dosage_qc.rds"))
cat("Wrote:", file.path(OUT, "G_matrix_hexaploid.rds"), "\n")
