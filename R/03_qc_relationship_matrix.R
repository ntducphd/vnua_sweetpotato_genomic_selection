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
# Missing-rate and MAF filtering (kept general/real-data-ready even though
# this simulator produces no missing calls): drop markers with >20% missing
# genotype calls, and markers with minor allele frequency < 0.05 among the
# progeny (uninformative for within-family prediction).
miss_rate <- rowMeans(is.na(est_dosage))
maf <- pmin(rowMeans(est_dosage, na.rm = TRUE) / PLOIDY,
            1 - rowMeans(est_dosage, na.rm = TRUE) / PLOIDY)

keep <- miss_rate <= 0.20 & maf >= 0.05
cat(sprintf("QC: %d / %d markers retained (missing<=20%%, MAF>=0.05)\n",
            sum(keep), length(keep)))

dosage_qc <- est_dosage[keep, , drop = FALSE]
# clip to the valid [0, ploidy] range (polyRAD's posterior mean can slightly
# undershoot/overshoot near the boundary for very confident calls)
dosage_qc[dosage_qc < 0] <- 0
dosage_qc[dosage_qc > PLOIDY] <- PLOIDY

# 05. Duplicate / near-identical individual check -----------------------------------
# Genomic analogue of the project's SSR-based identity check: flag progeny pairs
# whose marker dosage correlation exceeds 0.98 (near-clonal identity, which in a
# real dataset would indicate accidental replication or a labelling error; none
# is expected here, since the progeny are simulated independently).
cor_mat <- cor(dosage_qc, use = "pairwise.complete.obs")
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
