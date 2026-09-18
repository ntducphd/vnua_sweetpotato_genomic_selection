# Single-command reproducibility runner for the Stage-3 genomic-selection
# pipeline.
#
# Runs every numbered script (01-07) in a fresh Rscript subprocess each, in
# order, from this file's own directory. A fresh subprocess per script, rather
# than source()-ing them into one session, prevents variable and state leakage
# between scripts, the same pattern used by the authors' companion germplasm
# workflow (github.com/ntducphd/vnua_sweetpotato_improvement).
#
# Run on Windows via PowerShell, not Git Bash: Rscript under Git Bash/MSYS2
# is known to segfault on the PCA/SVD-heavy calls AGHmatrix's G-matrix
# construction and sommer's mixed-model fitting rely on internally.
#
# Usage:  Rscript run_all.R
# Exit status is non-zero if any script fails.

scripts <- c(
  "01_simulate_hexaploid_population.R",
  "02_dosage_genotype_calling.R",
  "03_qc_relationship_matrix.R",
  "04_gblup_cross_validation.R",
  "05_stage1_classical_baseline.R",
  "06_stage_gate_decision.R",
  "07_diagnostic_figure.R"
)

here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(here) == 0 || !nzchar(here)) here <- getwd()
# each script resolves its inputs and outputs relative to R/, so run from here
# whether this file was invoked as `Rscript run_all.R` or `Rscript R/run_all.R`
setwd(here)

failed <- character(0)
for (s in scripts) {
  cat(sprintf("\n===== %s =====\n", s))
  status <- system2("Rscript", shQuote(s), stdout = "", stderr = "", wait = TRUE)
  if (!identical(status, 0L)) {
    cat(sprintf("FAILED (exit %s): %s\n", status, s))
    failed <- c(failed, s)
  }
}

cat("\n============================================\n")
if (length(failed) == 0) {
  cat("run_all.R: ALL", length(scripts), "scripts completed successfully.\n")
} else {
  cat("run_all.R: FAILED scripts:\n")
  cat(paste(" -", failed, collapse = "\n"), "\n")
  quit(status = 1)
}
