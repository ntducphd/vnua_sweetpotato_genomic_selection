# Stage-3 Dosage-Aware Genomic Selection for Sweetpotato

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22822689.svg)](https://doi.org/10.5281/zenodo.22822689)

Reference R implementation of Stage 3 ("pilot dosage-aware genomic selection") from the staged,
resource-gated AI-adoption roadmap proposed in:

> Nguyen Trung Duc, Pham Quang Tuan, Nguyen Thi Thu, Doan Thu Thuy, Vu Thi Xuan Binh, Nguyen Van
> Loc. **AI-assisted genomic selection in sweetpotato: A staged roadmap for polyploid,
> clonal-propagated breeding.** *Ecological Genetics and Genomics* **41** (2026) 100530.
> https://doi.org/10.1016/j.egg.2026.100530

Sweetpotato (*Ipomoea batatas*, 2n = 6x = 90) is hexaploid, clonally propagated, and has almost no
genomic-selection track record. The article's Table 5 (Section 7) identifies one existing
feasibility study (Gemenet et al. 2020, *Theoretical and Applied Genetics* 133(12), 3345–3363) and
sets an explicit, statistically rigorous go/no-go gate for when a breeding programme should trust
genomic prediction over its existing classical multi-trait index. **No genotyping-by-sequencing data
exists for the authors' own sweetpotato panel** (only phenotypic MGIDI data and SSR fingerprinting
do). This repository is a complete, runnable pipeline validated end to end on a realistic simulated
hexaploid population, so the gate logic and every upstream method (dosage-uncertain genotype
calling, polyploid relationship-matrix construction, dosage-aware GBLUP, bootstrap-CI gate testing)
are ready to use as-is the moment real genotyping data becomes available; only the data-loading
step changes.

## What this workflow does

1. **Simulate** a 315-progeny biparental hexaploid population (matching the population size of the
   one existing feasibility study, Gemenet et al. 2020), with a realistic multi-trait (5-trait) QTL
   architecture.
2. **Call dosage genotypes** from simulated noisy GBS read counts using
   [polyRAD](https://github.com/lvclark/polyRAD) (Clark, Lipka & Sacks 2019, *G3* 9(3), 663-673),
   the precedent tool the article cites for this step.
3. **Build a hexaploid genomic relationship (G) matrix** with
   [AGHmatrix](https://github.com/prmunoz/AGHmatrix)'s polyploid-extended VanRaden method, after QC
   filtering and a duplicate/identity check, the genomic-marker analogue of the SSR-fingerprinting
   congruence check in the authors' companion germplasm workflow
   ([vnua_sweetpotato_improvement](https://github.com/ntducphd/vnua_sweetpotato_improvement)).
4. **Fit dosage-aware GBLUP across all five traits**, each with its own 5-fold cross-validation on
   one shared fold assignment, using [sommer](https://github.com/covaruber/sommer), then combine the
   five out-of-fold GEBVs per genotype into one composite predicted-merit score (standardise each
   trait, then average). All five traits are used, matched in scope to the Stage-1 index below.
5. **Compute the Stage-1 classical baseline** (MGIDI, via
   [metan](https://github.com/TiagoOlivoto/metan)) on the identical simulated phenotype data, the
   same MGIDI method the authors' companion germplasm workflow uses, and combine the same five
   observed traits into a composite *observed* score with the identical standardise-then-average
   formula.
6. **Validate two ways**: (A) an **oracle check**, both methods' scores correlated against the
   TRUE, noise-free simulated breeding value, possible only in silico, which confirms the pipeline
   recovers real genetic signal rather than an artifact of how the composite scores are built; and
   (B) **Table 5's literal operational gate**, both methods' scores correlated against the
   *observed* composite outcome (the only thing a real programme ever has), bootstrapped for a 95%
   CI on the accuracy difference, verdict = **GO** only if the CI excludes zero. See
   [Two-tier validation](#two-tier-validation-why-two-numbers) below for why both are reported.
7. **Render a diagnostic figure** summarising both validations
   (`figures/Fig_stage3_gate_diagnostic.png`). This figure belongs to the repository; it is not one
   of the article's figures.

The method, criteria and thresholds are set out in full in the article's Table 5 and its
surrounding Section 7 text; this repository is the executable version of that description.

## Two-tier validation: why two numbers

Stage-3 GBLUP and the Stage-1 MGIDI index are scored on the same five traits and compared against
the same shared outcome, in the comparison mode Section 7 specifies: "a specific target trait, or a
pre-specified composite performance score computed identically for every genotype". Two different
targets are available in a simulated population, and this repository reports both, because they
answer different questions.

MGIDI and a composite built by standardising each observed trait and averaging are both near-linear
combinations of the same five observed traits. A high correlation between them is therefore close
to a mathematical certainty and says little about predictive skill. Section 7 of the article
states this construction pitfall. The simulation makes a cleaner target available: the
true, noise-free breeding value, which no real breeding programme ever has.

- **(A) Oracle validation** (`06`, panel A of the figure): both methods' scores correlated against
  the TRUE simulated breeding value. This is the artifact-free check that the pipeline mechanics
  recover real genetic signal. **GBLUP 0.664 vs. MGIDI 0.441**, both below the **0.668** noise
  ceiling (the observed composite's own correlation with the true value); the 95% CI on the
  difference, **[0.117, 0.335], excludes zero**. With a negative genetic trade-off present
  (Yield/Carotenoid true breeding-value correlation r = −0.37), marker-based prediction
  outperforms a phenotype-only classical index on this population.
- **(B) Table 5's literal operational gate** (`06`, panel B): both methods' scores correlated
  against the *observed* composite outcome, bootstrapped for the go/no-go verdict, because this is
  the only comparison a real breeding programme, which never has a true breeding value to fall
  back on, can actually compute. **GBLUP 0.482 vs. MGIDI 0.561**, 95% CI on the difference
  **[−0.180, 0.022]**, spans zero, **NO-GO**. One limitation is structural: `metan::mgidi()` has no
  mechanism to score a genotype blind to its own data the way GBLUP's held-out-fold prediction
  can, so this number favours MGIDI beyond genuine predictive skill.

Read together, (A) and (B) are the practical lesson the article draws in Section 7: a real,
statistically detectable advantage for genomic prediction exists in these data, while the
operational test a resource-limited programme can actually run, on noisy observed phenotypes and
with MGIDI's in-sample advantage, cannot confirm it at this population size. That gap is the
Type-II-error risk Table 5's own text warns a small validation set could create.

The article states both comparisons qualitatively (Section 7 and the Fig. 1 caption); the numeric
values above are this repository's own output on the committed simulated population, reproducible
with the scripts in `R/`.

## What's real vs. simulated

| Component | Status |
|---|---|
| Population size (315 progeny), biparental hexaploid design | Real precedent: matches Gemenet et al. (2020), *Theor. Appl. Genet.* 133(12), the one existing sweetpotato genomic-selection feasibility study |
| Marker count (1,500), linkage-group count (15) | Matches sweetpotato's real basic chromosome number (x = 15); GBS-scale marker density |
| Progeny genotypes, GBS read counts, phenotypes | **Fully simulated**, see [Known limitations](#known-limitations) |
| polyRAD dosage-calling pipeline | Real tool, real API, run on simulated read counts |
| AGHmatrix G-matrix construction | Real tool/method (VanRaden, polyploid-extended), run on simulated dosage |
| sommer GBLUP + cross-validation | Real tool/method, run on simulated data |
| metan/MGIDI Stage-1 baseline | The same real method as the authors' companion germplasm workflow, applied here for a fair Stage-1-vs-Stage-3 comparison |
| Bootstrap gate decision | Real statistical procedure, implements Table 5's literal wording |

Every method and tool is used exactly as it would be on a real dataset; only the input
genotypes and phenotypes are synthetic. Reported numbers (dosage-calling accuracy, cross-validation
accuracy, gate verdict) are pipeline-validation figures, not claims about real sweetpotato genetics
or breeding-value predictability.

## Table 5 correspondence

| Table 5 (Stage-3 row) | Implemented as |
|---|---|
| Entry requirement: one biparental population, deep/accurate genotyping (benchmark: 315 progeny in the one existing feasibility study) | `R/01_simulate_hexaploid_population.R` |
| Precedent tool for dosage-uncertain genotype calling: polyRAD (Clark et al. 2019) | `R/02_dosage_genotype_calling.R` |
| Precedent tool family for dosage-aware prediction: polyGBLUP-family models | `R/03_qc_relationship_matrix.R` + `R/04_gblup_cross_validation.R` |
| Gate: prediction accuracy exceeds the Stage-1 classical-index ranking on a held-out validation subset, with the improvement's confidence interval excluding zero | `R/06_stage_gate_decision.R`: 10,000-resample bootstrap, 95% CI on (GBLUP − Stage-1) accuracy; verdict = GO only if the CI's lower bound > 0 |

Run on the committed simulated population (see [Two-tier validation](#two-tier-validation-why-two-numbers)
above for why two comparisons are reported):

- **Oracle validation** (vs. TRUE breeding value): GBLUP 0.664 vs. Stage-1 0.441, 95% CI on the
  difference **[0.117, 0.335]**, excludes zero; genomic prediction is detectably better here, both
  below the 0.668 noise ceiling.
- **Table 5's literal gate** (vs. observed composite): GBLUP 0.482 vs. Stage-1 0.561, 95% CI
  **[−0.180, 0.022]**, spans zero, **NO-GO (not distinguishable)**, read alongside the caveat that
  MGIDI is not held out in this comparison (above). The contrast with the oracle result is the
  point: a real, detectable advantage that the operational, noise-limited test still cannot confirm.

Both are visualised in `figures/Fig_stage3_gate_diagnostic.png`. Neither should be read as a claim
about real sweetpotato breeding-value predictability; see [Known limitations](#known-limitations).

## Repository structure

```
.
├── R/                    7 numbered scripts + run_all.R
├── data_simulated/       Simulated genotypes, phenotypes, dosage calls, G-matrix (01-03 output)
├── outputs/              Cross-validation and gate-decision tables (04-06 output)
├── figures/              Diagnostic figure, PNG + PDF (07 output)
├── environment/          R version and exact package versions used for these results
├── LICENSE               MIT
└── CITATION.cff          Machine-readable citation metadata
```

`data_simulated/`, `outputs/`, and `figures/` are included pre-generated so the results can be
inspected without running R, but every file in all three is fully reproducible via the scripts in
`R/`.

## Requirements

- **R ≥ 4.5** (developed and tested on R 4.5.1)
- Run on **Windows via PowerShell, not Git Bash**. `Rscript` under Git Bash/MSYS2 is confirmed to
  segfault on the SVD-family linear algebra that AGHmatrix's G-matrix construction and sommer's
  mixed-model fitting rely on internally. Not an issue on macOS/Linux or native Windows
  PowerShell/cmd.
- R packages (exact versions used for the results in this repository, see
  `environment/R_package_versions.csv`):

  | Package | Version | Package | Version |
  |---|---|---|---|
  | polyRAD | 2.0.1 | ggplot2 | 4.0.3 |
  | AGHmatrix | 3.0.1 | patchwork | 1.3.2 |
  | sommer | 4.4.6 | dplyr | 1.2.1 |
  | metan | 1.19.0 | tidyr | 1.3.2 |
  | purrr | 1.2.2 | easypackages | 0.1.0 |

  Install with:
  ```r
  install.packages(c("polyRAD", "AGHmatrix", "sommer", "metan", "dplyr", "tidyr",
                      "purrr", "ggplot2", "patchwork", "easypackages"))
  ```

  `install.packages()` fetches whatever CRAN currently serves, which may be newer than the versions
  above. To reproduce the committed results exactly, install the pinned versions of the four
  analysis packages:
  ```r
  remotes::install_version("polyRAD",   "2.0.1")
  remotes::install_version("AGHmatrix", "3.0.1")
  remotes::install_version("sommer",    "4.4.6")
  remotes::install_version("metan",     "1.19.0")
  ```

Full session details (platform, locale, complete dependency tree) are recorded in
`environment/R_sessionInfo.txt`.

## How to reproduce

From the `R/` directory, or from the repository root:

```powershell
Rscript run_all.R          # from R/
Rscript R/run_all.R        # from the repository root
```

This runs all 7 scripts in order, each in its own fresh `Rscript` subprocess, and reports which (if
any) fail. The whole pipeline takes under a minute on a desktop CPU.
The pipeline is deterministic: scripts 01, 02, 04 and 06 set a fixed seed and scripts 03, 05 and 07
draw no random numbers, so a clean run regenerates every committed CSV and `.rds` byte for byte. The
one exception is `figures/Fig_stage3_gate_diagnostic.pdf`, which embeds a creation timestamp in its
document metadata and therefore never byte-matches across runs even when every plotted value is
identical.

Scripts can also be run individually, in order (01 to 07); each reads the previous stage's saved
`.rds` output from `../data_simulated/` or `../outputs/`.

## Script → output map

| Script | Produces | Description |
|---|---|---|
| `01_simulate_hexaploid_population.R` | `data_simulated/progeny_dosage_matrix.rds`, `marker_map.rds`, `plot_level_phenotypes.csv`, `true_qtl_effects.rds`, `true_breeding_value_matrix.rds` | 315-progeny biparental hexaploid population, hexasomic segregation, 5-trait QTL-driven phenotypes, plus the pre-noise true breeding value (simulation ground truth, for oracle validation only) |
| `02_dosage_genotype_calling.R` | `data_simulated/estimated_dosage_matrix.rds` | polyRAD dosage calling from simulated GBS read counts |
| `03_qc_relationship_matrix.R` | `data_simulated/dosage_qc.rds`, `G_matrix_hexaploid.rds` | Missing-rate (≤ 20%) and MAF (≥ 0.05) filtering, duplicate check, AGHmatrix hexaploid G-matrix. On the committed simulated population this retains 683 of 1,500 markers, and the G-matrix is built on those 683 |
| `04_gblup_cross_validation.R` | `outputs/gblup_cv_results.rds`, `gblup_cv_accuracy_by_fold.csv` | sommer dosage-aware GBLUP, all 5 traits, one shared 5-fold assignment (25 fits, convergence checked and reported), out-of-fold composite score |
| `05_stage1_classical_baseline.R` | `outputs/stage1_mgidi_results.rds`, `stage1_mgidi_scores.csv` | metan/MGIDI baseline on the same simulated phenotypes, correlated against the shared observed composite |
| `06_stage_gate_decision.R` | `outputs/stage_gate_decision.rds`, `stage_gate_comparison_table.csv` | Two-tier validation: (A) oracle accuracy vs. TRUE breeding value; (B) Table 5's literal bootstrap-CI gate vs. observed composite |
| `07_diagnostic_figure.R` | `figures/Fig_stage3_gate_diagnostic.png` (+ `.pdf`) | Two-panel summary: (A) predicted vs. TRUE breeding value, both methods; (B) bootstrap CI for Table 5's literal gate |

## Swapping in real GBS data

Once real genotyping data exists for the sweetpotato panel, replace scripts 01–02 with a loader
that produces the same two inputs script 03 expects:

1. **Reference/alternate read-depth counts** (or already-called dosage) for polyRAD
   (script 02's `RADdata()` call). A real GBS pipeline (e.g. TASSEL-GBS, Stacks) produces these in
   the same taxa × allele shape. Skip script 01 entirely; skip script 02's synthetic read-count
   simulation but keep its `RADdata()` → `SetDonorParent()`/`SetRecurrentParent()` →
   `PipelineMapping2Parents()` → `GetWeightedMeanGenotypes()` structure, pointed at the real reads.
2. **`plot_level_phenotypes.csv`**: a real trial dataset with columns `GEN`, `REP`, and one column
   per trait, in the layout scripts 04 and 05 expect.

Scripts 03 to 07 need no changes; they operate on dosage matrices and phenotype tables in the same
shape regardless of whether the dosage came from simulation or real GBS calling. Re-check
`freqAllowedDeviation` (script 02) and the missing-rate/MAF QC thresholds (script 03) against the
real dataset's actual depth/missingness profile before trusting the defaults carried over from this
simulation.

## Known limitations

- **Inheritance model** (script 01): random hexasomic (polysomic) segregation, in which each parent's 6
  homologous copies at a linkage group are equally likely to be inherited, no preferential pairing,
  no within-linkage-group recombination. Real sweetpotato meiosis "does not fit auto- or
  allopolyploid models cleanly" (Gao et al. 2020, *PLoS ONE* 15(3), e0229624), so this simulator
  validates pipeline *mechanics*, not real sweetpotato transmission genetics.
- **Parental read counts** (script 02): true parental dosage is not carried over from script 01, so
  parent read counts are reconstructed from the population allele frequency rather than the
  parents' own actual genotypes. A real dataset would genotype the actual parents directly.
- **Trait architecture** (script 01): 40 QTL/trait (60 for Yield/Carotenoid and DryMatter/
  StorageScore, which additionally share 20 QTL each with a controlled sign relationship, realised
  at r = −0.37 and r = +0.32 respectively), additive-only, no dominance or epistasis. Heritability
  is targeted at h² = 0.40 on a genotype-mean basis, with the residual-noise scale accounting for
  the 2 field replicates being averaged (see script 01 for the formula), not derived from real
  sweetpotato trait-heritability estimates. The trade-off and association signs (Yield vs.
  Carotenoid negative, DryMatter vs. StorageScore positive) are a simplified illustration inspired
  by the β-carotene/starch trade-off reported by Gemenet et al. (2020, *Theor. Appl. Genet.* 133(1),
  23–36) and catalogued in Table 1 of the article, not a claim about the real magnitude or even
  direction of these relationships in sweetpotato.
- **Stage-1/Stage-3 in-sample asymmetry in comparison (B)** (see script 06's own comments and the
  [Two-tier validation](#two-tier-validation-why-two-numbers) section above): MGIDI is computed once
  on the full observed dataset, since it has no mechanism to extrapolate to a genotype it was not given
  phenotypes for, and `metan::mgidi()` exposes no built-in way to score a genotype "blind" to its
  own data the way GBLUP's held-out-fold prediction can. Comparison (A), the oracle validation
  against the true breeding value, does not have this problem (both methods' scores are computed the
  same way regardless; only the *target* they are compared against changes), which is why (A) is the
  more informative number, while (B) is what Table 5 literally asks a real deployment, which never
  has a true breeding value, to compute.
- **Independent per-trait GBLUP, not a joint multi-trait model** (script 04): the five traits are
  each fitted with their own univariate GBLUP (on one shared fold assignment) and combined
  mechanically into a composite score, rather than a single joint mixed model with an estimated
  cross-trait genetic covariance structure. A true multi-trait GBLUP would let traits borrow
  statistical strength from each other's genetic correlations and is a natural next refinement; the
  univariate-then-combine design used here keeps the two methods' predictive scope matched, which is
  what the Stage-1-vs-Stage-3 comparison requires.
- **GBLUP fitted on genotype-mean phenotypes, not raw plot-level records** (script 04): a standard
  two-step approach (collapse the 2 replicates to a genotype mean first, then fit the genomic
  model), matching the same simplification 05's MGIDI baseline makes by construction. A one-step
  model fitting variance components jointly from the raw plot data would be more statistically
  precise but is not implemented here.
- Numeric results (dosage-calling correlation, cross-validation accuracy, both gate verdicts)
  describe this simulated population only and should not be read as predictions about real
  sweetpotato breeding populations.

## Citation

If you use this code, please cite the article:

> Nguyen Trung Duc, Pham Quang Tuan, Nguyen Thi Thu, Doan Thu Thuy, Vu Thi Xuan Binh, Nguyen Van
> Loc. AI-assisted genomic selection in sweetpotato: A staged roadmap for polyploid,
> clonal-propagated breeding. *Ecological Genetics and Genomics* 41 (2026) 100530.
> https://doi.org/10.1016/j.egg.2026.100530

To cite this software release specifically:

> Nguyen Trung Duc, Pham Quang Tuan, Nguyen Thi Thu, Doan Thu Thuy, Vu Thi Xuan Binh, Nguyen Van
> Loc. Stage-3 Dosage-Aware Genomic Selection Pipeline for Hexaploid Sweetpotato: Simulation-Validated
> Reference Implementation (v1.0.0). Zenodo, 2026. https://doi.org/10.5281/zenodo.22822690

Machine-readable metadata is in `CITATION.cff`.

## License

Released under the [MIT License](LICENSE).

## Contact

Nguyen Trung Duc, Vietnam National University of Agriculture, ntduc@vnua.edu.vn
