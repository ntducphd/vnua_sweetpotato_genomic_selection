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
sets an explicit, provisional go/no-go gate for when a breeding programme should trust genomic
prediction over its existing classical multi-trait index. No genotyping-by-sequencing data yet
exists for the authors' own sweetpotato panel (only phenotypic MGIDI data and SSR fingerprinting
do). This repository validates the pipeline end to end on a realistic simulated hexaploid
population, so the gate logic and every upstream method (dosage-uncertain genotype calling,
polyploid relationship-matrix construction, dosage-aware GBLUP, bootstrap-CI gate testing) run
unchanged the moment real genotyping data becomes available; only the data-loading step changes.

## What this workflow does

1. **Simulate** a 315-progeny biparental hexaploid population (matching the population size of the
   one existing feasibility study, Gemenet et al. 2020), with a realistic multi-trait (5-trait) QTL
   architecture.
2. **Call dosage genotypes** from simulated noisy GBS read counts using
   [polyRAD](https://github.com/lvclark/polyRAD) (Clark, Lipka & Sacks 2019, *G3* 9(3), 663–673),
   the precedent tool the article cites for this step.
3. **Build a hexaploid genomic relationship (G) matrix** with
   [AGHmatrix](https://github.com/rramadeu/AGHmatrix)'s polyploid-extended VanRaden method, after QC
   filtering and a duplicate/identity check, the genomic-marker analogue of the SSR-fingerprinting
   congruence check in the authors' companion germplasm workflow
   ([vnua_sweetpotato_improvement](https://github.com/ntducphd/vnua_sweetpotato_improvement)).
4. **Fit dosage-aware GBLUP across all five traits**, each with its own 5-fold cross-validation on
   one shared fold assignment, using [sommer](https://github.com/covaruber/sommer), then combine the
   five out-of-fold GEBVs per genotype into one composite predicted-merit score (standardise each
   trait, then average). All five traits are used, matched in scope to the Stage-1 index below.
5. **Compute the Stage-1 classical baseline** (MGIDI, via
   [metan](https://github.com/nepem-ufsc/metan)) on the identical simulated phenotype data, the
   same MGIDI method the authors' companion germplasm workflow uses, and combine the same five
   observed traits into a composite *observed* score with the identical standardise-then-average
   formula.
6. **Validate two ways**: (A) an **oracle check**, both methods' scores correlated against the
   TRUE, noise-free simulated breeding value, possible only in silico, which confirms the pipeline
   recovers real genetic signal rather than an artifact of how the composite scores are built; and
   (B) **Table 5's literal operational gate**, both methods' scores correlated against the
   *observed* composite outcome (the only thing a real programme ever has), bootstrapped for a 95%
   CI on the accuracy difference, verdict = **GO** only if the CI excludes zero. Both comparisons
   are defined under [Two-tier validation](#two-tier-validation).
7. **Render a diagnostic figure** summarising both validations
   (`figures/Fig_stage3_gate_diagnostic.png`), a repository figure rather than one of the
   article's.

The method, criteria and thresholds are set out in full in the article's Table 5 and its
surrounding Section 7 text; this repository is the executable version of that description.

## Two-tier validation

Stage-3 GBLUP and the Stage-1 MGIDI index are scored on the same five traits and compared against
the same shared outcome, in the comparison mode Section 7 specifies: "a specific target trait, or a
pre-specified composite performance score computed identically for every genotype". Two targets are
available in a simulated population, and both are reported.

MGIDI and a composite built by standardising each observed trait and averaging are both functions
of the same five observed traits, so a correlation between them partly reflects that shared
construction rather than predictive skill. Section 7 of the article describes this pitfall. The
simulation makes a cleaner target available: the true, noise-free breeding value, which no real
breeding programme ever has.

- **(A) Oracle validation** (`06`, panel A of the figure): both methods' scores correlated against
  the TRUE simulated breeding value, a target free of the construction pitfall above. **GBLUP 0.664
  vs. MGIDI 0.441**, both below the **0.668** noise ceiling. That ceiling is the observed
  composite's own correlation with the true value, a measure of how much true signal survives this
  population's phenotypic noise, not a theoretical bound: the equal-weight observed composite in
  fact edges out the out-of-fold GBLUP composite here. The 95% CI on the difference between the two
  methods, **[0.117, 0.335], excludes zero**. With a negative genetic trade-off present
  (Yield/Carotenoid true breeding-value correlation r = −0.37), marker-based prediction outperforms
  a phenotype-only classical index on this population.
- **(B) Table 5's literal operational gate** (`06`, panel B): both methods' scores correlated
  against the *observed* composite outcome, bootstrapped for the go/no-go verdict, because this is
  the only comparison a real breeding programme, which never has a true breeding value to fall
  back on, can actually compute. **GBLUP 0.482 vs. MGIDI 0.561**, 95% CI on the difference
  **[−0.180, 0.022]**, spans zero, **NO-GO**. One limitation is structural: `metan::mgidi()` has no
  mechanism to score a genotype blind to its own data the way GBLUP's held-out-fold prediction
  can, so this number favours MGIDI beyond genuine predictive skill.

(A) and (B) together are the practical lesson the article draws in Section 7: a real,
statistically detectable advantage for genomic prediction exists in these data, while the
operational test a resource-limited programme can actually run, on noisy observed phenotypes and
with MGIDI's in-sample advantage, cannot confirm it at this population size. That gap is a concrete
case of the validation-set power limitation Table 5's own text warns about.

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
| metan/MGIDI Stage-1 baseline | The same method as the authors' companion germplasm workflow, matched in scope to the Stage-3 composite |
| Bootstrap gate decision | Real statistical procedure, implements Table 5's literal wording |

Every method and tool is the real implementation called through its public API; only the input
genotypes and phenotypes are synthetic. Three settings in script 02 are specific to this simulated
input and would be revisited on real data: `taxaPloidy = 2` alongside `possiblePloidies = list(6)`
to reach polyRAD's effective ploidy of 6, `useLinkage = FALSE` because the marker positions are
arbitrary indices rather than a real map, and parental read counts reconstructed from the
population allele frequency rather than from the parents' own genotypes.

Reported numbers (dosage-calling accuracy over the called cells, cross-validation accuracy, gate
verdict) are pipeline-validation figures, not claims about real sweetpotato genetics or
breeding-value predictability.

## Table 5 correspondence

| Table 5 (Stage-3 row) | Implemented as |
|---|---|
| Entry requirement: one biparental population, deep/accurate genotyping (benchmark: 315 progeny in the one existing feasibility study) | `R/01_simulate_hexaploid_population.R` |
| Precedent tool for dosage-uncertain genotype calling: polyRAD (Clark et al. 2019) | `R/02_dosage_genotype_calling.R` |
| Method: pilot dosage-aware genomic selection | `R/03_qc_relationship_matrix.R` + `R/04_gblup_cross_validation.R` |
| Precedent, third entry: linkage analysis and haplotype phasing for high-ploidy populations (Mollinari & Garcia 2019) | **Not implemented.** Script 02 runs `useLinkage = FALSE`, because the simulated marker positions are arbitrary indices rather than a genetic or physical map |
| Gate: prediction accuracy exceeds the Stage-1 classical-index ranking on a held-out validation subset, with the improvement's confidence interval excluding zero | `R/06_stage_gate_decision.R`: 10,000-resample bootstrap, 95% CI on (GBLUP − Stage-1) accuracy; verdict = GO only if the CI's lower bound > 0 |

Run on the committed simulated population (see [Two-tier validation](#two-tier-validation) above):

- **Oracle validation** (vs. TRUE breeding value): GBLUP 0.664 vs. Stage-1 0.441, 95% CI on the
  difference **[0.117, 0.335]**, excludes zero; genomic prediction is detectably better here, both
  below the 0.668 noise ceiling.
- **Table 5's literal gate** (vs. observed composite): GBLUP 0.482 vs. Stage-1 0.561, 95% CI
  **[−0.180, 0.022]**, spans zero, **NO-GO (not distinguishable)**, with the caveat that MGIDI is
  not held out in this comparison (above). The oracle comparison detects an advantage that this
  operational, noise-limited test does not.

Both are visualised in `figures/Fig_stage3_gate_diagnostic.png`. Neither is a claim about real
sweetpotato breeding-value predictability; see [Known limitations](#known-limitations).

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
- Run on **Windows via PowerShell, not Git Bash**. `Rscript` under Git Bash/MSYS2 segfaults on the
  SVD-family linear algebra that AGHmatrix's G-matrix construction and sommer's mixed-model fitting
  rely on internally. Not an issue on macOS/Linux or native Windows PowerShell/cmd.
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
  install.packages("remotes")
  remotes::install_version("polyRAD",   "2.0.1")
  remotes::install_version("AGHmatrix", "3.0.1")
  remotes::install_version("sommer",    "4.4.6")
  remotes::install_version("metan",     "1.19.0")
  ```

Full session details (platform, locale, complete dependency tree) are recorded in
`environment/R_sessionInfo.txt`, captured from a session that attaches every package the seven
scripts attach, in the order they attach them.

## How to reproduce

From the `R/` directory, or from the repository root:

```powershell
Rscript run_all.R          # from R/
Rscript R/run_all.R        # from the repository root
```

This runs all 7 scripts in order, each in its own fresh `Rscript` subprocess, and reports which (if
any) fail.

The pipeline is deterministic: scripts 01, 02, 04 and 06 set a fixed seed and scripts 03, 05 and 07
draw no random numbers, so a clean run regenerates every committed CSV and `.rds` byte for byte.
`figures/Fig_stage3_gate_diagnostic.pdf` embeds a creation timestamp in its document metadata, so
it does not byte-match across runs.

Scripts can also be run individually, in order (01 to 07), **from inside the `R/` directory**: each
reads what the earlier stages wrote to `../data_simulated/` or `../outputs/`. Most of those inputs
are `.rds`; scripts 04 and 05 both read `data_simulated/plot_level_phenotypes.csv`, which script 01
writes. Run them from `R/` and not from the repository root, because the paths are relative to the
working directory and from the root they resolve outside the clone. `run_all.R` is safe either way:
it sets its own working directory first.

## Script → output map

| Script | Produces | Description |
|---|---|---|
| `01_simulate_hexaploid_population.R` | `data_simulated/progeny_dosage_matrix.rds`, `marker_map.rds`, `plot_level_phenotypes.csv`, `true_qtl_effects.rds`, `true_breeding_value_matrix.rds` | 315-progeny biparental hexaploid population, hexasomic segregation, 5-trait QTL-driven phenotypes, plus the pre-noise true breeding value (simulation ground truth, for oracle validation only) |
| `02_dosage_genotype_calling.R` | `data_simulated/estimated_dosage_matrix.rds` | polyRAD dosage calling from simulated GBS read counts |
| `03_qc_relationship_matrix.R` | `data_simulated/dosage_qc.rds`, `G_matrix_hexaploid.rds` | Missing-rate (≤ 20%) and MAF (≥ 0.05) filtering, marker-centred duplicate/identity check, AGHmatrix hexaploid G-matrix. On the committed simulated population this retains 683 of 1,500 markers, and the G-matrix is built on those 683. The reduction is driven by polyRAD's call rate rather than by MAF filtering, see [Known limitations](#known-limitations) |
| `04_gblup_cross_validation.R` | `outputs/gblup_cv_results.rds`, `gblup_cv_accuracy_by_fold.csv` | sommer dosage-aware GBLUP, all 5 traits, one shared 5-fold assignment (25 fits, convergence checked and reported), out-of-fold composite score |
| `05_stage1_classical_baseline.R` | `outputs/stage1_mgidi_results.rds`, `stage1_mgidi_scores.csv` | metan/MGIDI baseline on the same simulated phenotypes, correlated against the shared observed composite |
| `06_stage_gate_decision.R` | `outputs/stage_gate_decision.rds`, `stage_gate_comparison_table.csv` | Two-tier validation: (A) oracle accuracy vs. TRUE breeding value; (B) Table 5's literal bootstrap-CI gate vs. observed composite |
| `07_diagnostic_figure.R` | `figures/Fig_stage3_gate_diagnostic.png` (+ `.pdf`) | Two-panel summary: (A) predicted vs. TRUE breeding value, both methods; (B) bootstrap CI for Table 5's literal gate |

## Swapping in real GBS data

Once real genotyping data exists for the sweetpotato panel, replace script 01's simulation with a
loader that supplies three inputs:

1. **Reference/alternate read-depth counts** (or already-called dosage) for polyRAD
   (script 02's `RADdata()` call). A real GBS pipeline (e.g. TASSEL-GBS, Stacks) produces these in
   the same taxa × allele shape. Keep script 02's `RADdata()` →
   `SetDonorParent()`/`SetRecurrentParent()` → `PipelineMapping2Parents()` →
   `GetWeightedMeanGenotypes()` structure, pointed at the real reads, and drop its synthetic
   read-count simulation.
2. **A marker map** in the shape script 01 writes (`marker_id`, `lg_id`, `allele_freq`), saved to
   `data_simulated/marker_map.rds`. Scripts 02 and 03 both read it, and with real map positions
   script 02 can run with `useLinkage = TRUE` instead of the `FALSE` this simulation requires.
3. **`plot_level_phenotypes.csv`**: a real trial dataset with columns `GEN`, `REP`, and one column
   per trait, in the layout scripts 04 and 05 expect.

Scripts 03 to 05 then carry over unchanged: they operate on a dosage matrix and a plot-level
phenotype table of the same shape whether the dosage came from simulation or from real GBS calling.
Script 06's oracle validation (A), and panel A of script 07 that plots it, are the exception, and
not for an implementation reason: they score each method against the pre-noise true breeding value,
which exists only because the population is simulated. On real data that target does not exist, so
(A) has no counterpart and only Table 5's literal gate (B) can be computed.

Re-check `freqAllowedDeviation` (script 02) and the missing-rate/MAF QC thresholds (script 03)
against the real dataset's actual depth and missingness profile before trusting the values carried
over from this simulation. On this simulated input those settings leave half the loci without a
dosage call, see [Known limitations](#known-limitations).

Two further changes are needed if the real data is a **trial series rather than one experiment**,
both required by the article's own Section 7.3 and neither implemented here, because this simulated
population has a single environment:

- **Compute the gate on across-environment genotype BLUPs, not per-environment means**, so that the
  value being gated is the one the programme would actually select on. The article points to a
  two-stage analysis developed for exactly the sparse, unevenly replicated design sweetpotato trial
  series usually have as the route from raw plot data to those BLUPs.
- **Partition the held-out set by environment or by year, never at random within the pooled
  dataset.** Script 04 assigns folds at random (`sample(rep(seq_len(N_FOLDS), ...))`), which is
  correct for one experiment but not for a pooled series: a random split leaves records from every
  environment on both sides, so a model can reach the threshold by learning environment means
  rather than genotype differences and the gate passes for the wrong reason.

Section 7.3 also asks for the pooled estimate to be the primary gate with the per-environment
values read as a dispersion check. None of this changes Table 5's thresholds; it changes what they
are computed on.

## Known limitations

- **Inheritance model** (script 01): random hexasomic (polysomic) segregation, in which each
  parent's 6 homologous copies at a linkage group are equally likely to be inherited, with no
  preferential pairing and no within-linkage-group recombination. The article records that real
  sweetpotato meiotic behaviour "does not fit auto- or allopolyploid models cleanly" (Table 2; Gao
  et al. 2020, *PLoS ONE* 15(3), e0229624), so this simulator validates pipeline *mechanics*, not
  real sweetpotato transmission genetics.
- **polyRAD leaves half the simulated loci without a dosage call** (script 02, carried into script
  03): on the committed population polyRAD returns no dosage for either allele at 675 of the 1,500
  loci, and at 75 more it calls only the reference allele, which script 02's alt-allele extraction
  then drops. So 750 markers reach script 03 entirely missing and are removed by its > 20%
  missing-rate filter; a further 67 go on MAF, leaving 683. The discarded markers are not the less
  polymorphic ones (mean MAF 0.318 against 0.313 for those retained), so this is a property of how
  the mapping pipeline is configured, chiefly the deliberately tight `freqAllowedDeviation = 0.01`,
  and not MAF filtering of a complete matrix. The setting is left as it is so that the committed
  results match the archived release; a real dataset should be re-tuned against its own depth and
  missingness profile.
- **Untruncated noise in the deposited phenotypes** (script 01): trait values are a trait mean plus
  a genotype effect plus a Gaussian residual, with no truncation, so
  `data_simulated/plot_level_phenotypes.csv` holds values outside any physically meaningful range:
  27 negative Yield plots (minimum −20.6), 82 negative Carotenoid values, and the two visual
  scores, centred on 5, running from −1.5 to 13.4. This does not bias the comparison, which is
  scale-free and gives both methods the same phenotypes, but the file should not be reused as a
  realistic sweetpotato phenotype set.
- **Parental read counts** (script 02): script 01 builds the two parents' dosage but does not save
  it, so 02 reconstructs parent read counts from the population allele frequency rather than from
  the parents' own genotypes. polyRAD's mapping pipeline is therefore given parents that are not
  the cross's actual parents, which is one reason so many loci end up without a call. A real
  dataset would genotype the actual parents directly, and a loader that supplies real parental
  reads removes this limitation.
- **Trait architecture** (script 01): 40 QTL/trait, plus 20 further shared QTL added to each of
  Yield/Carotenoid and DryMatter/StorageScore with a controlled sign relationship, realised at
  r = −0.37 and r = +0.32 respectively. That gives those four traits 60 effect entries each, and 58
  to 60 distinct loci, because a shared index can coincide with a trait's own draw. Additive-only,
  no dominance or epistasis. Heritability
  is targeted at h² = 0.40 on a genotype-mean basis, with the residual-noise scale accounting for
  the 2 field replicates being averaged (see script 01 for the formula), not derived from real
  sweetpotato trait-heritability estimates. The trade-off and association signs (Yield vs.
  Carotenoid negative, DryMatter vs. StorageScore positive) are a simplified illustration inspired
  by the β-carotene/starch trade-off reported by Gemenet et al. (2020, *Theor. Appl. Genet.* 133(1),
  23–36) and catalogued in Table 1 of the article, not a claim about the real magnitude or even
  direction of these relationships in sweetpotato.
- **Stage-1/Stage-3 in-sample asymmetry in comparison (B)** (see scripts 05 and 06, and the
  [Two-tier validation](#two-tier-validation) section above): MGIDI is computed once on the full
  observed dataset, since it has no mechanism to extrapolate to a genotype it was not given
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
  univariate-then-combine design used here keeps the two methods' predictive scope matched.
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
