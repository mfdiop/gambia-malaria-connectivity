# Snakemake pipeline for the Objective 3 analysis — design

**Date:** 2026-06-03
**Status:** Design approved — implementation plan pending
**Scope:** This spec covers the *first* Snakefile only: automation of the existing five-phase R analysis. The manuscript write-up and any additional analyses (selection scans, time-stratified runs, demographic inference, etc.) are out of scope for this spec; they are anticipated and the design supports them, but they will get their own specs when work begins.

---

## 1. Context

The repository contains a five-phase population-genetics analysis answering whether residual malaria in western Gambia is sustained by local transmission, importation from eastern Gambia, or cross-border importation from Senegal. The scientific framing lives in `notes/objective3_proposal.md`; the directional-IBD method spec lives in `notes/directional_connectivity_Amambua-Ngwamodel.md`.

Today the analysis runs as 19 R scripts under `codes/`, sourced interactively in RStudio. Outputs land under `results/phaseN_NN/`. The scripts work, but the pipeline has no orchestration, no DAG, and no parameter-sweep mechanism — re-running after a change requires manually tracking what is stale.

This spec defines the first Snakemake automation pass.

## 2. Drivers (ranked)

The user confirmed all three drivers matter, ranked:

1. **Reproducibility for the thesis / manuscript** — one-command rerun from raw VCF → all figures/tables, so reviewers (and future-you) can verify. Pushes toward script-level rules with pinned envs.
2. **Re-running cheaply after changes** — DAG-driven incremental reruns. Pushes toward finer-grained rules and accurate input/output declarations.
3. **Parameter sweeps / sensitivity** — sweep K for sNMF, m for TreeMix, cohort subsets. Pushes toward wildcards on the right axes.

The design resolves the tension between (1) and (2,3) by making **one rule per R script** (1) but giving rules **wildcards on K, m, subset** (2,3), with explicit input/output declarations so the DAG is accurate.

## 3. Key design decisions

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | One Snakemake rule per canonical R script | Matches existing `codes/phaseN_NN_*.R` granularity; keeps reproducibility unit ≈ script |
| D2 | Inputs (VCFs, metadata) stay hardcoded inside R scripts; Snakemake controls outputs only | Per user choice — minimises R-script edits; trade-off is the pipeline is not portable to a machine without the sibling `01_objective1/` repo, which is acceptable since this is a single-author thesis project |
| D3 | Wildcards expose three axes only: `{k}`, `{m}`, `{subset}` | Tier-2 trade-off chosen by user — captures the real sweeps without combinatorial explosion. Sensitivity scenarios are different analyses (named rules), not parameter sweeps |
| D4 | Conda envs managed by Snakemake for TreeMix, hmmIBD, ADMIXTURE only | These three are version-sensitive and reproducibility-critical. R, PLINK2, vcftools, bcftools assumed on PATH |
| D4b | R package versions pinned via `renv` (lockfile in repo) | Closes the R-side of driver A (thesis reproducibility). The `renv/library/` is excluded from git; the lockfile is committed; `.Rprofile` auto-activates the project library when any R or Rscript call starts in the project root, so Snakemake rules require no changes |
| D5 | R scripts accept Snakemake params via `commandArgs(trailingOnly=TRUE)` with a tiny header block; defaults preserved for interactive sourcing | Avoids the `snakemake@input` directive (which breaks interactive use); no new R package dependency; per-script and self-documenting |
| D6 | `phase4_02_treemix.R` split into `phase4_02a_treemix.R` (per-m run) + `phase4_02b_optm.R` (OptM aggregation across the m sweep) | OptM needs all `m{m}` runs to exist first — natural DAG cut |
| D7 | hmmIBD runs once over the full monoclonal cohort; downstream rules subset the IBD matrix in memory | Cheaper than re-running hmmIBD per `{subset}`; hmmIBD output is the expensive artefact |
| D8 | `phase2_02_pca_admixture.R` kept on disk but NOT wired into the DAG | LEA is the canonical Phase 2.2 path; ADMIXTURE re-enters supervised in Phase 3.1. The unsupervised ADMIXTURE script can be wired in later as an alternative-method rule for the supplement if needed |
| D9 | `clustering_approaches*.R` excluded from DAG (scratchpads per CLAUDE.md) | Same |

## 4. Repository layout (proposed)

```
gambia-malaria-connectivity/
├── Snakefile                          # top-level; imports phase rules; defines `all`
├── config/
│   └── config.yaml                    # k_range, m_range, subsets, sibling-repo paths
├── envs/
│   ├── treemix.yaml
│   ├── hmmibd.yaml
│   └── admixture.yaml
├── workflow/
│   └── rules/
│       ├── common.smk                 # shared wildcard constraints + helpers
│       ├── phase1.smk                 # QC, COI, annotation
│       ├── phase2.smk                 # FST/IBD/Mantel, PCA-LEA (sweep K), hmmIBD
│       ├── phase3.smk                 # supervised ADMIXTURE (sweep K), directional IBD, isolation
│       ├── phase4.smk                 # cross-border IBD, decision subset, TreeMix (sweep m), resistance
│       ├── phase5.smk                 # attribution model, sensitivity (named, not wildcarded)
│       └── manuscript.smk             # placeholder hook for future Rmd/Quarto report
├── codes/                             # R scripts; minor edits to add commandArgs header
└── results/phaseN_NN/                 # outputs; subdirectories carry {k}/{m}/{subset}
```

## 5. Wildcards

| Wildcard | Domain | Source | Used in phases |
|----------|--------|--------|----------------|
| `{k}`    | integers in `config.yaml: k_range` (default `2..10`) | sNMF / supervised ADMIXTURE | 2.2, 3.1 |
| `{m}`    | integers in `config.yaml: m_range` (default `0..5`)  | TreeMix migration edges | 4.2 |
| `{subset}` | strings in `config.yaml: subsets` (default `allgambia, westonly, withsenegal`) | cohort definitions | 4.1, 4.2, 4.3, 5.1 |

Wildcard constraints declared in `workflow/rules/common.smk` so they cannot ambiguously match other parts of a path (e.g. `k` constrained to `\d+`).

Output paths bake wildcards in:
- `results/phase2_02_lea/K{k}/qmatrix.tsv`
- `results/phase4_02/subset_{subset}/m{m}/treemix.treeout.gz`
- `results/phase5_01/subset_{subset}/attribution_posterior.rds`

A parameter change invalidates only its subtree.

## 6. R script parameter convention

Each canonical script gets a small header prepended once:

```r
# --- snakemake interface (no-op when sourced interactively) ---
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- parse_arg("output_dir", "results/phase2_02_lea/K3")
k          <- as.integer(parse_arg("k", 3))
# --- end snakemake interface ---
```

Snakemake invokes scripts as `Rscript codes/phase2_02_pca_lea.R --output_dir={output_dir} --k={wildcards.k}`. Defaults preserve the existing interactive RStudio workflow.

## 7. Rule inventory

### Phase 1 — QC, COI, annotation (no wildcards)

| Rule | Script | Snakemake inputs | Outputs | Env |
|------|--------|------------------|---------|-----|
| `phase1_qc` | `phase1_01_wgs_qc.R` | (sibling VCF read internally) | `results/phase1_01/{qc_passed.vcf.gz, qc_report.html}` | system R |
| `phase1_coi` | `phase1_02_coi_stratification.R` | `phase1_01/qc_passed.vcf.gz` | `results/phase1_02/{coi_table.tsv, monoclonal_samples.txt, gambia_2014_qc_monoclonal.vcf.gz}` | system R |
| `phase1_annotate` | `phase1_03_sample_annotation.R` | `phase1_02/coi_table.tsv` | `results/phase1_03/sample_master_metadata.rds` | system R |

### Phase 2 — population structure

| Rule | Script | Wildcards | Snakemake inputs | Outputs | Env |
|------|--------|-----------|------------------|---------|-----|
| `phase2_fst_ibd_mantel` | `phase2_01_fst_ibd_mantel.R` | — | `phase1_03/sample_master_metadata.rds` | `results/phase2_01/{fst.tsv, mantel.rds, figs/*.pdf}` | system R |
| `phase2_pca_lea` | `phase2_02_pca_lea.R` | `{k}` | annotation | `results/phase2_02_lea/K{k}/{qmatrix.tsv, ce.tsv, fig.pdf}` | system R |
| `phase2_pca_lea_extendK` | `phase2_02_pca_lea_extendK.R` | — (aggregates over K) | all `K{k}/ce.tsv` for k in `k_range` | `results/phase2_02_lea/{ce_vs_K.pdf, K_summary.tsv}` | system R |
| `phase2_hmmibd` | `phase2_03_hmmibd.R` | — | `phase1_02/monoclonal_samples.txt`, `phase1_02/gambia_2014_qc_monoclonal.vcf.gz` | `results/phase2_03/{hmm_fract.txt, ibd_matrix.rds, figs/*.pdf}` | **hmmibd env** |

### Phase 3 — within-Gambia directionality

| Rule | Script | Wildcards | Snakemake inputs | Outputs | Env |
|------|--------|-----------|------------------|---------|-----|
| `phase3_supervised_admixture` | `phase3_01_supervised_admixture.R` | `{k}` | annotation, `phase2_03/ibd_matrix.rds` | `results/phase3_01/K{k}/{qmatrix.tsv, fig.pdf}` | **admixture env** |
| `phase3_directional_ibd` | `phase3_02_directional_ibd.R` | — | `phase2_03/ibd_matrix.rds`, annotation | `results/phase3_02/{P_trans.tsv, asymmetry.tsv, figs/*.pdf}` | system R |
| `phase3_isolation` | `phase3_03_isolation.R` | — | `phase3_02/P_trans.tsv` | `results/phase3_03/{isolation_index.tsv, figs/*.pdf}` | system R |

### Phase 4 — cross-border + selection

| Rule | Script | Wildcards | Snakemake inputs | Outputs | Env |
|------|--------|-----------|------------------|---------|-----|
| `phase4_decision_subset` | `phase4_01_decision_subset.R` | `{subset}` | annotation | `results/phase4_01_decision_subset/{subset}/{sample_list.txt, manifest.tsv}` | system R |
| `phase4_crossborder_ibd` | `phase4_01_crossborder_ibd.R` | `{subset}` | decision_subset sample list, `phase2_03/ibd_matrix.rds` | `results/phase4_01/{subset}/{ibd_pairs.tsv, figs/*.pdf}` | system R |
| `phase4_treemix` | `phase4_02a_treemix.R` (split out) | `{subset}, {m}` | decision_subset, phase1 VCF | `results/phase4_02/{subset}/m{m}/{treemix.treeout.gz, treemix.edges, treemix.cov.gz}` | **treemix env** |
| `phase4_treemix_optm` | `phase4_02b_optm.R` (split out) | `{subset}` | all `m{m}/treemix.*` for m in `m_range` | `results/phase4_02/{subset}/{optm_summary.tsv, optm_plot.pdf}` | system R |
| `phase4_resistance` | `phase4_03_resistance_haplotype_ancestry.R` | `{subset}` | decision_subset, annotation | `results/phase4_03/{subset}/{haplotype_freq.tsv, ancestry.tsv, figs/*.pdf}` | system R |

### Phase 5 — attribution + sensitivity

| Rule | Script | Wildcards | Snakemake inputs | Outputs | Env |
|------|--------|-----------|------------------|---------|-----|
| `phase5_attribution` | `phase5_01_attribution_model.R` | `{subset}` | `phase3_02/P_trans.tsv`, `phase4_01/{subset}/ibd_pairs.tsv`, `phase4_03/{subset}/ancestry.tsv` | `results/phase5_01/{subset}/{posterior.rds, attribution.tsv, figs/*.pdf}` | system R |
| `phase5_sensitivity` | `phase5_02_sensitivity.R` | — (named scenarios inside script) | all `phase5_01/{subset}/posterior.rds` | `results/phase5_02/{scenario_grid.tsv, figs/*.pdf}` | system R |

## 8. Top-level targets

| Target | What it builds | Use case |
|--------|----------------|----------|
| `all` | Canonical-defaults result set across all phases (chosen K and m only, all subsets) | Full pipeline rerun (driver 1) |
| `phase1` … `phase5` | Terminal outputs of that phase | Iterating on one phase |
| `phase2_lea_sweep` | All `K{k}` runs + `ce_vs_K.pdf` aggregation | Picking K for sNMF |
| `phase4_treemix_sweep` | All `m{m}` runs for default subset + OptM plot | Picking m for TreeMix |
| `phase5_sensitivity` | The sensitivity grid (depends on `phase5_01` for all subsets) | Sensitivity analysis |
| `manuscript` | Placeholder; errors with friendly message until implemented | Future Quarto/Rmd report |
| `clean_phase{N}` | Removes `results/phase{N}_*` only | Forcing a phase rerun |

`all` intentionally targets the **canonical-defaults** result set, not the full Cartesian sweep. The intended workflow is:

1. Run the sweep targets (`phase2_lea_sweep`, `phase4_treemix_sweep`) to explore K and m.
2. Inspect cross-entropy / OptM outputs and choose the canonical K and m.
3. Pin those choices as `canonical_k` / `canonical_m` in `config.yaml`.
4. Run `snakemake all` — which uses the pinned canonical values, not the sweep ranges.

This means `all` is deterministic given a fixed `config.yaml` and never accidentally launches the full Cartesian sweep just from typing `snakemake`. Downstream rules (Phase 3.1, Phase 4.2 OptM, Phase 5) read `canonical_k` / `canonical_m` from config when building their `all`-target dependencies.

## 9. Extensibility model

| Slot | Action | Cost |
|------|--------|------|
| New rule in an existing phase | Add R script under `codes/phaseN_NN_*.R`; add rule block to `workflow/rules/phaseN.smk` | ~15 lines of Snakefile |
| New phase | Create `codes/phase6_*.R`; create `workflow/rules/phase6.smk`; add one `include:` line to `Snakefile`; add `phase6` to `all` | one new file + two existing-file edits |
| New sweep axis | Add wildcard's value list to `config.yaml`; reference `{newaxis}` in relevant output paths; update affected R scripts to parse `--newaxis=` | Localised — change scales with number of affected rules |

## 10. Anticipated future analyses (not implemented in this pass)

Listed for traceability so the spec's scope decisions can be re-evaluated when they arrive. Per user direction, no placeholder rules are added now — they'd be added via the slot-1 path above when the work starts.

- **Selection scans** (iHS / XP-EHH / nSL on the cross-border subset) — would slot in as Phase 4.4.
- **Time-stratified comparison** (split monoclonal samples by collection year) — would introduce a `{year}` wildcard and meta-rules.
- **Demographic inference / Ne estimation** (IBDNe, SMC++) — would slot in as Phase 4.5 consuming the hmmIBD matrix.
- **Manuscript write-up** (Quarto/Rmd render) — `manuscript.smk` placeholder is the hook.

## 11. Out of scope for this spec

- Refactoring R script *internals* (the `Sys.sleep(0.1)` in `get_ibs_matrix()`, restructuring scratchpad logic, etc.). The R scripts are touched only to add the `commandArgs` header block.
- Migrating inputs to local `data/raw/` (decision D2 keeps them hardcoded).
- Containerisation (Docker / Singularity) — conda envs only.
- HPC submission wrappers (e.g. SLURM profile) — Snakefile will be cluster-compatible by virtue of using only standard rule mechanics, but no profile is included.
- The `manuscript` rule body itself.

## 12. Open decisions for the implementation plan

These are deliberately deferred to the writing-plans step, not the design:

- Concrete `config.yaml` values (current default-K, default-m, exact subset names).
- Whether the `commandArgs` header is added by hand to each script or via a single edit pass.
- Order of implementation (likely: Snakefile skeleton + Phase 1 first, then add phases bottom-up).
- Whether `phase4_02_treemix.R` split happens as a separate commit before the Snakefile work or together with it.
