# Snakemake Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Automate the 5-phase Objective 3 population-genetics analysis as a single Snakemake pipeline with parameter sweeps for K, m, and cohort subset, while keeping R scripts interactively runnable.

**Architecture:** One rule per canonical R script; wildcards `{k}`, `{m}`, `{subset}` for sweeps; conda envs only for TreeMix / hmmIBD / ADMIXTURE; R scripts gain a tiny `commandArgs` header so Snakemake controls outputs while interactive sourcing still works. Inputs (sibling-repo VCFs) stay hardcoded in R per design decision D2.

**Tech Stack:** Snakemake ≥7, R + LEA/hierfstat/vcfR/tidyverse, conda (miniforge3) for tool envs, TreeMix v1.13, hmmIBD, ADMIXTURE v1.3, PLINK2.

**Design source:** `docs/superpowers/specs/2026-06-03-snakemake-pipeline-design.md`. Read it before starting Task 1.

**Testing model:** This is an analysis project with no test suite. The "test" for each Snakemake rule is:
1. `snakemake --lint` passes.
2. `snakemake -n <target>` produces a sensible DAG (correct dependencies, no missing inputs).
3. Running the rule produces the expected output files (non-empty, parseable by the next rule).
4. Where possible, compare key outputs against the existing manual runs under `results/phaseN_NN/`.

Skip step 4 if no prior manual output exists for that rule.

---

## File Structure

**Create:**
- `Snakefile` — top-level; `include:`s all phase rule files; defines `all` and convenience targets
- `config/config.yaml` — k_range, m_range, subsets, canonical_k, canonical_m, paths
- `envs/treemix.yaml`, `envs/hmmibd.yaml`, `envs/admixture.yaml` — conda specs
- `workflow/rules/common.smk` — wildcard constraints, config-load helpers
- `workflow/rules/phase1.smk`, `phase2.smk`, `phase3.smk`, `phase4.smk`, `phase5.smk`, `manuscript.smk`
- `codes/phase4_02a_treemix.R`, `codes/phase4_02b_optm.R` — split from existing `phase4_02_treemix.R`

**Modify:** (add ~10-line `commandArgs` header — see header template below)
- `codes/phase1_01_wgs_qc.R`
- `codes/phase1_02_coi_stratification.R`
- `codes/phase1_03_sample_annotation.R`
- `codes/phase2_01_fst_ibd_mantel.R`
- `codes/phase2_02_pca_lea.R`
- `codes/phase2_02_pca_lea_extendK.R`
- `codes/phase2_03_hmmibd.R`
- `codes/phase3_01_supervised_admixture.R`
- `codes/phase3_02_directional_ibd.R`
- `codes/phase3_03_isolation.R`
- `codes/phase4_01_decision_subset.R`
- `codes/phase4_01_crossborder_ibd.R`
- `codes/phase4_02a_treemix.R` (newly split)
- `codes/phase4_02b_optm.R` (newly split)
- `codes/phase4_03_resistance_haplotype_ancestry.R`
- `codes/phase5_01_attribution_model.R`
- `codes/phase5_02_sensitivity.R`

**Untouched (kept as scratchpads):** `codes/phase2_02_pca_admixture.R`, `codes/clustering_approaches.R`, `codes/clustering_approaches_v1.R`.

---

## Header Template (reused across all canonical R scripts)

Paste this block immediately after the existing top-of-file comment block, before any `library()` call:

```r
# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
# Each script reads only the args it needs; e.g.:
#   output_dir <- .parse_arg("output_dir", "results/phaseX_YY")
#   k          <- as.integer(.parse_arg("k", 3))
# Existing hardcoded output paths inside the script should be replaced by
# `output_dir` (and `file.path(output_dir, ...)` for sub-files).
# --- end snakemake interface ------------------------------------------------
```

For each script-modification step the plan will note which args to parse and which hardcoded output strings to replace.

---

## Task 0: Initialise `renv` for R package reproducibility

The project currently uses the system R library — running `library(LEA)` picks up whatever version of LEA `install.packages()` last wrote. For driver A (thesis-grade reproducibility) we want a per-project library snapshotted into `renv.lock`. `renv` is the standard R tool for this; it auto-discovers package usage by scanning the codebase, writes a lockfile, and reactivates the project library every time R starts in the project directory.

**Files:**
- Create: `renv.lock` (written by `renv::init()`)
- Create: `.Rprofile` (written by `renv::init()` — auto-activates renv on R startup)
- Create: `renv/activate.R`, `renv/settings.json` (written by `renv::init()` — committed)
- Modify: `.gitignore` (renv writes its own `renv/.gitignore` excluding `renv/library/` — verify)

- [ ] **Step 1: Verify R and find version**

Run: `R --version`
Expected: prints R version (e.g., `R version 4.4.x`). Note the version — `renv.lock` will pin it.

- [ ] **Step 2: Install `renv` into the user library (one-off)**

Run: `R --quiet -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")'`
Expected: either silently does nothing (renv already installed) or installs it. No error.

- [ ] **Step 3: Initialise renv in the project**

Run: `R --quiet -e 'renv::init(bare = FALSE)'`

Expected: renv scans `codes/` for `library()` / `require()` calls, prompts (or auto-decides — depends on renv version) to copy your currently-installed package versions into a new `renv/library/` symlink farm, and writes:
- `renv.lock` — package versions
- `.Rprofile` — auto-activates renv
- `renv/activate.R`
- `renv/settings.json`
- `renv/.gitignore` (excludes `renv/library/`, `renv/staging/`, etc.)

Watch for "the following package(s) are used in your project, but are not installed" warnings — those flag packages your scripts reference but aren't on your system yet (e.g., `OptM`, `dcifer`, `isoRelate`, `LEA`). Note them; install in Step 4.

- [ ] **Step 4: Install any missing packages renv flagged**

For each missing package `X` from Step 3:

Run: `R --quiet -e 'renv::install("X")'`

For Bioconductor packages (e.g. `LEA`, `SNPRelate`): `R --quiet -e 'renv::install("bioc::LEA")'`.
For GitHub-only packages (e.g. `dcifer`, `isoRelate`): `R --quiet -e 'renv::install("OJWatson/isoRelate")'` (adjust org/repo as appropriate; check each package's README for the canonical install instruction).

After all missing packages are installed:

Run: `R --quiet -e 'renv::snapshot()'`
Expected: prompts (or auto-decides) to update `renv.lock` to include the newly installed packages. Answer yes.

- [ ] **Step 5: Verify state is consistent**

Run: `R --quiet -e 'renv::status()'`
Expected: prints `* The project is already synchronized with the lockfile.` (or equivalent). If it reports discrepancies, run `renv::snapshot()` again and re-check.

- [ ] **Step 6: Sanity-check a script still parses under renv**

Run: `Rscript -e 'parse(file = "codes/phase2_02_pca_lea.R"); cat("OK\n")'`
Expected: prints `OK`. (We're not running the script — just confirming the project's `.Rprofile` doesn't break Rscript invocation. Snakemake will invoke `Rscript codes/...` from the project root, which picks up `.Rprofile` automatically and activates the project library.)

- [ ] **Step 7: Confirm `.gitignore` covers `renv/library`**

Run: `cat renv/.gitignore`
Expected: includes `library/`, `local/`, `cellar/`, `lock/`, `python/`, `sandbox/`, `staging/`. This is renv's own gitignore — it correctly excludes the library symlinks while keeping `activate.R`, `settings.json`, and `renv.lock` trackable.

If you also want renv-related entries in the top-level `.gitignore`, no change needed — renv's own file is sufficient.

- [ ] **Step 8: Commit**

```bash
git add renv.lock .Rprofile renv/activate.R renv/settings.json renv/.gitignore
git commit -m "Initialise renv for R package version pinning"
```

**Notes for subsequent tasks:**
- Whenever a later task installs a new R package, also run `renv::snapshot()` and amend the lockfile in that task's commit (or as a follow-up commit).
- The Snakemake rules in Tasks 4-8 invoke `Rscript codes/...` from the project root. `.Rprofile` activates the renv library before the script runs, so no Snakemake-side changes are needed.
- If you ever need to give a fresh machine a clean restore, the recipe is: `git clone`, `cd <repo>`, `R --quiet -e 'renv::restore()'`. That's the reproducibility payoff.

---

## Task 1: Repository skeleton

**Files:**
- Create: `Snakefile`
- Create: `config/config.yaml`
- Create: `workflow/rules/common.smk`
- Create: `workflow/rules/phase1.smk`, `phase2.smk`, `phase3.smk`, `phase4.smk`, `phase5.smk`, `manuscript.smk`

- [ ] **Step 1: Verify Snakemake is installed**

Run: `snakemake --version`
Expected: prints a version string ≥ 7.0. If not installed: `conda install -n base -c bioconda snakemake-minimal` (or use existing env).

- [ ] **Step 2: Create directory tree**

Run: `mkdir -p config envs workflow/rules`

- [ ] **Step 3: Write `config/config.yaml`**

```yaml
# Sweep ranges
k_range: [2, 3, 4, 5, 6, 7, 8, 9, 10]
m_range: [0, 1, 2, 3, 4, 5]

# Canonical defaults (picked after sweep + inspection; updated when chosen)
canonical_k: 3
canonical_m: 1

# Cohort subsets
subsets: ["allgambia", "westonly", "withsenegal"]

# Cohort used for canonical `all` target (typically the one Phase 5 attributes on)
canonical_subset: "withsenegal"

# Paths to sibling-repo inputs (read internally by R scripts; declared here
# for documentation, not used as Snakemake inputs per design decision D2)
sibling:
  raw_vcf_dir: "../01_objective1/01_data/raw_data"
  metadata_dir: "../01_objective1/01_data/meta"
```

- [ ] **Step 4: Write `workflow/rules/common.smk`**

```python
# Shared wildcard constraints and config helpers.
configfile: "config/config.yaml"

wildcard_constraints:
    k       = r"\d+",
    m       = r"\d+",
    subset  = r"[a-zA-Z0-9_]+",

def canonical_k():
    return config["canonical_k"]

def canonical_m():
    return config["canonical_m"]

def canonical_subset():
    return config["canonical_subset"]

def all_subsets():
    return config["subsets"]
```

- [ ] **Step 5: Write empty phase rule files**

Create each of `workflow/rules/phase{1..5}.smk` and `workflow/rules/manuscript.smk` with a one-line comment:

```python
# Rules for Phase N — populated by later tasks.
```

For `manuscript.smk`, use:

```python
# Manuscript render — placeholder; implemented after analysis is complete.
rule manuscript:
    output:
        "manuscript/objective3.pdf"
    shell:
        "echo 'Manuscript rule not implemented yet — add Quarto/Rmd render command here.' >&2 && exit 1"
```

- [ ] **Step 6: Write top-level `Snakefile`**

```python
include: "workflow/rules/common.smk"
include: "workflow/rules/phase1.smk"
include: "workflow/rules/phase2.smk"
include: "workflow/rules/phase3.smk"
include: "workflow/rules/phase4.smk"
include: "workflow/rules/phase5.smk"
include: "workflow/rules/manuscript.smk"

# Convenience targets (populated as rules are added in Tasks 4-8)
rule all:
    input:
        []  # populated in Task 9
```

- [ ] **Step 7: Verify Snakemake parses the skeleton**

Run: `snakemake --lint`
Expected: completes without error (may emit zero-rule warnings — those are fine).

Run: `snakemake -n all`
Expected: "Nothing to be done" or "0 of 0 steps".

- [ ] **Step 8: Commit**

```bash
git add Snakefile config/config.yaml workflow/rules/
git commit -m "Add Snakemake skeleton (config, common, empty phase files)"
```

---

## Task 2: Conda environment specs

**Files:**
- Create: `envs/treemix.yaml`, `envs/hmmibd.yaml`, `envs/admixture.yaml`

- [ ] **Step 1: Write `envs/treemix.yaml`**

```yaml
name: treemix
channels:
  - bioconda
  - conda-forge
dependencies:
  - treemix=1.13
  - plink=1.9
  - bcftools
  - python>=3.9
```

- [ ] **Step 2: Write `envs/hmmibd.yaml`**

```yaml
name: hmmibd
channels:
  - bioconda
  - conda-forge
dependencies:
  - hmmibd
  - bcftools
  - python>=3.9
```

- [ ] **Step 3: Write `envs/admixture.yaml`**

```yaml
name: admixture
channels:
  - bioconda
  - conda-forge
dependencies:
  - admixture=1.3.0
  - plink=1.9
```

- [ ] **Step 4: Verify envs build**

Run: `snakemake --use-conda --conda-create-envs-only -n` (note: this still requires at least one rule with `conda:`. Skip verification here and defer to Task 6/7 where the first conda-using rules are added.)

Mark this step satisfied if the YAMLs are syntactically valid YAML (`python -c "import yaml; [yaml.safe_load(open(f)) for f in ['envs/treemix.yaml','envs/hmmibd.yaml','envs/admixture.yaml']]"`).

- [ ] **Step 5: Commit**

```bash
git add envs/
git commit -m "Add conda env specs for TreeMix, hmmIBD, ADMIXTURE"
```

---

## Task 3: Split `phase4_02_treemix.R` into per-m + OptM scripts

The existing `codes/phase4_02_treemix.R` runs the m-sweep, OptM aggregation, and bootstrap in one script. The Snakemake DAG needs these as separate rules so the m-sweep can fan out per `{m}` while OptM aggregates afterwards.

**Files:**
- Create: `codes/phase4_02a_treemix.R` (one m, one subset)
- Create: `codes/phase4_02b_optm.R` (aggregate over m; bootstrap + plots at best m)
- Modify: original `codes/phase4_02_treemix.R` — leave unchanged on disk for reference; the new scripts are derived from it

- [ ] **Step 1: Read and map the existing script**

Run: `wc -l codes/phase4_02_treemix.R` (expect ~600 lines per the header comment).
Open it and identify the section boundaries by comment headers:
- "Sample list per population; subset filtered VCF" → goes to **02a**
- "LD-prune (PLINK)" → goes to **02a**
- "Build TreeMix input" → goes to **02a**
- "Run TreeMix m = 0..M_MAX, each with N_REPS_PER_M" → **02a runs one m only**; the loop becomes parameterised by Snakemake `{m}` wildcard
- "Pick optimal m via OptM" → **02b**
- "Run N_BOOTSTRAP iterations at chosen m" → **02b**
- "Plot best tree, residuals, migration-edge bootstrap" → **02b**

- [ ] **Step 2: Create `codes/phase4_02a_treemix.R`**

Header (replaces the `Outputs:` block in the original):

```r
# =============================================================================
# Phase 4 — Step 4.2a: TreeMix run for one (subset, m) combination
# =============================================================================
# Split from phase4_02_treemix.R for Snakemake parameterisation. Runs N_REPS
# replicates at a single m value for a single cohort subset.
#
# Args: --output_dir, --subset, --m, --n_reps
# Outputs (under output_dir):
#   populations.txt, samples_kept.txt, subset.vcf.gz (+ .tbi)  [shared scaffold]
#   pruned.{prune.in,prune.out}                                [shared scaffold]
#   treemix_input.frq.gz                                       [shared scaffold]
#   runs/treemix_m{m}_rep{r}.{treeout.gz,llik,modelcov.gz,cov.gz,covse.gz}
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase4_02/subset_withsenegal/m1")
subset_name <- .parse_arg("subset", "withsenegal")
m_value     <- as.integer(.parse_arg("m", 1))
n_reps      <- as.integer(.parse_arg("n_reps", 10))
# --- end snakemake interface ------------------------------------------------
```

Then copy from the original everything from the first `library()` call through and including the per-(m, rep) TreeMix run loop, but:
- Restrict the m-loop to the single `m_value`.
- Replace hardcoded `M_MAX`, `N_REPS_PER_M` references with `m_value`, `n_reps`.
- Replace hardcoded `results/phase4_02/...` output prefixes with `file.path(output_dir, ...)`.
- Drop the OptM, bootstrap, and plotting sections (those go in 02b).

- [ ] **Step 3: Create `codes/phase4_02b_optm.R`**

Header:

```r
# =============================================================================
# Phase 4 — Step 4.2b: OptM aggregation + bootstrap at best m
# =============================================================================
# Consumes the per-m runs produced by phase4_02a_treemix.R across the full
# m_range for one cohort subset, picks best m via OptM, runs bootstrap at
# best m, produces final plots.
#
# Args: --output_dir, --subset, --runs_dirs (comma-separated list of per-m
#       output_dirs from 02a), --n_bootstrap
# Outputs (under output_dir):
#   optm/                          OptM-format inputs
#   optm_summary.tsv               OptM delta-m table
#   optm_curve.pdf                 log-lik / variance / delta-m plot
#   best_m.txt
#   bootstrap/treemix_boot{r}.*    per-rep bootstrap at best m
#   migration_edge_support.tsv     per-edge bootstrap support
#   tree_best.pdf
#   residuals_best.pdf
#   phase4_02b_summary.txt
#   phase4_02b_results.rds
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir   <- .parse_arg("output_dir", "results/phase4_02/subset_withsenegal")
subset_name  <- .parse_arg("subset", "withsenegal")
runs_dirs    <- strsplit(.parse_arg("runs_dirs", ""), ",", fixed = TRUE)[[1]]
n_bootstrap  <- as.integer(.parse_arg("n_bootstrap", 100))
# --- end snakemake interface ------------------------------------------------
```

Then copy from the original the OptM, bootstrap, and plotting sections. Adapt input-file discovery: instead of globbing `results/phase4_02/runs/`, iterate over `runs_dirs` (each entry is one m's `output_dir` from 02a) and glob `<dir>/runs/treemix_m*_rep*.*` from each.

- [ ] **Step 4: Sanity-check both new scripts parse**

Run: `Rscript -e 'parse(file="codes/phase4_02a_treemix.R"); parse(file="codes/phase4_02b_optm.R"); cat("OK\n")'`
Expected: prints `OK` (parse errors would show as syntax errors).

- [ ] **Step 5: Commit**

```bash
git add codes/phase4_02a_treemix.R codes/phase4_02b_optm.R
git commit -m "Split phase4_02_treemix.R into per-m run (02a) and OptM aggregation (02b)"
```

Note: the original `codes/phase4_02_treemix.R` stays in place (untouched, scratchpad).

---

## Task 4: Phase 1 rules

**Files:**
- Modify: `codes/phase1_01_wgs_qc.R`, `codes/phase1_02_coi_stratification.R`, `codes/phase1_03_sample_annotation.R`
- Modify: `workflow/rules/phase1.smk`

- [ ] **Step 1: Add `commandArgs` header to `codes/phase1_01_wgs_qc.R`**

Paste the header template (see top of plan) immediately after the existing top-of-file comment block. Add this line within the parsed-args section:

```r
output_dir <- .parse_arg("output_dir", "results/phase1_01")
```

Find every hardcoded `"results/phase1_01"` string in the script and replace with `output_dir`, and every `"results/phase1_01/<file>"` with `file.path(output_dir, "<file>")`. Use Grep / Edit to be exhaustive.

- [ ] **Step 2: Add `commandArgs` header to `codes/phase1_02_coi_stratification.R`**

Same procedure. Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase1_02")
qc_vcf     <- .parse_arg("qc_vcf", "results/phase1_01/qc_passed.vcf.gz")
```

Replace hardcoded paths with these.

- [ ] **Step 3: Add `commandArgs` header to `codes/phase1_03_sample_annotation.R`**

Same procedure. Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase1_03")
coi_table  <- .parse_arg("coi_table", "results/phase1_02/coi_table.tsv")
```

- [ ] **Step 4: Verify each script still runs interactively**

For each modified script, in R: `source("codes/phase1_01_wgs_qc.R")` — should behave exactly as before (defaults restore original paths). Skip if doing so would re-run an expensive step; instead just confirm the file parses: `Rscript -e 'parse(file="codes/phase1_01_wgs_qc.R")'`.

- [ ] **Step 5: Write Phase 1 rules in `workflow/rules/phase1.smk`**

```python
rule phase1_qc:
    output:
        vcf       = "results/phase1_01/qc_passed.vcf.gz",
        report    = "results/phase1_01/qc_report.html",
    params:
        outdir = "results/phase1_01",
    shell:
        "Rscript codes/phase1_01_wgs_qc.R --output_dir={params.outdir}"

rule phase1_coi:
    input:
        vcf = rules.phase1_qc.output.vcf,
    output:
        coi_table  = "results/phase1_02/coi_table.tsv",
        mono_list  = "results/phase1_02/monoclonal_samples.txt",
        mono_vcf   = "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz",
    params:
        outdir = "results/phase1_02",
    shell:
        "Rscript codes/phase1_02_coi_stratification.R "
        "--output_dir={params.outdir} --qc_vcf={input.vcf}"

rule phase1_annotate:
    input:
        coi = rules.phase1_coi.output.coi_table,
    output:
        meta = "results/phase1_03/sample_master_metadata.rds",
    params:
        outdir = "results/phase1_03",
    shell:
        "Rscript codes/phase1_03_sample_annotation.R "
        "--output_dir={params.outdir} --coi_table={input.coi}"
```

- [ ] **Step 6: Verify DAG builds**

Run: `snakemake -n phase1_annotate`
Expected: lists three jobs (qc → coi → annotate) in correct order. No "MissingInputException".

- [ ] **Step 7: Decide whether to actually run Phase 1**

Phase 1 outputs already exist under `results/phase1_01/`, `results/phase1_02/`, `results/phase1_03/` from manual runs. To avoid re-running, mark them as up-to-date:

Run: `snakemake --touch phase1_annotate`
Expected: timestamps refreshed; subsequent dry-runs show no work.

If you want a fresh end-to-end test, instead:
Run: `snakemake --cores 1 phase1_annotate` (will rebuild from scratch; takes ~minutes).

- [ ] **Step 8: Commit**

```bash
git add codes/phase1_01_wgs_qc.R codes/phase1_02_coi_stratification.R codes/phase1_03_sample_annotation.R workflow/rules/phase1.smk
git commit -m "Wire Phase 1 (QC, COI, annotation) into Snakemake"
```

---

## Task 5: Phase 2 rules

**Files:**
- Modify: `codes/phase2_01_fst_ibd_mantel.R`, `codes/phase2_02_pca_lea.R`, `codes/phase2_02_pca_lea_extendK.R`, `codes/phase2_03_hmmibd.R`
- Modify: `workflow/rules/phase2.smk`

- [ ] **Step 1: Add header to `phase2_01_fst_ibd_mantel.R`**

Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase2_01")
meta_in    <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
```

Replace hardcoded paths.

- [ ] **Step 2: Add header to `phase2_02_pca_lea.R` (per-K)**

Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase2_02_lea/K3")
k          <- as.integer(.parse_arg("k", 3))
meta_in    <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
mono_vcf   <- .parse_arg("mono_vcf", "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
```

The script currently sweeps K internally; restrict the K loop to `k:k` so each Snakemake job runs one K. Outputs go under `output_dir` (one K per dir).

**Note:** the LEA `snmfProject` requires writing all K runs to the same project file to enable cross-entropy aggregation. Two viable patterns:

- **Pattern A (recommended):** Run sNMF for the single K, write `qmatrix.tsv` + `ce.tsv` into `output_dir`. The aggregation step (Step 3 / `extendK`) reads the per-K `ce.tsv` files directly without touching the snmfProject file. This avoids contention on the project file across parallel jobs.
- **Pattern B:** Use Snakemake's `localrules:` directive and run all K serially in one job. Simpler but loses parallelism — fall back to this only if Pattern A doesn't reproduce the existing cross-entropy values.

Go with Pattern A. Persist `qmatrix.tsv` and `ce.tsv` per `K{k}/`.

- [ ] **Step 3: Add header to `phase2_02_pca_lea_extendK.R`**

This script becomes the **aggregation** rule, reading all per-K `ce.tsv` files produced by Step 2 and producing the cross-entropy curve and summary. Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase2_02_lea")
ce_files   <- strsplit(.parse_arg("ce_files", ""), ",", fixed = TRUE)[[1]]
```

Replace the internal `LEA::cross.entropy()` call with reading each path in `ce_files` and concatenating. Strip out the `snmfProject` "continue" logic — it's no longer needed in the aggregation context. Leave the plotting code as-is.

- [ ] **Step 4: Add header to `phase2_03_hmmibd.R`**

Add:

```r
output_dir <- .parse_arg("output_dir", "results/phase2_03")
mono_vcf   <- .parse_arg("mono_vcf", "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
mono_list  <- .parse_arg("mono_list", "results/phase1_02/monoclonal_samples.txt")
```

Replace hardcoded paths.

- [ ] **Step 5: Write Phase 2 rules in `workflow/rules/phase2.smk`**

```python
rule phase2_fst_ibd_mantel:
    input:
        meta = rules.phase1_annotate.output.meta,
    output:
        fst    = "results/phase2_01/fst.tsv",
        mantel = "results/phase2_01/mantel.rds",
    params:
        outdir = "results/phase2_01",
    shell:
        "Rscript codes/phase2_01_fst_ibd_mantel.R "
        "--output_dir={params.outdir} --meta={input.meta}"

rule phase2_pca_lea:
    input:
        meta     = rules.phase1_annotate.output.meta,
        mono_vcf = rules.phase1_coi.output.mono_vcf,
    output:
        qmatrix = "results/phase2_02_lea/K{k}/qmatrix.tsv",
        ce      = "results/phase2_02_lea/K{k}/ce.tsv",
    params:
        outdir = "results/phase2_02_lea/K{k}",
    shell:
        "Rscript codes/phase2_02_pca_lea.R "
        "--output_dir={params.outdir} --k={wildcards.k} "
        "--meta={input.meta} --mono_vcf={input.mono_vcf}"

rule phase2_pca_lea_extendK:
    input:
        ce_files = expand("results/phase2_02_lea/K{k}/ce.tsv", k=config["k_range"]),
    output:
        plot    = "results/phase2_02_lea/ce_vs_K.pdf",
        summary = "results/phase2_02_lea/K_summary.tsv",
    params:
        outdir = "results/phase2_02_lea",
        ce_files_csv = lambda w, input: ",".join(input.ce_files),
    shell:
        "Rscript codes/phase2_02_pca_lea_extendK.R "
        "--output_dir={params.outdir} --ce_files={params.ce_files_csv}"

rule phase2_hmmibd:
    input:
        mono_vcf  = rules.phase1_coi.output.mono_vcf,
        mono_list = rules.phase1_coi.output.mono_list,
    output:
        hmm_fract  = "results/phase2_03/hmm_fract.txt",
        ibd_matrix = "results/phase2_03/ibd_matrix.rds",
    params:
        outdir = "results/phase2_03",
    conda:
        "../../envs/hmmibd.yaml"
    shell:
        "Rscript codes/phase2_03_hmmibd.R "
        "--output_dir={params.outdir} "
        "--mono_vcf={input.mono_vcf} --mono_list={input.mono_list}"
```

- [ ] **Step 6: Verify DAG**

Run: `snakemake -n phase2_hmmibd phase2_fst_ibd_mantel "results/phase2_02_lea/ce_vs_K.pdf"`
Expected: lists hmmibd, fst, all K runs, and the extendK aggregation. Correct dependency order.

- [ ] **Step 7: Touch existing outputs (skip rerun) or actually run**

If existing outputs are still valid: `snakemake --touch phase2_hmmibd phase2_fst_ibd_mantel`

If running fresh: `snakemake --use-conda --cores 4 phase2_hmmibd` (note: hmmibd rule needs the conda env, which will be built on first run — expect ~minutes).

- [ ] **Step 8: Commit**

```bash
git add codes/phase2_*.R workflow/rules/phase2.smk
git commit -m "Wire Phase 2 (FST/IBD/Mantel, PCA-LEA sweep K, hmmIBD) into Snakemake"
```

---

## Task 6: Phase 3 rules

**Files:**
- Modify: `codes/phase3_01_supervised_admixture.R`, `codes/phase3_02_directional_ibd.R`, `codes/phase3_03_isolation.R`
- Modify: `workflow/rules/phase3.smk`

- [ ] **Step 1: Add header to each Phase 3 script**

For `phase3_01_supervised_admixture.R`:

```r
output_dir <- .parse_arg("output_dir", "results/phase3_01/K3")
k          <- as.integer(.parse_arg("k", 3))
meta_in    <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
ibd_in     <- .parse_arg("ibd_matrix", "results/phase2_03/ibd_matrix.rds")
```

For `phase3_02_directional_ibd.R`:

```r
output_dir <- .parse_arg("output_dir", "results/phase3_02")
meta_in    <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
ibd_in     <- .parse_arg("ibd_matrix", "results/phase2_03/ibd_matrix.rds")
```

For `phase3_03_isolation.R`:

```r
output_dir <- .parse_arg("output_dir", "results/phase3_03")
ptrans_in  <- .parse_arg("p_trans", "results/phase3_02/P_trans.tsv")
```

Replace hardcoded paths in each script.

- [ ] **Step 2: Write Phase 3 rules in `workflow/rules/phase3.smk`**

```python
rule phase3_supervised_admixture:
    input:
        meta       = rules.phase1_annotate.output.meta,
        ibd_matrix = rules.phase2_hmmibd.output.ibd_matrix,
    output:
        qmatrix = "results/phase3_01/K{k}/qmatrix.tsv",
    params:
        outdir = "results/phase3_01/K{k}",
    conda:
        "../../envs/admixture.yaml"
    shell:
        "Rscript codes/phase3_01_supervised_admixture.R "
        "--output_dir={params.outdir} --k={wildcards.k} "
        "--meta={input.meta} --ibd_matrix={input.ibd_matrix}"

rule phase3_directional_ibd:
    input:
        meta       = rules.phase1_annotate.output.meta,
        ibd_matrix = rules.phase2_hmmibd.output.ibd_matrix,
    output:
        p_trans   = "results/phase3_02/P_trans.tsv",
        asymmetry = "results/phase3_02/asymmetry.tsv",
    params:
        outdir = "results/phase3_02",
    shell:
        "Rscript codes/phase3_02_directional_ibd.R "
        "--output_dir={params.outdir} "
        "--meta={input.meta} --ibd_matrix={input.ibd_matrix}"

rule phase3_isolation:
    input:
        p_trans = rules.phase3_directional_ibd.output.p_trans,
    output:
        idx = "results/phase3_03/isolation_index.tsv",
    params:
        outdir = "results/phase3_03",
    shell:
        "Rscript codes/phase3_03_isolation.R "
        "--output_dir={params.outdir} --p_trans={input.p_trans}"
```

- [ ] **Step 3: Verify DAG**

Run: `snakemake -n phase3_isolation "results/phase3_01/K3/qmatrix.tsv"`
Expected: jobs include hmmIBD → directional_ibd → isolation, plus the K=3 ADMIXTURE rule.

- [ ] **Step 4: Touch or run**

Outputs likely exist from manual runs — use `snakemake --touch` if so.

- [ ] **Step 5: Commit**

```bash
git add codes/phase3_*.R workflow/rules/phase3.smk
git commit -m "Wire Phase 3 (supervised ADMIXTURE, directional IBD, isolation) into Snakemake"
```

---

## Task 7: Phase 4 rules (including `{subset}` and `{m}` sweep)

**Files:**
- Modify: `codes/phase4_01_decision_subset.R`, `codes/phase4_01_crossborder_ibd.R`, `codes/phase4_03_resistance_haplotype_ancestry.R`
- `codes/phase4_02a_treemix.R` and `codes/phase4_02b_optm.R` already have headers from Task 3
- Modify: `workflow/rules/phase4.smk`

- [ ] **Step 1: Add header to `phase4_01_decision_subset.R`**

```r
output_dir <- .parse_arg("output_dir", "results/phase4_01_decision_subset/withsenegal")
subset_name <- .parse_arg("subset", "withsenegal")
meta_in    <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
```

The script likely has a hardcoded "which subset to build" choice — replace it with `subset_name` driving the conditional logic that selects which samples go into the list.

- [ ] **Step 2: Add header to `phase4_01_crossborder_ibd.R`**

```r
output_dir   <- .parse_arg("output_dir", "results/phase4_01/withsenegal")
subset_name  <- .parse_arg("subset", "withsenegal")
sample_list  <- .parse_arg("sample_list", "results/phase4_01_decision_subset/withsenegal/sample_list.txt")
ibd_in       <- .parse_arg("ibd_matrix", "results/phase2_03/ibd_matrix.rds")
```

- [ ] **Step 3: Add header to `phase4_03_resistance_haplotype_ancestry.R`**

```r
output_dir   <- .parse_arg("output_dir", "results/phase4_03/withsenegal")
subset_name  <- .parse_arg("subset", "withsenegal")
sample_list  <- .parse_arg("sample_list", "results/phase4_01_decision_subset/withsenegal/sample_list.txt")
meta_in      <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
```

- [ ] **Step 4: Write Phase 4 rules in `workflow/rules/phase4.smk`**

```python
rule phase4_decision_subset:
    input:
        meta = rules.phase1_annotate.output.meta,
    output:
        sample_list = "results/phase4_01_decision_subset/{subset}/sample_list.txt",
        manifest    = "results/phase4_01_decision_subset/{subset}/manifest.tsv",
    params:
        outdir = "results/phase4_01_decision_subset/{subset}",
    shell:
        "Rscript codes/phase4_01_decision_subset.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--meta={input.meta}"

rule phase4_crossborder_ibd:
    input:
        sample_list = rules.phase4_decision_subset.output.sample_list,
        ibd_matrix  = rules.phase2_hmmibd.output.ibd_matrix,
    output:
        ibd_pairs = "results/phase4_01/{subset}/ibd_pairs.tsv",
    params:
        outdir = "results/phase4_01/{subset}",
    shell:
        "Rscript codes/phase4_01_crossborder_ibd.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--sample_list={input.sample_list} --ibd_matrix={input.ibd_matrix}"

rule phase4_treemix:
    input:
        sample_list = rules.phase4_decision_subset.output.sample_list,
    output:
        # The script writes many files; declare a sentinel to anchor the rule.
        sentinel = "results/phase4_02/{subset}/m{m}/runs/.done",
    params:
        outdir = "results/phase4_02/{subset}/m{m}",
        n_reps = 10,
    conda:
        "../../envs/treemix.yaml"
    shell:
        "Rscript codes/phase4_02a_treemix.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--m={wildcards.m} --n_reps={params.n_reps} && "
        "touch {output.sentinel}"

def _treemix_run_dirs(wildcards):
    return [f"results/phase4_02/{wildcards.subset}/m{m}/runs/.done" for m in config["m_range"]]

rule phase4_treemix_optm:
    input:
        sentinels = _treemix_run_dirs,
    output:
        summary = "results/phase4_02/{subset}/optm_summary.tsv",
        plot    = "results/phase4_02/{subset}/optm_curve.pdf",
        best_m  = "results/phase4_02/{subset}/best_m.txt",
    params:
        outdir      = "results/phase4_02/{subset}",
        runs_dirs   = lambda w: ",".join(
            f"results/phase4_02/{w.subset}/m{m}" for m in config["m_range"]
        ),
        n_bootstrap = 100,
    conda:
        "../../envs/treemix.yaml"
    shell:
        "Rscript codes/phase4_02b_optm.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--runs_dirs={params.runs_dirs} --n_bootstrap={params.n_bootstrap}"

rule phase4_resistance:
    input:
        sample_list = rules.phase4_decision_subset.output.sample_list,
        meta        = rules.phase1_annotate.output.meta,
    output:
        haplotype_freq = "results/phase4_03/{subset}/haplotype_freq.tsv",
        ancestry       = "results/phase4_03/{subset}/ancestry.tsv",
    params:
        outdir = "results/phase4_03/{subset}",
    shell:
        "Rscript codes/phase4_03_resistance_haplotype_ancestry.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--sample_list={input.sample_list} --meta={input.meta}"
```

- [ ] **Step 5: Verify DAG**

Run: `snakemake -n "results/phase4_02/withsenegal/optm_curve.pdf"`
Expected: lists decision_subset, all m runs (m=0..5), then optm aggregation. ~7 jobs.

Run: `snakemake -n phase4_resistance --config subsets="['withsenegal']"` (sanity-check no cross-subset leakage).

- [ ] **Step 6: Touch or run**

Treemix and resistance outputs for `withsenegal` exist under `results/phase4_*/` — use `snakemake --touch` to mark them current.

- [ ] **Step 7: Commit**

```bash
git add codes/phase4_*.R workflow/rules/phase4.smk
git commit -m "Wire Phase 4 (decision subset, cross-border IBD, TreeMix m sweep, OptM, resistance) into Snakemake"
```

---

## Task 8: Phase 5 rules

**Files:**
- Modify: `codes/phase5_01_attribution_model.R`, `codes/phase5_02_sensitivity.R`
- Modify: `workflow/rules/phase5.smk`

- [ ] **Step 1: Add header to `phase5_01_attribution_model.R`**

```r
output_dir   <- .parse_arg("output_dir", "results/phase5_01/withsenegal")
subset_name  <- .parse_arg("subset", "withsenegal")
p_trans_in   <- .parse_arg("p_trans", "results/phase3_02/P_trans.tsv")
ibd_pairs_in <- .parse_arg("ibd_pairs", "results/phase4_01/withsenegal/ibd_pairs.tsv")
ancestry_in  <- .parse_arg("ancestry", "results/phase4_03/withsenegal/ancestry.tsv")
```

- [ ] **Step 2: Add header to `phase5_02_sensitivity.R`**

```r
output_dir       <- .parse_arg("output_dir", "results/phase5_02")
posterior_files  <- strsplit(.parse_arg("posterior_files", ""), ",", fixed = TRUE)[[1]]
```

Replace the internal subset-iteration over hardcoded paths with iteration over `posterior_files`.

- [ ] **Step 3: Write Phase 5 rules in `workflow/rules/phase5.smk`**

```python
rule phase5_attribution:
    input:
        p_trans   = rules.phase3_directional_ibd.output.p_trans,
        ibd_pairs = rules.phase4_crossborder_ibd.output.ibd_pairs,
        ancestry  = rules.phase4_resistance.output.ancestry,
    output:
        posterior   = "results/phase5_01/{subset}/posterior.rds",
        attribution = "results/phase5_01/{subset}/attribution.tsv",
    params:
        outdir = "results/phase5_01/{subset}",
    shell:
        "Rscript codes/phase5_01_attribution_model.R "
        "--output_dir={params.outdir} --subset={wildcards.subset} "
        "--p_trans={input.p_trans} --ibd_pairs={input.ibd_pairs} "
        "--ancestry={input.ancestry}"

rule phase5_sensitivity:
    input:
        posteriors = expand("results/phase5_01/{subset}/posterior.rds", subset=config["subsets"]),
    output:
        grid = "results/phase5_02/scenario_grid.tsv",
    params:
        outdir = "results/phase5_02",
        posterior_files_csv = lambda w, input: ",".join(input.posteriors),
    shell:
        "Rscript codes/phase5_02_sensitivity.R "
        "--output_dir={params.outdir} --posterior_files={params.posterior_files_csv}"
```

- [ ] **Step 4: Verify DAG**

Run: `snakemake -n phase5_sensitivity`
Expected: lists attribution rule expanded over all subsets, then sensitivity aggregation.

- [ ] **Step 5: Touch or run**

- [ ] **Step 6: Commit**

```bash
git add codes/phase5_*.R workflow/rules/phase5.smk
git commit -m "Wire Phase 5 (attribution model, sensitivity grid) into Snakemake"
```

---

## Task 9: Top-level convenience targets

**Files:**
- Modify: `Snakefile`

- [ ] **Step 1: Populate `rule all` and add convenience targets**

Replace the placeholder `rule all` in `Snakefile` with:

```python
# === Convenience targets ===================================================

rule all:
    input:
        # Phase 1
        rules.phase1_annotate.output.meta,
        # Phase 2 (canonical K only — sweep targets handle the rest)
        rules.phase2_fst_ibd_mantel.output.fst,
        f"results/phase2_02_lea/K{canonical_k()}/qmatrix.tsv",
        rules.phase2_hmmibd.output.ibd_matrix,
        # Phase 3 (canonical K only)
        f"results/phase3_01/K{canonical_k()}/qmatrix.tsv",
        rules.phase3_directional_ibd.output.p_trans,
        rules.phase3_isolation.output.idx,
        # Phase 4 (canonical subset only for cross-border + TreeMix + resistance)
        f"results/phase4_01/{canonical_subset()}/ibd_pairs.tsv",
        f"results/phase4_02/{canonical_subset()}/optm_curve.pdf",
        f"results/phase4_03/{canonical_subset()}/ancestry.tsv",
        # Phase 5 (all subsets — sensitivity aggregates over them)
        rules.phase5_sensitivity.output.grid,

rule phase1:
    input:
        rules.phase1_annotate.output.meta,

rule phase2:
    input:
        rules.phase2_fst_ibd_mantel.output.fst,
        f"results/phase2_02_lea/K{canonical_k()}/qmatrix.tsv",
        rules.phase2_hmmibd.output.ibd_matrix,

rule phase3:
    input:
        f"results/phase3_01/K{canonical_k()}/qmatrix.tsv",
        rules.phase3_directional_ibd.output.p_trans,
        rules.phase3_isolation.output.idx,

rule phase4:
    input:
        f"results/phase4_01/{canonical_subset()}/ibd_pairs.tsv",
        f"results/phase4_02/{canonical_subset()}/optm_curve.pdf",
        f"results/phase4_03/{canonical_subset()}/ancestry.tsv",

rule phase5:
    input:
        rules.phase5_sensitivity.output.grid,

# === Sweep targets =========================================================

rule phase2_lea_sweep:
    input:
        "results/phase2_02_lea/ce_vs_K.pdf",
        expand("results/phase2_02_lea/K{k}/qmatrix.tsv", k=config["k_range"]),

rule phase4_treemix_sweep:
    input:
        f"results/phase4_02/{canonical_subset()}/optm_curve.pdf",
        expand("results/phase4_02/{subset}/m{m}/runs/.done",
               subset=[canonical_subset()], m=config["m_range"]),

# === Cleanup ===============================================================

rule clean_phase1: shell: "rm -rf results/phase1_*"
rule clean_phase2: shell: "rm -rf results/phase2_*"
rule clean_phase3: shell: "rm -rf results/phase3_*"
rule clean_phase4: shell: "rm -rf results/phase4_*"
rule clean_phase5: shell: "rm -rf results/phase5_*"
```

- [ ] **Step 2: Verify all convenience targets**

For each of `all`, `phase1`, `phase2`, `phase3`, `phase4`, `phase5`, `phase2_lea_sweep`, `phase4_treemix_sweep`:

Run: `snakemake -n <target>`
Expected: builds a DAG (may show "nothing to be done" if outputs already touched). No "MissingInputException" or "AmbiguousRuleException".

- [ ] **Step 3: Commit**

```bash
git add Snakefile
git commit -m "Add top-level targets (all, phase1-5, sweeps, clean)"
```

---

## Task 10: End-to-end smoke test

- [ ] **Step 1: Full DAG dry-run**

Run: `snakemake -n all --quiet`
Expected: prints final job-count summary, no errors. If existing outputs are touched, may show 0 jobs.

- [ ] **Step 2: Force a small re-run to test live execution**

Pick one cheap rule (e.g. `phase3_isolation`) and force-rerun:
Run: `snakemake --cores 1 --force phase3_isolation`
Expected: rule executes, writes `results/phase3_03/isolation_index.tsv`, exit 0.

Inspect: `head results/phase3_03/isolation_index.tsv` — should match the columns/structure of the prior manual output (if any).

- [ ] **Step 3: Confirm conda envs build on demand**

Run: `snakemake --use-conda --conda-create-envs-only --cores 1 phase2_hmmibd phase3_supervised_admixture "results/phase4_02/withsenegal/m0/runs/.done"`
Expected: each conda env (hmmibd, admixture, treemix) builds successfully. May take several minutes on first run.

- [ ] **Step 4: Lint pass**

Run: `snakemake --lint`
Expected: no errors. Warnings about unused params or missing `threads:` directives are acceptable.

- [ ] **Step 5: Final commit (if any fixups needed)**

If smoke testing turned up any bugs, fix them and:

```bash
git add <fixed files>
git commit -m "Fix issues found in end-to-end smoke test"
```

If nothing changed: no commit needed.

---

## Notes for the implementer

- **Existing outputs are precious.** Before any `snakemake --cores N <target>` that would rebuild, prefer `snakemake --touch <target>` if the existing files under `results/phaseN_NN/` are still scientifically valid. CLAUDE.md flags that hmmIBD inputs were a previous footgun (memory: `feedback_hmmibd_input_header.md`) — don't re-run hmmIBD casually.
- **The `Sys.sleep(0.1)` in `get_ibs_matrix()` (clustering scratchpad)** is out of scope — that file is not in the DAG.
- **The `phase4_02_treemix.R` original is left in place** as a scratchpad. The new `02a` / `02b` scripts are the canonical ones.
- **R script path-replacements must be exhaustive.** Use `Grep` for each script to find every `"results/phaseN_NN"` literal before declaring the header-adding step done.
- **macOS-BSD awk quirk** (memory: `feedback_macos_awk.md`) — if any Snakemake `shell:` command uses awk, verify it works on BSD awk; otherwise use `gawk` only after confirming it's installed.
- **Treemix env path** (memory: `reference_conda_envs.md`) — the user's conda installation is `miniforge3`. Snakemake should pick it up automatically via `which conda` if the shell is configured; otherwise pass `--conda-frontend mamba` or `--conda-prefix /Users/mdiop/miniforge3/envs`.
