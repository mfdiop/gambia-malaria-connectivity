# Phase 6 — Reviewer Extensions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new Phase 6 to the Snakemake pipeline delivering five analyses that close the four code gaps and two documentation gaps raised in the supervisor-style review of the Objective 3 proposal — plus an exploratory EEMS migration surface.

**Architecture:** Five new `codes/phase6_NN_*.R` scripts, each parameterised by the existing `commandArgs` + `.parse_arg` pattern. A new `workflow/rules/phase6.smk` registers one rule per script and wires inputs from Phases 1–3. A single surgical edit to `phase2.smk` exposes hmmIBD's already-written per-segment file. Two new conda envs (`chromopainter.yaml`, `eems.yaml`) cover the only non-R tools.

**Tech Stack:** R + tidyverse / vcfR / igraph / ggraph / circlize, ChromoPainter v2 + fineSTRUCTURE (bioconda), EEMS C++ binary + rEEMSplots (bioconda + GitHub), Snakemake ≥7, miniforge3 conda.

**Design source:** `docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md`. Read it before starting Task 0.

**Testing model.** This is an analysis project with no test suite. Per CLAUDE.md, the "test" for each rule is:

1. `snakemake --lint` passes.
2. `snakemake -n <target>` builds a sensible DAG (correct dependencies, no missing inputs).
3. The script's interactive default args resolve to existing files under `results/`.
4. Where a prior manual run exists (e.g., the IBD network in `codes/clustering_approaches_v1.R`), spot-check the new output is in the same ballpark.

**Commit cadence.** Per the user's standing instruction, **do not auto-commit**. Each task ends at a clean save point; the user will batch-commit at their discretion.

---

## File Structure

**Create:**
- `codes/phase6_01_haplotype_coancestry.R` — ChromoPainter + fineSTRUCTURE on monoclonal Gambia
- `codes/phase6_02_ibd_tracts.R` — parse hmmIBD per-segment output, tract-length structure
- `codes/phase6_03_network_topology.R` — IBD network communities, betweenness, bridges (+ directed bonus)
- `codes/phase6_04_symmetric_migration_null.R` — stratified label-permutation null for Phase 3.2 asymmetry
- `codes/phase6_05_eems.R` — EEMS migration / diversity surface
- `workflow/rules/phase6.smk` — five rules
- `envs/chromopainter.yaml` — bioconda chromopainter + finestructure
- `envs/eems.yaml` — bioconda eems
- `docs/superpowers/plans/2026-06-03-phase6-reviewer-extensions.md` — this plan

**Modify:**
- `CLAUDE.md` — add three-level importation table + source-sink framing
- `notes/objective3_proposal.md` — add Phase 6 addendum + deferred-methods note
- `workflow/rules/phase2.smk` — add `hmm_segments` to `rule phase2_hmmibd.output`
- `Snakefile` — `include: phase6.smk`, extend `rule all`, add `rule phase6`, add `rule clean_phase6`

**Untouched (do not edit):** all phase1–phase5 R scripts, all other rule files, the existing spec.

---

## Header Template (reused across all Phase 6 R scripts)

Every Phase 6 R script starts with this block, identical to Phases 1–5:

```r
# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
# Script-specific args parsed below.
# --- end snakemake interface ------------------------------------------------
```

Each task lists the exact `.parse_arg(...)` lines for that script.

---

## Task 0: Documentation — three-level importation table + source-sink framing in CLAUDE.md

**Files:**
- Modify: `CLAUDE.md` (append a new section after the existing "Conventions worth respecting" section)

- [ ] **Step 1: Read CLAUDE.md end to find insertion point**

Run: read the bottom of `CLAUDE.md` and note the last section. Insert after it.

- [ ] **Step 2: Append the three-level importation table and source-sink framing**

Append this block verbatim to `CLAUDE.md`:

```markdown
## Conceptual framework: importation vs local transmission

The chapter operationalises three distinct levels of "importation". Conflating them is the single largest interpretation hazard in this work (and the central reviewer concern). Use this table whenever you reach for the word "importation" in code comments, summaries, or figure captions:

| Level | What it measures | Timescale | Method in this pipeline |
|---|---|---|---|
| Historical connectivity | Long-term shared ancestry / drift | Years–decades | Phase 2.1 FST + Mantel; Phase 2.2 PCA + sNMF |
| Recent parasite migration | Recently transmitted shared ancestry | Months–years | Phase 2.3 hmmIBD pairwise IBD; Phase 6.1 haplotype co-ancestry; Phase 6.2 IBD tract lengths |
| Active transmission importation | Source-to-recipient event within the parasite generation time | Weeks–months | Phase 3.2 P_trans; Phase 4.1 cross-border decision-era subset; Phase 5.1 Bayesian attribution |

ADMIXTURE / sNMF clusters live at level 1 and **must not** be cited as evidence of level 3.

### Source–sink framing

The chapter's central question is operationalised as:

- **Self-sustaining local transmission:** R_local > 1 in western Gambia — persistence maintained internally.
- **Importation-maintained sink:** R_local < 1 in western Gambia — persistence sustained by parasite influx from eastern Gambia and Senegal.

Phase 5.1 Bayesian attribution is the inferential bridge: P(H_local), P(H_east), P(H_senegal) is the source–sink decomposition.
```

- [ ] **Step 3: Verify the edit landed**

Run: open `CLAUDE.md` and confirm the new section is at the bottom, with the table rendering correctly (markdown pipes aligned).

---

## Task 1: Documentation — Phase 6 addendum in the proposal

**Files:**
- Modify: `notes/objective3_proposal.md` (append a new section after the existing Phase 5 description, before the tools/figures appendices)

- [ ] **Step 1: Locate the insertion point**

Find the end of the Phase 5 description in `notes/objective3_proposal.md`. The Phase 6 addendum goes immediately after, before the "Expected figures" / "Tools" appendices.

- [ ] **Step 2: Insert the Phase 6 addendum**

Insert this block verbatim:

````markdown
## Phase 6 — Reviewer-driven extensions

A supervisor-style review raised four code gaps and one out-of-scope exploration. Phase 6 closes them, all on the monoclonal Gambia subset built in Phase 1:

| Step | Analysis | Reviewer concern addressed |
|---|---|---|
| 6.1 | ChromoPainter + fineSTRUCTURE haplotype co-ancestry | sNMF / ADMIXTURE alone over-interprets long-term ancestry as recent migration |
| 6.2 | IBD tract-length structure (parsed from hmmIBD `.hmm.txt`) | Total IBD conflates recent and ancient sharing |
| 6.3 | IBD network topology (Louvain communities, betweenness, bridges) | Pair-level IBD misses hubs / bridges / community structure |
| 6.4 | Stratified label-permutation null for Phase 3.2 directional asymmetry | No explicit symmetric-migration null exists for the asymmetry claim |
| 6.5 | EEMS spatial migration / diversity surface (exploratory) | Where, on the map, does effective migration deviate from isolation by distance? |

### Deferred methods (post-manuscript)

The reviewer also flagged MASCOT/BASTA (structured coalescent), Relate (genome-wide genealogies), and tsinfer/tsdate (ARG inference). All three are explicitly deferred to post-manuscript exploration. EEMS is the only one of the four advanced spatial-genomic methods admitted into the pipeline for this chapter, on the grounds that it extends the Phase 2 IBD-by-distance story most cleanly and is the only one that finishes in PhD-realistic time without Pf-phasing pain.
````

- [ ] **Step 3: Verify**

Open the file, confirm the table renders and the new section sits before the appendices.

---

## Task 2: Conda env files

**Files:**
- Create: `envs/chromopainter.yaml`
- Create: `envs/eems.yaml`

- [ ] **Step 1: Write `envs/chromopainter.yaml`**

```yaml
name: chromopainter
channels:
  - bioconda
  - conda-forge
dependencies:
  - chromopainter
  - finestructure
  - bcftools
  - python>=3.9
```

- [ ] **Step 2: Write `envs/eems.yaml`**

```yaml
name: eems
channels:
  - bioconda
  - conda-forge
dependencies:
  - eems
  - bcftools
  - python>=3.9
  - cxx-compiler   # source-build fallback if the binary fails to link on macOS
```

- [ ] **Step 3: Sanity check the YAML syntax**

Run: `python -c "import yaml; yaml.safe_load(open('envs/chromopainter.yaml')); yaml.safe_load(open('envs/eems.yaml'))"`
Expected: no output (clean exit). If yaml errors, fix indentation.

Note: do **not** create the conda envs yet (per the user's "skip env builds" preference from the Phase 5 smoke-test choice). Env creation happens at first real run.

---

## Task 3: Expose hmmIBD per-segment output in `phase2.smk`

The hmmIBD binary already writes `results/phase2_03/hmmibd.hmm.txt` on every run (see `codes/phase2_03_hmmibd.R:141-150` — `system2(hmmibd_bin, args = c("-o", out_prefix, ...))` produces both `.hmm.txt` and `.hmm_fract.txt`). The file just isn't declared as a Snakemake output. Phase 6.2 needs it as a tracked input, so we add the declaration.

**Files:**
- Modify: `workflow/rules/phase2.smk` (add one line to `rule phase2_hmmibd.output`)

- [ ] **Step 1: Edit `workflow/rules/phase2.smk`**

In `rule phase2_hmmibd`, change the `output:` block from:

```python
    output:
        hmm_fract  = "results/phase2_03/hmmibd.hmm_fract.txt",
        hap_input  = "results/phase2_03/hmmibd_input.txt",
        pair_class = "results/phase2_03/pair_ibd_classified.tsv",
        ibd_matrix = "results/phase2_03/ibd_matrix.rds",
        summary    = "results/phase2_03/phase2_03_summary.txt",
```

to:

```python
    output:
        hmm_fract    = "results/phase2_03/hmmibd.hmm_fract.txt",
        hmm_segments = "results/phase2_03/hmmibd.hmm.txt",
        hap_input    = "results/phase2_03/hmmibd_input.txt",
        pair_class   = "results/phase2_03/pair_ibd_classified.tsv",
        ibd_matrix   = "results/phase2_03/ibd_matrix.rds",
        summary      = "results/phase2_03/phase2_03_summary.txt",
```

- [ ] **Step 2: Verify the DAG still resolves**

Run: `snakemake -n results/phase2_03/hmmibd.hmm.txt`
Expected: the DAG lists `phase2_hmmibd` as the producing rule with no errors. (If hmmIBD has already been run, snakemake will report "Nothing to be done" — that's fine; the output is also satisfied.)

No script edit needed — the file is already written by hmmIBD itself.

---

## Task 4: Phase 6.1 — Haplotype co-ancestry script

**Files:**
- Create: `codes/phase6_01_haplotype_coancestry.R`

- [ ] **Step 1: Create the script with the standard header**

```r
# =============================================================================
# Phase 6 — Step 6.1: Haplotype co-ancestry (ChromoPainter v2 + fineSTRUCTURE)
# =============================================================================
# Per docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md §5.1:
#   - Resolve recent co-ancestry between Gambian monoclonal P. falciparum
#     samples without the long-term-ancestry confound of sNMF/ADMIXTURE.
#   - Pf monoclonal samples are haploid, so no SHAPEIT phasing is required;
#     ChromoPainter is run in haploid mode (-j 1).
#   - fineSTRUCTURE MCMC clusters the chunkcounts matrix.
#
# Sample set : results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz
#              (filtered further to COI == 1 samples via meta_rds$coi).
#
# Output : results/phase6_01/
#            chromopainter_inputs/  .phase, .recombfile, .idfile
#            em_results/            EM-estimated Ne / mu per chromosome
#            chunkcounts.tsv        N x N chunkcount matrix
#            coancestry_matrix.tsv  long-format N x N matrix
#            chunkcounts_heatmap.pdf
#            fs_mcmc.xml            fineSTRUCTURE MCMC trace
#            fs_tree.newick         fineSTRUCTURE tree
#            phase6_01_summary.txt
#            phase6_01_results.rds
#
# Tools  : bioconda chromopainter + finestructure (see envs/chromopainter.yaml).
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase6_01")
mono_vcf   <- .parse_arg("mono_vcf",   "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
meta_arg   <- .parse_arg("meta",       "results/phase1_03/sample_master_metadata.rds")
em_chrom   <- as.integer(.parse_arg("em_chrom",   "11"))
mcmc_burn  <- as.integer(.parse_arg("mcmc_burn",  "1000000"))
mcmc_post  <- as.integer(.parse_arg("mcmc_post",  "1000000"))
mcmc_thin  <- as.integer(.parse_arg("mcmc_thin",  "1000"))
tree_iter  <- as.integer(.parse_arg("tree_iter",  "100000"))
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vcfR)
  library(pheatmap)
  library(ape)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cp_dir <- file.path(output_dir, "chromopainter_inputs"); dir.create(cp_dir, showWarnings = FALSE)
em_dir <- file.path(output_dir, "em_results");           dir.create(em_dir, showWarnings = FALSE)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }
log_msg("=== Phase 6.1: haplotype co-ancestry (ChromoPainter + fineSTRUCTURE) ===")
```

- [ ] **Step 2: Restrict samples to COI == 1 and write ChromoPainter inputs**

Append:

```r
# ============================================================================
# 1. Load VCF + metadata; restrict to true monoclonal (COI == 1)
# ============================================================================
vcf  <- read.vcfR(mono_vcf, verbose = FALSE)
meta <- readRDS(meta_arg)

mono_ids <- meta %>% filter(coi == 1) %>% pull(SampleID)
vcf_keep <- intersect(colnames(vcf@gt)[-1], mono_ids)
log_msg(sprintf("VCF samples: %d ; meta COI==1: %d ; intersect: %d",
                ncol(vcf@gt) - 1, length(mono_ids), length(vcf_keep)))
stopifnot(length(vcf_keep) >= 50)

vcf <- vcf[, c("FORMAT", vcf_keep)]

# ============================================================================
# 2. Convert to ChromoPainter .phase / .recombfile / .idfile (haploid mode)
# ============================================================================
# Per-chromosome files. ChromoPainter expects one .phase per chromosome:
#   line 1: ploidy (1 for haploid)
#   line 2: n_haps
#   line 3: n_snps
#   line 4: P <space-separated bp positions>
#   lines 5..: one haplotype per line, 0/1 string
PF_CHRS <- sprintf("Pf3D7_%02d_v3", 1:14)

gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
encode_hap <- function(x) {
  out <- rep(NA_integer_, length(x))
  out[x %in% c("0", "0/0", "0|0")] <- 0L
  out[x %in% c("1", "1/1", "1|1")] <- 1L
  out
}
geno <- matrix(encode_hap(as.vector(gt_chr)),
               nrow = nrow(gt_chr), ncol = ncol(gt_chr),
               dimnames = dimnames(gt_chr))

# ChromoPainter cannot handle missing data; drop SNPs with any NA
keep_snp <- complete.cases(geno)
log_msg(sprintf("Complete-case SNPs: %d / %d", sum(keep_snp), nrow(geno)))
geno <- geno[keep_snp, , drop = FALSE]
fix  <- vcf@fix[keep_snp, , drop = FALSE]

writeLines(vcf_keep, file.path(cp_dir, "samples.idfile"))

for (ch_idx in seq_along(PF_CHRS)) {
  ch_name <- PF_CHRS[ch_idx]
  rows <- fix[, "CHROM"] == ch_name
  if (sum(rows) < 10) next
  pos  <- as.integer(fix[rows, "POS"])
  g_ch <- geno[rows, , drop = FALSE]

  phase_path <- file.path(cp_dir, sprintf("chr%02d.phase", ch_idx))
  con <- file(phase_path, "w")
  writeLines(c("1", as.character(ncol(g_ch)), as.character(sum(rows)),
               paste("P", paste(pos, collapse = " "))), con)
  # One haplotype per column → one row per sample in the output
  for (j in seq_len(ncol(g_ch))) {
    writeLines(paste(g_ch[, j], collapse = ""), con)
  }
  close(con)

  # Uniform recombination rate placeholder: 1e-8 per bp between adjacent SNPs.
  # (ChromoPainter EM step will re-estimate Ne; a flat map is the standard
  #  starting point for Pf.)
  recomb_path <- file.path(cp_dir, sprintf("chr%02d.recombfile", ch_idx))
  dgap <- c(diff(pos), 0)
  rec  <- data.frame(start_position = pos, recom_rate_perbp = 1e-8)
  write.table(rec, recomb_path, sep = " ", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}
log_msg(sprintf("ChromoPainter inputs written to %s", cp_dir))
```

- [ ] **Step 3: Run the EM step on one chromosome, then full painting**

Append:

```r
# ============================================================================
# 3. EM step on em_chrom to estimate Ne and mu, then genome-wide painting
# ============================================================================
cp_bin <- Sys.getenv("CHROMOPAINTER", "ChromoPainterv2")

run_em <- function() {
  phase_file <- file.path(cp_dir, sprintf("chr%02d.phase",   em_chrom))
  recomb_file <- file.path(cp_dir, sprintf("chr%02d.recombfile", em_chrom))
  out_prefix <- file.path(em_dir, sprintf("em_chr%02d", em_chrom))
  args <- c("-g", phase_file,
            "-r", recomb_file,
            "-t", file.path(cp_dir, "samples.idfile"),
            "-j",                       # haploid mode
            "-i", "10",                 # 10 EM iterations
            "-in", "-iM",               # estimate Ne and mu
            "-o", out_prefix)
  res <- system2(cp_bin, args = args, stdout = TRUE, stderr = TRUE)
  writeLines(res, paste0(out_prefix, ".log"))
  # Parse Ne and mu from the EM log
  ne  <- as.numeric(sub(".*Ne final estimate ", "",
                        grep("Ne final estimate", res, value = TRUE)[1]))
  mu  <- as.numeric(sub(".*mu final estimate ", "",
                        grep("mu final estimate", res, value = TRUE)[1]))
  list(ne = ne, mu = mu)
}
em <- run_em()
log_msg(sprintf("EM on chr%02d: Ne = %.3g, mu = %.3g", em_chrom, em$ne, em$mu))

# Full painting genome-wide with the EM-fitted Ne / mu
paint_dir <- file.path(output_dir, "paint"); dir.create(paint_dir, showWarnings = FALSE)
for (ch_idx in seq_along(PF_CHRS)) {
  phase_file <- file.path(cp_dir, sprintf("chr%02d.phase", ch_idx))
  if (!file.exists(phase_file)) next
  recomb_file <- file.path(cp_dir, sprintf("chr%02d.recombfile", ch_idx))
  out_prefix  <- file.path(paint_dir, sprintf("chr%02d", ch_idx))
  args <- c("-g", phase_file,
            "-r", recomb_file,
            "-t", file.path(cp_dir, "samples.idfile"),
            "-j", "-a", "0", "0",                      # all donors x all recipients
            "-n", sprintf("%.6f", em$ne),
            "-M", sprintf("%.6e", em$mu),
            "-o", out_prefix)
  system2(cp_bin, args = args, stdout = FALSE, stderr = FALSE)
}
log_msg("Painting complete; combining chunkcounts across chromosomes.")
```

- [ ] **Step 4: Combine chunkcounts and run fineSTRUCTURE**

Append:

```r
# ============================================================================
# 4. Sum chunkcounts across chromosomes -> N x N matrix
# ============================================================================
chunk_files <- list.files(paint_dir, pattern = "\\.chunkcounts\\.out$",
                          full.names = TRUE)
stopifnot(length(chunk_files) > 0)
read_chunk <- function(f) {
  m <- as.matrix(read.table(f, header = TRUE, row.names = 1, check.names = FALSE))
  m
}
mats <- lapply(chunk_files, read_chunk)
cc   <- Reduce(`+`, mats)
# Align row/col order
ord <- sort(rownames(cc))
cc  <- cc[ord, ord]
write.table(as.data.frame(cc) %>% rownames_to_column("SampleID"),
            file.path(output_dir, "chunkcounts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cc, file.path(output_dir, "coancestry_matrix.rds"))

# ============================================================================
# 5. fineSTRUCTURE MCMC + tree
# ============================================================================
fs_bin   <- Sys.getenv("FINESTRUCTURE", "fs")
fs_input <- file.path(output_dir, "chunkcounts.tsv")
fs_mcmc  <- file.path(output_dir, "fs_mcmc.xml")
fs_tree  <- file.path(output_dir, "fs_tree.xml")

# fineSTRUCTURE expects a chunkcounts file with the cfg header; build it.
fs_cc_path <- file.path(output_dir, "chunkcounts.cfg")
con <- file(fs_cc_path, "w")
writeLines(sprintf("#Cfactor 1"), con)
writeLines(c(paste(c("Recipient", colnames(cc)), collapse = " ")), con)
for (i in seq_len(nrow(cc))) {
  writeLines(paste(c(rownames(cc)[i], format(cc[i, ], scientific = FALSE)),
                   collapse = " "), con)
}
close(con)

system2(fs_bin, c("mcmc", "-x", mcmc_burn, "-y", mcmc_post, "-z", mcmc_thin,
                  fs_cc_path, fs_mcmc),
        stdout = FALSE, stderr = FALSE)
system2(fs_bin, c("tree", "-x", tree_iter, fs_cc_path, fs_mcmc, fs_tree),
        stdout = FALSE, stderr = FALSE)
log_msg("fineSTRUCTURE MCMC + tree complete.")
```

- [ ] **Step 5: Heatmap, summary, persist**

Append:

```r
# ============================================================================
# 6. Heatmap ordered by fineSTRUCTURE tree leaves
# ============================================================================
# Extract leaf order from fs_tree.xml; the fineSTRUCTURE tree is embedded as
# Newick inside the XML's <Tree> element.
tree_lines <- readLines(fs_tree)
nwk_line   <- grep("<Tree>", tree_lines, value = TRUE)[1]
nwk        <- sub(".*<Tree>", "", sub("</Tree>.*", "", nwk_line))
writeLines(nwk, file.path(output_dir, "fs_tree.newick"))
fs_tree_phy <- ape::read.tree(text = nwk)
leaf_order  <- fs_tree_phy$tip.label
# Restrict to samples in the matrix (in case fineSTRUCTURE relabels)
leaf_order  <- intersect(leaf_order, rownames(cc))
cc_ord      <- cc[leaf_order, leaf_order]

ann <- meta %>%
  filter(SampleID %in% leaf_order) %>%
  select(SampleID, region) %>%
  column_to_rownames("SampleID")

pdf(file.path(output_dir, "chunkcounts_heatmap.pdf"), width = 10, height = 9)
pheatmap(log10(cc_ord + 1),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = ann, annotation_col = ann,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Haplotype co-ancestry (log10 chunkcounts), Gambia monoclonal")
dev.off()

# ============================================================================
# 7. Summary + persist
# ============================================================================
saveRDS(list(coancestry = cc, leaf_order = leaf_order,
             em = em, n_samples = length(vcf_keep)),
        file.path(output_dir, "phase6_01_results.rds"))

log_msg(sprintf("Samples: %d ; chunkcount matrix: %d x %d",
                length(vcf_keep), nrow(cc), ncol(cc)))
log_msg(sprintf("EM-estimated Ne = %.3g, mu = %.3g", em$ne, em$mu))
writeLines(log_lines, file.path(output_dir, "phase6_01_summary.txt"))
message("Done. Outputs in ", output_dir)
```

- [ ] **Step 6: Smoke-check the script parses**

Run: `Rscript --vanilla -e 'parse("codes/phase6_01_haplotype_coancestry.R")'`
Expected: no output (clean parse). If syntax errors, fix them.

---

## Task 5: Phase 6.2 — IBD tract-length script

**Files:**
- Create: `codes/phase6_02_ibd_tracts.R`

- [ ] **Step 1: Write the full script**

```r
# =============================================================================
# Phase 6 — Step 6.2: IBD tract-length structure
# =============================================================================
# Per docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md §5.2:
#   - Parse hmmIBD's per-segment output (hmmibd.hmm.txt).
#   - Long tracts = recent shared ancestry; short tracts = ancient.
#   - Compare within-region vs cross-region tract length distributions.
#   - Focus contrast: west-west vs west-east.
#
# Output : results/phase6_02/
#            tract_table.tsv                long-form per-segment table
#            tract_summary_by_regionpair.tsv
#            wilcoxon_within_vs_cross.tsv
#            tract_density_facet.pdf
#            tract_length_vs_total_ibd.pdf
#            phase6_02_summary.txt
#            phase6_02_results.rds
#
# Tools  : pure R (tidyverse, ggplot2).
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir   <- .parse_arg("output_dir",   "results/phase6_02")
segments_arg <- .parse_arg("segments",     "results/phase2_03/hmmibd.hmm.txt")
fract_arg    <- .parse_arg("fract",        "results/phase2_03/hmmibd.hmm_fract.txt")
meta_arg     <- .parse_arg("meta",         "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load hmmIBD per-segment table + metadata
# ============================================================================
# Schema check: hmmIBD .hmm.txt columns are
#   sample1  sample2  chr  start  end  different  Nsnp
# Fail loudly if the schema drifts.
seg <- read_tsv(segments_arg, show_col_types = FALSE)
expected_cols <- c("sample1", "sample2", "chr", "start", "end", "different", "Nsnp")
stopifnot(all(expected_cols %in% colnames(seg)))

meta <- readRDS(meta_arg)
meta_min <- meta %>% select(SampleID, region)

# `different == 0` means IBD; `different == 1` means non-IBD. Keep IBD only.
seg_ibd <- seg %>%
  filter(different == 0) %>%
  mutate(length_bp = end - start + 1L) %>%
  rename(SampleID1 = sample1, SampleID2 = sample2)

log_msg(sprintf("Per-segment rows total: %d ; IBD segments: %d",
                nrow(seg), nrow(seg_ibd)))

# ============================================================================
# 2. Annotate region pair
# ============================================================================
seg_ibd <- seg_ibd %>%
  left_join(meta_min %>% rename(region1 = region), by = c("SampleID1" = "SampleID")) %>%
  left_join(meta_min %>% rename(region2 = region), by = c("SampleID2" = "SampleID")) %>%
  mutate(
    r1 = pmin(as.character(region1), as.character(region2)),
    r2 = pmax(as.character(region1), as.character(region2)),
    region_pair = paste(r1, r2, sep = "-"),
    pair_type   = ifelse(r1 == r2, "within-region", "cross-region"))

write_tsv(seg_ibd, file.path(output_dir, "tract_table.tsv"))

# ============================================================================
# 3. Summary by region pair
# ============================================================================
summary_rp <- seg_ibd %>%
  group_by(region_pair, pair_type) %>%
  summarise(n_segments  = n(),
            n_pairs     = n_distinct(paste(SampleID1, SampleID2)),
            mean_length = mean(length_bp),
            median_length = median(length_bp),
            p90_length  = quantile(length_bp, 0.9),
            .groups = "drop")
write_tsv(summary_rp, file.path(output_dir, "tract_summary_by_regionpair.tsv"))
log_msg("Tract summary by region pair written.")
print(summary_rp)

# ============================================================================
# 4. Wilcoxon tests
# ============================================================================
wt_all <- wilcox.test(length_bp ~ pair_type, data = seg_ibd)

ww <- seg_ibd %>% filter(region_pair == "west-west")    %>% pull(length_bp)
we <- seg_ibd %>% filter(region_pair == "east-west")    %>% pull(length_bp)
wt_we <- if (length(we) > 0 && length(ww) > 0) {
  wilcox.test(ww, we)
} else { list(statistic = NA, p.value = NA) }

wt_tbl <- tibble(
  test = c("within_vs_cross_all", "west_west_vs_west_east"),
  W    = c(unname(wt_all$statistic), unname(wt_we$statistic)),
  p    = c(wt_all$p.value, wt_we$p.value))
write_tsv(wt_tbl, file.path(output_dir, "wilcoxon_within_vs_cross.tsv"))

# ============================================================================
# 5. Plots
# ============================================================================
p_density <- ggplot(seg_ibd, aes(x = length_bp, fill = pair_type)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = scales::label_number(scale = 1e-3, suffix = " kb")) +
  facet_wrap(~ region_pair, scales = "free_y") +
  geom_vline(data = summary_rp,
             aes(xintercept = median_length),
             linetype = "dashed") +
  labs(title = "IBD tract length density by region pair",
       x = "Tract length (log10 bp)", y = "Density") +
  theme_minimal()
ggsave(file.path(output_dir, "tract_density_facet.pdf"), p_density,
       width = 11, height = 8, dpi = 600)

# Scatter: per-pair mean tract length vs total fract_sites_IBD
fract <- read_tsv(fract_arg, show_col_types = FALSE) %>%
  rename(SampleID1 = sample1, SampleID2 = sample2)
pair_mean <- seg_ibd %>%
  group_by(SampleID1, SampleID2, region_pair) %>%
  summarise(mean_tract = mean(length_bp), .groups = "drop") %>%
  left_join(fract %>% select(SampleID1, SampleID2, fract_sites_IBD),
            by = c("SampleID1", "SampleID2"))
p_scatter <- ggplot(pair_mean, aes(x = fract_sites_IBD, y = mean_tract,
                                    colour = region_pair)) +
  geom_point(alpha = 0.6) +
  scale_y_log10(labels = scales::label_number(scale = 1e-3, suffix = " kb")) +
  labs(title = "Per-pair mean tract length vs total IBD",
       x = "fract_sites_IBD", y = "Mean tract length (log10 bp)") +
  theme_minimal()
ggsave(file.path(output_dir, "tract_length_vs_total_ibd.pdf"), p_scatter,
       width = 8, height = 6, dpi = 600)

# ============================================================================
# 6. Persist
# ============================================================================
saveRDS(list(seg_ibd = seg_ibd, summary_rp = summary_rp,
             wilcoxon = wt_tbl, pair_mean = pair_mean),
        file.path(output_dir, "phase6_02_results.rds"))
writeLines(log_lines, file.path(output_dir, "phase6_02_summary.txt"))
message("Done. Outputs in ", output_dir)
```

- [ ] **Step 2: Parse-check the script**

Run: `Rscript --vanilla -e 'parse("codes/phase6_02_ibd_tracts.R")'`
Expected: clean parse.

---

## Task 6: Phase 6.3 — Network topology script

**Files:**
- Create: `codes/phase6_03_network_topology.R`

- [ ] **Step 1: Write the script header + graph construction**

```r
# =============================================================================
# Phase 6 — Step 6.3: IBD network topology
# =============================================================================
# Per docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md §5.3:
#   - Build undirected weighted graphs at two thresholds (>=0.25, >=0.50).
#   - Louvain communities; betweenness centrality; bridge identification.
#   - Bonus: directed graph from Phase 3.2 transmission_pairs.tsv with P_trans.
#
# Output : results/phase6_03/
#            communities_0.25.tsv
#            communities_0.50.tsv
#            betweenness.tsv
#            bridges.tsv
#            network_0.25.pdf
#            network_0.50.pdf
#            directed_indeg_outdeg.tsv
#            phase6_03_summary.txt
#            phase6_03_results.rds
#
# Tools  : pure R (igraph, ggraph, tidyverse).
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase6_03")
pair_arg   <- .parse_arg("pair_class", "results/phase2_03/pair_ibd_classified.tsv")
trans_arg  <- .parse_arg("trans_pairs", "results/phase3_02/transmission_pairs.tsv")
meta_arg   <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

pairs <- read_tsv(pair_arg, show_col_types = FALSE)
meta  <- readRDS(meta_arg) %>% select(SampleID, region, Location, VillageCode)

build_graph <- function(df, threshold) {
  edges <- df %>%
    filter(fract_sites_IBD >= threshold) %>%
    select(SampleID1, SampleID2, weight = fract_sites_IBD)
  nodes <- meta %>%
    filter(SampleID %in% unique(c(edges$SampleID1, edges$SampleID2)))
  graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
}

g25 <- build_graph(pairs, 0.25)
g50 <- build_graph(pairs, 0.50)
log_msg(sprintf("Network @0.25: %d nodes, %d edges", vcount(g25), ecount(g25)))
log_msg(sprintf("Network @0.50: %d nodes, %d edges", vcount(g50), ecount(g50)))
```

- [ ] **Step 2: Append community detection, betweenness, bridges**

```r
# ============================================================================
# Communities + betweenness for each threshold
# ============================================================================
analyse_graph <- function(g, label) {
  if (ecount(g) == 0) return(NULL)
  comm <- cluster_louvain(g, weights = E(g)$weight)
  btw  <- betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE)
  tibble(SampleID = V(g)$name,
         community = membership(comm),
         betweenness = btw,
         region = V(g)$region,
         Location = V(g)$Location,
         threshold = label)
}

a25 <- analyse_graph(g25, "0.25")
a50 <- analyse_graph(g50, "0.50")

if (!is.null(a25))
  write_tsv(a25, file.path(output_dir, "communities_0.25.tsv"))
if (!is.null(a50))
  write_tsv(a50, file.path(output_dir, "communities_0.50.tsv"))

btw_all <- bind_rows(a25, a50) %>%
  select(SampleID, threshold, betweenness, region)
write_tsv(btw_all, file.path(output_dir, "betweenness.tsv"))

# Bridges: top-decile betweenness AND edges span >= 2 regions
find_bridges <- function(g, a) {
  if (is.null(a) || nrow(a) == 0) return(tibble())
  top_q <- quantile(a$betweenness, 0.9, na.rm = TRUE)
  candidates <- a %>% filter(betweenness >= top_q) %>% pull(SampleID)
  bridge_rows <- map_dfr(candidates, function(s) {
    nbrs <- neighbors(g, s)$name
    n_regions <- meta %>%
      filter(SampleID %in% nbrs) %>%
      pull(region) %>% unique() %>% length()
    tibble(SampleID = s, n_regions_touched = n_regions,
           betweenness = a$betweenness[a$SampleID == s])
  }) %>%
    filter(n_regions_touched >= 2)
  bridge_rows
}
br25 <- find_bridges(g25, a25) %>% mutate(threshold = "0.25")
br50 <- find_bridges(g50, a50) %>% mutate(threshold = "0.50")
write_tsv(bind_rows(br25, br50), file.path(output_dir, "bridges.tsv"))
log_msg(sprintf("Bridges @0.25: %d ; @0.50: %d", nrow(br25), nrow(br50)))
```

- [ ] **Step 3: Append force-directed plots + directed bonus + persist**

```r
# ============================================================================
# Force-directed plots
# ============================================================================
plot_network <- function(g, a, label, outfile) {
  if (is.null(a)) { file.create(outfile); return(invisible()) }
  V(g)$community   <- a$community[match(V(g)$name, a$SampleID)]
  V(g)$betweenness <- a$betweenness[match(V(g)$name, a$SampleID)]
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = weight), colour = "grey40") +
    geom_node_point(aes(size = betweenness, fill = factor(community)),
                    shape = 21, colour = "black") +
    scale_size_continuous(range = c(1.5, 8)) +
    labs(title = sprintf("IBD network (threshold = %s)", label),
         fill = "community") +
    theme_void()
  ggsave(outfile, p, width = 10, height = 8, dpi = 600)
}
plot_network(g25, a25, "0.25", file.path(output_dir, "network_0.25.pdf"))
plot_network(g50, a50, "0.50", file.path(output_dir, "network_0.50.pdf"))

# ============================================================================
# Directed bonus: in / out degree from Phase 3.2 transmission_pairs.tsv
# ============================================================================
trans <- if (file.exists(trans_arg)) {
  read_tsv(trans_arg, show_col_types = FALSE)
} else {
  log_msg("No transmission_pairs.tsv found; skipping directed analysis.")
  NULL
}
if (!is.null(trans) && nrow(trans) > 0) {
  trans_edges <- trans %>%
    filter(P_trans >= 0.5) %>%
    transmute(from = source_id, to = recipient_id, weight = P_trans)
  gd <- graph_from_data_frame(trans_edges, directed = TRUE)
  deg_tbl <- tibble(
    SampleID = V(gd)$name,
    in_degree  = degree(gd, mode = "in"),
    out_degree = degree(gd, mode = "out"))
  write_tsv(deg_tbl, file.path(output_dir, "directed_indeg_outdeg.tsv"))
  log_msg(sprintf("Directed network (P_trans >= 0.5): %d nodes, %d edges",
                  vcount(gd), ecount(gd)))
}

# ============================================================================
# Persist
# ============================================================================
saveRDS(list(g25 = g25, g50 = g50, a25 = a25, a50 = a50,
             bridges_25 = br25, bridges_50 = br50),
        file.path(output_dir, "phase6_03_results.rds"))
writeLines(log_lines, file.path(output_dir, "phase6_03_summary.txt"))
message("Done. Outputs in ", output_dir)
```

- [ ] **Step 4: Parse-check**

Run: `Rscript --vanilla -e 'parse("codes/phase6_03_network_topology.R")'`
Expected: clean parse.

---

## Task 7: Phase 6.4 — Symmetric-migration permutation null

**Files:**
- Create: `codes/phase6_04_symmetric_migration_null.R`

- [ ] **Step 1: Write the full script**

```r
# =============================================================================
# Phase 6 — Step 6.4: Symmetric-migration permutation null
# =============================================================================
# Per docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md §5.4:
#   - Observed statistic A_obs = (n_e2w - n_w2e) / (n_e2w + n_w2e),
#     over high-P_trans pairs (>= 0.5).
#   - Null: permute `region` WITHIN year + transmission-season strata
#     (preserves temporal structure; symmetric-migration hypothesis).
#   - 1000 permutations; empirical two-sided p-value.
#
# Output : results/phase6_04/
#            null_distribution.tsv
#            observed_vs_null.tsv
#            observed_vs_null.pdf
#            phase6_04_summary.txt
#            phase6_04_results.rds
#
# Tools  : pure R (tidyverse).
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase6_04")
trans_arg  <- .parse_arg("trans_pairs", "results/phase3_02/transmission_pairs.tsv")
meta_arg   <- .parse_arg("meta", "results/phase1_03/sample_master_metadata.rds")
n_perm     <- as.integer(.parse_arg("n_perm", "1000"))
ptrans_thr <- as.numeric(.parse_arg("ptrans_threshold", "0.5"))
seed       <- as.integer(.parse_arg("seed", "42"))
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load + filter
# ============================================================================
trans <- read_tsv(trans_arg, show_col_types = FALSE)
meta  <- readRDS(meta_arg)

# Sanity: expected columns
stopifnot(all(c("source_id", "recipient_id", "P_trans",
                "source_region", "recipient_region") %in% colnames(trans)))

# Restrict to high-P_trans pairs (consistent with Phase 3.2)
high <- trans %>% filter(P_trans >= ptrans_thr)
log_msg(sprintf("High-P_trans pairs (P_trans >= %.2f): %d", ptrans_thr, nrow(high)))

# ============================================================================
# 2. Observed asymmetry: east -> west vs west -> east
# ============================================================================
count_dir <- function(df) {
  n_e2w <- sum(df$source_region == "east" & df$recipient_region == "west")
  n_w2e <- sum(df$source_region == "west" & df$recipient_region == "east")
  tibble(n_e2w = n_e2w, n_w2e = n_w2e,
         A = if (n_e2w + n_w2e == 0) NA_real_
              else (n_e2w - n_w2e) / (n_e2w + n_w2e))
}
obs <- count_dir(high)
A_obs <- obs$A
log_msg(sprintf("Observed: n_e2w = %d, n_w2e = %d, A_obs = %.4f",
                obs$n_e2w, obs$n_w2e, A_obs))

# ============================================================================
# 3. Stratified permutation: region labels exchangeable within year + season
# ============================================================================
# Build a sample-level region label table with stratum keys (year + season).
# `meta` must carry `year` and `transmission_season` (or equivalent). If the
# season column has a different name, alias it.
season_col <- intersect(c("transmission_season", "season"), colnames(meta))[1]
stopifnot(!is.na(season_col))
year_col   <- intersect(c("Year", "year"), colnames(meta))[1]
stopifnot(!is.na(year_col))

sample_labels <- meta %>%
  select(SampleID,
         region,
         year   = all_of(year_col),
         season = all_of(season_col)) %>%
  mutate(stratum = paste(year, season, sep = "_"))

permute_labels <- function(labels) {
  labels %>%
    group_by(stratum) %>%
    mutate(region = sample(region)) %>%
    ungroup()
}

apply_labels <- function(df, labels) {
  df %>%
    select(source_id, recipient_id, P_trans) %>%
    left_join(labels %>% select(SampleID, region) %>%
                rename(source_region = region),
              by = c("source_id" = "SampleID")) %>%
    left_join(labels %>% select(SampleID, region) %>%
                rename(recipient_region = region),
              by = c("recipient_id" = "SampleID"))
}

A_perm <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  perm_lab <- permute_labels(sample_labels)
  perm_df  <- apply_labels(high, perm_lab)
  A_perm[i] <- count_dir(perm_df)$A
  if (i %% 100 == 0) log_msg(sprintf("Permutation %d/%d", i, n_perm))
}

null_tbl <- tibble(perm = seq_len(n_perm), A_perm = A_perm)
write_tsv(null_tbl, file.path(output_dir, "null_distribution.tsv"))

# ============================================================================
# 4. Two-sided empirical p-value
# ============================================================================
p_right <- mean(A_perm >= A_obs, na.rm = TRUE)
p_left  <- mean(A_perm <= A_obs, na.rm = TRUE)
p_two   <- 2 * min(p_right, p_left)
p_two   <- min(p_two, 1)

obs_tbl <- tibble(
  n_e2w = obs$n_e2w, n_w2e = obs$n_w2e,
  A_obs = A_obs, n_perm = n_perm,
  p_right = p_right, p_left = p_left, p_two_sided = p_two)
write_tsv(obs_tbl, file.path(output_dir, "observed_vs_null.tsv"))
log_msg(sprintf("A_obs = %.4f ; p_two = %.4g", A_obs, p_two))

# ============================================================================
# 5. Plot
# ============================================================================
p_null <- ggplot(null_tbl, aes(x = A_perm)) +
  geom_histogram(bins = 60, fill = "steelblue", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = A_obs, colour = "tomato", linewidth = 1) +
  annotate("text", x = A_obs, y = Inf,
           label = sprintf("A_obs = %.3f\np = %.3g", A_obs, p_two),
           vjust = 1.5, hjust = -0.05, colour = "tomato") +
  labs(title = "Symmetric-migration null vs observed asymmetry",
       subtitle = "Region labels permuted within year + transmission-season strata",
       x = "A = (n_e2w - n_w2e) / (n_e2w + n_w2e)",
       y = "Permutations") +
  theme_minimal()
ggsave(file.path(output_dir, "observed_vs_null.pdf"), p_null,
       width = 8, height = 5, dpi = 600)

# ============================================================================
# 6. Persist
# ============================================================================
saveRDS(list(obs = obs_tbl, null = null_tbl,
             n_perm = n_perm, ptrans_threshold = ptrans_thr,
             seed = seed),
        file.path(output_dir, "phase6_04_results.rds"))
writeLines(log_lines, file.path(output_dir, "phase6_04_summary.txt"))
message("Done. Outputs in ", output_dir)
```

- [ ] **Step 2: Parse-check**

Run: `Rscript --vanilla -e 'parse("codes/phase6_04_symmetric_migration_null.R")'`
Expected: clean parse.

---

## Task 8: Phase 6.5 — EEMS script

**Files:**
- Create: `codes/phase6_05_eems.R`

- [ ] **Step 1: Header + inputs**

```r
# =============================================================================
# Phase 6 — Step 6.5: EEMS effective migration surface (EXPLORATORY)
# =============================================================================
# Per docs/superpowers/specs/2026-06-03-phase6-reviewer-extensions-design.md §5.5:
#   - Compute pairwise allele-difference (Dij) matrix from the monoclonal VCF.
#   - Assign samples to village centroid coordinates.
#   - Habitat outer polygon = convex hull of centroids, expanded by 0.2 deg.
#   - Run 3 EEMS MCMC chains (5M iter, 2M burn, thin 1k).
#   - Post-process with rEEMSplots: posterior mean log(m), log(q), diagnostics.
#
# Output : results/phase6_05/
#            eems_inputs/{gambia.coord, gambia.outer, gambia.diffs}
#            eems_mcmc/chain{1,2,3}/*
#            eems_mig_surface.pdf
#            eems_div_surface.pdf
#            eems_diagnostics.pdf
#            phase6_05_summary.txt
#            phase6_05_results.rds
#
# Tools  : bioconda eems (binary `runeems_snps`) + rEEMSplots (GitHub).
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir   <- .parse_arg("output_dir",   "results/phase6_05")
mono_vcf     <- .parse_arg("mono_vcf",     "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
meta_arg     <- .parse_arg("meta",         "results/phase1_03/sample_master_metadata.rds")
centroid_arg <- .parse_arg("centroid",     "results/phase1_03/village_centroids.tsv")
n_iter       <- as.integer(.parse_arg("n_iter", "5000000"))
n_burn       <- as.integer(.parse_arg("n_burn", "2000000"))
n_thin       <- as.integer(.parse_arg("n_thin", "1000"))
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(vcfR)
  library(rEEMSplots)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
in_dir   <- file.path(output_dir, "eems_inputs");  dir.create(in_dir, showWarnings = FALSE)
mcmc_dir <- file.path(output_dir, "eems_mcmc");    dir.create(mcmc_dir, showWarnings = FALSE)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }
log_msg("=== Phase 6.5: EEMS migration surface (exploratory) ===")
```

- [ ] **Step 2: Build the diffs matrix + EEMS inputs**

Append:

```r
# ============================================================================
# 1. Pairwise allele-difference matrix from the monoclonal VCF
# ============================================================================
vcf  <- read.vcfR(mono_vcf, verbose = FALSE)
meta <- readRDS(meta_arg)
cent <- read_tsv(centroid_arg, show_col_types = FALSE)

samples <- colnames(vcf@gt)[-1]
log_msg(sprintf("Samples in VCF: %d", length(samples)))

# Haploid encoding: 0 or 1; NA elsewhere.
gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
encode_hap <- function(x) {
  out <- rep(NA_real_, length(x))
  out[x %in% c("0", "0/0", "0|0")] <- 0
  out[x %in% c("1", "1/1", "1|1")] <- 1
  out
}
geno <- matrix(encode_hap(as.vector(gt_chr)),
               nrow = nrow(gt_chr), ncol = ncol(gt_chr),
               dimnames = dimnames(gt_chr))

# Pairwise mean(|g_i - g_j|) over jointly-called sites
n <- ncol(geno)
diffs <- matrix(0, n, n, dimnames = list(samples, samples))
for (i in seq_len(n - 1)) {
  for (j in seq.int(i + 1, n)) {
    both <- !is.na(geno[, i]) & !is.na(geno[, j])
    diffs[i, j] <- diffs[j, i] <- mean(abs(geno[both, i] - geno[both, j]))
  }
  if (i %% 50 == 0) log_msg(sprintf("Diffs row %d/%d", i, n))
}

# ============================================================================
# 2. Sample coordinates from village centroids
# ============================================================================
# Each sample inherits its village centroid; samples sharing a village share
# coordinates (EEMS handles this).
sample_loc <- meta %>%
  filter(SampleID %in% samples) %>%
  select(SampleID, VillageCode) %>%
  left_join(cent %>% select(VillageCode, lon, lat), by = "VillageCode") %>%
  filter(!is.na(lon))

# Align matrix to samples that have coordinates
keep_samp <- sample_loc$SampleID
diffs <- diffs[keep_samp, keep_samp]
log_msg(sprintf("Samples with centroid coords: %d", length(keep_samp)))

# ============================================================================
# 3. Write EEMS inputs
# ============================================================================
write.table(sample_loc %>% select(lon, lat),
            file.path(in_dir, "gambia.coord"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(diffs,
            file.path(in_dir, "gambia.diffs"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Outer polygon: convex hull of village centroids, expanded by 0.2 deg
hull_idx <- chull(sample_loc$lon, sample_loc$lat)
hull <- sample_loc[hull_idx, c("lon", "lat")]
centroid_xy <- colMeans(hull)
hull_exp <- hull %>%
  mutate(lon = lon + 0.2 * sign(lon - centroid_xy["lon"]),
         lat = lat + 0.2 * sign(lat - centroid_xy["lat"]))
# Close the polygon
hull_exp <- rbind(hull_exp, hull_exp[1, ])
write.table(hull_exp,
            file.path(in_dir, "gambia.outer"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
log_msg("EEMS inputs written.")
```

- [ ] **Step 3: Run three MCMC chains and post-process**

Append:

```r
# ============================================================================
# 4. Run three MCMC chains (different seeds)
# ============================================================================
eems_bin <- Sys.getenv("EEMS", "runeems_snps")
n_demes  <- 200L

run_chain <- function(chain_id, seed) {
  ch_dir <- file.path(mcmc_dir, sprintf("chain%d", chain_id))
  dir.create(ch_dir, showWarnings = FALSE)
  params <- c(
    sprintf("datapath = %s/gambia", in_dir),
    sprintf("mcmcpath = %s", ch_dir),
    sprintf("nIndiv = %d", length(keep_samp)),
    sprintf("nSites = %d", nrow(geno)),
    sprintf("nDemes = %d", n_demes),
    "diploid = false",
    sprintf("numMCMCIter = %d", n_iter),
    sprintf("numBurnIter = %d", n_burn),
    sprintf("numThinIter = %d", n_thin),
    sprintf("seed = %d", seed))
  param_file <- file.path(ch_dir, "params.ini")
  writeLines(params, param_file)
  log_msg(sprintf("Launching EEMS chain %d (seed %d)", chain_id, seed))
  system2(eems_bin, c("--params", param_file),
          stdout = file.path(ch_dir, "stdout.log"),
          stderr = file.path(ch_dir, "stderr.log"))
}
for (i in 1:3) run_chain(i, 100 * i + 42)

# ============================================================================
# 5. Post-process with rEEMSplots
# ============================================================================
chain_dirs <- file.path(mcmc_dir, sprintf("chain%d", 1:3))

eems.plots(mcmcpath = chain_dirs,
           plotpath = file.path(output_dir, "eems_mig_surface"),
           longlat  = TRUE,
           plot.height = 6, plot.width = 8,
           res = 600,
           add.grid = TRUE, add.outline = TRUE, add.demes = TRUE)
# rEEMSplots writes both mig + diversity panels with that prefix.

# Diagnostics: trace + chain mixing
eems.posterior.draws(mcmcpath = chain_dirs,
                     plotpath = file.path(output_dir, "eems_diagnostics"))

# ============================================================================
# 6. Persist
# ============================================================================
saveRDS(list(diffs = diffs, sample_loc = sample_loc, hull = hull_exp,
             n_iter = n_iter, n_burn = n_burn, n_thin = n_thin),
        file.path(output_dir, "phase6_05_results.rds"))
writeLines(log_lines, file.path(output_dir, "phase6_05_summary.txt"))
message("Done. Outputs in ", output_dir)
```

- [ ] **Step 4: Parse-check**

Run: `Rscript --vanilla -e 'parse("codes/phase6_05_eems.R")'`
Expected: clean parse.

---

## Task 9: `workflow/rules/phase6.smk`

**Files:**
- Create: `workflow/rules/phase6.smk`

- [ ] **Step 1: Write the rule file**

```python
# =============================================================================
# Phase 6 — Reviewer-driven extensions
# =============================================================================
# All five rules are single-pass (no wildcards). All work on the Gambia-only
# monoclonal subset; no {subset} wildcard, no cross-border duplication.
# Conda envs only where non-R tools are involved.
# =============================================================================

rule phase6_haplotype_coancestry:
    input:
        script   = "codes/phase6_01_haplotype_coancestry.R",
        mono_vcf = rules.phase1_coi.output.mono,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        chunkcounts  = "results/phase6_01/chunkcounts.tsv",
        coancestry   = "results/phase6_01/coancestry_matrix.rds",
        heatmap      = "results/phase6_01/chunkcounts_heatmap.pdf",
        tree         = "results/phase6_01/fs_tree.newick",
        rds          = "results/phase6_01/phase6_01_results.rds",
        summary      = "results/phase6_01/phase6_01_summary.txt",
    params:
        out_dir = "results/phase6_01",
    conda:
        "../../envs/chromopainter.yaml"
    log:
        "logs/phase6_haplotype_coancestry.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} > {log} 2>&1"


rule phase6_ibd_tracts:
    input:
        script   = "codes/phase6_02_ibd_tracts.R",
        segments = rules.phase2_hmmibd.output.hmm_segments,
        fract    = rules.phase2_hmmibd.output.hmm_fract,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        tract_table = "results/phase6_02/tract_table.tsv",
        summary_rp  = "results/phase6_02/tract_summary_by_regionpair.tsv",
        wilcoxon    = "results/phase6_02/wilcoxon_within_vs_cross.tsv",
        density_pdf = "results/phase6_02/tract_density_facet.pdf",
        scatter_pdf = "results/phase6_02/tract_length_vs_total_ibd.pdf",
        rds         = "results/phase6_02/phase6_02_results.rds",
        summary     = "results/phase6_02/phase6_02_summary.txt",
    params:
        out_dir = "results/phase6_02",
    log:
        "logs/phase6_ibd_tracts.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--segments={input.segments} --fract={input.fract} "
        "--meta={input.meta} > {log} 2>&1"


rule phase6_network_topology:
    input:
        script     = "codes/phase6_03_network_topology.R",
        pair_class = rules.phase2_hmmibd.output.pair_class,
        trans      = rules.phase3_directional_ibd.output.trans_pairs,
        meta       = rules.phase1_annotate.output.meta_rds,
    output:
        comm_25     = "results/phase6_03/communities_0.25.tsv",
        comm_50     = "results/phase6_03/communities_0.50.tsv",
        betweenness = "results/phase6_03/betweenness.tsv",
        bridges     = "results/phase6_03/bridges.tsv",
        net_25_pdf  = "results/phase6_03/network_0.25.pdf",
        net_50_pdf  = "results/phase6_03/network_0.50.pdf",
        directed    = "results/phase6_03/directed_indeg_outdeg.tsv",
        rds         = "results/phase6_03/phase6_03_results.rds",
        summary     = "results/phase6_03/phase6_03_summary.txt",
    params:
        out_dir = "results/phase6_03",
    log:
        "logs/phase6_network_topology.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--pair_class={input.pair_class} --trans_pairs={input.trans} "
        "--meta={input.meta} > {log} 2>&1"


rule phase6_symmetric_null:
    input:
        script = "codes/phase6_04_symmetric_migration_null.R",
        trans  = rules.phase3_directional_ibd.output.trans_pairs,
        meta   = rules.phase1_annotate.output.meta_rds,
    output:
        null_dist    = "results/phase6_04/null_distribution.tsv",
        obs_vs_null  = "results/phase6_04/observed_vs_null.tsv",
        obs_pdf      = "results/phase6_04/observed_vs_null.pdf",
        rds          = "results/phase6_04/phase6_04_results.rds",
        summary      = "results/phase6_04/phase6_04_summary.txt",
    params:
        out_dir = "results/phase6_04",
    log:
        "logs/phase6_symmetric_null.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--trans_pairs={input.trans} --meta={input.meta} > {log} 2>&1"


rule phase6_eems:
    input:
        script    = "codes/phase6_05_eems.R",
        mono_vcf  = rules.phase1_coi.output.mono,
        meta      = rules.phase1_annotate.output.meta_rds,
        centroids = rules.phase1_annotate.output.centroids,
    output:
        coord       = "results/phase6_05/eems_inputs/gambia.coord",
        outer       = "results/phase6_05/eems_inputs/gambia.outer",
        diffs       = "results/phase6_05/eems_inputs/gambia.diffs",
        mig_pdf     = "results/phase6_05/eems_mig_surface-mrates01.pdf",
        diag_pdf    = "results/phase6_05/eems_diagnostics-pilogl01.pdf",
        rds         = "results/phase6_05/phase6_05_results.rds",
        summary     = "results/phase6_05/phase6_05_summary.txt",
    params:
        out_dir = "results/phase6_05",
    conda:
        "../../envs/eems.yaml"
    log:
        "logs/phase6_eems.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} "
        "--centroid={input.centroids} > {log} 2>&1"
```

- [ ] **Step 2: Syntax sanity-check**

Run: `python -c "import ast, pathlib; ast.parse(pathlib.Path('workflow/rules/phase6.smk').read_text())"`
Expected: no output. (Snakemake rule files are Python-syntactic; this catches indentation / bracket errors before snakemake itself runs.)

---

## Task 10: Snakefile integration

**Files:**
- Modify: `Snakefile`

- [ ] **Step 1: Add the `include:` line**

Change line 6 area of `Snakefile` from:

```python
include: "workflow/rules/phase5.smk"
include: "workflow/rules/manuscript.smk"
```

to:

```python
include: "workflow/rules/phase5.smk"
include: "workflow/rules/phase6.smk"
include: "workflow/rules/manuscript.smk"
```

- [ ] **Step 2: Extend `rule all`**

After the existing Phase 5 entries:

```python
        # Phase 5
        rules.phase5_attribution.output.summary,
        rules.phase5_sensitivity.output.summary,
```

append:

```python
        # Phase 6 (reviewer extensions; Gambia-only, no wildcards)
        rules.phase6_haplotype_coancestry.output.summary,
        rules.phase6_ibd_tracts.output.summary,
        rules.phase6_network_topology.output.summary,
        rules.phase6_symmetric_null.output.summary,
        rules.phase6_eems.output.summary,
```

- [ ] **Step 3: Add the `rule phase6` convenience target**

After the existing `rule phase5:` block, append:

```python
rule phase6:
    input:
        rules.phase6_haplotype_coancestry.output.summary,
        rules.phase6_ibd_tracts.output.summary,
        rules.phase6_network_topology.output.summary,
        rules.phase6_symmetric_null.output.summary,
        rules.phase6_eems.output.summary,
```

- [ ] **Step 4: Add `clean_phase6`**

After the existing `clean_phase5` rule, append:

```python
rule clean_phase6:
    shell:
        "rm -rf results/phase6_*"
```

---

## Task 11: Verification

This task replaces the TDD pass/fail pattern. No test suite exists; we verify the DAG and the lints.

- [ ] **Step 1: Lint the Snakefile**

Run: `snakemake --lint 2>&1 | tail -30`
Expected: no errors related to the new Phase 6 rules. (Pre-existing lint warnings unrelated to Phase 6 are fine.)

- [ ] **Step 2: DAG dry-run for the new convenience target**

Run: `snakemake -n phase6 2>&1 | tail -40`
Expected: five jobs listed (`phase6_haplotype_coancestry`, `phase6_ibd_tracts`, `phase6_network_topology`, `phase6_symmetric_null`, `phase6_eems`), each with their declared inputs/outputs. No "MissingInputException" or "AmbiguousRuleException".

- [ ] **Step 3: DAG dry-run for the full pipeline**

Run: `snakemake -n all 2>&1 | tail -50`
Expected: Phase 6 jobs appear in the listing alongside Phases 1–5. If Phase 6 jobs say "Nothing to be done" that's fine (means existing outputs already satisfy the rules); we want no errors, not necessarily fresh runs.

- [ ] **Step 4: Verify `hmm_segments` is reachable as a target**

Run: `snakemake -n results/phase2_03/hmmibd.hmm.txt 2>&1 | tail -10`
Expected: snakemake either reports "Nothing to be done" (if the file already exists from previous runs) or schedules `phase2_hmmibd`. No "MissingInputException".

- [ ] **Step 5: List all rules and confirm Phase 6 rules are visible**

Run: `snakemake --list 2>&1 | grep -E '^(phase6|clean_phase6)'`
Expected output:

```
phase6
phase6_eems
phase6_haplotype_coancestry
phase6_ibd_tracts
phase6_network_topology
phase6_symmetric_null
clean_phase6
```

- [ ] **Step 6: Stop here. Do not run real jobs.**

Per the user's standing "skip env builds" instruction (Phase 5 smoke-test choice), do not invoke `snakemake --use-conda --conda-create-envs-only` and do not force-rerun real Phase 6 jobs. Real execution waits for the user to schedule it.

---

## Summary of what's done after Task 11

- Three-level importation framework and source-sink framing documented in `CLAUDE.md`.
- Phase 6 addendum and deferred-methods note added to `notes/objective3_proposal.md`.
- Two conda env specs (`chromopainter.yaml`, `eems.yaml`) on disk, not yet built.
- hmmIBD `.hmm.txt` per-segment file now a declared Snakemake output.
- Five Phase 6 R scripts on disk, each parse-checked.
- `workflow/rules/phase6.smk` with five rules.
- `Snakefile` includes Phase 6 in `rule all`, exposes `rule phase6`, exposes `rule clean_phase6`.
- DAG dry-run + lint + rule list all clean.

Real execution and conda env materialisation happen later, on the user's schedule. The deferred methods (MASCOT / Relate / tsinfer) remain on the post-manuscript reminder list in `memory/project_deferred_methods.md`.
