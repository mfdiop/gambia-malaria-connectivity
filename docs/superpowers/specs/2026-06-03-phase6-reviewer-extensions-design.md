# Phase 6 — Reviewer-driven extensions: design spec

**Author:** Mouhamadou Fadel Diop
**Date:** 2026-06-03
**Status:** Design approved, plan pending

---

## 1. Motivation

A supervisor-style review of the Objective 3 proposal (recorded in `notes/proposal_reviewed.md` and `notes/objective3_proposal_reviewed.md`) raised six critical issues. The first five are already addressed implicitly by Phases 1–5 of our Snakemake pipeline. Four are not, and one is a documentation-only fix:

1. **"Importation" is not formally defined.** Documentation gap only — fixed by adding a three-level definition table to `CLAUDE.md` and the proposal.
2. **Source–sink framing (R_local < 1 vs > 1) is buried.** Documentation gap only — fixed by elevating it to the central question of the chapter.
3. **ADMIXTURE / sNMF is overinterpreted.** Code gap — we need a haplotype-based alternative (ChromoPainter + fineSTRUCTURE) and to demote sNMF to a supportive role.
4. **No IBD tract-length analysis.** Code gap — tract length is the recency signal; we currently only use total IBD.
5. **No IBD network topology.** Code gap — we have pair-level IBD but no community / hub / bridge analysis.
6. **No explicit symmetric-migration null** for the Phase 3.2 directional asymmetry. Code gap — strengthens the asymmetry claim against the obvious reviewer attack.

Separately, the user wants one of the reviewer's "out-of-scope" advanced spatial-genomic methods explored. After tradeoff discussion (recorded in `memory/project_deferred_methods.md`), **EEMS** was chosen because it extends the Phase 2 isolation-by-distance story cleanly and is the only one that finishes in PhD-realistic time without phasing pain.

The remaining three (MASCOT structured coalescent, Relate, tsinfer/tsdate) are deferred to post-manuscript and will be revisited once the Obj3 manuscript draft is complete.

## 2. Goal

Add a new **Phase 6** to the pipeline that delivers five analyses, each addressing one of the gaps above, plus the two documentation updates. Maintain the existing Phase 1–5 invariants:

- Same script architecture (R workbooks under `codes/`, parameterised via `--key=value` args parsed by `.parse_arg`, runnable interactively *or* by Snakemake).
- Same output convention (`results/phase6_NN/`, one `summary.txt` and one `.rds` per rule).
- Same Snakemake conventions (rule per analysis, `--vanilla` Rscript, log to `logs/`, conda only where needed).

## 3. Conceptual additions (documentation)

### 3.1 Three-level definition of importation

Add this table to `CLAUDE.md` (and the proposal addendum) to head off reviewer confusion about ADMIXTURE/PCA being misread as evidence of recent transmission importation:

| Level | What it measures | Timescale | Method in this pipeline |
|---|---|---|---|
| Historical connectivity | Long-term shared ancestry / drift | Years–decades | Phase 2.1 FST + Mantel; Phase 2.2 PCA + sNMF |
| Recent parasite migration | Recently transmitted shared ancestry | Months–years | Phase 2.3 hmmIBD pairwise IBD; **Phase 6.1 haplotype co-ancestry**; **Phase 6.2 IBD tract lengths** |
| Active transmission importation | Source-to-recipient event within the parasite generation time | Weeks–months | Phase 3.2 P_trans; Phase 4.1 cross-border decision-era subset; Phase 5.1 Bayesian attribution |

### 3.2 Source–sink framing

Add a paragraph to the proposal making explicit that the chapter operationalises:

- **Self-sustaining local transmission:** R_local > 1 in western Gambia, with persistence maintained internally.
- **Importation-maintained sink:** R_local < 1 in western Gambia, with persistence sustained by parasite influx from eastern Gambia and Senegal.

Phase 5.1 Bayesian attribution is the inferential bridge: P(H_local), P(H_east), P(H_senegal) is the source-sink decomposition. No new code; this is a framing edit.

## 4. Phase 6 architecture

```
codes/
  phase6_01_haplotype_coancestry.R
  phase6_02_ibd_tracts.R
  phase6_03_network_topology.R
  phase6_04_symmetric_migration_null.R
  phase6_05_eems.R

workflow/rules/phase6.smk

envs/
  chromopainter.yaml   # ChromoPainter v2 + fineSTRUCTURE
  eems.yaml            # EEMS C++ binary + rEEMSplots

results/
  phase6_01/ ... phase6_05/
```

Snakefile additions:

- `include: "workflow/rules/phase6.smk"`
- Each phase6 rule's summary added to `rule all`
- New convenience target `rule phase6` aggregating all five summaries

All five scripts follow the same arg-parsing pattern used throughout Phases 1–5 (defaults map to results paths so the scripts remain runnable interactively in RStudio).

### 4.1 Sample scope (decided)

All five analyses are **Gambia-only**, keyed off:

- `rules.phase1_coi.output.mono` — the QC'd monoclonal Gambia VCF (input to Phase 2)
- `rules.phase1_annotate.output.meta_rds` — sample metadata with region, year, season, location
- `rules.phase1_annotate.output.centroids` — village centroids (needed only by EEMS)

Cross-border analyses already live in Phase 4; Phase 6 does not duplicate that.

### 4.2 Cross-rule extension

Phase 6.2 requires hmmIBD's per-segment output (`.hmm.txt`), which Phase 2.3 currently runs but does not save. The Phase 6 work includes a small surgical edit to Phase 2.3:

- Add `hmm_segments = "results/phase2_03/hmmibd.hmm.txt"` to `rule phase2_hmmibd.output`
- Add a corresponding `cp` step in `codes/phase2_03_hmmibd.R` to persist that file alongside the existing `hmm_fract.txt`

No other Phase 1–5 code is touched.

## 5. Per-analysis specifications

### 5.1 Phase 6.1 — Haplotype co-ancestry (ChromoPainter + fineSTRUCTURE)

**Scientific role.** Direct, haplotype-resolution alternative to sNMF/ADMIXTURE. Reviewer's central objection is that ADMIXTURE clusters reflect long-term ancestry and are over-interpreted as recent migration evidence; haplotype-sharing tools resolve recent co-ancestry without that confound.

**Inputs.**
- `phase1_coi.mono` (monoclonal-only Gambia VCF)
- `phase1_annotate.meta_rds` for region/year/location annotations

**Method.**
1. Restrict to samples with COI = 1 (true monoclonal). Pf monoclonal samples are effectively haploid; this avoids phasing pain.
2. Convert VCF → ChromoPainter `.phase` and `.recombfile` formats. Use ChromoPainter haploid mode (`-j` flag, ploidy = 1) — no SHAPEIT required.
3. Estimate switch rate (Ne) and mutation rate (μ) via 10 EM iterations on one chromosome (typically chr 11, mid-size, well-behaved).
4. Run full painting genome-wide with the EM-estimated Ne/μ.
5. fineSTRUCTURE MCMC on the chunkcounts matrix: 1M burn-in + 1M post-burn, thin 1k. Tree inference with 100k MCMC tree-building steps.
6. Extract the co-ancestry matrix (chunkcounts) and the fineSTRUCTURE tree. Order the heatmap by fineSTRUCTURE tree leaves.

**Outputs.**
- `coancestry_matrix.tsv` — N × N chunkcount matrix
- `chunkcounts_heatmap.pdf` — ordered heatmap with region annotation bars
- `fs_tree.newick` — fineSTRUCTURE tree
- `phase6_01_summary.txt` — sample count, EM-estimated Ne/μ, top-level fineSTRUCTURE clusters by region composition
- `phase6_01_results.rds` — persisted state

**Manuscript figure target.** Main figure or supplementary, replacing the central claim from sNMF ancestry. Caption: "Haplotype co-ancestry of Gambian monoclonal P. falciparum infections, ordered by fineSTRUCTURE clustering."

**Conda env.** `envs/chromopainter.yaml` — pulls `chromopainter` and `finestructure` from bioconda.

### 5.2 Phase 6.2 — IBD tract length structure

**Scientific role.** Total IBD says "how much of the genome is shared." Tract lengths say "how recently was it shared." Long cross-regional tracts = recent east↔west migration; short tracts = ancient shared ancestry. This separates the recent-migration signal from the historical-connectivity signal within a single dataset.

**Inputs.**
- `phase2_hmmibd.hmm_segments` (newly exposed — see §4.2)
- `phase1_annotate.meta_rds`

**Method.**
1. Parse `hmmibd.hmm.txt` into a long table: (sample1, sample2, chr, start_bp, end_bp, length_bp, ibd_state).
2. Filter to `ibd_state = 1` (true IBD segments, not non-IBD).
3. Annotate each pair with region pair (west–west, west–central, west–east, central–central, central–east, east–east) from `meta_rds`.
4. Compute summary stats per region pair: median tract length, mean tract length, p90 tract length, n_segments, n_pairs contributing.
5. Two-sided Wilcoxon test: within-region tract lengths vs cross-region tract lengths (pooled). Then specifically: west–west tracts vs west–east tracts.
6. Plot 1: density of tract lengths, faceted by region pair, log-x scale, vertical line at median.
7. Plot 2: scatter of mean tract length per pair vs total fract_sites_IBD per pair, coloured by region pair.

**Outputs.**
- `tract_table.tsv` — long-form per-segment table
- `tract_summary_by_regionpair.tsv` — summary stats by region pair
- `wilcoxon_within_vs_cross.tsv` — test results
- `tract_density_facet.pdf`
- `tract_length_vs_total_ibd.pdf`
- `phase6_02_summary.txt`
- `phase6_02_results.rds`

**Manuscript figure target.** Supplementary or part of the IBD figure. Caption: "IBD tract lengths reveal recent (long) vs ancient (short) connectivity within Gambia."

**Conda env.** None — pure R.

### 5.3 Phase 6.3 — IBD network topology

**Scientific role.** Pair-level IBD becomes a network when you cross a threshold. Network topology answers: are there transmission hubs? Bridge samples connecting otherwise separate communities? Distinct community structure that aligns with geography? Reviewer's recommendation pulled from population-genetics standard practice — already done in scratchpad `clustering_approaches_v1.R` but not in the live pipeline.

**Inputs.**
- `phase2_hmmibd.pair_class` (pair_ibd_classified.tsv with fract_sites_IBD)
- `phase3_directional_ibd.trans_pairs` (transmission_pairs.tsv with P_trans) — optional input for the directed bonus
- `phase1_annotate.meta_rds`

**Method.**
1. Build two undirected weighted graphs from pair_ibd_classified.tsv:
   - **lenient network:** edges with fract_sites_IBD ≥ 0.25
   - **strict network:** edges with fract_sites_IBD ≥ 0.50
2. For each graph: Louvain community detection (`igraph::cluster_louvain`).
3. For each node: betweenness centrality (`igraph::betweenness`).
4. Bridge identification: nodes in the top decile of betweenness whose edges span ≥ 2 regions.
5. Force-directed layout plot (`ggraph` Fruchterman–Reingold) for each network, coloured by community, node size scaled by betweenness, region annotation legend.
6. **Bonus directed analysis:** weighted directed graph from transmission_pairs.tsv (source_id → recipient_id, weight = P_trans). Per-node in-degree (sources reaching this sample) and out-degree (samples this one putatively seeded). Report top-10 sources, top-10 sinks.

**Outputs.**
- `communities_0.25.tsv`, `communities_0.50.tsv`
- `betweenness.tsv`
- `bridges.tsv`
- `network_0.25.pdf`, `network_0.50.pdf`
- `directed_indeg_outdeg.tsv` (bonus)
- `phase6_03_summary.txt`
- `phase6_03_results.rds`

**Manuscript figure target.** Supplementary figure. Caption: "Network topology of Gambian P. falciparum IBD reveals communities and bridge samples."

**Conda env.** None — pure R (`igraph`, `ggraph`).

### 5.4 Phase 6.4 — Symmetric-migration permutation null

**Scientific role.** Phase 3.2 currently reports a directional asymmetry index but without an explicit null. The obvious reviewer attack: "Maybe the asymmetry comes from sampling imbalance, not biology." A label-permutation test exposes how big the asymmetry would look under a symmetric-migration null.

**Inputs.**
- `phase3_directional_ibd.trans_pairs` (transmission pairs with directional labels)
- `phase3_directional_ibd.asymmetry` (the observed asymmetry index table)
- `phase1_annotate.meta_rds`

**Method.**
1. Observed statistic: A_obs = (n_east→west − n_west→east) / (n_east→west + n_west→east), restricted to high-P_trans pairs (≥ 0.5). Same threshold used in Phase 3.2.
2. Null: under the symmetric-migration hypothesis, region labels are exchangeable across samples **within year + transmission-season strata** (preserving temporal structure so the null doesn't artificially break seasonality). Permute `region` within strata 1000 times. For each permutation:
   - Re-derive each transmission pair's source_region and recipient_region from the permuted labels.
   - Recompute A_perm.
3. Empirical two-sided p-value: 2 × min(P(A_perm ≥ A_obs), P(A_perm ≤ A_obs)).
4. Plot: histogram of A_perm under the null, vertical line at A_obs, p-value annotated.

**Outputs.**
- `null_distribution.tsv` — A_perm for each permutation
- `observed_vs_null.tsv` — A_obs, p_value, n_perm
- `observed_vs_null.pdf`
- `phase6_04_summary.txt`
- `phase6_04_results.rds`

**Manuscript figure target.** Inset to the Phase 3.2 directional figure, or supplementary. Caption: "Observed east→west directional asymmetry exceeds the symmetric-migration null (p_perm)."

**Conda env.** None — pure R.

### 5.5 Phase 6.5 — EEMS (exploratory)

**Scientific role.** Spatial migration surface across Gambia. Where are the migration corridors and barriers between west, central, east? Extends the Phase 2.1 IBD-by-distance story from "is there an IBD-distance correlation" to "where on the map does the effective migration rate diverge from the isolation-by-distance expectation." Marked as exploratory: included in the supplement and the thesis, **not** the primary manuscript figure set.

**Inputs.**
- `phase1_coi.mono` (monoclonal Gambia VCF — sample list and genotypes for the diffs matrix)
- `phase1_annotate.meta_rds` (sample → village mapping)
- `phase1_annotate.centroids` (village → lat/lon)

**Method.**
1. Pairwise allele-difference matrix: compute fresh from the VCF using `1 − mean(|g_i − g_j|)` over called sites, where g is the haploid genotype call (0 or 1). EEMS wants raw allelic diffs, not IBD; do not reuse `phase2_03/ibd_matrix.rds`.
2. Assign each sample to its village centroid coordinates. (EEMS allows samples to share coordinates.)
3. Habitat polygon = convex hull of all village centroids, expanded outward by 0.2°.
4. Write the three EEMS input files: `gambia.coord` (lon lat per sample), `gambia.outer` (polygon vertices), `gambia.diffs` (the N × N diffs matrix).
5. Run EEMS MCMC via the bioconda binary: 3 chains, 5M iterations, 2M burn-in, thin 1k. Different `--seed` per chain.
6. Post-process with `rEEMSplots`: posterior mean log(m) migration surface, posterior mean log(q) diversity surface, MCMC trace diagnostics across chains.

**Outputs.**
- `eems_inputs/{gambia.coord, gambia.outer, gambia.diffs}` — persisted EEMS inputs
- `eems_mcmc/chain{1,2,3}/` — raw MCMC outputs
- `eems_mig_surface.pdf` — migration surface
- `eems_div_surface.pdf` — diversity surface
- `eems_diagnostics.pdf` — MCMC convergence + chain mixing
- `phase6_05_summary.txt` — runtime, posterior summaries
- `phase6_05_results.rds`

**Manuscript figure target.** Supplementary or thesis-only figure, explicitly flagged as exploratory. Caption: "Estimated effective migration surface across The Gambia (EEMS)."

**Conda env.** `envs/eems.yaml` — `bioconda::eems`. R script post-processing uses `rEEMSplots` (CRAN-not-available; installed via `remotes::install_github('dipetkov/eems/plotting/rEEMSplots')` in renv).

## 6. Top-level integration

### 6.1 Snakefile additions

- `include: "workflow/rules/phase6.smk"`
- `rule all` input list extended with each phase6 summary
- `rule phase6` convenience target aggregating all five summaries

### 6.2 New convenience pattern

The clean rules already include `clean_phase5`; add `clean_phase6` for symmetry.

### 6.3 Verification

- `snakemake -n phase6` builds DAG with the five rules.
- `snakemake -n all` builds DAG including phase6 outputs.
- `snakemake --lint` produces no new errors.

## 7. Out of scope (explicitly)

The reviewer suggested several additional advanced methods. After tradeoff discussion, three are explicitly **deferred to post-manuscript**, recorded in `memory/project_deferred_methods.md`:

- **MASCOT / BASTA (structured coalescent)** — would directly estimate m_east→west vs m_west→east in a Bayesian framework. Most aligned with the source-sink question. Postponed because compute time is days–weeks and convergence requires careful babysitting.
- **Relate** — genome-wide genealogies with branch dates. Phasing tooling for Pf monoclonal is awkward.
- **tsinfer / tsdate** — modern ARG inference. Cool but only loosely tied to the source-sink question.

The user has asked to be reminded about these once the Obj3 manuscript draft is complete.

Mobility modelling is not in scope — no mobility data is available (per the reviewer's own follow-up).

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| ChromoPainter haploid mode is poorly documented | Test the conversion on a small chr 11 subset first; sanity-check chunkcount column sums against sample counts |
| fineSTRUCTURE MCMC very slow on 800+ samples | Down-sample to monoclonal-only (already in spec) reduces this materially; document runtime in summary |
| hmmIBD `.hmm.txt` format differs across versions | Add a one-line schema check in `phase6_02` (expect columns: sample1, sample2, chr, start, end, different, Nsnp); fail loudly if not |
| EEMS conda binary linked against incompatible libc on macOS | Fallback documented: build from source via `bioconda::cxx-compiler` and the EEMS source tarball, both in the conda env |
| Permutation null in 6.4 misinterprets sampling design | Stratified permutation (within year + transmission-season) is built into §5.4 step 2 |

## 9. Deliverables checklist

- [ ] `codes/phase6_01_haplotype_coancestry.R`
- [ ] `codes/phase6_02_ibd_tracts.R`
- [ ] `codes/phase6_03_network_topology.R`
- [ ] `codes/phase6_04_symmetric_migration_null.R`
- [ ] `codes/phase6_05_eems.R`
- [ ] `workflow/rules/phase6.smk`
- [ ] `envs/chromopainter.yaml`
- [ ] `envs/eems.yaml`
- [ ] Edit `workflow/rules/phase2.smk` to expose `hmm_segments`
- [ ] Edit `codes/phase2_03_hmmibd.R` to persist the segment file
- [ ] Edit `Snakefile` (`rule all`, `rule phase6`, `rule clean_phase6`)
- [ ] Edit `CLAUDE.md` (three-level importation table + R_local framing)
- [ ] Edit `notes/objective3_proposal.md` (Phase 6 addendum + deferred-methods note)
