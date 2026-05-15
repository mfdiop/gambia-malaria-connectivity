# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

PhD analysis objective 3: a population-genetics study answering whether residual malaria in western Gambia is sustained by **local transmission**, **importation from eastern Gambia**, or **cross-border importation from Senegal** — and whether that balance shifts as The Gambia approaches elimination. The full scientific framing, five-phase analysis plan, expected figures and tool stack live in `notes/objective3_proposal.md`; treat that document as the source of truth for scope and methodology.

This is an analysis project, not a software package: there is no build system, no test suite, no linter, no CI. Work happens inside the RStudio project `gambia-malaria-connectivity.Rproj` by sourcing R scripts interactively. Scripts under `codes/` are workbooks-in-progress, not modules — they contain commented-out experiments, scratch blocks, and dead code paths alongside live code.

## Repository layout

- `codes/` — R analysis scripts. `clustering_approaches.R` is the main scratchpad (IBS matrix construction, optimal-k selection, k-means/PAM/hierarchical/spectral clustering, validation metrics). `clustering_approaches_v1.R` is an earlier variant using `SNPRelate` and adding Louvain network clustering.
- `data/` — VCFs and metadata. `gambia.vcf` and `senegam_pf.vcf` are large (300–400 MB) and live only on disk, not in version control. `gambia_2014isolates/gambia.ann.vcf.gz` is the annotated 2014 working set. Metadata in `data/metadata/` as `.xlsx`.
- `notes/objective3_proposal.md` — full proposal with research objectives, phased analysis plan, tool table, expected figures. Read this first to understand any task framed around "Phase 2 step 2.1", "the FST analysis", "the attribution model", etc.
- `notes/directional_connectivity_Amambua-Ngwamodel.md` — self-contained implementation spec for the Amambua-Ngwa et al. 2019 directional-transmission model adapted for WGS. Lists the ten functions to implement (`compute_pairwise_distance`, `simulate_lifecycle_times`, … `compute_directional_asymmetry`). Follow this spec when implementing Phase 3.2 / Objective 4.
- `results/` — output directory (currently empty); analysis-script `ggsave()` calls should write here.

## Working with the analysis scripts

**Hardcoded paths reach outside this repo.** Scripts read VCFs from `../01_objective1/01_data/raw_data/` and metadata from `../01_objective1/01_data/meta/`, and write results to `../benchmarking_transmission_methods/results/`. These paths point to sibling project directories under `PhD/Analysis/`. Before running anything, either (a) confirm those sibling paths exist on the current machine, or (b) repoint them to the local `data/` and `results/` folders. Don't silently "fix" a path without checking — the sibling files may be the intended canonical inputs.

**The scripts are not idempotent and not parameterised.** Variables like `k` (cluster count) are assigned partway through the file and reused, sections assume earlier sections have been run in the same R session, and some chunks (`res.km <- eclust(wsaf_clean, "kmeans", k = k, …)`) will fail if you source the file from a fresh session without first running the k-selection section above. When editing, prefer keeping the existing flow over restructuring unless asked.

**IBS matrix construction is the bottleneck.** `get_ibs_matrix()` in `clustering_approaches.R` is an O(n²) nested loop with a `Sys.sleep(0.1)` in the progress-bar tick — for n ≈ 300 samples that alone adds ~30 s of artificial delay. Don't add more sleeps; if you touch this function, the `Sys.sleep` is a candidate for removal but check with the user first since it may be intentional pacing.

## Tooling expectations

The proposal table (`notes/objective3_proposal.md` §7) lists the canonical tools for each analysis. Honour these when extending the codebase:

- VCF processing: PLINK2, `vcfR`, `tidyverse`
- PCA / FST: PLINK2, `hierfstat`, `pegas`, `SNPRelate`
- Ancestry: ADMIXTURE v1.3, `LEA`, fineSTRUCTURE
- IBD: hmmIBD, dcifer, isoRelate
- Population graph: TreeMix v1.13
- Isolation-by-distance: `ade4` Mantel, `vegan`

If a task asks for "IBD relatedness" use hmmIBD (monoclonal) or dcifer (polyclonal) per the proposal — do not substitute `SNPRelate::snpgdsIBDMLE` even though it appears in the v1 script. The v1 approach is exploratory; the proposal commits to the hmmIBD/dcifer stack.

## Conventions worth respecting

- Heterozygous calls in WGS are handled via within-sample allele frequency (WSAF = altAD/DP) computed from `vcfR::extract.gt`, with optional rounding or NA-masking — see `get_ibs_matrix()` for the canonical pattern.
- Polyclonal samples: for distance/relatedness use the minor allele at mixed positions (per the Amambua-Ngwa spec, `notes/directional_connectivity_Amambua-Ngwamodel.md` Step 1).
- Sample-to-metadata joins are by `SampleID` (against `GamMetadata_2014.xlsx` / `GamMetadata_Final_imputemissingdate.xlsx`). Region/village/year columns are: `VillageCode`, `Year`, `Village`.
- Plots are saved with `ggsave(..., width = 12, height = 10, dpi = 600)` as PDFs. Match this style for new figures.
