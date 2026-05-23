# Gambia malaria connectivity — Importation versus local transmission

PhD analysis (Objective 3) resolving the source of residual *Plasmodium falciparum* transmission in western Gambia, using whole-genome parasite connectivity. We ask whether residual malaria in low-prevalence western regions is sustained by **local transmission**, **importation from the higher-transmission eastern Gambia**, or **cross-border importation from Senegal**, and whether that balance is shifting as The Gambia approaches elimination.

This work extends Amambua-Ngwa et al. (2019) — which used 54-SNP barcodes to detect 60 directional transmission events — into WGS resolution, removing four constraints simultaneously: barcode resolution, p-distance recombination artefacts, polyclonal infections, and absence of Senegalese comparators.

## Project structure

```
codes/                                   # Analysis scripts (phase-based + exploratory variants)
  phase1_01_wgs_qc.R
  phase1_02_coi_stratification.R
  phase1_03_sample_annotation.R
  phase2_01_fst_ibd_mantel.R
  phase2_02_pca_admixture.R
  phase2_02_pca_lea.R
  phase2_02_pca_lea_extendK.R
  phase2_03_hmmibd.R
  phase3_01_supervised_admixture.R
  phase3_02_directional_ibd.R
  phase3_03_isolation.R
  clustering_approaches.R
  clustering_approaches_v1.R
data/                                    # Primary genomics inputs and metadata
  gambia.vcf
  senegam_pf.vcf
  gambia_2014isolates/
  metadata/
notes/                                   # Study plan and modelling notes
  objective3_proposal.md
  objective3_proposal.docx
  directional_connectivity_Amambua-Ngwamodel.md
results/                                 # Script outputs organized by phase
  phase1_01/
  phase1_02/
  phase1_03/
  phase2_01/
  phase2_02_lea/
  phase2_03/
  phase3_01/
  phase3_02/
  phase3_03/
  phase4_01/
tools/
  hmmIBD/                                # Local helper assets/utilities for hmmIBD workflows
gambia-malaria-connectivity.Rproj        # RStudio project entry point
CLAUDE.md                                # Local coding-agent guidance
SampleMetadata.txt                       # Auxiliary metadata text file
```

## Data

| Dataset | Source | Use |
|---|---|---|
| Gambian WGS 2014 | This cohort | Primary set for all analyses |
| 54-SNP barcodes (Amambua-Ngwa 2019) | Published | Temporal validation, barcode-vs-WGS comparison |
| Senegalese WGS | MalariaGEN Pf7 | Cross-border connectivity |

Reference: *P. falciparum* 3D7 (`Pf3D7_*_v3` contigs). VCF is already annotated with `RegionType` (Core, Centromere, Subtelomeric*, InternalHypervariable), so hypervariable antigen loci (*var*, *rifin*, *stevor*) are filtered by retaining `RegionType=Core`.

## Analysis pipeline

Current implemented phases map to numbered scripts in `codes/`:

| Phase | Step | Script |
|---|---|---|
| 1. Data prep | 1.1 WGS QC | `phase1_01_wgs_qc.R` |
|  | 1.2 COI stratification | `phase1_02_coi_stratification.R` |
|  | 1.3 Regional/temporal annotation | `phase1_03_sample_annotation.R` |
| 2. Pop-genetic background | 2.1 FST + isolation-by-distance | `phase2_01_fst_ibd_mantel.R` |
|  | 2.2 PCA + ADMIXTURE | `phase2_02_pca_admixture.R` |
|  | 2.3 Baseline IBD landscape | `phase2_03_hmmibd.R` |
| 3. Importation | 3.1 Supervised admixture | `phase3_01_supervised_admixture.R` |
|  | 3.2 Directional IBD/transmission chains | `phase3_02_directional_ibd.R` |
|  | 3.3 Genomic isolation tests | `phase3_03_isolation.R` |

Phases 4–5 in `notes/objective3_proposal.md` are planned but do not yet have corresponding scripts in `codes/`.

Each script is self-contained and re-runnable: reads VCF + metadata, writes outputs to `results/`, and persists intermediate objects (RDS) for downstream phases.

## Reproducibility

- **R version:** 4.5.x (tested with R 4.5.3)
- **Working directory:** open `gambia-malaria-connectivity.Rproj` in RStudio; all paths are relative to the project root
- **External tools:** `bcftools`, `vcftools`, `plink`, `hmmIBD`, `ADMIXTURE`, `TreeMix` — install separately and ensure they are on `PATH`
- **Run order:** scripts must run in numerical order. Phase *n* depends on RDS outputs of phase *n−1*.

## Status

In progress. See `notes/objective3_proposal.md` for the full plan and `CLAUDE.md` for development conventions.

## Citation

Diop MF, Amambua-Ngwa A. Importation versus local transmission: resolving the source of residual malaria in western Gambia using whole-genome parasite connectivity. *In preparation*.

## References

1. Amambua-Ngwa A, et al. Major subpopulations of *Plasmodium falciparum* in sub-Saharan Africa. *Science* (2019).
2. Okafor UA, et al. Projection of future malaria prevalence in the upper river region of The Gambia. *Malar J* **24**, 108 (2025).
