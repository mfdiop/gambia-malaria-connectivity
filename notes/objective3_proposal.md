---
output:
  word_document: default
  pdf_document: default
---

# Importation versus local transmission: Resolving the source of residual malaria in Western Gambia Using Whole-genome parasite connectivity

**Author:** Mouhamadou Fadel Diop  
**Supervisory team:** Alfred Amambua-Ngwa  
**Date:** 15 May 2026  
**Datasets:** Gambian WGS from 2014 (primary), 54-SNP barcodes, WGS from Senegal, MalariaGEN Pf8 


## 1. Central Question

> **Is the residual malaria transmission persisting in low-prevalence western Gambia sustained by local mosquito-mediated transmission, by parasite importation from the higher-transmission eastern Gambia, or by cross-border importation from Senegal and has this balance shifted as the Gambia approaches elimination?**


## 2. Scientific Background and Rationale

### 2.1 The elimination problem in western Gambia

Malaria transmission in The Gambia is geographically polarised or presents a contrasting picture [1]. Western regions, particularly around the Atlantic coast and the Greater Banjul Area, have achieved prevalence rates below 2% in recent surveys. Eastern regions, clustered around the North and South Banks in the Upper River Division, maintain prevalence rates of 10–40% during peak season [2]. Despite decades of intensive intervention in the west, high LLIN coverage, indoor residual spraying, ACT availability, and seasonal chemoprevention, the western transmission hotspots have proven resistant to elimination.

Two explanations have been proposed but not yet resolved with high-resolution data. The first is that residual local transmission is maintained by a small, hard-to-reach infectious reservoir, asymptomatic carriers, or untreated children that continuously regenerates local incidence independent of importation. The second is that the western Gambia parasites are repeatedly re-seeded by importation from the east or from Senegal, such that local transmission is interrupted but reintroduction prevents elimination. 
The practical consequence of choosing the wrong model is substantial: the first explanation calls for intensified local detection and treatment; the second calls for border-level surveillance and coordinated regional intervention.

### 2.2 What the barcode study showed and its limits

Amambua-Ngwa et al. (2019), the foundational paper this objective builds on, demonstrated with 54-SNP barcode data that the most likely directional flow of parasites at peak transmission was east-to-west, with three cross-regional transmission paths from eastern to western Gambia at the October peak. However, this analysis was constrained in several important ways: only 60 transmission events could be detected under the travel-constrained model; barcodes could not resolve parasite genetic backgrounds at fine scale; COI above 1 rendered many samples uninformative; and no Senegalese comparative data were available to evaluate cross-border importation. We assume that WGS data will remove all four constraints simultaneously.


## 3. What this project try to achieve

**Answer a question.** The importation-versus-local-transmission question has a factual answer that real data can reveal. This is distinct from simulation studies, which generate plausible scenarios but cannot resolve what is actually happening in The Gambia right now.

**Immediate programmatic impact.** The Gambia National Malaria Control Programme and WHO AFRO have both identified resolution of the importation question as a strategic priority for 2030 elimination roadmap. A genomics-based answer, validated with the temporal depth of your WGS archive, would directly inform national strategy.

**Extends rather than replaces existing published work.** Rather than starting from zero, this proposal extends the 2019 Alfred et al. barcode paper into WGS resolution.

**Classical and well-established methods.** IBD-based relatedness analysis, ADMIXTURE ancestry estimation, and isolation-by-distance modelling are population genetics. They require no novel algorithmic development and are immediately interpretable.


## 4. Research Objectives

1. Determine whether highly related parasite pairs in western Gambia are genetically more similar to eastern Gambian parasites or to local western parasites, quantifying the proportion of western Gambia infections attributable to eastern and Senegalese sources.
2. Characterise the seasonal and year-to-year dynamics of east-to-west parasite flow, testing whether importation peaks at high-transmission season and falls in the dry season as predicted by the mosquito lifecycle and human mobility patterns.
3. Evaluate the contribution of cross-border importation from Senegal relative to domestic east-to-west flow, using the first WGS-resolution comparison of Gambian and Senegalese parasite genomes.
4. Apply the life-cycle-constrained directional transmission model from Amambua-Ngwa et al. (2019) to WGS data to identify directional transmission paths with substantially greater resolution and confidence than was possible with 54-SNP barcodes.


## 5. Data

### 5.1 Gambian WGS dataset

- WGS samples collected across seasonal malaria surveys, The Gambia, 2014 (primary), 1966–2022 (secondary)
- Geographic coverage: 3 main regions (west, central, east)
- Temporal coverage: to enable longitudinal comparison

### 5.2 Barcode dataset

- 355 isolates from Amambua-Ngwa et al. (2019): 12 villages, 2013 season, monthly sampling
- Used for: temporal validation, household-level ground truth for method calibration, comparison of barcode vs WGS resolution

### 5.3 Senegalese comparator

- MalariaGEN Pf7 release: Senegalese samples, particularly from border zones (Ziguinchor, Kolda, Tambacounda)
- Used for: cross-border connectivity analysis, ancestry deconvolution of Gambian samples


## 6. Analysis Plan

### Phase 1 — Data preparation

**Step 1.1 — WGS quality control**

Apply uniform quality filters:
- Minimum mean coverage: 5×
- Maximum missing genotype rate per sample: 20%
- Minor allele frequency > 0.01 across the full dataset
- Exclude hypervariable antigen loci (*var*, *rifin*, *stevor*)
- Retain biallelic SNPs only for relatedness analyses

**Step 1.2 — COI estimation and sample stratification**

Use `THE REAL McCOIL` or `moimix` to estimate COI per sample. Stratify samples into:
- Monoclonal (COI = 1): primary set for haplotype-level analyses
- Low polyclonal (COI 2–3): include in IBD analyses using dcifer
- High polyclonal (COI ≥ 4): flag for sensitivity analysis; exclude from directional inference

**Step 1.3 — Regional and temporal annotation**

Assign each sample to: region (west/central/east), village, household (if available), collection month, transmission season, and intervention era (if available). This metadata table is the backbone of all downstream analyses.


### Phase 2 — Population genetic background characterisation

**Step 2.1 — Genetic differentiation between regions**

Compute pairwise FST (Weir–Cockerham) between all village pairs using R or PLINK2. Construct an FST distance matrix and visualise using:
- Hierarchical clustering (to reveal natural regional groupings)
- Geographic overlay (test whether FST increases monotonically with geographic distance)

Test for isolation-by-distance (IBD) using a Mantel test (genetic distance vs geographic distance, 10000 permutations). A significant IBD relationship supports the expectation that gene flow decreases with geographic separation, the baseline model against which importation events are detected as deviations.

**Step 2.2 — Population structure: PCA and ADMIXTURE**

Run PCA using PLINK2 on LD-pruned SNPs (r² < 0.1, 50kb windows). Run ADMIXTURE for K = 2–8 with 10-fold cross-validation to select optimal K.

Key question: Do western and eastern Gambia form distinct genetic clusters, or does the analysis show a continuous gradient? A discrete western cluster would suggest local adaptation or prolonged isolation; a gradient would suggest ongoing gene flow.

**Step 2.3 — Baseline relatedness landscape**

Run hmmIBD on all pairwise sample combinations to produce a genome-wide pairwise IBD matrix. Classify pairs as:
- Highly related (IBD > 0.5): candidate recent transmission pairs
- Moderately related (IBD 0.125–0.5): probable 2nd–4th degree genealogical connections
- Unrelated (IBD < 0.125): background

Visualise the spatial distribution of highly related pairs as a chord diagram and geographic network: are high-IBD pairs spatially clustered within regions, or do cross-regional high-IBD pairs form a directional pattern?


### Phase 3 — Quantifying importation: three complementary analyses

**Step 3.1 — Ancestry deconvolution approach**

Fit an ADMIXTURE model treating the three Gambian regions (west, central, east) as three source populations. For each western Gambia sample, estimate what proportion of its genome traces to the eastern Gambian gene pool versus the local western gene pool. This proportion is the individual-level importation score.

Extension: Include Senegalese samples (Pf7) as a fourth source population. For each Gambian sample, compute: % ancestry from west Gambia, % from east Gambia, % from Senegal. Western Gambia samples with >20% eastern or Senegalese ancestry are flagged as probable imports.

Temporal analysis: Has the proportion of western Gambia samples with substantial eastern ancestry increased or decreased across seasons? This tests whether importation pressure is growing or declining as transmission diverges.

**Step 3.2 — IBD-based transmission chain reconstruction**

Apply the life-cycle-constrained probability model from Amambua-Ngwa et al. (2019) to the WGS dataset. The model constrains transmission paths by:
1. Time window: recipient must be sampled 32–67 days after donor (the *P. falciparum* life cycle constraint)
2. Genetic similarity: pairs must share sufficient genomic similarity to be within the detectable relationship range
3. Travel history: if available, cross-regional paths constrained to individuals with documented inter-regional travel

WGS offers three major improvements over the original barcode implementation:
- **Resolution:** > 10.000 biallelic SNPs versus 54 SNPs. Genetic distance estimates are far more precise.
- **IBD over IBS:** WGS enables hmmIBD estimation rather than p-distance, removing the recombination-confounded artefacts that forced the 30% distance cutoff in the barcode study.
- **Polyclonal handling:** dcifer and isoRelate allow relatedness estimation for polyclonal samples that were uninformative with barcodes.

For each detected transmission path, classify as:
- Within-village
- Within-region (between villages, same region)
- Cross-regional (west ↔ east)
- Cross-border (Gambia ↔ Senegal, if Senegalese samples are temporally aligned)

Compute: number of paths in each category, directional asymmetry (east→west vs west→east), seasonal timing of cross-regional paths.

**Step 3.3 — Genomic isolation analysis: are western Gambia parasites genetically distinct?**

If western Gambia is primarily maintained by local transmission, the local parasite population should show signs of genetic isolation: elevated within-region IBD, lower diversity than the east, more pronounced LD at short distances. Test each of these predictions:

- Compare mean within-region IBD (west vs east, using hmmIBD)
- Compare π per region per season
- Compare LD decay per region

If the west shows patterns consistent with local isolation (high IBD, low π, slow LD decay), this supports the local-reservoir model. If the west shows genomic diversity inconsistent with its low transmission intensity — i.e., more diversity than expected given its prevalence — this supports continuous re-seeding from the east.


### Phase 4 — Cross-border importation from Senegal

**Step 4.1 — WGS-level comparison with Senegalese isolates**

Using MalariaGEN Pf7 Senegalese samples (matched for collection period where possible), compute pairwise IBD between all Gambian and Senegalese sample pairs. Identify cross-border highly related pairs (IBD > 0.25) and characterise:
- Their geographic origin within Senegal (border zone vs interior)
- Their seasonal distribution (peak at October–November, consistent with peak-season movement)
- Their directional assignment using the life-cycle constraint model

**Step 4.2 — TreeMix population graph**

Apply TreeMix (Pickrell and Pritchard 2012) to model the population graph including western Gambia, central Gambia, eastern Gambia, and Senegal nodes. Allow up to 5 migration edges and test whether a migration edge from Senegal to eastern Gambia is consistently detected across bootstrap replicates. This provides a graph-based test of the cross-border importation hypothesis independent of the pairwise IBD analysis.

**Step 4.3 — Senegalese ancestry in Gambian drug resistance haplotypes**

Compare the genetic background of drug resistance haplotypes in Gambian samples to their Senegalese counterparts. If Gambian resistance alleles are embedded in genetic backgrounds more closely related to Senegal than to local Gambia at flanking neutral loci, this suggests resistance was imported rather than locally evolved — a finding of direct relevance to SMC programme planning.


### Phase 5 — Synthesis and policy implications

**Step 5.1 — Attribution model**

Integrate evidence from Steps 3.1–3.3 and Phase 4 into a coherent attribution statement. Using a simple Bayesian framework, estimate the posterior probability that a randomly selected western Gambia infection in a given season was:
- Locally transmitted (P_local)
- Imported from eastern Gambia (P_east)  
- Imported from Senegal (P_senegal)

Priors derived from: ADMIXTURE ancestry proportions (informative prior), travel history data (informative likelihood), IBD-based transmission path counts (likelihood). This produces a quantitative, uncertainty-bounded answer to the central question.

**Step 5.2 — Sensitivity analysis**

Repeat the attribution analysis under: (a) different IBD thresholds, (b) excluding travel-constrained samples, (c) using barcode data as a validation dataset for the 2013 season where both barcode and WGS are available.


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


## 7. Statistical Framework

| Analysis | Tool / Package |
|---|---|
| VCF processing and filtering | PLINK2, `vcfR`, `tidyverse` |
| PCA and FST | PLINK2, `hierfstat`, `pegas`, `SNPRelate` |
| ADMIXTURE ancestry | ADMIXTURE v1.3, `LEA`, fineSTRUCTURE |
| IBD estimation | hmmIBD, dcifer, isoRelate |
| Transmission path model | Custom R (extends Amambua-Ngwa 2019 model) |
| Population graph | TreeMix v1.13 |
| Isolation-by-distance | `ade4` (Mantel test), `vegan` |
| LD analysis | PLINK2, `genetics` |
| Visualisation | `ggplot2`, `patchwork`, `leaflet`, `circlize` |


## 8. Expected Results and Key Figures

**Figure 1 — Genetic landscape of Gambian malaria**
- Panel A: Geographic map with pie charts showing ADMIXTURE ancestry proportions per village
- Panel B: IBD network plot (chord diagram, nodes = villages, edges = high-IBD pairs)
- Panel C: FST matrix heatmap with hierarchical clustering
- *Message:* Gambian parasites show detectable regional structure, with the east and west forming partially distinct genetic clusters while central Gambia acts as a genetic bridge.

**Figure 2 — Importation signal in western Gambia**
- Panel A: Distribution of eastern ancestry proportions in western Gambia samples (violin by season)
- Panel B: Comparison of within-region IBD (west vs east vs central)
- Panel C: π per region per season — does western Gambia show excess diversity relative to its prevalence?
- *Message:* Western Gambia infections carry significant eastern-ancestry genomic signal, particularly at peak season, inconsistent with primarily local transmission.

**Figure 3 — Directional transmission paths (WGS resolution)**
- Panel A: Transmission path map (extends Fig 5 of Amambua-Ngwa 2019 at higher resolution)
- Panel B: Directional asymmetry: east→west vs west→east path counts per season
- Panel C: Seasonal timing of cross-regional paths
- *Message:* East-to-west parasite flow is confirmed at WGS resolution with substantially greater statistical confidence and a greater number of detected paths than was possible with barcodes.

**Figure 4 — Cross-border importation from Senegal**
- Panel A: Pairwise IBD heatmap, Gambia vs Senegal
- Panel B: TreeMix graph
- Panel C: Ancestry proportions — Senegalese component in Gambian samples per region
- *Message:* Cross-border importation from Senegal is detectable and contributes disproportionately to eastern Gambia, with a further trickle to the west.

**Figure 5 — Attribution and policy implications**
- Panel A: Posterior attribution for western Gambia infection sources (local vs east vs Senegal)
- Panel B: Temporal trend in attribution (has importation fraction changed over studied years?)
- Panel C: Intervention map — geographic targets implied by the attribution model
- *Message:* Elimination of malaria in western Gambia will require coordinated intervention in eastern Gambia and border-zone surveillance; local intervention alone is insufficient.


## 9. Potential Challenges and Mitigations

| Challenge | Mitigation |
|---|---|
| Temporal mismatch between Gambian and Senegalese samples | Use Pf7 data with overlapping collection periods; treat non-contemporaneous comparisons as sensitivity analyses |
| Limited cross-border travel metadata | Use IBD-based pairs as primary evidence; treat travel constraint as a supplementary filter rather than a prerequisite |
| Polyclonal infections diluting the importation signal | Use dcifer for polyclonal-aware relatedness; sensitivity analysis using monoclonal samples only |
| East–west genetic differentiation may be confounded by distance | Explicitly model isolation-by-distance baseline; test for excess cross-regional relatedness above the IBD expectation |
| Definition of "importation" is ambiguous at the parasite level | Distinguish parasite-level importation (a genome that originated outside the region) from human importation (an infected traveller). State explicitly which is being measured. |


## 10. Specific Scientific Contribution

1. **The first WGS-resolution quantification of parasite importation contributing to residual malaria transmission in a near-elimination African setting**, providing a genomic answer to a question of direct programmatic importance.

2. **A validated, high-resolution extension of the Amambua-Ngwa et al. (2019) transmission path model** demonstrating the improvement in detection power achievable with WGS over barcode data, a methodological contribution applicable to any setting with a genomic surveillance archive.

3. **The first cross-border Gambia–Senegal parasite connectivity analysis at WGS resolution**, characterising the role of Senegalese importation in sustaining Gambian transmission and providing a template for genomic border surveillance applicable across the sub-region.

4. **A quantitative attribution framework** — local versus imported transmission that can be updated as new data arrive and used directly by national programmes to allocate intervention resources.


## 12. Target Journals

**Primary:** *Nature Communications* | *The Lancet Infectious Diseases* | *eLife*  
**Secondary:** *Malaria Journal* | *PLOS Medicine*  


# References

1. Okafor, U.A., Kakou, PC.K., D’Alessandro, U. et al. Projection of future malaria prevalence in the upper river region of The Gambia. Malar J 24, 108 (2025). https://doi.org/10.1186/s12936-025-05348-z

2. Khan O, Ajadi JO, Hossain MP. Predicting malaria outbreak in The Gambia using machine learning techniques. PLoS One. 2024 May 16;19(5):e0299386. doi: 10.1371/journal.pone.0299386. PMID: 38753678; PMCID: PMC11098333.
