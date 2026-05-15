
## Implementation Guide for the Amambua-Ngwa et al. 2019 Method

This is a complete, self-contained guide written so that Claude Code can implement it accurately from scratch without needing to read the original paper.

---

### What the Method Does

The model estimates the probability that a pair of parasite isolates represents a true person-to-person transmission event, given:
1. Their genetic distance
2. The time separating their sampling dates
3. Whether at least one person in the pair travelled between their sampling locations

It does this by simulating the expected genetic distance between transmission-linked pairs under the *P. falciparum* life cycle timing constraints, then asking: given the observed genetic distance and time gap, what is the probability this is a true transmission pair?

---

### Step 1 — Compute pairwise genetic distance

For each pair of isolates $(i, j)$, compute the p-distance:

$$d_{ij} = \frac{\text{number of SNP positions where alleles differ}}{\text{number of positions with non-missing calls in both samples}}$$

For WGS data, use only biallelic SNPs passing quality filters. For polyclonal samples, use the minor (rarer) allele at mixed positions — this maximises sensitivity for detecting shared rare variants between transmission-linked pairs, as described in the original paper.

Restrict analysis to pairs where $d_{ij} \leq 0.30$. Beyond 0.30, the relationship between p-distance and evolutionary distance becomes non-linear due to recombination and homoplasy, making probability assignment unreliable.

---

### Step 2 — Simulate the life cycle time constraint

The model constrains which pairs can be transmission-linked by requiring that the sampling time difference $\Delta t_{ij} = |t_j - t_i|$ (in days) is consistent with at least one full *P. falciparum* person-to-vector-to-person transmission cycle.

Simulate the total cycle duration by drawing 100,000 realisations of the following stochastic components:

| Component | Symbol | Distribution | Parameters |
|---|---|---|---|
| Time from bite to blood-stage detectable | $T_i$ | Triangular | min=5, max=15, mode=10 days |
| Time for vector to become infectious (sporogony) | $T_v$ | Triangular | min=8, max=15, mode=11.5 days |
| Time from infectious bite to next sampling | $T_d$ | Derived | $T_d = \text{vector lifespan} - T_v$ |
| Vector lifespan | — | Triangular | min=18, max=24, mode=21 days |
| Time since last sampling for source | $T_A$ | Uniform | min=0, max=30 days |
| Time since last sampling for recipient | $T_B$ | Uniform | min=0, max=30 days |
| Bite-to-sampling lag for source | $b$ | Uniform | min=0, max=$T_A$ |

For a **direct transmission path** ($A \rightarrow B$), the total evolutionary time between the sampled genotypes is:

$$\Delta t_{\text{evo}} = T_v + T_d + T_i + T_B - b$$

For a **common source path** (both $A$ and $B$ infected by the same mosquito event), the evolutionary time is:

$$\Delta t_{\text{evo}} = T_d + T_i^B + T_B - T_i^A + T_A$$

Run 100,000 simulations of each formula. The resulting distribution of $\Delta t_{\text{evo}}$ values defines the plausible range of observation time differences consistent with a single transmission event. Restrict to pairs whose observed $\Delta t_{ij}$ falls to the right of the 97.5th percentile of the same-source distribution — this ensures you are detecting directional transmission rather than common-source co-infection. In the original paper this produced the window of **32 to 67 days**.

---

### Step 3 — Build the bivariate probability density

Using the 100,000 simulated evolutionary times from Step 2, construct a joint probability density of genetic distance and observation time difference for transmission-linked pairs:

1. For each simulated $\Delta t_{\text{evo}}$, evolve 10,000 random barcodes of 54 SNPs (or for WGS: of length $L$ SNPs) for that duration using the maximum-likelihood substitution model (TVM+I+G4 for the original barcode study; select the best model using IQ-TREE2 for your WGS data). This generates a simulated genetic distance $d_{\text{sim}}$.

2. Collect all $(d_{\text{sim}}, \Delta t_{\text{evo}})$ pairs. Apply 2D kernel density estimation to produce the joint density $f(d, \Delta t \mid \text{transmission})$.

3. For pairs with $d_{ij} > 0.30$, assign probability zero.

4. For each observed pair $(i, j)$ with $d_{ij} \leq 0.30$ and $\Delta t_{ij}$ in the valid window, compute:

$$P(\text{transmission} \mid d_{ij}, \Delta t_{ij}) = \frac{f(d_{ij}, \Delta t_{ij} \mid \text{transmission})}{f(d_{ij}, \Delta t_{ij} \mid \text{transmission}) + f(d_{ij}, \Delta t_{ij} \mid \text{background})}$$

The background density $f(d, \Delta t \mid \text{background})$ is estimated from all pairs in your dataset — it is the empirical distribution of genetic distances and time differences across all observed pairs regardless of transmission status.

---

### Step 4 — Assign maximum-likelihood transmission paths

Build a transmission probability matrix $M$ where entry $M_{ij}$ is the transmission probability computed in Step 3 for the pair $(i, j)$.

Apply a bilinear maximisation to resolve ambiguous paths:
- For each source $i$, find the recipient $j^*$ that maximises $M_{ij}$
- For each recipient $j$, find the source $i^*$ that maximises $M_{ij}$
- A path is accepted only if the maximisation is consistent in both directions: $j^* = \arg\max_j M_{ij}$ AND $i^* = \arg\max_i M_{ij^*}$

This produces a set of non-overlapping, maximum-likelihood transmission paths where each isolate is the source of at most one detected transmission and each recipient has at most one detected source.

---

### Step 5 — Apply the travel constraint

For cross-regional paths (east-to-west or west-to-east), the original paper requires that at least one member of the pair has documented inter-regional travel within the relevant time window:
- Source: travel within 30 days before sampling
- Recipient: travel within 60 days before sampling (the full life cycle window)

For your 2014 dataset, apply this constraint using whatever travel or movement data is available from the cohort. For pairs where no travel data exists, report results both with and without the constraint as a sensitivity analysis.

---

### Step 6 — Extension for WGS data

The original paper used 54 SNPs and p-distance. For WGS, extend as follows:

**Replace p-distance with IBD estimation.** Use hmmIBD to estimate pairwise IBD proportions for monoclonal pairs. For polyclonal pairs, use dcifer. The IBD score replaces p-distance as the genetic relatedness metric. Adjust the threshold: instead of $d \leq 0.30$, use $\text{IBD} \geq 0.10$ as the lower bound for candidate transmission pairs (pairs below this IBD are almost certainly unrelated).

**Replace the substitution model simulation with tree-sequence IBD.** If your WGS data includes tree sequences from the simulation pipeline (which yours does), you can compute exact pairwise IBD from the true genealogy rather than simulating it. For real WGS data, use the empirical IBD distribution from household pairs (which you know are co-resident and therefore likely transmission-linked) as the transmission distribution $f(d \mid \text{transmission})$ rather than simulating it from a substitution model.

**Adjust the time window.** The 32–67 day window was derived for a 54-SNP barcode with the TVM+I+G4 model. For WGS, the window should be re-derived using the WGS substitution rate rather than the barcode substitution rate. Re-run the Step 2 simulation using your WGS-calibrated substitution rate.

---

### Deliverable for Claude Code

The full method requires these code components, which Claude Code should implement as separate, testable functions:

```
1. compute_pairwise_distance(genotype_matrix, method = "pdist" | "IBD")
   → returns n×n distance matrix

2. simulate_lifecycle_times(n_sims = 100000)
   → returns vector of valid delta_t_evo values and the 32–67 day window

3. fit_substitution_model(alignment, model = "TVM+I+G4")
   → returns rate parameters for genome evolution simulation

4. simulate_genetic_distances(delta_t_vec, rate_params, n_snps, n_sims = 10000)
   → returns simulated (d_sim, delta_t) pairs

5. build_bivariate_density(sim_pairs, bandwidth = "silverman")
   → returns 2D KDE object for f(d, delta_t | transmission)

6. compute_transmission_probabilities(dist_matrix, time_matrix, density_obj, threshold = 0.30)
   → returns n×n probability matrix M

7. apply_travel_constraint(M, travel_df, source_window = 30, recipient_window = 60)
   → returns filtered M

8. maximise_transmission_paths(M)
   → returns list of accepted (source_id, recipient_id, probability) triples

9. classify_paths_by_region(paths_df, sample_metadata)
   → returns paths annotated as within-west | within-east | east-to-west | west-to-east | cross-border

10. compute_directional_asymmetry(classified_paths)
    → returns counts and bootstrap CI for each direction
```

Each function should be independently testable on synthetic data before being connected into the pipeline. Start with function 1 and verify it reproduces the pairwise distance distribution from the original paper's Supplementary Fig. 7 before proceeding.