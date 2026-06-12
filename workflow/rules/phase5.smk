# =============================================================================
# Phase 5 — Bayesian attribution model + sensitivity analysis
# =============================================================================
# Both scripts are hardcoded to TARGET_REGION = "west" and
# DECISION_YEARS = c(2014, 2015), so neither rule uses a {subset} wildcard.
# Phase 5.1 has an internal sensitivity grid (ALPHA x SENEGAL_BASELINE x
# P_TRANS_FLOOR x EVIDENCE); Phase 5.2 explores the orthogonal axes
# (within-Gambia P_trans floor x cross-border IBD threshold, travel filter
# hook, barcode comparator hook).
# =============================================================================

rule phase5_attribution:
    input:
        script    = "codes/phase5_01_attribution_model.R",
        ancestry  = rules.phase3_supervised_admixture.output.ancestry_q,
        pairs     = rules.phase3_directional_ibd.output.trans_pairs,
        crossb    = rules.phase4_decision_subset.output.high_pairs,
        meta      = rules.phase1_annotate.output.meta_rds,
    output:
        per_sample   = "results/phase5_01/per_sample_posteriors.tsv",
        season       = "results/phase5_01/season_posteriors.tsv",
        sensitivity  = "results/phase5_01/sensitivity_grid.tsv",
        evidence     = "results/phase5_01/evidence_per_sample.tsv",
        rds          = "results/phase5_01/phase5_01_results.rds",
        summary      = "results/phase5_01/phase5_01_summary.txt",
    params:
        out_dir = "results/phase5_01",
    log:
        "logs/phase5_attribution.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--ancestry={input.ancestry} --pairs={input.pairs} "
        "--crossb={input.crossb} --meta={input.meta} > {log} 2>&1"


rule phase5_sensitivity:
    input:
        script     = "codes/phase5_02_sensitivity.R",
        ancestry   = rules.phase3_supervised_admixture.output.ancestry_q,
        pairs      = rules.phase3_directional_ibd.output.trans_pairs,
        crossb_all = rules.phase4_decision_subset.output.pair_ibd,
        meta       = rules.phase1_annotate.output.meta_rds,
    output:
        threshold_per_sample = "results/phase5_02/threshold_grid_per_sample.tsv",
        threshold_season     = "results/phase5_02/threshold_grid_season.tsv",
        rds                  = "results/phase5_02/phase5_02_results.rds",
        summary              = "results/phase5_02/phase5_02_summary.txt",
    params:
        out_dir     = "results/phase5_02",
        barcode_dir = "data/barcode_amambua",
    log:
        "logs/phase5_sensitivity.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--ancestry={input.ancestry} --pairs={input.pairs} "
        "--crossb_all={input.crossb_all} --meta={input.meta} "
        "--barcode_dir={params.barcode_dir} > {log} 2>&1"
