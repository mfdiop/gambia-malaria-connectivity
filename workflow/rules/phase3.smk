# =============================================================================
# Phase 3 — supervised ancestry, directional IBD, genomic isolation
# =============================================================================
# phase3_01 is supervised admixture with K = 3 fixed by the regional design,
# so no {k} wildcard here — it consumes the canonical K from Phase 2.2 only
# for the per-sample geno matrix and sample order.
# =============================================================================

CANONICAL_K = config["canonical_k"]

rule phase3_supervised_admixture:
    input:
        script    = "codes/phase3_01_supervised_admixture.R",
        geno      = f"results/phase2_02_lea/K{CANONICAL_K}/snmf/gambia_mono_pruned.geno",
        phase22   = f"results/phase2_02_lea/K{CANONICAL_K}/phase2_02_lea_results.rds",
        meta      = rules.phase1_annotate.output.meta_rds,
        centroids = rules.phase1_annotate.output.centroids,
    output:
        ancestry_q   = "results/phase3_01/ancestry_Q.tsv",
        import_flags = "results/phase3_01/importation_flags.tsv",
        region_summ  = "results/phase3_01/ancestry_region_summary.tsv",
        rds          = "results/phase3_01/phase3_01_results.rds",
        summary      = "results/phase3_01/phase3_01_summary.txt",
    params:
        out_dir = "results/phase3_01",
    log:
        "logs/phase3_supervised_admixture.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--geno={input.geno} --phase22={input.phase22} "
        "--meta={input.meta} --centroid={input.centroids} > {log} 2>&1"


rule phase3_directional_ibd:
    input:
        script     = "codes/phase3_02_directional_ibd.R",
        pair_class = rules.phase2_hmmibd.output.pair_class,
        meta       = rules.phase1_annotate.output.meta_rds,
    output:
        trans_pairs = "results/phase3_02/transmission_pairs.tsv",
        accepted    = "results/phase3_02/accepted_paths.tsv",
        asymmetry   = "results/phase3_02/directional_asymmetry.tsv",
        rds         = "results/phase3_02/phase3_02_results.rds",
        summary     = "results/phase3_02/phase3_02_summary.txt",
    params:
        out_dir = "results/phase3_02",
    log:
        "logs/phase3_directional_ibd.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--pair_class={input.pair_class} --meta={input.meta} > {log} 2>&1"


rule phase3_isolation:
    input:
        script     = "codes/phase3_03_isolation.R",
        hap        = rules.phase2_hmmibd.output.hap_input,
        pair_class = rules.phase2_hmmibd.output.pair_class,
        meta       = rules.phase1_annotate.output.meta_rds,
    output:
        ibd_within = "results/phase3_03/ibd_within_region.tsv",
        pi_table   = "results/phase3_03/pi_region_season.tsv",
        ld_decay   = "results/phase3_03/ld_decay.tsv",
        rds        = "results/phase3_03/phase3_03_results.rds",
        summary    = "results/phase3_03/phase3_03_summary.txt",
    params:
        out_dir = "results/phase3_03",
    log:
        "logs/phase3_isolation.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--hap={input.hap} --pair_class={input.pair_class} "
        "--meta={input.meta} > {log} 2>&1"
