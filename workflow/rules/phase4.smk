# =============================================================================
# Phase 4 — cross-border IBD, TreeMix m-sweep + OptM, resistance haplotype
# =============================================================================
# Only the TreeMix rules use the {subset} wildcard (per design D6). The
# crossborder-IBD and resistance scripts have fixed sample-set logic and run
# once. The {m} wildcard sweeps TreeMix migration-edge counts; phase4_02b
# aggregates across m per subset.
# =============================================================================

rule phase4_crossborder_ibd:
    input:
        script      = "codes/phase4_01_crossborder_ibd.R",
        phase1_meta = rules.phase1_annotate.output.meta_rds,
    output:
        pair_ibd   = "results/phase4_01/crossborder_pair_ibd.tsv",
        high_pairs = "results/phase4_01/crossborder_high_ibd_pairs.tsv",
        senegam_vcf = "results/phase4_01/senegam_core_biSNP_maf.vcf.gz",
        summary    = "results/phase4_01/phase4_01_summary.txt",
    params:
        out_dir = "results/phase4_01",
    log:
        "logs/phase4_crossborder_ibd.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--phase1_meta={input.phase1_meta} > {log} 2>&1"


rule phase4_decision_subset:
    input:
        script  = "codes/phase4_01_decision_subset.R",
        pair_in = rules.phase4_crossborder_ibd.output.pair_ibd,
    output:
        pair_ibd   = "results/phase4_01_decision_subset/pair_ibd.tsv",
        high_pairs = "results/phase4_01_decision_subset/high_ibd_pairs.tsv",
        summary    = "results/phase4_01_decision_subset/summary.txt",
    params:
        out_dir = "results/phase4_01_decision_subset",
    log:
        "logs/phase4_decision_subset.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--pair_in={input.pair_in} > {log} 2>&1"


rule phase4_treemix:
    input:
        script      = "codes/phase4_02a_treemix.R",
        senegam_vcf = rules.phase4_crossborder_ibd.output.senegam_vcf,
    output:
        # phase4_02a_summary.txt is written last; serves as sentinel for the
        # variable-arity runs/ dir.
        summary = "results/phase4_02/{subset}/m{m}/phase4_02a_summary.txt",
    params:
        out_dir = "results/phase4_02/{subset}/m{m}",
        n_reps  = 10,
    conda:
        "../../envs/treemix.yaml"
    log:
        "logs/phase4_treemix_{subset}_m{m}.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} --subset={wildcards.subset} "
        "--m={wildcards.m} --n_reps={params.n_reps} > {log} 2>&1"


def _treemix_m_sentinels(wildcards):
    return [
        f"results/phase4_02/{wildcards.subset}/m{m}/phase4_02a_summary.txt"
        for m in config["m_range"]
    ]


rule phase4_treemix_optm:
    input:
        script    = "codes/phase4_02b_optm.R",
        sentinels = _treemix_m_sentinels,
    output:
        optm_summary = "results/phase4_02/{subset}/optm_summary.tsv",
        optm_curve   = "results/phase4_02/{subset}/optm_curve.pdf",
        best_m       = "results/phase4_02/{subset}/best_m.txt",
        summary      = "results/phase4_02/{subset}/phase4_02b_summary.txt",
    params:
        out_dir     = "results/phase4_02/{subset}",
        runs_dirs   = lambda w: ",".join(
            f"results/phase4_02/{w.subset}/m{m}" for m in config["m_range"]
        ),
        n_bootstrap = 100,
    conda:
        "../../envs/treemix.yaml"
    log:
        "logs/phase4_treemix_optm_{subset}.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} --subset={wildcards.subset} "
        "--runs_dirs={params.runs_dirs} "
        "--n_bootstrap={params.n_bootstrap} > {log} 2>&1"


rule phase4_resistance:
    input:
        script      = "codes/phase4_03_resistance_haplotype_ancestry.R",
        senegam_vcf = rules.phase4_crossborder_ibd.output.senegam_vcf,
    output:
        carriers = "results/phase4_03/carriers_per_locus.tsv",
        locus_ibs = "results/phase4_03/locus_flank_ibs.tsv",
        stats    = "results/phase4_03/locus_test_stats.tsv",
        summary  = "results/phase4_03/phase4_03_summary.txt",
    params:
        out_dir = "results/phase4_03",
    log:
        "logs/phase4_resistance.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--vcf={input.senegam_vcf} > {log} 2>&1"
