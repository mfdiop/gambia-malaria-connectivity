# =============================================================================
# Phase 2 — FST + IBD-Mantel, PCA + sNMF sweep, hmmIBD
# =============================================================================
# The sNMF K-sweep is split into per-K jobs (phase2_pca_lea) plus a cross-K
# aggregator (phase2_pca_lea_extendK). See plan Task 5 Pattern A.
# =============================================================================

rule phase2_fst_ibd_mantel:
    input:
        script   = "codes/phase2_01_fst_ibd_mantel.R",
        mono_vcf = rules.phase1_coi.output.mono,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        fst_rds  = "results/phase2_01/fst_results.rds",
        mantel   = "results/phase2_01/mantel_results.tsv",
        summary  = "results/phase2_01/fst_summary.txt",
    params:
        out_dir = "results/phase2_01",
    log:
        "logs/phase2_fst_ibd_mantel.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} > {log} 2>&1"


rule phase2_pca_lea:
    input:
        script   = "codes/phase2_02_pca_lea.R",
        mono_vcf = rules.phase1_coi.output.mono,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        qmatrix = "results/phase2_02_lea/K{k}/qmatrix.tsv",
        ce      = "results/phase2_02_lea/K{k}/ce.tsv",
        geno    = "results/phase2_02_lea/K{k}/snmf/gambia_mono_pruned.geno",
        rds     = "results/phase2_02_lea/K{k}/phase2_02_lea_results.rds",
    params:
        out_dir = "results/phase2_02_lea/K{k}",
    log:
        "logs/phase2_pca_lea_K{k}.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} --k={wildcards.k} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} > {log} 2>&1"


rule phase2_pca_lea_extendK:
    input:
        script   = "codes/phase2_02_pca_lea_extendK.R",
        ce_files = expand("results/phase2_02_lea/K{k}/ce.tsv", k=config["k_range"]),
    output:
        plot    = "results/phase2_02_lea/ce_vs_K.pdf",
        summary = "results/phase2_02_lea/K_summary.tsv",
    params:
        out_dir      = "results/phase2_02_lea",
        ce_files_csv = lambda w, input: ",".join(input.ce_files),
    log:
        "logs/phase2_pca_lea_extendK.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--ce_files={params.ce_files_csv} > {log} 2>&1"


rule phase2_hmmibd:
    input:
        script    = "codes/phase2_03_hmmibd.R",
        mono_vcf  = rules.phase1_coi.output.mono,
        meta      = rules.phase1_annotate.output.meta_rds,
        centroids = rules.phase1_annotate.output.centroids,
    output:
        hmm_fract    = "results/phase2_03/hmmibd.hmm_fract.txt",
        hmm_segments = "results/phase2_03/hmmibd.hmm.txt",
        hap_input    = "results/phase2_03/hmmibd_input.txt",
        pair_class   = "results/phase2_03/pair_ibd_classified.tsv",
        ibd_matrix   = "results/phase2_03/ibd_matrix.rds",
        summary      = "results/phase2_03/phase2_03_summary.txt",
    params:
        out_dir = "results/phase2_03",
    log:
        "logs/phase2_hmmibd.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} "
        "--centroid={input.centroids} > {log} 2>&1"
