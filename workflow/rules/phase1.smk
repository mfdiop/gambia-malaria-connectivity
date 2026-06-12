# =============================================================================
# Phase 1 — QC, COI stratification, sample annotation
# =============================================================================
# Each rule wraps one R workbook in codes/. Scripts read sibling-repo VCFs and
# metadata via their own hardcoded paths (D2 in the design spec); Snakemake
# controls only the per-rule output directory via --output_dir=.
# =============================================================================

rule phase1_qc:
    input:
        script = "codes/phase1_01_wgs_qc.R",
    output:
        vcf      = "results/phase1_01/gambia_2014_qc.vcf.gz",
        keep     = "results/phase1_01/samples_keep.txt",
        metrics  = "results/phase1_01/sample_qc_metrics.tsv",
        dropped  = "results/phase1_01/qc_dropped_samples.tsv",
        summary  = "results/phase1_01/qc_summary.txt",
    params:
        out = "results/phase1_01",
    log:
        "logs/phase1_qc.log",
    shell:
        "Rscript --vanilla {input.script} --output_dir={params.out} > {log} 2>&1"


rule phase1_coi:
    input:
        script  = "codes/phase1_02_coi_stratification.R",
        vcf     = rules.phase1_qc.output.vcf,
        keep    = rules.phase1_qc.output.keep,
    output:
        strat   = "results/phase1_02/sample_coi_stratification.tsv",
        mono    = "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz",
        polylow = "results/phase1_02/gambia_2014_qc_polyclonal_low.vcf.gz",
        summary = "results/phase1_02/coi_summary.txt",
    params:
        in_dir  = "results/phase1_01",
        out_dir = "results/phase1_02",
    log:
        "logs/phase1_coi.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--input_dir={params.in_dir} --output_dir={params.out_dir} > {log} 2>&1"


rule phase1_annotate:
    input:
        script  = "codes/phase1_03_sample_annotation.R",
        keep    = rules.phase1_qc.output.keep,
        strat   = rules.phase1_coi.output.strat,
    output:
        meta_tsv  = "results/phase1_03/sample_master_metadata.tsv",
        meta_rds  = "results/phase1_03/sample_master_metadata.rds",
        centroids = "results/phase1_03/village_centroids.tsv",
        summary   = "results/phase1_03/annotation_summary.txt",
    params:
        qc_dir  = "results/phase1_01",
        coi_dir = "results/phase1_02",
        out_dir = "results/phase1_03",
    log:
        "logs/phase1_annotate.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--qc_dir={params.qc_dir} --coi_dir={params.coi_dir} "
        "--output_dir={params.out_dir} > {log} 2>&1"
