# =============================================================================
# Phase 6 — Reviewer-driven extensions
# =============================================================================
# All five rules are single-pass (no wildcards). All work on the Gambia-only
# monoclonal subset; no {subset} wildcard, no cross-border duplication.
# Conda envs only where non-R tools are involved.
# =============================================================================

rule phase6_haplotype_coancestry:
    input:
        script   = "codes/phase6_01_haplotype_coancestry.R",
        mono_vcf = rules.phase1_coi.output.mono,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        chunkcounts  = "results/phase6_01/chunkcounts.tsv",
        coancestry   = "results/phase6_01/coancestry_matrix.rds",
        heatmap      = "results/phase6_01/chunkcounts_heatmap.pdf",
        tree         = "results/phase6_01/fs_tree.newick",
        rds          = "results/phase6_01/phase6_01_results.rds",
        summary      = "results/phase6_01/phase6_01_summary.txt",
    params:
        out_dir = "results/phase6_01",
    conda:
        "../../envs/chromopainter.yaml"
    log:
        "logs/phase6_haplotype_coancestry.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} > {log} 2>&1"


rule phase6_ibd_tracts:
    input:
        script   = "codes/phase6_02_ibd_tracts.R",
        segments = rules.phase2_hmmibd.output.hmm_segments,
        fract    = rules.phase2_hmmibd.output.hmm_fract,
        meta     = rules.phase1_annotate.output.meta_rds,
    output:
        tract_table = "results/phase6_02/tract_table.tsv",
        summary_rp  = "results/phase6_02/tract_summary_by_regionpair.tsv",
        wilcoxon    = "results/phase6_02/wilcoxon_within_vs_cross.tsv",
        density_pdf = "results/phase6_02/tract_density_facet.pdf",
        scatter_pdf = "results/phase6_02/tract_length_vs_total_ibd.pdf",
        rds         = "results/phase6_02/phase6_02_results.rds",
        summary     = "results/phase6_02/phase6_02_summary.txt",
    params:
        out_dir = "results/phase6_02",
    log:
        "logs/phase6_ibd_tracts.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--segments={input.segments} --fract={input.fract} "
        "--meta={input.meta} > {log} 2>&1"


rule phase6_network_topology:
    input:
        script     = "codes/phase6_03_network_topology.R",
        pair_class = rules.phase2_hmmibd.output.pair_class,
        trans      = rules.phase3_directional_ibd.output.trans_pairs,
        meta       = rules.phase1_annotate.output.meta_rds,
    output:
        comm_25     = "results/phase6_03/communities_0.25.tsv",
        comm_50     = "results/phase6_03/communities_0.50.tsv",
        betweenness = "results/phase6_03/betweenness.tsv",
        bridges     = "results/phase6_03/bridges.tsv",
        net_25_pdf  = "results/phase6_03/network_0.25.pdf",
        net_50_pdf  = "results/phase6_03/network_0.50.pdf",
        directed    = "results/phase6_03/directed_indeg_outdeg.tsv",
        rds         = "results/phase6_03/phase6_03_results.rds",
        summary     = "results/phase6_03/phase6_03_summary.txt",
    params:
        out_dir = "results/phase6_03",
    log:
        "logs/phase6_network_topology.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--pair_class={input.pair_class} --trans_pairs={input.trans} "
        "--meta={input.meta} > {log} 2>&1"


rule phase6_symmetric_null:
    input:
        script = "codes/phase6_04_symmetric_migration_null.R",
        trans  = rules.phase3_directional_ibd.output.trans_pairs,
        meta   = rules.phase1_annotate.output.meta_rds,
    output:
        null_dist    = "results/phase6_04/null_distribution.tsv",
        obs_vs_null  = "results/phase6_04/observed_vs_null.tsv",
        obs_pdf      = "results/phase6_04/observed_vs_null.pdf",
        rds          = "results/phase6_04/phase6_04_results.rds",
        summary      = "results/phase6_04/phase6_04_summary.txt",
    params:
        out_dir = "results/phase6_04",
    log:
        "logs/phase6_symmetric_null.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--trans_pairs={input.trans} --meta={input.meta} > {log} 2>&1"


rule phase6_eems:
    input:
        script    = "codes/phase6_05_eems.R",
        mono_vcf  = rules.phase1_coi.output.mono,
        meta      = rules.phase1_annotate.output.meta_rds,
        centroids = rules.phase1_annotate.output.centroids,
    output:
        coord       = "results/phase6_05/eems_inputs/gambia.coord",
        outer       = "results/phase6_05/eems_inputs/gambia.outer",
        diffs       = "results/phase6_05/eems_inputs/gambia.diffs",
        mig_pdf     = "results/phase6_05/eems_mig_surface-mrates01.pdf",
        diag_pdf    = "results/phase6_05/eems_diagnostics-pilogl01.pdf",
        rds         = "results/phase6_05/phase6_05_results.rds",
        summary     = "results/phase6_05/phase6_05_summary.txt",
    params:
        out_dir = "results/phase6_05",
    conda:
        "../../envs/eems.yaml"
    log:
        "logs/phase6_eems.log",
    shell:
        "Rscript --vanilla {input.script} "
        "--output_dir={params.out_dir} "
        "--mono_vcf={input.mono_vcf} --meta={input.meta} "
        "--centroid={input.centroids} > {log} 2>&1"
