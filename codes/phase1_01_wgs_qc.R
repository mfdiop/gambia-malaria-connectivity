# =============================================================================
# Phase 1 — Step 1.1: WGS quality control
# =============================================================================
# Implements the QC filters specified in notes/objective3_proposal.md §6 Phase 1.
#
# Input  : data/gambia_2014isolates/gambia.ann.vcf.gz   (160 samples, 16,390 biallelic SNPs)
#          data/metadata/GamMetadata_2014.xlsx          (180 metadata rows)
#
# Filters (per proposal):
#   - Mean per-sample coverage >= 5x
#   - Per-sample missing genotype rate <= 20%
#   - Minor allele frequency > 0.01 across full dataset
#   - Exclude hypervariable antigen loci (var, rifin, stevor) — enforced by
#     RegionType == "Core" (subtelomeric/internal hypervariable regions are
#     pre-annotated in the VCF and contain var/rifin/stevor families)
#   - Retain biallelic SNPs only (already enforced upstream; re-checked)
#
# Output : results/phase1_01/
#            gambia_2014_qc.vcf.gz          filtered VCF (bgzipped + indexed)
#            sample_qc_metrics.tsv          per-sample coverage + missingness
#            site_qc_summary.tsv            per-step site counts
#            qc_dropped_samples.tsv         samples removed and reason
#            qc_summary.txt                 human-readable summary log
#            sample_missingness.pdf         diagnostic plot
#            sample_meandepth.pdf           diagnostic plot
#
# Tools  : bcftools (system), vcftools (system), R (vcfR, tidyverse, readxl)
# =============================================================================

suppressPackageStartupMessages({
   library(tidyverse)
   library(readxl)
   library(vcfR)
})

# ----- Paths ----------------------------------------------------------------
vcf_in      <- "data/gambia_2014isolates/gambia.ann.vcf.gz"
meta_in     <- "data/metadata/GamMetadata_2014.xlsx"
out_dir     <- "results/phase1_01"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vcf_core    <- file.path(out_dir, "gambia_2014_core.vcf.gz")           # after RegionType==Core
vcf_biSNP   <- file.path(out_dir, "gambia_2014_core_biSNP.vcf.gz")     # after biallelic SNP filter
vcf_maf     <- file.path(out_dir, "gambia_2014_core_biSNP_maf.vcf.gz") # after MAF filter
vcf_final   <- file.path(out_dir, "gambia_2014_qc.vcf.gz")             # after sample drop

# ----- Thresholds (proposal §6 Phase 1, Step 1.1) ---------------------------
MIN_MEAN_DP    <- 5      # mean coverage per sample
MAX_MISSING    <- 0.20   # max per-sample missing-genotype rate
MIN_MAF        <- 0.01   # site-level MAF

# ----- Helper to log --------------------------------------------------------
log_lines <- c()
log_msg <- function(...) {
   msg <- paste0(...)
   log_lines <<- c(log_lines, msg)
   message(msg)
}

count_sites <- function(vcf_path) {
   as.integer(system(sprintf("bcftools view -H %s | wc -l", vcf_path), intern = TRUE))
}

# ============================================================================
# 1. Restrict to Core regions (excludes var/rifin/stevor)
# ============================================================================
log_msg("=== Phase 1 Step 1.1: WGS QC ===")
log_msg("Input VCF: ", vcf_in)

n0 <- count_sites(vcf_in)
log_msg("Sites in raw VCF: ", n0)

system(sprintf(
   "bcftools view -i 'INFO/RegionType==\"Core\"' -Oz -o %s %s && bcftools index -t -f %s",
   vcf_core, vcf_in, vcf_core
))
n_core <- count_sites(vcf_core)
log_msg("Sites after RegionType==Core (excludes var/rifin/stevor): ", n_core,
        " (dropped ", n0 - n_core, ")")

# ============================================================================
# 2. Biallelic SNPs only
# ============================================================================
system(sprintf(
   "bcftools view -m2 -M2 -v snps -Oz -o %s %s && bcftools index -t -f %s",
   vcf_biSNP, vcf_core, vcf_biSNP
))
n_bi <- count_sites(vcf_biSNP)
log_msg("Sites after biallelic-SNP filter: ", n_bi, " (dropped ", n_core - n_bi, ")")

# ============================================================================
# 3. Site-level MAF > 0.01
# ============================================================================
system(sprintf(
   "bcftools view -q %f:minor -Oz -o %s %s && bcftools index -t -f %s",
   MIN_MAF, vcf_maf, vcf_biSNP, vcf_maf
))
n_maf <- count_sites(vcf_maf)
log_msg("Sites after MAF > ", MIN_MAF, ": ", n_maf, " (dropped ", n_bi - n_maf, ")")

# ============================================================================
# 4. Per-sample coverage + missingness
# ============================================================================
log_msg("Computing per-sample mean depth and missingness …")

# Mean depth per sample (vcftools)
system(sprintf("vcftools --gzvcf %s --depth --out %s 2> /dev/null",
               vcf_maf, file.path(out_dir, "sample_depth")))
dep <- read_tsv(file.path(out_dir, "sample_depth.idepth"),
                show_col_types = FALSE) %>%
   rename(SampleID = INDV, mean_depth = MEAN_DEPTH, n_sites = N_SITES)

# Missingness per sample (vcftools)
system(sprintf("vcftools --gzvcf %s --missing-indv --out %s 2> /dev/null",
               vcf_maf, file.path(out_dir, "sample_missing")))
mis <- read_tsv(file.path(out_dir, "sample_missing.imiss"),
                show_col_types = FALSE) %>%
   rename(SampleID = INDV, f_missing = F_MISS) %>%
   select(SampleID, f_missing)

sample_qc <- dep %>%
   inner_join(mis, by = "SampleID") %>%
   mutate(pass_depth   = mean_depth >= MIN_MEAN_DP,
          pass_missing = f_missing  <= MAX_MISSING,
          pass         = pass_depth & pass_missing)

write_tsv(sample_qc, file.path(out_dir, "sample_qc_metrics.tsv"))

n_drop_dp  <- sum(!sample_qc$pass_depth)
n_drop_mis <- sum(!sample_qc$pass_missing)
n_pass     <- sum(sample_qc$pass)
log_msg(sprintf("Samples failing depth (<%dx): %d", MIN_MEAN_DP, n_drop_dp))
log_msg(sprintf("Samples failing missingness (>%.0f%%): %d",
                100 * MAX_MISSING, n_drop_mis))
log_msg(sprintf("Samples passing all filters: %d / %d", n_pass, nrow(sample_qc)))

# Record dropped samples
dropped <- sample_qc %>%
   filter(!pass) %>%
   mutate(reason = case_when(
      !pass_depth   & !pass_missing ~ "low_depth+high_missing",
      !pass_depth                   ~ "low_depth",
      !pass_missing                 ~ "high_missing"))
write_tsv(dropped, file.path(out_dir, "qc_dropped_samples.tsv"))

# Write keep-list for bcftools
keep_file <- file.path(out_dir, "samples_keep.txt")
writeLines(sample_qc$SampleID[sample_qc$pass], keep_file)

# ============================================================================
# 5. Apply sample filter and re-check MAF (alleles can shift after dropping samples)
# ============================================================================
system(sprintf(
   "bcftools view -S %s --force-samples -Oz %s | bcftools view -q %f:minor -Oz -o %s - && bcftools index -t -f %s",
   keep_file, vcf_maf, MIN_MAF, vcf_final, vcf_final
))
n_final <- count_sites(vcf_final)
log_msg("Final QC VCF: ", vcf_final)
log_msg("Final sites: ", n_final)
log_msg("Final samples: ", n_pass)

# ============================================================================
# 6. Per-step site count summary
# ============================================================================
site_summary <- tibble(
   step  = c("raw", "core_only", "biallelic_snp", "maf_filter", "after_sample_drop_remaf"),
   sites = c(n0, n_core, n_bi, n_maf, n_final)) %>%
   mutate(dropped = lag(sites) - sites)
write_tsv(site_summary, file.path(out_dir, "site_qc_summary.tsv"))

# ============================================================================
# 7. Cross-check QC-passed samples against metadata
# ============================================================================
meta <- read_xlsx(meta_in)
meta_2014 <- meta %>% filter(Year == 2014)

vcf_samples <- sample_qc$SampleID[sample_qc$pass]
n_in_meta   <- sum(vcf_samples %in% meta_2014$SampleID)

log_msg(sprintf("QC-passed samples present in 2014 metadata: %d / %d",
                n_in_meta, length(vcf_samples)))
missing_meta <- setdiff(vcf_samples, meta_2014$SampleID)
if (length(missing_meta) > 0) {
   log_msg("WARNING — samples passing QC but absent from 2014 metadata: ",
           paste(missing_meta, collapse = ", "))
}

# ============================================================================
# 8. Diagnostic plots
# ============================================================================
p_miss <- ggplot(sample_qc, aes(x = f_missing)) +
   geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
   geom_vline(xintercept = MAX_MISSING, linetype = 2, colour = "red") +
   labs(title = "Per-sample missingness",
        subtitle = sprintf("Threshold: F_MISS <= %.2f (red)", MAX_MISSING),
        x = "Missing-genotype fraction", y = "Samples") +
   theme_minimal()

p_dp <- ggplot(sample_qc, aes(x = mean_depth)) +
   geom_histogram(bins = 40, fill = "darkgreen", colour = "white") +
   geom_vline(xintercept = MIN_MEAN_DP, linetype = 2, colour = "red") +
   labs(title = "Per-sample mean coverage",
        subtitle = sprintf("Threshold: mean DP >= %dx (red)", MIN_MEAN_DP),
        x = "Mean depth (x)", y = "Samples") +
   theme_minimal()

ggsave(file.path(out_dir, "sample_missingness.pdf"), p_miss,
       width = 6, height = 4, dpi = 600)
ggsave(file.path(out_dir, "sample_meandepth.pdf"),   p_dp,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 9. Write human-readable summary
# ============================================================================
writeLines(log_lines, file.path(out_dir, "qc_summary.txt"))

# Clean intermediate VCFs to save space; keep final + indices
unlink(c(vcf_core, paste0(vcf_core, ".tbi"),
         vcf_biSNP, paste0(vcf_biSNP, ".tbi"),
         vcf_maf, paste0(vcf_maf, ".tbi")))

message("\nDone. Outputs in ", out_dir)
