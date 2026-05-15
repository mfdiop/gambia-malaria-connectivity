# =============================================================================
# Phase 1 — Step 1.2: COI stratification
# =============================================================================
# Reads pre-computed COI estimates (COIL and McCOIL columns in the 2014
# metadata) and stratifies QC-passed samples into the categories defined in
# notes/objective3_proposal.md §6 Phase 1, Step 1.2:
#
#   - Monoclonal       (COI = 1) — primary set for haplotype-level analyses
#   - Low polyclonal   (COI 2–3) — include in IBD analyses using dcifer
#   - High polyclonal  (COI >= 4) — flag for sensitivity analysis; exclude
#                                   from directional inference
#
# Primary classifier  : McCOIL (THE REAL McCOIL) — the proposal's named tool
# Secondary classifier: COIL — kept for sensitivity analyses
#
# Input  : results/phase1_01/gambia_2014_qc.vcf.gz
#          results/phase1_01/samples_keep.txt
#          data/metadata/GamMetadata_2014.xlsx
#
# Output : results/phase1_02/
#            sample_coi_stratification.tsv   per-sample stratum (both COIs)
#            coi_concordance.tsv             COIL vs McCOIL cross-tab
#            coi_summary.txt                 human-readable log
#            coi_distribution.pdf            stacked bar of strata
#            gambia_2014_qc_monoclonal.vcf.gz       VCF subset, monoclonal
#            gambia_2014_qc_polyclonal_low.vcf.gz   VCF subset, COI 2-3
#            gambia_2014_qc_polyclonal_high.vcf.gz  VCF subset, COI >= 4 (if any)
# =============================================================================

suppressPackageStartupMessages({
   library(tidyverse)
   library(readxl)
})

# ----- Paths ----------------------------------------------------------------
vcf_qc   <- "results/phase1_01/gambia_2014_qc.vcf.gz"
keep_in  <- "results/phase1_01/samples_keep.txt"
meta_in  <- "data/metadata/GamMetadata_2014.xlsx"

out_dir  <- "results/phase1_02"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Log helper -----------------------------------------------------------
log_lines <- c()
log_msg <- function(...) {
   m <- paste0(...); log_lines <<- c(log_lines, m); message(m)
}

# ============================================================================
# 1. Load QC-passed samples and metadata
# ============================================================================
log_msg("=== Phase 1 Step 1.2: COI stratification ===")

qc_samples <- readLines(keep_in)
log_msg("QC-passed samples: ", length(qc_samples))

meta <- read_xlsx(meta_in) %>%
   filter(Year == 2014, SampleID %in% qc_samples)
log_msg("2014 metadata rows for QC-passed samples: ", nrow(meta))

missing <- setdiff(qc_samples, meta$SampleID)
if (length(missing) > 0) {
   log_msg("WARNING — QC samples without 2014 metadata: ", paste(missing, collapse=","))
}

# ============================================================================
# 2. Stratify samples
# ============================================================================
# Note: in this cohort COIL/McCOIL are stored as character "1"/"2".
# Coerce to integer; downstream stratification supports COI >= 4 even though
# no such samples are present here — keeps the pipeline general.

stratify <- function(coi) {
   case_when(
      is.na(coi)   ~ NA_character_,
      coi == 1     ~ "monoclonal",
      coi >= 2 & coi <= 3 ~ "polyclonal_low",
      coi >= 4     ~ "polyclonal_high",
      TRUE         ~ NA_character_)
}

strat <- meta %>%
   transmute(SampleID,
             COIL_n   = suppressWarnings(as.integer(COIL)),
             McCOIL_n = suppressWarnings(as.integer(McCOIL)),
             stratum_McCOIL = stratify(McCOIL_n),
             stratum_COIL   = stratify(COIL_n))

write_tsv(strat, file.path(out_dir, "sample_coi_stratification.tsv"))

# ============================================================================
# 3. Distributions and concordance
# ============================================================================
log_msg("\n-- McCOIL (primary) strata --")
tab_mc <- table(strat$stratum_McCOIL, useNA = "ifany")
for (nm in names(tab_mc)) log_msg(sprintf("  %-18s %d", nm, tab_mc[[nm]]))

log_msg("\n-- COIL (sensitivity) strata --")
tab_cl <- table(strat$stratum_COIL, useNA = "ifany")
for (nm in names(tab_cl)) log_msg(sprintf("  %-18s %d", nm, tab_cl[[nm]]))

concord <- strat %>%
   count(stratum_McCOIL, stratum_COIL, name = "n") %>%
   arrange(desc(n))
write_tsv(concord, file.path(out_dir, "coi_concordance.tsv"))

n_agree <- sum(strat$stratum_McCOIL == strat$stratum_COIL, na.rm = TRUE)
n_compared <- sum(!is.na(strat$stratum_McCOIL) & !is.na(strat$stratum_COIL))
log_msg(sprintf("\nCOIL/McCOIL agreement: %d / %d (%.1f%%)",
                n_agree, n_compared, 100 * n_agree / n_compared))

# ============================================================================
# 4. Write per-stratum keep-lists and subset VCFs (primary = McCOIL)
# ============================================================================
strata <- c("monoclonal", "polyclonal_low", "polyclonal_high")
strat_files <- c(monoclonal      = file.path(out_dir, "gambia_2014_qc_monoclonal.vcf.gz"),
                 polyclonal_low  = file.path(out_dir, "gambia_2014_qc_polyclonal_low.vcf.gz"),
                 polyclonal_high = file.path(out_dir, "gambia_2014_qc_polyclonal_high.vcf.gz"))

log_msg("\nWriting per-stratum VCF subsets (primary = McCOIL):")
for (s in strata) {
   ids <- strat$SampleID[which(strat$stratum_McCOIL == s)]
   keep <- file.path(out_dir, sprintf("samples_%s.txt", s))
   writeLines(ids, keep)
   if (length(ids) == 0) {
      log_msg(sprintf("  %-18s n=0  (no VCF written)", s))
      next
   }
   system(sprintf(
      "bcftools view -S %s --force-samples -Oz %s | bcftools view -q 0.01:minor -Oz -o %s - && bcftools index -t -f %s",
      keep, vcf_qc, strat_files[[s]], strat_files[[s]]))
   nsites <- as.integer(system(sprintf("bcftools view -H %s | wc -l",
                                       strat_files[[s]]), intern = TRUE))
   log_msg(sprintf("  %-18s n=%d samples, %d sites after re-MAF",
                   s, length(ids), nsites))
}

# Samples with NA McCOIL: flag separately
na_ids <- strat$SampleID[is.na(strat$stratum_McCOIL)]
if (length(na_ids) > 0) {
   writeLines(na_ids, file.path(out_dir, "samples_coi_NA.txt"))
   log_msg(sprintf("  %-18s n=%d (no VCF; flagged for review)", "NA_McCOIL", length(na_ids)))
}

# ============================================================================
# 5. Diagnostic plot
# ============================================================================
plot_df <- bind_rows(
   strat %>% transmute(method = "McCOIL", stratum = stratum_McCOIL),
   strat %>% transmute(method = "COIL",   stratum = stratum_COIL)) %>%
   mutate(stratum = factor(replace_na(stratum, "NA"),
                           levels = c("monoclonal", "polyclonal_low",
                                      "polyclonal_high", "NA")))

p <- ggplot(plot_df, aes(x = method, fill = stratum)) +
   geom_bar(position = "stack", colour = "white") +
   scale_fill_manual(values = c(monoclonal      = "#2c7fb8",
                                polyclonal_low  = "#f7b733",
                                polyclonal_high = "#d7191c",
                                `NA`            = "grey70")) +
   labs(title = "COI stratification — McCOIL (primary) vs COIL (sensitivity)",
        subtitle = sprintf("QC-passed samples (n=%d)", nrow(strat)),
        x = NULL, y = "Samples", fill = "Stratum") +
   theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "coi_distribution.pdf"), p,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 6. Save summary
# ============================================================================
writeLines(log_lines, file.path(out_dir, "coi_summary.txt"))
message("\nDone. Outputs in ", out_dir)
