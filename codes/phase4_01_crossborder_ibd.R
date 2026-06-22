# =============================================================================
# Phase 4 — Step 4.1: Cross-border importation (Gambia vs Senegal)
# =============================================================================
# Proposal §6 Phase 4.1:
#   - Compute pairwise IBD between Gambian and Senegalese isolates.
#   - Identify cross-border highly related pairs (IBD > 0.25).
#   - Characterise origin within Senegal (border vs interior proxy),
#     seasonal pattern, and lifecycle-window-compatible directionality.
#
# Inputs : data/senegam_pf.vcf
#          data/metadata/senegam_pf.txt
#          results/phase1_03/sample_master_metadata.rds (optional, for season/date)
#
# Outputs: results/phase4_01/
#          senegam_core_biSNP_maf.vcf.gz
#          hmmibd_input.txt
#          hmmibd.hmm.txt
#          hmmibd.hmm_fract.txt
#          crossborder_pair_ibd.tsv
#          crossborder_high_ibd_pairs.tsv
#          senegal_origin_summary.tsv
#          crossborder_seasonality.tsv
#          crossborder_directionality_lifecycle.tsv
#          phase4_01_crossborder_ibd_distribution.pdf
#          phase4_01_senegal_origin_barplot.pdf
#          phase4_01_seasonality_barplot.pdf
#          phase4_01_summary.txt
#          phase4_01_results.rds
#
# Notes:
#   - This script uses hmmIBD on the combined Senegal+Gambia VCF.
#   - Lifecycle directionality requires exact dates for both samples.
#     The provided senegam metadata has Year-level timing; therefore
#     directional assignment may be mostly/entirely unassigned unless richer
#     Senegal dates are supplied in future metadata.
# =============================================================================
rm(list = ls())

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir    <- .parse_arg("output_dir",    "results/phase4_01")
phase1_meta_arg <- .parse_arg("phase1_meta", "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
})

# ----- Paths ----------------------------------------------------------------
vcf_in         <- "data/senegam_pf.vcf"
meta_in        <- "data/metadata/senegam_pf.txt"
phase1_meta_in <- phase1_meta_arg

out_dir        <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bcftools_bin   <- Sys.which("bcftools")
hmmibd_bin     <- Sys.getenv("HMMIBD", "tools/hmmIBD/hmmIBD")

vcf_filt       <- file.path(out_dir, "senegam_core_biSNP_maf.vcf.gz")
hmm_input      <- file.path(out_dir, "hmmibd_input.txt")
out_prefix     <- file.path(out_dir, "hmmibd")

# ----- Parameters -----------------------------------------------------------
MIN_MAF            <- 0.01
IBD_HIGH_THRESHOLD <- 0.25
MAX_ITER           <- 5
LIFECYCLE_MIN_DAYS <- 32
LIFECYCLE_MAX_DAYS <- 67
SEED               <- 42L
set.seed(SEED)

# Border-zone proxy for Senegal locations (proposal mentions border emphasis).
BORDER_LOCATIONS <- c("Bounkiling", "Kolda", "Ziguinchor", "Tambacounda")

log_lines <- c()
log_msg <- function(...) {
   m <- paste0(...)
   log_lines <<- c(log_lines, m)
   message(m)
}

run_cmd <- function(cmd) {
   log_msg("$ ", cmd)
   status <- system(cmd)
   if (status != 0)
      stop("Command failed with exit code ", status, ": ", cmd)
}

empty_tbl <- function(...) tibble(...)

# ============================================================================
# 1. Preconditions + metadata
# ============================================================================
log_msg("=== Phase 4.1: Cross-border IBD (Gambia vs Senegal) ===")
stopifnot(file.exists(vcf_in))
stopifnot(file.exists(meta_in))
stopifnot(file.exists(hmmibd_bin))
if (!nzchar(bcftools_bin))
   stop("bcftools not found on PATH.")

meta <- read_tsv(meta_in, show_col_types = FALSE) %>%
   distinct(SampleID, .keep_all = TRUE) %>%
   mutate(
      Country = str_to_title(Country),
      country_lc = str_to_lower(Country),
      Year = suppressWarnings(as.integer(Year)),
      region = case_when(
         str_to_lower(Region) %in% c("center", "central") ~ "central",
         str_to_lower(Region) == "west"                   ~ "west",
         str_to_lower(Region) == "east"                   ~ "east",
         str_to_lower(Region) == "sen"                    ~ "sen",
         TRUE                                             ~ str_to_lower(Region)
      ),
      senegal_zone = case_when(
         country_lc != "senegal" ~ NA_character_,
         str_to_lower(Location) %in% str_to_lower(BORDER_LOCATIONS) ~ "border",
         TRUE ~ "interior"
      )
   )

log_msg(sprintf("Metadata rows (deduplicated): %d", nrow(meta)))
log_msg("Metadata country counts:")
log_msg(paste(capture.output(print(table(meta$Country))), collapse = "\n"))

# Optional Phase 1 metadata to recover exact dates/seasons for 2014 Gambia set
if (file.exists(phase1_meta_in)) {
   p1 <- readRDS(phase1_meta_in) %>%
      transmute(SampleID,
                gambia_visit_date = as.Date(VisitDate),
                gambia_transmission_season = as.character(transmission_season)) %>%
      distinct(SampleID, .keep_all = TRUE)
   log_msg(sprintf("Loaded phase1 metadata for seasonal/date joins: %d rows", nrow(p1)))
} else {
   p1 <- tibble(SampleID = character(0),
                gambia_visit_date = as.Date(character(0)),
                gambia_transmission_season = character(0))
   log_msg("Phase1 metadata not found; season/date enrichment disabled.")
}

# ============================================================================
# 2. Filter VCF to core biallelic SNPs + MAF
# ============================================================================
cmd_filter <- paste0(
   shQuote(bcftools_bin),
   " view -f PASS -i 'INFO/RegionType==\"Core\"' -m2 -M2 -v snps ",
   shQuote(vcf_in),
   " -Ou | ",
   shQuote(bcftools_bin),
   " view -q ", MIN_MAF, ":minor -Oz -o ", shQuote(vcf_filt), " -"
)
run_cmd(cmd_filter)
run_cmd(paste(shQuote(bcftools_bin), "index -t -f", shQuote(vcf_filt)))

vcf_samples <- system2(bcftools_bin, c("query", "-l", vcf_filt), stdout = TRUE)
log_msg(sprintf("Samples in filtered VCF: %d", length(vcf_samples)))

meta_vcf <- meta %>% filter(SampleID %in% vcf_samples)
missing_meta <- setdiff(vcf_samples, meta_vcf$SampleID)
if (length(missing_meta) > 0) {
   writeLines(sort(missing_meta), file.path(out_dir, "samples_missing_metadata.txt"))
   log_msg(sprintf("WARNING: %d VCF samples missing metadata (written to samples_missing_metadata.txt)",
                   length(missing_meta)))
}

country_tab <- meta_vcf %>% count(Country, name = "n") %>% arrange(desc(n))
log_msg("VCF+metadata country counts:")
log_msg(paste(capture.output(print(country_tab)), collapse = "\n"))

if (!all(c("Gambia", "Senegal") %in% country_tab$Country)) {
   stop("Expected both Gambia and Senegal samples in filtered VCF+metadata join.")
}

# ============================================================================
# 3. Build hmmIBD input and run hmmIBD
# ============================================================================
# hmmIBD requires a header row: chrom \t pos \t sample1 \t sample2 \t ...
# `bcftools query -l` returns samples in the same order as their genotype
# columns in `bcftools query -f '...[\t%GT]'`, so we use it for the header.
header_line <- paste(c("chrom", "pos", vcf_samples), collapse = "\t")
writeLines(header_line, hmm_input)

cmd_hmm_input <- paste(
   shQuote(bcftools_bin), "query -f '%CHROM\\t%POS[\\t%GT]\\n'", shQuote(vcf_filt),
   "| awk 'BEGIN{OFS=\"\\t\"}",
   "{if ($1 ~ /^Pf3D7_[0-9][0-9]_v3$/) chr = substr($1, 7, 2) + 0; else next;",
   "printf \"%d\\t%s\", chr, $2;",
   "for (i=3;i<=NF;i++){gt=$i; val=-1;",
   "if (gt==\"0/0\" || gt==\"0|0\" || gt==\"0\") val=0;",
   "else if (gt==\"1/1\" || gt==\"1|1\" || gt==\"1\") val=1;",
   "printf \"\\t%s\", val}",
   "printf \"\\n\"}' >>",
   shQuote(hmm_input)
)
run_cmd(cmd_hmm_input)

log_msg("Running hmmIBD ...")
hmm_log <- system2(hmmibd_bin,
                   args = c("-i", hmm_input, "-o", out_prefix, "-m", MAX_ITER),
                   stdout = TRUE, stderr = TRUE)
writeLines(hmm_log, file.path(out_dir, "hmmibd_run.log"))
status <- attr(hmm_log, "status")
if (!is.null(status) && status != 0) {
   stop("hmmIBD failed with exit code ", status,
        ". See results/phase4_01/hmmibd_run.log")
}

fract_path <- paste0(out_prefix, ".hmm_fract.txt")
if (!file.exists(fract_path))
   stop("hmmIBD output not found: ", fract_path)

# ============================================================================
# 4. Parse pairwise IBD and isolate cross-border pairs
# ============================================================================
pairs <- read_tsv(fract_path,
                  col_types = cols(sample1 = "c", sample2 = "c",
                                   .default = col_guess())) %>%
   rename(SampleID1 = sample1, SampleID2 = sample2)

meta_min <- meta_vcf %>%
   select(SampleID, Location, Country, Year, region, senegal_zone)

pairs <- pairs %>%
   left_join(meta_min %>%
                rename(SampleID1 = SampleID,
                       Location1 = Location, Country1 = Country, Year1 = Year,
                       region1 = region, senegal_zone1 = senegal_zone),
             by = "SampleID1") %>%
   left_join(meta_min %>%
                rename(SampleID2 = SampleID,
                       Location2 = Location, Country2 = Country, Year2 = Year,
                       region2 = region, senegal_zone2 = senegal_zone),
             by = "SampleID2") %>%
   mutate(
      country1_lc = str_to_lower(Country1),
      country2_lc = str_to_lower(Country2),
      cross_border =
         (country1_lc == "gambia"  & country2_lc == "senegal") |
            (country1_lc == "senegal" & country2_lc == "gambia")
   )

cross <- pairs %>%
   filter(cross_border) %>%
   transmute(
      SampleID1, SampleID2, fract_sites_IBD,
      gambia_id       = if_else(country1_lc == "gambia",  SampleID1, SampleID2),
      senegal_id      = if_else(country1_lc == "senegal", SampleID1, SampleID2),
      gambia_location = if_else(country1_lc == "gambia",  Location1, Location2),
      senegal_location = if_else(country1_lc == "senegal", Location1, Location2),
      gambia_region   = if_else(country1_lc == "gambia",  region1, region2),
      senegal_zone    = if_else(country1_lc == "senegal", senegal_zone1, senegal_zone2),
      gambia_year     = if_else(country1_lc == "gambia",  Year1, Year2),
      senegal_year    = if_else(country1_lc == "senegal", Year1, Year2)
   ) %>%
   left_join(p1, by = c("gambia_id" = "SampleID")) %>%
   mutate(
      senegal_visit_date = as.Date(NA_character_),  # placeholder unless richer metadata is added
      gambia_transmission_season = replace_na(gambia_transmission_season, "unknown")
   ) %>%
   arrange(desc(fract_sites_IBD))

high <- cross %>%
   filter(fract_sites_IBD >= IBD_HIGH_THRESHOLD) %>%
   mutate(
      has_pair_dates = !is.na(gambia_visit_date) & !is.na(senegal_visit_date),
      dt_days = if_else(has_pair_dates,
                        as.numeric(abs(difftime(gambia_visit_date, senegal_visit_date,
                                                units = "days"))),
                        NA_real_),
      in_lifecycle_window = has_pair_dates &
         dt_days >= LIFECYCLE_MIN_DAYS & dt_days <= LIFECYCLE_MAX_DAYS,
      inferred_direction = case_when(
         !has_pair_dates ~ "unassigned_missing_dates",
         senegal_visit_date < gambia_visit_date ~ "Senegal -> Gambia",
         senegal_visit_date > gambia_visit_date ~ "Gambia -> Senegal",
         TRUE ~ "same_day"
      )
   )

write_tsv(cross, file.path(out_dir, "crossborder_pair_ibd.tsv"))
write_tsv(high,  file.path(out_dir, "crossborder_high_ibd_pairs.tsv"))

log_msg(sprintf("Cross-border pairs: %d", nrow(cross)))
log_msg(sprintf("High cross-border pairs (IBD >= %.2f): %d",
                IBD_HIGH_THRESHOLD, nrow(high)))

# ============================================================================
# 5. Summaries for origin, seasonality, directionality
# ============================================================================
if (nrow(high) > 0) {
   origin_summary <- high %>%
      count(senegal_location, senegal_zone, name = "n_high_pairs") %>%
      arrange(desc(n_high_pairs))

   season_summary <- high %>%
      count(gambia_transmission_season, name = "n_high_pairs") %>%
      arrange(desc(n_high_pairs))

   direction_summary <- high %>%
      count(inferred_direction, name = "n_high_pairs") %>%
      arrange(desc(n_high_pairs))
} else {
   origin_summary <- empty_tbl(
      senegal_location = character(0),
      senegal_zone = character(0),
      n_high_pairs = integer(0))
   season_summary <- empty_tbl(
      gambia_transmission_season = character(0),
      n_high_pairs = integer(0))
   direction_summary <- empty_tbl(
      inferred_direction = character(0),
      n_high_pairs = integer(0))
}

write_tsv(origin_summary, file.path(out_dir, "senegal_origin_summary.tsv"))
write_tsv(season_summary, file.path(out_dir, "crossborder_seasonality.tsv"))
write_tsv(direction_summary, file.path(out_dir, "crossborder_directionality_lifecycle.tsv"))

log_msg("Senegal origin among high-IBD cross-border pairs:")
log_msg(paste(capture.output(print(origin_summary)), collapse = "\n"))
log_msg("Seasonality among high-IBD cross-border pairs (Gambia sample season):")
log_msg(paste(capture.output(print(season_summary)), collapse = "\n"))
log_msg("Lifecycle directionality assignment summary:")
log_msg(paste(capture.output(print(direction_summary)), collapse = "\n"))

# ============================================================================
# 6. Plots
# ============================================================================
p_ibd <- ggplot(cross, aes(x = fract_sites_IBD)) +
   geom_histogram(bins = 80, fill = "#2E86AB", colour = "white") +
   geom_vline(xintercept = IBD_HIGH_THRESHOLD, linetype = "dashed", colour = "red") +
   labs(title = "Cross-border pairwise IBD (Gambia–Senegal)",
        subtitle = sprintf("High-related threshold = %.2f", IBD_HIGH_THRESHOLD),
        x = "fract_sites_IBD (hmmIBD)", y = "Number of cross-border pairs") +
   theme_minimal()
ggsave(file.path(out_dir, "phase4_01_crossborder_ibd_distribution.pdf"),
       p_ibd, width = 8, height = 5, dpi = 600)

if (nrow(origin_summary) > 0) {
   p_origin <- ggplot(origin_summary,
                      aes(x = reorder(senegal_location, n_high_pairs),
                          y = n_high_pairs, fill = senegal_zone)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c(border = "#D7263D", interior = "#2E86AB"),
                        na.value = "grey70") +
      labs(title = "Origin of high-IBD cross-border pairs in Senegal",
           x = "Senegal location", y = "High-IBD pairs", fill = "Zone proxy") +
      theme_minimal()
} else {
   p_origin <- ggplot() +
      annotate("text", x = 1, y = 1,
               label = "No high-IBD cross-border pairs at current threshold") +
      theme_void()
}
ggsave(file.path(out_dir, "phase4_01_senegal_origin_barplot.pdf"),
       p_origin, width = 8, height = 5, dpi = 600)

if (nrow(season_summary) > 0) {
   p_season <- ggplot(season_summary,
                      aes(x = gambia_transmission_season, y = n_high_pairs)) +
      geom_col(fill = "#F46036") +
      labs(title = "Seasonality of high-IBD cross-border pairs",
           subtitle = "Season assigned from Gambian sample metadata where available",
           x = "Gambia transmission season", y = "High-IBD pairs") +
      theme_minimal()
} else {
   p_season <- ggplot() +
      annotate("text", x = 1, y = 1,
               label = "No high-IBD cross-border pairs at current threshold") +
      theme_void()
}
ggsave(file.path(out_dir, "phase4_01_seasonality_barplot.pdf"),
       p_season, width = 8, height = 5, dpi = 600)

# ============================================================================
# 7. Persist + summary
# ============================================================================
saveRDS(
   list(
      params = list(
         MIN_MAF = MIN_MAF,
         IBD_HIGH_THRESHOLD = IBD_HIGH_THRESHOLD,
         LIFECYCLE_MIN_DAYS = LIFECYCLE_MIN_DAYS,
         LIFECYCLE_MAX_DAYS = LIFECYCLE_MAX_DAYS,
         MAX_ITER = MAX_ITER
      ),
      metadata_counts = country_tab,
      cross_pairs = cross,
      high_pairs = high,
      origin_summary = origin_summary,
      season_summary = season_summary,
      direction_summary = direction_summary
   ),
   file.path(out_dir, "phase4_01_results.rds")
)

log_msg(sprintf("Pairs with assignable lifecycle direction (both exact dates): %d",
                sum(high$has_pair_dates, na.rm = TRUE)))
if (sum(high$has_pair_dates, na.rm = TRUE) == 0) {
   log_msg("NOTE: Senegal sample dates are absent in current metadata; lifecycle")
   log_msg("      directionality is therefore unassigned for high cross-border pairs.")
}

writeLines(log_lines, file.path(out_dir, "phase4_01_summary.txt"))
message("\nDone. Outputs in ", out_dir)
