# =============================================================================
# Phase 1 — Step 1.3: Regional and temporal sample annotation
# =============================================================================
# Builds the master metadata table that drives all downstream analyses.
# Each QC-passed, COI-stratified sample is annotated with:
#
#   - region (west / central / east)        derived from village longitude
#   - village name and VillageCode          from metadata
#   - household (HHCode) and compound       from metadata
#   - collection month (numeric)            from VisitDate
#   - transmission_season                   early / peak / late (Gambia 2014)
#   - intervention_era                      placeholder — single 2014 cohort
#   - COI stratum (McCOIL primary)          from Step 1.2
#
# Regions follow the natural longitude clustering of the 10 villages in this
# cohort (see results/phase1_03/village_centroids.tsv):
#   west    longitude <= -15.5   Besse, Chogen, Yallal Ba
#   central -15.5 < longitude <= -14.5   Sinchu, Dongoro Ba, Sare Seedy
#   east    longitude >  -14.5   SareWuro, Njayel, Madina Samako, Gunjur Kunta
#
# Seasons follow the Gambian malaria transmission calendar:
#   early   Jul–Aug      (pre-peak, start of wet season)
#   peak    Sep–Nov      (high-transmission window; Amambua-Ngwa 2019 peak)
#   late    Dec onward   (late wet / early dry transition)
#
# Input  : data/metadata/GamMetadata_2014.xlsx
#          results/phase1_01/samples_keep.txt
#          results/phase1_02/sample_coi_stratification.tsv
#
# Output : results/phase1_03/
#            sample_master_metadata.tsv   the canonical metadata table
#            sample_master_metadata.rds   same, as RDS for downstream scripts
#            region_counts.tsv            samples per region x season
#            village_centroids.tsv        village -> lon/lat/region map
#            sampling_map.pdf             geographic overview
#            sampling_calendar.pdf        month-by-region sampling figure
#            annotation_summary.txt       log
# =============================================================================

theme_fig <- function() {
   theme_bw() +
      theme(
         legend.title = element_text(size = 15, color = "black", face = "bold", hjust = 0.5),
         legend.text = element_text(size = 12, color = "black", face = "bold"),
         plot.title        = element_text(size = 15, color = "black", face = "bold"),
         axis.title        = element_text(size = 14, color = "black", face = "bold"),
         axis.text         = element_text(size = 12, color = "black"),
         axis.line         = element_line(linewidth = 1, colour = "black", lineend = "square"),
         axis.ticks        = element_line(color = "black", linewidth = 0.7),
         axis.ticks.length = unit(0.22, "cm")
      )
}

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
qc_dir     <- .parse_arg("qc_dir",     "results/phase1_01")
coi_dir    <- .parse_arg("coi_dir",    "results/phase1_02")
output_dir <- .parse_arg("output_dir", "results/phase1_03")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(readxl)
   library(lubridate)
})

# ----- Paths ----------------------------------------------------------------
meta_in    <- "data/metadata/GamMetadata_2014.xlsx"
keep_in    <- file.path(qc_dir,  "samples_keep.txt")
coi_in     <- file.path(coi_dir, "sample_coi_stratification.tsv")

out_dir    <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Cutoffs --------------------------------------------------------------
LON_WEST_MAX    <- -15.5
LON_CENTRAL_MAX <- -14.5
PEAK_MONTHS     <- c(9, 10, 11)
EARLY_MONTHS    <- c(6, 7, 8)
LATE_MONTHS     <- c(12, 1, 2)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load and join inputs
# ============================================================================
log_msg("=== Phase 1 Step 1.3: Sample annotation ===")

qc_samples <- readLines(keep_in)
coi <- read_tsv(coi_in, show_col_types = FALSE) %>%
   select(SampleID, COIL_n, McCOIL_n, stratum_McCOIL, stratum_COIL)

meta <- read_xlsx(meta_in) %>%
   filter(SampleID %in% qc_samples)

log_msg("QC-passed samples:  ", length(qc_samples))
log_msg("Metadata rows kept: ", nrow(meta))

# ============================================================================
# 2. Derive region, month, season, intervention era
# ============================================================================
master <- meta %>%
   left_join(coi, by = "SampleID") %>%
   mutate(
      Location = case_when(
         VillageCode == 'A' ~ 'Bessi (A)',      VillageCode == 'C' ~ 'Chogen (C)',
         VillageCode == 'D' ~ 'Yallal Ba (D)',  VillageCode == 'E' ~ 'Sinchu (E)',
         VillageCode == 'F' ~ 'Dongoro Ba (F)', VillageCode == 'G' ~ 'Sare Seedy (G)',
         VillageCode == 'J' ~ 'Njayel (J)',     VillageCode == 'K' ~ 'Madina Samako (K)',
         VillageCode == 'L' ~ 'Sare Wuro (L)',  VillageCode == 'M' ~ 'Gunjur Kunta (M)'
      ),
      region = case_when(
         longitude <= LON_WEST_MAX     ~ "Western",
         longitude <= LON_CENTRAL_MAX  ~ "Central",
         TRUE                          ~ "Eastern"),
      region = factor(region, levels = c("Western", "Central", "Eastern")),
      collection_month = month(VisitDate),
      collection_year  = year(VisitDate),
      transmission_season = case_when(
         collection_month %in% PEAK_MONTHS  ~ "peak",
         collection_month %in% EARLY_MONTHS ~ "early",
         collection_month %in% LATE_MONTHS  ~ "late",
         TRUE                               ~ NA_character_),
      transmission_season = factor(transmission_season,
                                   levels = c("early", "peak", "late")),
      intervention_era = "Gambia_2014_cohort") %>%
   select(SampleID, region, Location, VillageCode,
          HHCode, CompoundCode,
          VisitDate, collection_year, collection_month, transmission_season,
          latitude, longitude,
          COIL_n, McCOIL_n, stratum_McCOIL, stratum_COIL,
          intervention_era,
          GeoPop, YearPop, MatchID2)

# ============================================================================
# 3. Sanity checks
# ============================================================================
n_per_region <- master %>% count(region)
log_msg("\nSamples per region:")
for (i in seq_len(nrow(n_per_region)))
   log_msg(sprintf("  %-8s %d", as.character(n_per_region$region[i]), n_per_region$n[i]))

n_per_season <- master %>% count(transmission_season)
log_msg("\nSamples per transmission season:")
for (i in seq_len(nrow(n_per_season)))
   log_msg(sprintf("  %-8s %d", as.character(n_per_season$transmission_season[i]),
                   n_per_season$n[i]))

# Region x season
rs <- master %>% count(region, transmission_season) %>%
   pivot_wider(names_from = transmission_season, values_from = n, values_fill = 0)

write_tsv(rs, file.path(out_dir, "region_counts.tsv"))

log_msg("\nRegion x season cross-tab:")
print(rs); log_lines <<- c(log_lines, capture.output(print(rs)))

# Village centroids -> region mapping (for documentation)
vc <- master %>%
   group_by(Location, VillageCode, region) %>%
   summarise(lon = mean(longitude, na.rm = TRUE),
             lat = mean(latitude,  na.rm = TRUE),
             n   = n(), .groups = "drop") %>%
   arrange(lon)

write_tsv(vc, file.path(out_dir, "village_centroids.tsv"))

# Sanity: no village should straddle two regions
straddle <- master %>%
   group_by(VillageCode) %>%
   summarise(n_regions = n_distinct(region), .groups = "drop") %>%
   filter(n_regions > 1)

if (nrow(straddle) > 0) {
   log_msg("WARNING — villages spanning multiple regions: ",
           paste(straddle$VillageCode, collapse = ", "))
} else {
   log_msg("\nOK: every village assigned to exactly one region.")
}

# Sanity: HH/Compound counts
log_msg(sprintf("\nUnique households (HHCode): %d", n_distinct(master$HHCode)))
log_msg(sprintf("Unique compounds (CompoundCode): %d", n_distinct(master$CompoundCode)))
log_msg(sprintf("Unique villages: %d", n_distinct(master$VillageCode)))

# ============================================================================
# 4. Write master metadata
# ============================================================================
write_tsv(master, file.path(out_dir, "sample_master_metadata.tsv"))
saveRDS(master, file.path(out_dir, "sample_master_metadata.rds"))

log_msg(sprintf("\nMaster metadata written: %d samples x %d columns",
                nrow(master), ncol(master)))

# ============================================================================
# 5. Diagnostic figures
# ============================================================================

# --- 5a. Sampling map (lon/lat with region colour, size = N)
village_pts <- master %>%
   group_by(Location, VillageCode, region) %>%
   summarise(lon = mean(longitude, na.rm = TRUE),
             lat = mean(latitude,  na.rm = TRUE),
             N   = n(), .groups = "drop")

p_map <- ggplot(village_pts, aes(x = lon, y = lat,
                                 colour = region, size = N)) +
   geom_point(alpha = 0.85) +
   ggrepel::geom_text_repel(aes(label = paste0(Location, " (", VillageCode, ")")),
                            size = 3, max.overlaps = 20, show.legend = FALSE) +
   geom_vline(xintercept = c(LON_WEST_MAX, LON_CENTRAL_MAX),
              linetype = 2, colour = "grey40") +
   scale_colour_manual(values = c(Western = "#2c7fb8", Central = "#f7b733",
                                  Eastern = "#d7191c")) +
   labs(title = "Gambia 2014 cohort — sampling sites",
        # subtitle = "Vertical lines: region cutoffs (-15.5, -14.5 longitude)",
        x = "Longitude", y = "Latitude", colour = "Regions") +
   theme_minimal() + theme_fig()

if (!requireNamespace("ggrepel", quietly = TRUE)) {
   # fallback without ggrepel
   p_map <- ggplot(village_pts, aes(x = lon, y = lat, colour = region, size = N)) +
      geom_point(alpha = 0.85) +
      geom_text(aes(label = paste0(Location, " (", VillageCode, ")")),
                vjust = -1, size = 3) +
      geom_vline(xintercept = c(LON_WEST_MAX, LON_CENTRAL_MAX),
                 linetype = 2, colour = "grey50") +
      scale_colour_manual(values = c(Western = "#2c7fb8", Central = "#f7b733",
                                     Eastern = "#d7191c")) +
      labs(title = "Gambia 2014 cohort — sampling sites",
           x = "Longitude", y = "Latitude", colour = "Regions") +
      theme_fig()
}

print(p_map)

ggsave(file.path(out_dir, "sampling_map.pdf"), plot = p_map,
       width = 9, height = 5, dpi = 600)

# Add also a MAP as alternative or main figure

# =============================================
# --- 5b. Sampling calendar — month by region
# =============================================
p_cal <- ggplot(master, aes(x = factor(collection_month), fill = region)) +
   geom_bar(position = "stack") +
   scale_fill_manual(values = c(Western = "#2c7fb8", 
                                Central = "#f7b733",
                                Eastern = "#d7191c")) +
   labs(title = "Sampling calendar by region (2014)",
        x = "Collection month", y = "Samples", fill = "Regions") +
   theme_fig()

print(p_cal)

ggsave(file.path(out_dir, "sampling_calendar.pdf"), plot = p_cal,
       width = 7, height = 4, dpi = 600)

# ============================================================================
# 6. Save summary
# ============================================================================
writeLines(log_lines, file.path(out_dir, "annotation_summary.txt"))
message("\nDone. Outputs in ", out_dir)
