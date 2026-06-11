# =============================================================================
# Phase 2 — Step 2.1: Pairwise FST and isolation-by-distance
# =============================================================================
# Per notes/objective3_proposal.md §6 Phase 2 Step 2.1:
#   - Compute pairwise FST (Weir-Cockerham) between all village pairs
#   - Construct an FST distance matrix and visualise with hierarchical
#     clustering
#   - Test for isolation-by-distance (IBD) with a Mantel test
#     (genetic distance vs geographic distance, 10,000 permutations)
#
# Sample set : monoclonal subset (results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz)
#              FST assumes unambiguous diploid/haploid genotype calls;
#              polyclonal samples violate that assumption and are excluded
#              here. They re-enter for IBD-based phases via dcifer.
#
# Grouping  : three nested levels run in parallel
#   - village  : 10 VillageCodes (small ones flagged: n < MIN_VILLAGE_N)
#   - location : 6 named villages aggregating adjacent VillageCodes
#   - region   : 3 regions (west / central / east)
#
# Output : results/phase2_01/
#            fst_matrix_<level>.tsv           pairwise FST (Weir-Cockerham)
#            fst_long_<level>.tsv             long form with geographic distance
#            mantel_results.tsv               Mantel test summaries
#            fst_heatmap_<level>.pdf          heatmap + hierarchical clustering
#            fst_vs_geodist_<level>.pdf       IBD scatter with Mantel stats
#            fst_dendrogram_village.pdf       Ward dendrogram of village FST
#            fst_summary.txt                  human-readable log
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase2_01")
mono_vcf   <- .parse_arg("mono_vcf",   "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
meta_arg   <- .parse_arg("meta",       "results/phase1_03/sample_master_metadata.rds")

# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(vcfR)
   library(adegenet)
   library(hierfstat)
   library(vegan)
   library(geosphere)
   library(pheatmap)
})

# ----- Paths ----------------------------------------------------------------
vcf_in   <- mono_vcf
meta_in  <- meta_arg

out_dir  <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Parameters -----------------------------------------------------------
N_PERM         <- 10000      # Mantel permutations (proposal)
MIN_VILLAGE_N  <- 5          # flag villages with fewer samples
SEED           <- 42
set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load VCF + metadata, align
# ============================================================================
log_msg("=== Phase 2 Step 2.1: FST + IBD Mantel ===")
vcf  <- read.vcfR(vcf_in, verbose = FALSE)
meta <- readRDS(meta_in)
vcf_samples <- colnames(vcf@gt)[-1]

log_msg(sprintf("Monoclonal VCF samples: %d", length(vcf_samples)))
log_msg(sprintf("Sites: %d", nrow(vcf@fix)))

meta_aln <- meta %>%
   filter(SampleID %in% vcf_samples) %>%
   arrange(match(SampleID, vcf_samples))

stopifnot(all(meta_aln$SampleID == vcf_samples))
log_msg(sprintf("Metadata aligned: %d samples", nrow(meta_aln)))

# ============================================================================
# 2. Convert to genind once; reuse for each grouping
# ============================================================================
log_msg("Converting VCF to genind ...")
gi <- vcfR2genind(vcf)
indNames(gi) <- vcf_samples
log_msg(sprintf("genind: %d individuals, %d loci", nInd(gi), nLoc(gi)))

# ============================================================================
# 3. Helper: compute pairwise FST + village centroids + Mantel test
# ============================================================================
compute_fst_and_ibd <- function(level_label) {
   log_msg(sprintf("\n--- Level: %s ---", level_label))

   # Assign populations
   pop_vec <- meta_aln[[ if (level_label == "village")  "VillageCode"
                         else if (level_label == "location") "Location"
                         else "region" ]]
   pop_vec <- as.character(pop_vec)

   gi_lvl <- gi
   pop(gi_lvl) <- factor(pop_vec)

   # Group sizes
   ns <- table(pop_vec)
   log_msg(paste0("Group sizes: ",
                  paste(sprintf("%s=%d", names(ns), as.integer(ns)),
                        collapse = ", ")))

   if (level_label == "village") {
      small <- names(ns)[ns < MIN_VILLAGE_N]
      if (length(small) > 0)
         log_msg(sprintf("  Note: groups with n < %d (FST noisy): %s",
                         MIN_VILLAGE_N, paste(small, collapse = ", ")))
   }

   # Pairwise Weir-Cockerham FST
   log_msg("Computing pairwise Weir-Cockerham FST ...")
   fst_mat <- genet.dist(gi_lvl, method = "WC84")     # returns dist
   fst_mat <- as.matrix(fst_mat)
   pops <- rownames(fst_mat)

   write_tsv(as_tibble(fst_mat, rownames = "pop"),
             file.path(out_dir, sprintf("fst_matrix_%s.tsv", level_label)))

   # Group centroids (lon/lat) for geographic distance
   cent <- meta_aln %>%
      mutate(grp = pop_vec) %>%
      group_by(grp) %>%
      summarise(lon = mean(longitude, na.rm = TRUE),
                lat = mean(latitude,  na.rm = TRUE), .groups = "drop") %>%
      filter(grp %in% pops) %>%
      arrange(match(grp, pops))
   stopifnot(all(cent$grp == pops))

   # Haversine distance (km)
   geo_km <- distm(cbind(cent$lon, cent$lat), fun = distHaversine) / 1000
   rownames(geo_km) <- colnames(geo_km) <- pops

   # Long form
   pairs <- t(combn(pops, 2))
   long_df <- tibble(
      pop1 = pairs[, 1],
      pop2 = pairs[, 2],
      fst       = fst_mat[pairs],
      geo_dist_km = geo_km[pairs])
   
   write_tsv(long_df, file.path(out_dir, sprintf("fst_long_%s.tsv", level_label)))

   # Mantel test (only meaningful if >=4 pops, else skip with note)
   mantel_res <- NULL
   if (length(pops) >= 4) {
      m <- mantel(as.dist(fst_mat), as.dist(geo_km),
                  permutations = N_PERM, method = "pearson")
      mantel_res <- tibble(
         level       = level_label,
         n_pops      = length(pops),
         n_pairs     = choose(length(pops), 2),
         mantel_r    = m$statistic,
         p_value     = m$signif,
         n_perm      = N_PERM)
      log_msg(sprintf("Mantel (FST vs geo, %d perms): r = %.3f, p = %.4f",
                      N_PERM, m$statistic, m$signif))
   } else {
      log_msg(sprintf("Mantel skipped (n_pops = %d < 4)", length(pops)))
   }

   # IBD scatter
   p_ibd <- ggplot(long_df, aes(x = geo_dist_km, y = fst)) +
      geom_point(size = 2.5, alpha = 0.85, colour = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
      labs(title = sprintf("Isolation-by-distance (%s)", level_label),
           subtitle = if (!is.null(mantel_res))
              sprintf("Mantel r = %.2f, p = %.3f (%d permutations)",
                      mantel_res$mantel_r, mantel_res$p_value, N_PERM) else NULL,
           x = "Geographic distance (km)",
           y = "Pairwise FST (Weir-Cockerham)") +
      theme_minimal()
   
   ggsave(file.path(out_dir, sprintf("fst_vs_geodist_%s.pdf", level_label)),
          p_ibd, width = 6, height = 4.5, dpi = 600)

   # FST heatmap with hierarchical clustering
   diag(fst_mat) <- 0
   fst_plot <- pmax(fst_mat, 0)   # negative FST values clipped at 0 for display
   pheatmap(fst_plot,
            clustering_distance_rows = as.dist(fst_plot),
            clustering_distance_cols = as.dist(fst_plot),
            clustering_method = "ward.D2",
            display_numbers = TRUE, number_format = "%.3f",
            main = sprintf("Pairwise FST — %s level (hierarchical clustering)",
                           level_label),
            filename = file.path(out_dir,
                                 sprintf("fst_heatmap_%s.pdf", level_label)),
            width = max(6, 0.6 * length(pops) + 2),
            height = max(5, 0.55 * length(pops) + 1.5))

   list(fst_mat = fst_mat, geo_km = geo_km,
        long = long_df, mantel = mantel_res, pops = pops)
}

# ============================================================================
# 4. Run for the three nested levels
# ============================================================================
res_village  <- compute_fst_and_ibd("village")
res_location <- compute_fst_and_ibd("location")
res_region   <- compute_fst_and_ibd("region")

# ============================================================================
# 5. Save combined Mantel summary
# ============================================================================
mantel_all <- bind_rows(res_village$mantel,
                        res_location$mantel,
                        res_region$mantel)
write_tsv(mantel_all, file.path(out_dir, "mantel_results.tsv"))

# ============================================================================
# 6. Village-level Ward dendrogram (separate figure for the proposal Fig 1C)
# ============================================================================
fm <- res_village$fst_mat
diag(fm) <- 0
hc <- hclust(as.dist(pmax(fm, 0)), method = "ward.D2")
pdf(file.path(out_dir, "fst_dendrogram_village.pdf"), width = 7, height = 4.5)
plot(hc, hang = -1, main = "Village-level FST dendrogram (Ward.D2)",
     xlab = "", sub = "", cex = 0.9)
dev.off()

# ============================================================================
# 7. Persist objects for downstream phases
# ============================================================================
saveRDS(list(village = res_village, location = res_location, region = res_region),
        file.path(out_dir, "fst_results.rds"))

writeLines(log_lines, file.path(out_dir, "fst_summary.txt"))
message("\nDone. Outputs in ", out_dir)
