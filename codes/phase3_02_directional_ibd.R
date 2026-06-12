# =============================================================================
# Phase 3 — Step 3.2: Directional transmission inference (Amambua-Ngwa 2019,
# WGS adaptation)
# =============================================================================
# Implements the Amambua-Ngwa et al. 2019 lifecycle-constrained transmission
# model (notes/directional_connectivity_Amambua-Ngwamodel.md), with the
# WGS-specific substitutions called out in Step 6 of that spec:
#
#   - p-distance  --> hmmIBD fract_sites_IBD (Phase 2.3 output)
#   - simulated f(d | transmission) --> empirical density from co-resident
#     (within-household) pairs, which are presumed transmission-linked
#   - travel constraint --> SKIPPED (no travel data for this cohort);
#     reported in the writeup as a sensitivity caveat
#
# Model summary:
#   For each pair (i, j) with observed IBD x_ij and time gap dt_ij = |t_j - t_i|,
#   compute P(transmission | x, dt) = f_T(x, dt) / (f_T(x, dt) + f_B(x, dt))
#   where f_T is built from within-household pairs and f_B from all pairs.
#   Restrict to pairs with x >= IBD_THRESH and dt in the lifecycle-derived
#   directional window (lower = 97.5%-ile of simulated common-source dt_evo,
#   upper = 97.5%-ile of simulated direct dt_evo).
#   Apply bilinear (mutual-argmax) maximisation to extract non-overlapping
#   accepted transmission paths. Direction: earlier sample = source.
#   Classify accepted paths by region pair and compute asymmetry with
#   bootstrap 95% CI.
#
# Inputs : results/phase2_03/pair_ibd_classified.tsv  (10,296 pairs)
#          results/phase1_03/sample_master_metadata.rds  (VisitDate, HHCode, region)
#
# Outputs: results/phase3_02/
#          lifecycle_window.tsv          simulated dt_evo summary
#          transmission_pairs.tsv        all candidate pairs with P_trans
#          accepted_paths.tsv            mutual-argmax surviving paths
#          directional_asymmetry.tsv     counts + bootstrap CI per direction
#          phase3_02_lifecycle_sim.pdf   common-source vs direct dt_evo
#          phase3_02_bivariate.pdf       (IBD, dt) for HH / all / accepted
#          phase3_02_paths_map.pdf       transmission paths on geographic map
#          phase3_02_asymmetry.pdf       barplot of direction counts + CI
#          phase3_02_summary.txt         human-readable log
#          phase3_02_results.rds         persisted state
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase3_02")
pairs_arg  <- .parse_arg("pair_class",  "results/phase2_03/pair_ibd_classified.tsv")
meta_arg   <- .parse_arg("meta",        "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(MASS)        # kde2d
   library(ggrepel)
   library(scales)
})

# ----- Paths ----------------------------------------------------------------
pairs_in <- pairs_arg
meta_in  <- meta_arg
out_dir  <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Parameters -----------------------------------------------------------
IBD_THRESH    <- 0.10        # spec §6: candidate transmission pairs floor
N_SIM         <- 100000L     # lifecycle simulations
KDE_N_GRID    <- 200L        # KDE grid resolution
BOOT_REPS     <- 5000L
SEED          <- 42L
set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load IBD pairs (Phase 2.3) and metadata
# ============================================================================
log_msg("=== Phase 3.2: Directional transmission (Amambua-Ngwa, WGS-adapted) ===")
pairs <- read_tsv(pairs_in, show_col_types = FALSE)
meta  <- readRDS(meta_in) %>%
   dplyr::select(SampleID, region, Location, VillageCode, HHCode, CompoundCode,
                 VisitDate, latitude, longitude)

# Attach dates and household codes for both samples
pairs <- pairs %>%
   left_join(meta %>% dplyr::select(SampleID, VisitDate1 = VisitDate,
                                    HHCode1 = HHCode, CompCode1 = CompoundCode),
             by = c("SampleID1" = "SampleID")) %>%
   left_join(meta %>% dplyr::select(SampleID, VisitDate2 = VisitDate,
                                    HHCode2 = HHCode, CompCode2 = CompoundCode),
             by = c("SampleID2" = "SampleID")) %>%
   mutate(dt = as.numeric(abs(difftime(VisitDate2, VisitDate1, units = "days"))),
          same_HH       = HHCode1 == HHCode2,
          same_compound = CompCode1 == CompCode2,
          ibd = fract_sites_IBD)

log_msg(sprintf("Loaded %d pairs; date range: %s .. %s",
                nrow(pairs),
                format(min(c(pairs$VisitDate1, pairs$VisitDate2))),
                format(max(c(pairs$VisitDate1, pairs$VisitDate2)))))

log_msg(sprintf("Within-household pairs: %d", sum(pairs$same_HH, na.rm = TRUE)))
log_msg(sprintf("Within-compound pairs:  %d", sum(pairs$same_compound, na.rm = TRUE)))

# ============================================================================
# 2. Lifecycle simulation -> directional dt window
# ============================================================================
rtri <- function(n, a, b, c) {     # min=a, max=b, mode=c
   u  <- runif(n)
   Fc <- (c - a) / (b - a)
   ifelse(u < Fc,
          a + sqrt(u * (b - a) * (c - a)),
          b - sqrt((1 - u) * (b - a) * (b - c)))
}

T_i_dir  <- rtri(N_SIM, 5, 15, 10)            # bite -> detectable (recipient)
T_v_dir  <- rtri(N_SIM, 8, 15, 11.5)          # sporogony
T_v_life <- rtri(N_SIM, 18, 24, 21)           # vector lifespan
T_d_dir  <- T_v_life - T_v_dir                # remaining vector window
T_A_dir  <- runif(N_SIM, 0, 30)               # last-sampling lag, source
b_dir    <- runif(N_SIM, 0, T_A_dir)          # bite-to-sampling lag, source
T_B_dir  <- runif(N_SIM, 0, 30)               # last-sampling lag, recipient
dt_direct <- T_v_dir + T_d_dir + T_i_dir + T_B_dir - b_dir

T_i_A <- rtri(N_SIM, 5, 15, 10)
T_i_B <- rtri(N_SIM, 5, 15, 10)
T_v_cs   <- rtri(N_SIM, 8, 15, 11.5)
T_vlife2 <- rtri(N_SIM, 18, 24, 21)
T_d_cs   <- T_vlife2 - T_v_cs
T_A_cs   <- runif(N_SIM, 0, 30)
T_B_cs   <- runif(N_SIM, 0, 30)
dt_csource <- T_d_cs + T_i_B + T_B_cs - T_i_A + T_A_cs

# Spec window interpretation (Amambua-Ngwa 2019, paper Fig.):
#   Plausibility window = 2.5%..97.5% of *direct* dt_evo.
# We *also* compute the 97.5%-ile of common-source dt_evo so we can flag the
# common-source overlap as a caveat. With these parameters (T_v + T_d + T_i +
# T_B - b vs T_d + T_iB + T_B - T_iA + T_A) the two distributions sit on top
# of one another and the 97.5%-of-common-source-as-lower-bound recipe inverts
# the window. We default to the direct-only plausibility band; flag the
# common-source mass inside it as the caveat in the summary.
dt_lower <- quantile(dt_direct,  0.025)
dt_upper <- quantile(dt_direct,  0.975)
cs_97    <- quantile(dt_csource, 0.975)
overlap_frac <- mean(dt_csource >= dt_lower & dt_csource <= dt_upper)

log_msg(sprintf("Common-source dt_evo: median=%.1f, 97.5%%=%.1f",
                median(dt_csource), cs_97))
log_msg(sprintf("Direct dt_evo:        median=%.1f, 2.5%%=%.1f, 97.5%%=%.1f",
                median(dt_direct), dt_lower, dt_upper))
log_msg(sprintf("Directional window (2.5-97.5%% of direct): %.1f .. %.1f days",
                dt_lower, dt_upper))
log_msg(sprintf("Caveat: %.1f%% of common-source mass falls inside this window",
                100 * overlap_frac))

write_tsv(tibble(component = c("common_source", "direct"),
                 median = c(median(dt_csource), median(dt_direct)),
                 q025   = c(quantile(dt_csource, 0.025), quantile(dt_direct, 0.025)),
                 q975   = c(quantile(dt_csource, 0.975), quantile(dt_direct, 0.975))),
          file.path(out_dir, "lifecycle_window.tsv"))

# Lifecycle plot
p_life <- tibble(value = c(dt_csource, dt_direct),
                 source = rep(c("common_source", "direct"), each = N_SIM)) %>%
   ggplot(aes(x = value, fill = source)) +
   geom_density(alpha = 0.4) +
   geom_vline(xintercept = c(dt_lower, dt_upper), linetype = "dashed", colour = "red") +
   annotate("text", x = dt_lower, y = 0, label = sprintf("%.1f", dt_lower),
            vjust = -0.5, hjust = 1.05, colour = "red", size = 3) +
   annotate("text", x = dt_upper, y = 0, label = sprintf("%.1f", dt_upper),
            vjust = -0.5, hjust = -0.1, colour = "red", size = 3) +
   scale_fill_manual(values = c(common_source = "#F46036", direct = "#2E86AB")) +
   labs(title = sprintf("Lifecycle dt_evo simulation (n=%d each)", N_SIM),
        subtitle = sprintf("Directional window: %.1f .. %.1f days", dt_lower, dt_upper),
        x = "dt_evo (days)", y = "density") +
   theme_minimal() +
   theme(
      legend.title  = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
      legend.text   = element_text(size = 11, color = "black"),
      plot.title    = element_text(size = 14, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text     = element_text(size = 11, color = "black"),
      axis.title    = element_text(size = 13, color = "black", face = "bold"),
      axis.line     = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks    = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm")
   )

ggsave(file.path(out_dir, "phase3_02_lifecycle_sim.pdf"), p_life,
       width = 7, height = 4, dpi = 600)

# ============================================================================
# 3. Empirical bivariate densities
# ============================================================================
# Filter set: pairs with dt > 0 (drop same-day pairs from KDE)
pairs_kde <- pairs %>% filter(!is.na(dt), !is.na(ibd))
hh_set    <- pairs_kde %>% filter(same_HH)
log_msg(sprintf("KDE pairs: all=%d, within-household=%d", nrow(pairs_kde), nrow(hh_set)))

# Domain for KDE: IBD in [0,1], dt in [0, max observed]
dt_max <- max(pairs_kde$dt) + 1
fit_kde <- function(df) {
   if (nrow(df) < 5) return(NULL)
   MASS::kde2d(df$ibd, df$dt,
               n = KDE_N_GRID,
               lims = c(0, 1, 0, dt_max))
}
kde_T <- fit_kde(hh_set)         # transmission (within-household)
kde_B <- fit_kde(pairs_kde)      # background (all pairs)
if (is.null(kde_T))
   stop("Too few within-household pairs to fit transmission KDE.")

# Evaluate a KDE at arbitrary (x, y) by bilinear interpolation on its grid
eval_kde <- function(kde, x, y) {
   gx <- kde$x; gy <- kde$y; z <- kde$z
   ix <- findInterval(x, gx, all.inside = TRUE)
   iy <- findInterval(y, gy, all.inside = TRUE)
   x1 <- gx[ix]; x2 <- gx[ix + 1]
   y1 <- gy[iy]; y2 <- gy[iy + 1]
   tx <- (x - x1) / (x2 - x1); ty <- (y - y1) / (y2 - y1)
   z11 <- z[cbind(ix,     iy    )]
   z12 <- z[cbind(ix,     iy + 1)]
   z21 <- z[cbind(ix + 1, iy    )]
   z22 <- z[cbind(ix + 1, iy + 1)]
   (1 - tx) * ((1 - ty) * z11 + ty * z12) +
        tx  * ((1 - ty) * z21 + ty * z22)
}

# ============================================================================
# 4. Candidate pairs + transmission probability
# ============================================================================
cand <- pairs_kde %>%
   mutate(in_window = dt >= dt_lower & dt <= dt_upper,
          above_ibd = ibd >= IBD_THRESH,
          candidate = in_window & above_ibd)
log_msg(sprintf("Candidate pairs (IBD>=%.2f AND dt in [%.1f, %.1f]): %d",
                IBD_THRESH, dt_lower, dt_upper, sum(cand$candidate)))

if (sum(cand$candidate) == 0) {
   warning("No candidate pairs in directional window above IBD threshold.")
}

cand_pos <- cand %>% filter(candidate) %>%
   mutate(f_T  = pmax(eval_kde(kde_T, ibd, dt), 0),
          f_B  = pmax(eval_kde(kde_B, ibd, dt), 0),
          P_trans = ifelse(f_T + f_B > 0, f_T / (f_T + f_B), 0))

# Direction: earlier visit = source. Tie-break by SampleID (stable).
cand_pos <- cand_pos %>%
   mutate(source_id    = ifelse(VisitDate1 <= VisitDate2, SampleID1, SampleID2),
          recipient_id = ifelse(VisitDate1 <= VisitDate2, SampleID2, SampleID1),
          source_region    = ifelse(VisitDate1 <= VisitDate2, region1, region2),
          recipient_region = ifelse(VisitDate1 <= VisitDate2, region2, region1),
          source_loc    = ifelse(VisitDate1 <= VisitDate2, Location1, Location2),
          recipient_loc = ifelse(VisitDate1 <= VisitDate2, Location2, Location1),
          source_lon    = ifelse(VisitDate1 <= VisitDate2, longitude1, longitude2),
          source_lat    = ifelse(VisitDate1 <= VisitDate2, latitude1,  latitude2),
          recipient_lon = ifelse(VisitDate1 <= VisitDate2, longitude2, longitude1),
          recipient_lat = ifelse(VisitDate1 <= VisitDate2, latitude2,  latitude1))

write_tsv(cand_pos %>%
             dplyr::select(source_id, recipient_id, source_region, recipient_region,
                           source_loc, recipient_loc, dt, ibd, P_trans),
          file.path(out_dir, "transmission_pairs.tsv"))

# ============================================================================
# 5. Bilinear (mutual-argmax) maximisation -> accepted paths
# ============================================================================
# For each candidate source i, best recipient j*(i) = argmax_j M_ij
# For each candidate recipient j, best source  i*(j) = argmax_i M_ij
# Accept (i, j) iff j*(i) = j AND i*(j) = i.
build_M <- function(df) {
   sids <- sort(unique(c(df$source_id, df$recipient_id)))
   if (length(sids) == 0)
      return(matrix(0, 0, 0, dimnames = list(character(0), character(0))))
   M <- matrix(0, length(sids), length(sids), dimnames = list(sids, sids))
   for (k in seq_len(nrow(df))) {
      i <- df$source_id[k]; j <- df$recipient_id[k]
      M[i, j] <- max(M[i, j], df$P_trans[k])
   }
   M
}
empty_paths <- function() {
   tibble(source_id    = character(0),
          recipient_id = character(0),
          P_trans      = numeric(0))
}
mutual_argmax <- function(M) {
   if (nrow(M) == 0) return(empty_paths())
   row_best <- apply(M, 1, function(r) if (max(r) == 0) NA_integer_ else which.max(r))
   col_best <- apply(M, 2, function(c) if (max(c) == 0) NA_integer_ else which.max(c))
   accepted <- list()
   for (i in seq_len(nrow(M))) {
      j <- row_best[i]
      if (is.na(j)) next
      if (!is.na(col_best[j]) && col_best[j] == i) {
         accepted[[length(accepted) + 1]] <-
            tibble(source_id    = rownames(M)[i],
                   recipient_id = colnames(M)[j],
                   P_trans      = M[i, j])
      }
   }
   if (length(accepted) == 0) empty_paths() else bind_rows(accepted)
}

M <- build_M(cand_pos)
accepted_raw <- mutual_argmax(M)
if (nrow(accepted_raw) > 0) {
   accepted <- accepted_raw %>%
      left_join(cand_pos %>%
                   dplyr::select(source_id, recipient_id, source_region, recipient_region,
                                 source_loc, recipient_loc, dt, ibd,
                                 source_lon, source_lat, recipient_lon, recipient_lat),
                by = c("source_id", "recipient_id"))
} else {
   accepted <- accepted_raw
}

log_msg(sprintf("Accepted transmission paths (mutual-argmax): %d", nrow(accepted)))
write_tsv(accepted, file.path(out_dir, "accepted_paths.tsv"))

# ============================================================================
# 6. Directional asymmetry + bootstrap CI
# ============================================================================
classify_path <- function(sr, rr) {
   if (sr == rr) paste0("within-", sr)
   else paste(sr, "->", rr)
}
if (nrow(accepted) > 0) {
   accepted <- accepted %>%
      rowwise() %>%
      mutate(direction = classify_path(source_region, recipient_region)) %>%
      ungroup()

   dir_counts <- accepted %>% count(direction, name = "n") %>% arrange(desc(n))
   log_msg("Direction counts:")
   log_msg(paste(capture.output(print(dir_counts)), collapse = "\n"))

   # Bootstrap CI: resample accepted paths with replacement
   boot_mat <- replicate(BOOT_REPS, {
      idx <- sample.int(nrow(accepted), replace = TRUE)
      tbl <- table(accepted$direction[idx])
      tbl
   }, simplify = FALSE)
   all_dirs <- sort(unique(unlist(lapply(boot_mat, names))))
   boot_df <- do.call(rbind, lapply(boot_mat, function(t) {
      v <- setNames(rep(0L, length(all_dirs)), all_dirs)
      v[names(t)] <- as.integer(t); v
   }))
   ci <- apply(boot_df, 2, quantile, probs = c(0.025, 0.5, 0.975))
   asym_df <- tibble(direction = colnames(boot_df),
                     n         = dir_counts$n[match(colnames(boot_df), dir_counts$direction)],
                     median_boot = ci["50%", ],
                     lo95        = ci["2.5%", ],
                     hi95        = ci["97.5%", ]) %>%
      mutate(n = replace_na(n, 0L))
   write_tsv(asym_df, file.path(out_dir, "directional_asymmetry.tsv"))

   p_asym <- ggplot(asym_df, aes(x = reorder(direction, n), y = n)) +
      geom_col(fill = "#2E86AB") +
      geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.25) +
      coord_flip() +
      labs(title = "Inferred directional transmission paths",
           subtitle = sprintf("%d accepted paths; bars = 95%% bootstrap CI (n=%d reps)",
                              nrow(accepted), BOOT_REPS),
           x = NULL, y = "Number of paths") +
      theme_minimal() +
      theme(
         plot.title    = element_text(size = 14, color = "black", face = "bold"),
         plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
         axis.text     = element_text(size = 11, color = "black"),
         axis.title    = element_text(size = 13, color = "black", face = "bold"),
         axis.line     = element_line(linewidth = 1, colour = "black", lineend = "square"),
         axis.ticks    = element_line(color = "black", linewidth = 0.7),
         axis.ticks.length = unit(0.22, "cm")
      )
   
   ggsave(file.path(out_dir, "phase3_02_asymmetry.pdf"), p_asym,
          width = 6.5, height = 4, dpi = 600)
} else {
   asym_df <- tibble()
   log_msg("No accepted paths -- skipping asymmetry plot.")
}

# ============================================================================
# 7. Diagnostic plots
# ============================================================================
# 7a. Bivariate (IBD, dt): HH vs all, with candidate window + accepted overlay
pkde <- pairs_kde %>%
   mutate(class = if_else(same_HH, "within-household", "other"))

p_biv <- ggplot(pkde, aes(x = ibd, y = dt)) +
   geom_point(aes(colour = class), alpha = 0.5, size = 1) +
   geom_hline(yintercept = c(dt_lower, dt_upper), linetype = "dashed", colour = "red") +
   geom_vline(xintercept = IBD_THRESH, linetype = "dashed", colour = "red") +
   { if (nrow(accepted) > 0)
        geom_point(data = accepted, aes(x = ibd, y = dt),
                   shape = 21, colour = "black", fill = "yellow",
                   size = 2.5, stroke = 0.6) } +
   scale_colour_manual(values = c(`within-household` = "#D7263D", other = "grey60")) +
   labs(title = "IBD vs time-between-samples",
        subtitle = sprintf("Red lines: candidate region (IBD >= %.2f and dt in window); yellow = accepted paths",
                           IBD_THRESH),
        x = "fraction sites IBD (hmmIBD)", y = "Days between samples") +
   theme_minimal() +
   theme(
      legend.title  = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
      legend.text   = element_text(size = 11, color = "black"),
      plot.title    = element_text(size = 14, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text     = element_text(size = 11, color = "black"),
      axis.title    = element_text(size = 13, color = "black", face = "bold"),
      axis.line     = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks    = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm")
   )

ggsave(file.path(out_dir, "phase3_02_bivariate.pdf"), p_biv,
       width = 7.5, height = 5.5, dpi = 600)

# 7b. Geographic transmission paths map
if (nrow(accepted) > 0) {
   cent <- meta %>% group_by(Location, region) %>%
      summarise(lon = mean(longitude, na.rm = TRUE),
                lat = mean(latitude,  na.rm = TRUE),
                .groups = "drop")
   
   p_map <- ggplot() +
      geom_point(data = cent,
                 aes(x = lon, y = lat, colour = region), size = 5) +
      geom_segment(data = accepted,
                   aes(x = source_lon, y = source_lat,
                       xend = recipient_lon, yend = recipient_lat,
                       alpha = P_trans),
                   arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                   colour = "darkred") +
      # scale_colour_viridis_c() +  # viridis for continuous
      geom_text_repel(data = cent, 
                      aes(x = lon, y = lat, label = Location),
                      size = 3, max.overlaps = 30, color = "black") +
      scale_colour_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                     Eastern = "#2E86AB"), name = "Region") +
      scale_alpha_continuous(name = "P(transmission)", range = c(0.1, 1), ) +
      labs(title = sprintf("Inferred directional transmission paths (n=%d)",
                           nrow(accepted)),
           subtitle = "Earlier-sampled = source; arrowhead = recipient",
           x = "Longitude", y = "Latitude") +
      coord_quickmap() +
      theme_minimal() +
      theme(
         legend.title  = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
         legend.text   = element_text(size = 11, color = "black"),
         plot.title    = element_text(size = 14, color = "black", face = "bold"),
         plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
         axis.text     = element_text(size = 11, color = "black"),
         axis.title    = element_text(size = 13, color = "black", face = "bold"),
         axis.ticks    = element_line(color = "black", linewidth = 0.7),
         axis.ticks.length = unit(0.22, "cm")
      )
   
   ggsave(file.path(out_dir, "phase3_02_paths_map.pdf"), p_map,
          width = 9, height = 6, dpi = 600)
}

# ============================================================================
# 8. Persist + log
# ============================================================================
saveRDS(list(pairs       = cand,
             cand_pos    = cand_pos,
             accepted    = accepted,
             asymmetry   = asym_df,
             lifecycle   = list(dt_csource = dt_csource, dt_direct = dt_direct,
                                dt_lower = dt_lower, dt_upper = dt_upper),
             kde_T = kde_T, kde_B = kde_B),
        file.path(out_dir, "phase3_02_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase3_02_summary.txt"))
message("\nDone. Outputs in ", out_dir)
