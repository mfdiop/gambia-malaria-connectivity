# =============================================================================
# Phase 3 — Step 3.3: Genomic isolation analysis
# =============================================================================
# Proposal §6 Phase 3 Step 3.3:
#   If western Gambia is sustained by *local* transmission, the Western
#   parasite population should look genetically isolated — elevated within-
#   region IBD, lower nucleotide diversity (pi), slower LD decay than the
#   East. If instead Western is continuously re-seeded from elsewhere, its
#   diversity should be inconsistent with low transmission intensity (i.e.,
#   similar to or higher than East).
#
#   Three tests:
#     (1) Within-region IBD: density + Wilcoxon (West vs Central vs East)
#     (2) Per-region per-season nucleotide diversity (pi)
#     (3) Per-region LD decay (r^2 vs physical distance)
#
# Inputs : results/phase2_03/hmmibd_input.txt        (chrom, pos, 144 haploid {0,1,-1})
#          results/phase2_03/pair_ibd_classified.tsv (pair IBD + region annotations)
#          results/phase1_03/sample_master_metadata.rds
#
# Outputs: results/phase3_03/
#          ibd_within_region.{tsv,pdf}    within-region IBD comparison
#          pi_region_season.{tsv,pdf}     per-region per-season pi (bootstrap CI)
#          ld_decay.{tsv,pdf}             r^2 vs distance per region
#          phase3_03_summary.txt
#          phase3_03_results.rds
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir <- .parse_arg("output_dir", "results/phase3_03")
hap_arg    <- .parse_arg("hap",        "results/phase2_03/hmmibd_input.txt")
pairs_arg  <- .parse_arg("pair_class", "results/phase2_03/pair_ibd_classified.tsv")
meta_arg   <- .parse_arg("meta",       "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(scales)
})

# ----- Paths ----------------------------------------------------------------
hap_in   <- hap_arg
pairs_in <- pairs_arg
meta_in  <- meta_arg
out_dir  <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Parameters -----------------------------------------------------------
REGIONS      <- c("Western", "Central", "Eastern")
SEASONS      <- c("Early", "Peak", "Late")
LD_MAX_BP    <- 50000L
LD_BIN_BP    <- 500L          # distance bin width
LD_MIN_MAF   <- 0.05
PI_BOOT_REPS <- 10000L         # bootstrap SNPs for pi CIs
SEED         <- 42L
set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load haploid matrix + metadata
# ============================================================================
log_msg("=== Phase 3.3: Genomic isolation analysis ===")
hap_raw <- read_tsv(hap_in, show_col_types = FALSE)
chrom <- as.integer(hap_raw$chrom)
pos   <- as.integer(hap_raw$pos)
sids  <- colnames(hap_raw)[-(1:2)]
H     <- as.matrix(hap_raw[, -(1:2)])
storage.mode(H) <- "integer"
# H values: 0/1 = allele dose, -1 = missing
H[H == -1L] <- NA_integer_
L <- nrow(H); n <- ncol(H)
log_msg(sprintf("Loaded haploid matrix: %d SNPs x %d samples", L, n))

meta <- readRDS(meta_in) %>% 
   filter(SampleID %in% sids) %>%
   arrange(match(SampleID, sids)) %>% 
   mutate(
      transmission_season = recode(transmission_season, "early"= "Early", "peak" = "Peak", "late" = "Late"),
      transmission_season = factor(transmission_season, 
                                   levels = c("Early", "Peak", "Late")))
stopifnot(all(meta$SampleID == sids))
sample_region <- factor(meta$region, levels = REGIONS)
sample_season <- factor(meta$transmission_season, levels = SEASONS)

log_msg(paste0("Region sample sizes: ",
               paste(sprintf("%s=%d", REGIONS, table(sample_region)[REGIONS]),
                     collapse = ", ")))
log_msg(paste0("Season sample sizes: ",
               paste(sprintf("%s=%d", SEASONS, table(sample_season)[SEASONS]),
                     collapse = ", ")))

# ============================================================================
# 2. Within-region IBD
# ============================================================================
pairs <- read_tsv(pairs_in, show_col_types = FALSE)
within_pairs <- pairs %>% filter(region1 == region2) %>%
   mutate(region = factor(region1, levels = REGIONS),
          ibd    = fract_sites_IBD)

within_summary <- within_pairs %>%
   group_by(region) %>%
   summarise(n_pairs   = n(),
             mean_ibd  = mean(ibd, na.rm = TRUE),
             median_ibd = median(ibd, na.rm = TRUE),
             frac_high = mean(ibd >= 0.50, na.rm = TRUE),
             .groups = "drop")
log_msg("Within-region IBD summary:")
log_msg(paste(capture.output(print(within_summary)), collapse = "\n"))
write_tsv(within_summary, file.path(out_dir, "ibd_within_region.tsv"))

# Pairwise Wilcoxon w/ Bonferroni
wpair_test <- pairwise.wilcox.test(within_pairs$ibd, within_pairs$region,
                                    p.adjust.method = "bonferroni")
log_msg("Pairwise Wilcoxon (Bonferroni) on within-region IBD:")
log_msg(paste(capture.output(print(wpair_test)), collapse = "\n"))

p_ibd <- ggplot(within_pairs, aes(x = region, y = ibd, fill = region)) +
   geom_violin(alpha = 0.4, scale = "width") +
   geom_boxplot(width = 0.15, outlier.shape = NA) +
   scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.02),
                      breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1.0)) +
   scale_fill_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                Eastern = "#2E86AB"), guide = "none") +
   labs(title = "Within-region pairwise IBD (hmmIBD)",
        subtitle = "Pseudo-log y-axis; box = IQR; line = median",
        x = NULL, y = "fraction sites IBD") +
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

ggsave(file.path(out_dir, "ibd_within_region.pdf"), p_ibd,
       width = 6.5, height = 4.5, dpi = 600)

# ============================================================================
# 3. Per-region per-season nucleotide diversity (pi)
# ============================================================================
# Haploid pi at SNP l in pop P: pi_l = 2 * p * (1 - p) * n / (n - 1)
# where p is alt freq among non-missing samples in P, n is # non-missing.
# Genome-wide pi: mean over SNPs (we report per kept-SNP, not per accessible
# basepair, since our SNP set is filtered upstream).
pi_per_snp <- function(rows_idx, cols_idx) {
   sub <- H[rows_idx, cols_idx, drop = FALSE]
   k   <- rowSums(!is.na(sub))
   alt <- rowSums(sub == 1L, na.rm = TRUE)
   p   <- ifelse(k > 0, alt / k, NA_real_)
   ok  <- k >= 2 & !is.na(p)
   pi  <- ifelse(ok, 2 * p * (1 - p) * k / (k - 1), NA_real_)
   pi
}
pi_mean_ci <- function(pi_vec, reps = PI_BOOT_REPS) {
   x <- pi_vec[!is.na(pi_vec)]
   if (length(x) < 10) return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_,
                                 n_snp = length(x)))
   m <- mean(x)
   boot <- replicate(reps, mean(sample(x, replace = TRUE)))
   c(mean = m, lo = quantile(boot, 0.025), hi = quantile(boot, 0.975),
     n_snp = length(x))
}

pi_rows <- list()
for (r in REGIONS) {
   cols_r <- which(sample_region == r)
   if (length(cols_r) < 2) next
   pi_overall <- pi_per_snp(seq_len(L), cols_r)
   ci <- pi_mean_ci(pi_overall)
   pi_rows[[length(pi_rows) + 1]] <- tibble(region = r, season = "all",
                                             pi = ci["mean"], lo95 = ci["lo.2.5%"],
                                             hi95 = ci["hi.97.5%"], n_samp = length(cols_r),
                                             n_snp = ci["n_snp"])
   for (s in SEASONS) {
      cols_rs <- which(sample_region == r & sample_season == s)
      if (length(cols_rs) < 2) next
      pi_rs <- pi_per_snp(seq_len(L), cols_rs)
      ci    <- pi_mean_ci(pi_rs)
      pi_rows[[length(pi_rows) + 1]] <- tibble(region = r, season = s,
                                                pi = ci["mean"], lo95 = ci["lo.2.5%"],
                                                hi95 = ci["hi.97.5%"], n_samp = length(cols_rs),
                                                n_snp = ci["n_snp"])
   }
}
pi_df <- bind_rows(pi_rows) %>%
   mutate(region = factor(region, levels = REGIONS),
          season = factor(season, levels = c("all", SEASONS)))
write_tsv(pi_df, file.path(out_dir, "pi_region_season.tsv"))
log_msg("Pi per region per season:")
log_msg(paste(capture.output(print(pi_df)), collapse = "\n"))

p_pi <- ggplot(pi_df, aes(x = season, y = pi, fill = region)) +
   geom_col(position = position_dodge(width = 0.85), width = 0.8) +
   geom_errorbar(aes(ymin = lo95, ymax = hi95),
                 position = position_dodge(width = 0.85), width = 0.2) +
   geom_text(aes(label = sprintf("n=%d", n_samp), y = hi95 + 0.005),
             position = position_dodge(width = 0.85), size = 2.8) +
   scale_fill_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                Eastern = "#2E86AB")) +
   labs(title = "Nucleotide diversity (pi) per region per season",
        subtitle = sprintf("Per-SNP pi, mean +/- 95%% bootstrap CI (%d reps); %d SNPs",
                           PI_BOOT_REPS, L),
        x = "Season", y = "pi (per SNP)", fill = "Regions") +
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
ggsave(file.path(out_dir, "pi_region_season.pdf"), p_pi,
       width = 7.5, height = 4.5, dpi = 600)

# ============================================================================
# 4. LD decay per region
# ============================================================================
# Pairwise r^2 between SNPs on same chromosome with pos diff <= LD_MAX_BP.
# Filter to SNPs with MAF >= LD_MIN_MAF in the region (else r^2 is unstable).
ld_one_region <- function(cols_r, max_bp = LD_MAX_BP, bin_bp = LD_BIN_BP) {
   Hr <- H[, cols_r, drop = FALSE]
   k  <- rowSums(!is.na(Hr))
   alt <- rowSums(Hr == 1L, na.rm = TRUE)
   p  <- alt / pmax(k, 1L)
   maf <- pmin(p, 1 - p)
   keep <- !is.na(p) & k >= 5 & maf >= LD_MIN_MAF
   idx  <- which(keep)
   if (length(idx) < 2) return(tibble(dist_bin = integer(0), n = integer(0),
                                       mean_r2 = numeric(0)))

   chr <- chrom[idx]; ps <- pos[idx]; Hk <- Hr[idx, , drop = FALSE]
   r2_list <- list()
   for (c in unique(chr)) {
      i_c <- which(chr == c)
      if (length(i_c) < 2) next
      for (a in seq_along(i_c)) {
         ia <- i_c[a]
         # SNPs to the right within max_bp
         j_c <- i_c[(a + 1):length(i_c)]
         if (length(j_c) < 1) next
         within <- j_c[ps[j_c] - ps[ia] <= max_bp]
         if (length(within) == 0) next
         g1 <- Hk[ia, ]
         for (jb in within) {
            g2 <- Hk[jb, ]
            both <- !is.na(g1) & !is.na(g2)
            if (sum(both) < 5) next
            x <- g1[both]; y <- g2[both]
            pA <- mean(x); pB <- mean(y)
            if (pA == 0 || pA == 1 || pB == 0 || pB == 1) next
            pAB <- mean(x * y)
            D   <- pAB - pA * pB
            r2  <- D^2 / (pA * (1 - pA) * pB * (1 - pB))
            r2_list[[length(r2_list) + 1]] <-
               c(dist = ps[jb] - ps[ia], r2 = r2)
         }
      }
   }
   if (length(r2_list) == 0)
      return(tibble(dist_bin = integer(0), n = integer(0), mean_r2 = numeric(0)))
   m <- do.call(rbind, r2_list)
   tibble(dist = as.integer(m[, "dist"]), r2 = m[, "r2"]) %>%
      mutate(dist_bin = (dist %/% bin_bp) * bin_bp + bin_bp / 2L) %>%
      group_by(dist_bin) %>%
      summarise(n = n(), mean_r2 = mean(r2), .groups = "drop") %>%
      filter(n >= 5)
}

ld_rows <- list()
for (r in REGIONS) {
   cols_r <- which(sample_region == r)
   if (length(cols_r) < 5) next
   log_msg(sprintf("LD decay: region=%s (n=%d)...", r, length(cols_r)))
   ld_r <- ld_one_region(cols_r)
   if (nrow(ld_r) > 0)
      ld_rows[[length(ld_rows) + 1]] <- ld_r %>% mutate(region = r)
}
ld_df <- bind_rows(ld_rows) %>%
   mutate(region = factor(region, levels = REGIONS))
write_tsv(ld_df, file.path(out_dir, "ld_decay.tsv"))

# Half-decay distance (r^2 falls below half of the near-zero baseline)
ld_half <- ld_df %>% group_by(region) %>%
   summarise(r2_at_short = mean(mean_r2[dist_bin < 1000]),
             r2_at_far   = mean(mean_r2[dist_bin > 40000]),
             half_target = (mean(mean_r2[dist_bin < 1000]) +
                            mean(mean_r2[dist_bin > 40000])) / 2,
             half_dist   = {
                tgt <- (mean(mean_r2[dist_bin < 1000]) +
                        mean(mean_r2[dist_bin > 40000])) / 2
                first_below <- dist_bin[mean_r2 <= tgt][1]
                first_below
             },
             .groups = "drop")
log_msg("LD half-decay distance per region:")
log_msg(paste(capture.output(print(ld_half)), collapse = "\n"))
write_tsv(ld_half, file.path(out_dir, "ld_half_decay.tsv"))

p_ld <- ggplot(ld_df, aes(x = dist_bin, y = mean_r2, colour = region)) +
   geom_smooth(se = TRUE, method = "loess", span = 0.3, linewidth = 0.7) +
   geom_point(alpha = 0.5, size = 1) +
   scale_x_continuous(labels = label_number(scale = 1e-3, suffix = " kb")) +
   scale_colour_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                  Eastern = "#2E86AB")) +
   labs(title = "LD decay per region",
        subtitle = sprintf("Mean r^2 in %.0f-bp bins; pairs within %.0f kb; MAF >= %.2f",
                           LD_BIN_BP, LD_MAX_BP / 1000, LD_MIN_MAF),
        x = "Distance between SNPs", y = "Mean r^2", colour = "Regions") +
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

ggsave(file.path(out_dir, "ld_decay.pdf"), p_ld,
       width = 7.5, height = 4.5, dpi = 600)

# ============================================================================
# 5. Persist + log
# ============================================================================
saveRDS(list(within_summary = within_summary,
             wpair_test    = wpair_test,
             pi_df         = pi_df,
             ld_df         = ld_df,
             ld_half       = ld_half),
        file.path(out_dir, "phase3_03_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase3_03_summary.txt"))
message("\nDone. Outputs in ", out_dir)
