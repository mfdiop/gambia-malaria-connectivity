# =============================================================================
# Phase 3 — Step 3.1: Supervised ancestry deconvolution (Gambia-only)
# =============================================================================
# Per notes/objective3_proposal.md §6 Phase 3 Step 3.1:
#   - Fit an ADMIXTURE model treating the three Gambian regions (west, central,
#     east) as three source populations.
#   - For each sample, estimate Q = (Q_west, Q_central, Q_east). For western
#     samples, Q_east is the individual-level "importation score" from east.
#   - Flag western samples with Q_east > 0.20 (proposal threshold) as probable
#     imports.
#   - Temporal: does the fraction of west samples with substantial eastern
#     ancestry shift across early / peak / late season?
#
# Why supervised, not unsupervised:
#   The Phase 2.2 unsupervised K = 3 run did not produce three components that
#   aligned with west / central / east — one component dominated. To answer
#   the proposal's question, we need each component anchored to a region's
#   allele frequencies. ADMIXTURE's --supervised flag forces labeled samples
#   to Q = 1 on their own label and so cannot return their "non-local"
#   ancestry; the standard fix is leave-one-out (LOO): when scoring sample i,
#   its region's allele frequencies are recomputed *without* sample i.
#
# Model:
#   Pf monoclonal, haploid: g_il in {0, 1} (allele dose of alt), with NA for
#   missing / het-in-monoclonal sites. Reference frequencies F_rl per region.
#   Per-sample Q_i = (Q_west, Q_central, Q_east), simplex-constrained.
#   Likelihood: a_il = g_il ~ Bernoulli(p_il), p_il = sum_r Q_ir * F_rl.
#   EM (per sample):
#     E: gamma_irl = Q_ir * P(a_il | r) / sum_r' Q_ir' * P(a_il | r')
#     M: Q_ir_new = mean_l gamma_irl
#   Converges to MLE on the simplex.
#
# Inputs : results/phase2_02_lea/snmf/gambia_mono_pruned.geno
#          results/phase2_02_lea/phase2_02_lea_results.rds  (sample order)
#          results/phase1_03/sample_master_metadata.rds
#          results/phase1_03/village_centroids.tsv
#
# Output : results/phase3_01/
#            ancestry_Q.tsv                 per-sample Q + region + season
#            importation_flags.tsv          west samples with Q_east > 0.20
#            ancestry_barplot.pdf           stacked bar, ordered by region/Q_east
#            ancestry_Q_east_by_region.pdf  density / violin of Q_east per region
#            ancestry_temporal_west.pdf     Q_east vs season for west samples
#            ancestry_village_means.pdf     per-village mean ancestry barplot
#            phase3_01_summary.txt          human-readable log
#            phase3_01_results.rds          persisted state
#
# Tools  : R only (tidyverse, scales, ggrepel)
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
# K is fixed to 3 by the supervised design (west / central / east), so no --k
# wildcard here. Geno/phase22 default to canonical K=3 dir from Phase 2.2.
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir  <- .parse_arg("output_dir", "results/phase3_01")
geno_arg    <- .parse_arg("geno",       "results/phase2_02_lea/K3/snmf/gambia_mono_pruned.geno")
phase22_arg <- .parse_arg("phase22",    "results/phase2_02_lea/K3/phase2_02_lea_results.rds")
meta_arg    <- .parse_arg("meta",       "results/phase1_03/sample_master_metadata.rds")
cent_arg    <- .parse_arg("centroid",   "results/phase1_03/village_centroids.tsv")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(scales)
   library(ggrepel)
})

# ----- Paths ----------------------------------------------------------------
geno_in    <- geno_arg
phase22_in <- phase22_arg
meta_in    <- meta_arg
cent_in    <- cent_arg

out_dir   <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Parameters -----------------------------------------------------------
REGIONS         <- c("Western", "Central", "Eastern")
IMPORT_THRESH   <- 0.20         # proposal: >20% eastern ancestry -> import flag
EM_MAX_ITER     <- 5000
EM_TOL_Q        <- 1e-8         # max-abs step in Q (loose; fallback)
EM_TOL_LL       <- 1e-8         # relative log-likelihood change (primary)
FREQ_EPS        <- 1e-6         # clip F to (eps, 1-eps) to avoid log(0)
MIN_REGION_N    <- 3            # require >= N samples per region for LOO

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load genotypes (.geno from Phase 2.2) + sample IDs + metadata
# ============================================================================
log_msg("=== Phase 3 Step 3.1: Supervised ancestry deconvolution ===")
phase22 <- readRDS(phase22_in)
sample_ids <- phase22$eigenvec$SampleID
log_msg(sprintf("Samples: %d (order from phase2_02 .geno)", length(sample_ids)))

# .geno is rows = SNPs, cols = samples, no separator, single digit per sample
geno_lines <- readLines(geno_in)
L <- length(geno_lines)
n <- length(sample_ids)
geno_pseudo <- matrix(0L, nrow = L, ncol = n)

for (i in seq_len(L)) {
   chars <- strsplit(geno_lines[i], "", fixed = TRUE)[[1]]
   stopifnot(length(chars) == n)
   geno_pseudo[i, ] <- as.integer(chars)
}
log_msg(sprintf("Loaded .geno: %d SNPs x %d samples", L, n))

# Convert pseudo-haploid {0, 2, 9} -> haploid {0, 1, NA}
# (Phase 2.2 encoded haploid 0->0, 1->2, missing/het->9 for diploid-model snmf.)
hap <- geno_pseudo
hap[hap == 0L] <- 0L
hap[hap == 2L] <- 1L
hap[geno_pseudo == 9L | geno_pseudo == 1L] <- NA_integer_
storage.mode(hap) <- "integer"
n_missing <- sum(is.na(hap))
log_msg(sprintf("Haploid matrix: %d observed, %d missing (%.2f%%)",
                sum(!is.na(hap)), n_missing, 100 * n_missing / length(hap)))

meta <- readRDS(meta_in)
meta_aln <- meta %>% 
   filter(SampleID %in% sample_ids) %>%
   arrange(match(SampleID, sample_ids)) %>% 
   mutate(
      transmission_season = recode(transmission_season, "early"= "Early", "peak" = "Peak", "late" = "Late"),
      transmission_season = factor(transmission_season, 
                                   levels = c("Early", "Peak", "Late"))
      )

stopifnot(all(meta_aln$SampleID == sample_ids))
sample_region <- as.character(meta_aln$region)
stopifnot(all(sample_region %in% REGIONS))

region_n <- table(factor(sample_region, levels = REGIONS))
log_msg(paste0("Region sample sizes: ",
               paste(sprintf("%s=%d", names(region_n), region_n),
                     collapse = ", ")))

if (any(region_n < MIN_REGION_N))
   stop(sprintf("At least one region has < %d samples; LOO unsafe.",
                MIN_REGION_N))

# ============================================================================
# 2. Per-region allele counts + reference frequencies (full, not LOO)
# ============================================================================
# sum_rl  : sum of haploid alleles in region r at SNP l (observed only)
# n_rl    : observed-allele count in region r at SNP l
# F_rl    : sum_rl / n_rl, clipped to (eps, 1-eps) for log-stability
sum_rl <- matrix(0L, nrow = L, ncol = length(REGIONS))
n_rl   <- matrix(0L, nrow = L, ncol = length(REGIONS))

for (r_idx in seq_along(REGIONS)) {
   sub <- hap[, sample_region == REGIONS[r_idx], drop = FALSE]
   sum_rl[, r_idx] <- rowSums(sub, na.rm = TRUE)
   n_rl[, r_idx]   <- rowSums(!is.na(sub))
}

F_full <- sum_rl / pmax(n_rl, 1L)
F_full[n_rl == 0L] <- NA_real_
log_msg("Computed per-region allele frequencies (all samples).")

# Filter to SNPs with adequate coverage in every region
ok_snp <- rowSums(n_rl >= 2L) == length(REGIONS)
log_msg(sprintf("SNPs usable (>=2 observations per region): %d / %d",
                sum(ok_snp), L))
sum_rl <- sum_rl[ok_snp, , drop = FALSE]
n_rl   <- n_rl[ok_snp, , drop = FALSE]
F_full <- F_full[ok_snp, , drop = FALSE]
hap    <- hap[ok_snp, , drop = FALSE]

# ============================================================================
# 3. Per-sample EM with leave-one-out on own region
# ============================================================================
clip_freqs <- function(F) pmax(pmin(F, 1 - FREQ_EPS), FREQ_EPS)

em_fit <- function(a, F) {
   # a : integer vector length L, values in {0, 1, NA}
   # F : L x R reference freqs (LOO-adjusted, clipped, observation-aligned)
   obs <- !is.na(a)
   a_obs <- a[obs]
   F_obs <- F[obs, , drop = FALSE]
   R <- ncol(F_obs)
   Q <- rep(1 / R, R)

   alt <- a_obs == 1L
   F_alt <- F_obs[alt,  , drop = FALSE]
   F_ref <- 1 - F_obs[!alt, , drop = FALSE]
   nL <- nrow(F_obs)

   ll_prev <- -Inf
   converged <- FALSE
   for (it in seq_len(EM_MAX_ITER)) {
      pa <- as.vector(F_alt %*% Q)              # p(alt) at each alt site
      pr <- as.vector(F_ref %*% Q)              # p(ref) at each ref site
      ll <- sum(log(pa)) + sum(log(pr))         # observed-data log-likelihood
      gam_alt <- sweep(F_alt, 2, Q, `*`) / pa
      gam_ref <- sweep(F_ref, 2, Q, `*`) / pr
      Q_new <- (colSums(gam_alt) + colSums(gam_ref)) / nL
      Q_new <- Q_new / sum(Q_new)               # numerical safety
      step <- max(abs(Q_new - Q))
      ll_rel <- if (it == 1) Inf else abs(ll - ll_prev) / max(1, abs(ll_prev))
      Q <- Q_new
      ll_prev <- ll
      if (ll_rel < EM_TOL_LL || step < EM_TOL_Q) { converged <- TRUE; break }
   }
   list(Q = Q, iter = it, ll = ll_prev, step = step, converged = converged)
}

log_msg(sprintf("Running EM for %d samples (LOO on own region)...", n))

Q_all <- matrix(NA_real_, nrow = n, ncol = length(REGIONS),
                dimnames = list(sample_ids, REGIONS))
iters     <- integer(n)
final_step <- numeric(n)
converged <- logical(n)

for (i in seq_len(n)) {
   r_i <- match(sample_region[i], REGIONS)
   F_loo <- F_full
   # LOO adjustment: at sites where sample i is observed, remove its contribution
   mask <- !is.na(hap[, i])
   if (any(mask)) {
      new_sum <- sum_rl[mask, r_i] - hap[mask, i]
      new_n   <- n_rl[mask, r_i] - 1L
      # If new_n == 0 (singleton region at this site), fall back to full F (rare)
      safe <- new_n > 0L
      F_loo[mask, r_i][safe] <- new_sum[safe] / new_n[safe]
   }
   F_loo_clip <- clip_freqs(F_loo)
   res <- em_fit(hap[, i], F_loo_clip)
   Q_all[i, ]   <- res$Q
   iters[i]     <- res$iter
   final_step[i] <- res$step
   converged[i]  <- res$converged
}

log_msg(sprintf("EM done. Iter: median=%.0f max=%d; converged=%d/%d; final-step max=%.2e median=%.2e",
                median(iters), max(iters), sum(converged), n,
                max(final_step), median(final_step)))

# ============================================================================
# 4. Build ancestry table; flag imports
# ============================================================================
Q_df <- as_tibble(Q_all, rownames = "SampleID") %>%
   left_join(meta_aln %>% select(SampleID, region, Location, VillageCode,
                                  transmission_season, collection_month,
                                  latitude, longitude),
             by = "SampleID") %>%
   mutate(region = factor(region, levels = REGIONS),
          transmission_season = factor(transmission_season,
                                       levels = c("Early", "Peak", "Late")))

write_tsv(Q_df, file.path(out_dir, "ancestry_Q.tsv"))

import_west <- Q_df %>%
   filter(region == "Western", Eastern > IMPORT_THRESH) %>%
   arrange(desc(Eastern))
write_tsv(import_west, file.path(out_dir, "importation_flags.tsv"))

log_msg(sprintf("Western samples with Q_east > %.2f: %d / %d",
                IMPORT_THRESH, nrow(import_west),
                sum(Q_df$region == "Western")))

# Region-level Q summary
region_summary <- Q_df %>%
   group_by(region) %>%
   summarise(mean_Q_west    = mean(Western),
             mean_Q_central = mean(Central),
             mean_Q_east    = mean(Eastern),
             median_Q_east  = median(Eastern),
             frac_above_thresh = mean(Eastern > IMPORT_THRESH),
             .groups = "drop")

log_msg("Per-region mean Q:")
log_msg(paste(capture.output(print(region_summary)), collapse = "\n"))
write_tsv(region_summary, file.path(out_dir, "ancestry_region_summary.tsv"))

# ============================================================================
# 5. Plots
# ============================================================================
# 5a. Stacked-bar ancestry per sample, ordered by region then by Q_east
order_df <- Q_df %>% 
   arrange(region, Eastern) %>% 
   mutate(rank = row_number())

Q_long <- order_df %>%
   pivot_longer(all_of(REGIONS), names_to = "source", values_to = "prop") %>%
   mutate(source = factor(source, levels = REGIONS),
          SampleID = factor(SampleID, levels = order_df$SampleID))

p_bar <- ggplot(Q_long, aes(x = SampleID, y = prop, fill = source)) +
   geom_col(width = 1) +
   facet_grid(~ region, scales = "free_x", space = "free_x", switch = "x") +
   scale_y_continuous(expand = c(0, 0), labels = percent) +
   scale_fill_manual(values = c(Western = "#D7263D", 
                                Central = "#F46036",
                                Eastern = "#2E86AB"),
                     name = "Source region") +
   labs(title = "Supervised ancestry, K=3 regional sources (leave-one-out)",
        subtitle = sprintf("Samples within each region ordered by Q_east; n = %d", n),
        x = NULL, y = "Ancestry proportion") +
   theme_minimal(base_size = 13) +
   theme(
      legend.title  = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
      legend.text   = element_text(size = 11, color = "black"),
      plot.title    = element_text(size = 13, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text.x   = element_blank(),
      axis.text.y   = element_text(size = 12, color = "black"),
      axis.title    = element_text(size = 14, color = "black", face = "bold"),
      axis.ticks.x  = element_blank(),
      panel.spacing.x = unit(0.2, "lines"),
      strip.placement = "outside"
      )

ggsave(file.path(out_dir, "ancestry_barplot.pdf"), p_bar,
       width = 12, height = 4, dpi = 600)

# 5b. Q_east distribution per region (violin + jitter)
p_qe <- ggplot(Q_df, aes(x = region, y = Eastern, fill = region)) +
   geom_violin(alpha = 0.5) +
   geom_boxplot(width = 0.18, outlier.shape = NA) +
   geom_jitter(width = 0.12, alpha = 0.7, size = 1.4) +
   geom_hline(yintercept = IMPORT_THRESH, linetype = "dashed",
              colour = "red") +
   annotate("text", x = 0.55, y = IMPORT_THRESH + 0.02,
            label = sprintf("import threshold = %.2f", IMPORT_THRESH),
            hjust = 0, size = 3, colour = "red") +
   scale_fill_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                Eastern = "#2E86AB"), guide = "none") +
   labs(title = "Eastern-ancestry proportion (Q_east) by region",
        subtitle = "Supervised LOO admixture; western samples above red line = candidate imports",
        x = NULL, y = "Q_east") +
   theme_minimal() +
   theme(
      plot.title    = element_text(size = 13, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text     = element_text(size = 12, color = "black"),
      axis.title    = element_text(size = 14, color = "black", face = "bold"),
      axis.line         = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks        = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm")
   )

ggsave(file.path(out_dir, "ancestry_Q_east_by_region.pdf"), p_qe,
       width = 6.5, height = 4.5, dpi = 600)

# 5c. Temporal: Q_east per season among west samples
p_temp <- Q_df %>% filter(region == "Western") %>%
   ggplot(aes(x = transmission_season, y = Eastern, fill = transmission_season)) +
   geom_violin(alpha = 0.5) +
   geom_boxplot(width = 0.15, outlier.shape = NA) +
   geom_jitter(width = 0.12, alpha = 0.75, size = 1.5) +
   geom_hline(yintercept = IMPORT_THRESH, linetype = "dashed", colour = "red") +
   labs(title = "Western samples: Q_east by transmission season",
        subtitle = "Tests whether importation pressure shifts with season",
        x = "Transmission season", y = "Q_east") +
   theme_minimal() +
   theme(
      legend.position = "none",
      plot.title    = element_text(size = 13, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text     = element_text(size = 12, color = "black"),
      axis.title    = element_text(size = 14, color = "black", face = "bold"),
      axis.line         = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks        = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm")
   )

ggsave(file.path(out_dir, "ancestry_temporal_west.pdf"), p_temp,
       width = 6, height = 4.5, dpi = 600)

# Kruskal-Wallis on Q_east ~ season among west samples (>=2 seasons present)
west_qd <- Q_df %>% filter(region == "Western")
if (nlevels(droplevels(west_qd$transmission_season)) >= 2 &&
    nrow(west_qd) >= 6) {
   kw <- kruskal.test(Eastern ~ transmission_season, data = west_qd)
   log_msg(sprintf("Kruskal-Wallis Q_east ~ season (West): chi2 = %.3f, p = %.3g",
                   kw$statistic, kw$p.value))
} else {
   kw <- NULL
   log_msg("Skipped Kruskal-Wallis: insufficient seasonal variation among west samples.")
}

# 5d. Per-village mean ancestry, ordered by longitude (west -> east)
cent_raw <- read_tsv(cent_in, show_col_types = FALSE)
cent <- cent_raw %>%
   group_by(Location, region) %>%
   summarise(lon = weighted.mean(lon, n),
             lat = weighted.mean(lat, n),
             n   = sum(n), .groups = "drop")

village_means <- Q_df %>%
   group_by(Location) %>%
   summarise(across(all_of(REGIONS), mean), n_samp = n(), .groups = "drop") %>%
   left_join(cent %>% select(Location, lon), by = "Location") %>%
   arrange(lon) %>%
   mutate(Location = factor(Location, levels = Location))

vm_long <- village_means %>%
   pivot_longer(all_of(REGIONS), names_to = "source", values_to = "prop") %>%
   mutate(source = factor(source, levels = REGIONS))

p_vill <- ggplot(vm_long, aes(x = Location, y = prop, fill = source)) +
   geom_col(width = 0.85) +
   geom_text(data = village_means,
             aes(x = Location, y = 1.02, label = sprintf("n=%d", n_samp)),
             inherit.aes = FALSE, size = 3, vjust = 0) +
   scale_y_continuous(limits = c(0, 1.07), expand = c(0, 0), labels = percent) +
   scale_fill_manual(values = c(Western = "#D7263D", Central = "#F46036",
                                Eastern = "#2E86AB"), name = "Source") +
   labs(title = "Per-village mean ancestry composition",
        subtitle = "Villages ordered west -> east by longitude",
        x = NULL, y = "Mean ancestry proportion") +
   theme_minimal() + 
   theme(
      legend.title  = element_text(size = 14, color = "black", face = "bold", hjust = 0.5),
      legend.text   = element_text(size = 11, color = "black"),
      plot.title    = element_text(size = 13, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray80", face = "bold"),
      axis.text.y   = element_text(size = 12, color = "black"),
      axis.title    = element_text(size = 14, color = "black", face = "bold"),
      axis.text.x   = element_text(angle = 30, hjust = 1)
      )

ggsave(file.path(out_dir, "ancestry_village_means.pdf"), p_vill,
       width = 8, height = 4.5, dpi = 600)

# ============================================================================
# 6. Persist + log
# ============================================================================
saveRDS(list(Q             = Q_all,
             Q_df          = Q_df,
             region_summary = region_summary,
             import_west   = import_west,
             n_snps_used   = nrow(hap),
             iters         = iters,
             kruskal_west  = kw),
        file.path(out_dir, "phase3_01_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase3_01_summary.txt"))
message("\nDone. Outputs in ", out_dir)
