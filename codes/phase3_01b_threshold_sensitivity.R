# =============================================================================
# Phase 3 — Step 3.1b: Sensitivity analysis for the Q_east importation
#                       threshold and other fixed pipeline parameters
# =============================================================================
# Companion to phase3_01_supervised_admixture.R. Addresses the reviewer-
# facing question: "why Q_east > 0.20, and is it robust?" — plus a parallel
# question reviewers often ask next: "are MIN_REGION_N=3 and the implicit
# min-2-observations-per-region SNP filter also arbitrary?"
#
# Five analyses, run in order:
#   A. Threshold sweep        — does the main conclusion hold across a grid
#                                of thresholds? (cheap; reuses baseline fit)
#   B. Permutation null       — what Q_east would a NON-imported western
#                                sample show, just from estimation noise?
#                                (expensive; refits EM under shuffled labels)
#   C. Beta-mixture threshold — is there a natural two-component split in
#                                the western Q_east distribution, and where
#                                is its crossing point? (cheap)
#   D. Bootstrap-over-SNPs CI — per-sample uncertainty in Q_east from a
#                                finite SNP panel (expensive; refits EM per
#                                resample)
#   E. Parameter sweep        — how sensitive are the import calls to
#                                MIN_REGION_N and the per-SNP coverage
#                                filter? (moderate; refits EM per config)
#
# RUNTIME NOTE: B and D each refit the full LOO-EM pipeline many times.
# Defaults below (n_perm=199, B=200) are a reasonable first pass. Increase
# for final manuscript numbers (n_perm=999+) and use n_cores > 1 (Unix/Mac)
# to parallelise via parallel::mclapply.
#
# Output: results/phase3_01b_sensitivity/
#   threshold_sweep.tsv / .pdf
#   permutation_null.tsv / .pdf            (+ derived threshold)
#   mixture_threshold.tsv / .pdf           (+ derived threshold)
#   bootstrap_qeast_ci.tsv / .pdf
#   parameter_sweep.tsv / .pdf
#   threshold_comparison_summary.txt       (the headline reviewer-facing table)
#   phase3_01b_sensitivity_results.rds
# =============================================================================

.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
   hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
   if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
geno_arg    <- .parse_arg("geno",       "results/phase2_02_lea/K3/snmf/gambia_mono_pruned.geno")
phase22_arg <- .parse_arg("phase22",    "results/phase2_02_lea/K3/phase2_02_lea_results.rds")
meta_arg    <- .parse_arg("meta",       "results/phase1_03/sample_master_metadata.rds")
output_dir  <- .parse_arg("output_dir", "results/phase3_01b_sensitivity")
core_lib    <- .parse_arg("core_lib",   "codes/phase3_01b_sensitivity_core.R")

# Toggle the expensive analyses independently so you can iterate quickly on
# A/C/E first, then switch these on for the final run.
RUN_PERMUTATION <- as.logical(.parse_arg("run_permutation", "TRUE"))
RUN_BOOTSTRAP   <- as.logical(.parse_arg("run_bootstrap",   "TRUE"))
N_PERM   <- as.integer(.parse_arg("n_perm", "199"))
N_BOOT   <- as.integer(.parse_arg("n_boot", "200"))
N_CORES  <- as.integer(.parse_arg("n_cores", "1"))

suppressPackageStartupMessages({
   library(tidyverse)
   library(scales)
})

source(core_lib)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

REGIONS       <- c("Western", "Central", "Eastern")
WEST_LABEL    <- "Western"
EAST_LABEL    <- "Eastern"
IMPORT_THRESH <- 0.20   # the protocol value being stress-tested below

log_msg("=== Phase 3 Step 3.1b: Threshold & parameter sensitivity analysis ===")

# ============================================================================
# 0. Load raw inputs once; fit the BASELINE configuration (matches the
#    original script's defaults: min_obs_per_region_snp=2, min_region_n=3,
#    freq_eps=1e-6) for use by analyses A and C.
# ============================================================================
raw <- load_phase3_01_raw_inputs(geno_arg, phase22_arg, meta_arg, REGIONS)
log_msg(sprintf("Raw data: %d SNPs x %d samples", raw$L_total, length(raw$sample_ids)))

baseline <- fit_one_config(raw$hap, raw$sample_region, REGIONS,
                           min_obs_per_region_snp = 2, min_region_n = 3,
                           sample_ids = raw$sample_ids)

stopifnot(is.null(baseline$error) || is.na(baseline$error))
log_msg(sprintf("Baseline fit: %d SNPs used, %d/%d samples converged",
                baseline$n_snps_used, sum(baseline$converged), length(raw$sample_ids)))

Q_df <- build_Q_df(baseline$Q, raw$meta_aln, REGIONS)
west_qd <- Q_df %>% filter(region == WEST_LABEL)
log_msg(sprintf("Western samples in baseline fit: %d", nrow(west_qd)))

# ============================================================================
# A. THRESHOLD SWEEP — does the conclusion survive a grid of thresholds?
# ============================================================================
log_msg("--- A. Threshold sweep ---")

threshold_sweep <- function(Q_df, thresholds, west_label, east_label,
                            season_col = "transmission_season") {
   map_dfr(thresholds, function(thr) {
      west <- Q_df %>% filter(region == west_label)
      flagged <- west %>% filter(.data[[east_label]] > thr)
      kw_p <- NA_real_
      if (season_col %in% names(west) &&
          nlevels(droplevels(west[[season_col]])) >= 2 && nrow(west) >= 6) {
         kw <- suppressWarnings(
            kruskal.test(west[[east_label]] ~ droplevels(west[[season_col]])))
         kw_p <- kw$p.value
      }
      tibble(threshold = thr, n_flagged = nrow(flagged), n_west = nrow(west),
             frac_flagged = nrow(flagged) / nrow(west), kw_p_season = kw_p)
   })
}

thr_grid <- seq(0.05, 0.50, by = 0.025)
sweep_df <- threshold_sweep(Q_df, thr_grid, WEST_LABEL, EAST_LABEL)
write_tsv(sweep_df, file.path(output_dir, "threshold_sweep.tsv"))
log_msg(sprintf("Threshold sweep: at protocol threshold (%.2f), %d/%d (%.1f%%) flagged",
                IMPORT_THRESH,
                sweep_df$n_flagged[which.min(abs(sweep_df$threshold - IMPORT_THRESH))],
                sweep_df$n_west[1],
                100 * sweep_df$frac_flagged[which.min(abs(sweep_df$threshold - IMPORT_THRESH))]))

p_sweep <- ggplot(sweep_df, aes(x = threshold)) +
   geom_line(aes(y = frac_flagged), colour = "#2E86AB", linewidth = 1) +
   geom_point(aes(y = frac_flagged), colour = "#2E86AB", size = 2) +
   geom_vline(xintercept = IMPORT_THRESH, linetype = "dashed", colour = "red") +
   annotate("text", x = IMPORT_THRESH + 0.005, y = max(sweep_df$frac_flagged, na.rm = TRUE)-0.02,
            label = "protocol threshold", hjust = 0, size = 3, colour = "red") +
   scale_y_continuous(labels = percent, limits = c(0, NA)) +
   labs(title = "Sensitivity of importation calls to the Q_east threshold",
        subtitle = "Fraction of western samples flagged, across the candidate threshold range",
        x = "Q_east threshold", y = "Fraction of western samples flagged") +
   theme_minimal()

ggsave(file.path(output_dir, "threshold_sweep.pdf"), p_sweep, width = 6.5, height = 4.5, dpi = 600)

# ============================================================================
# B. PERMUTATION NULL — Q_east expected under NO true east/west structure
# ============================================================================
permutation_null_threshold <- function(hap, sample_region, REGIONS,
                                       west_label, east_label, n_perm = 199,
                                       min_obs_per_region_snp = 2,
                                       min_region_n = 3, seed = 1,
                                       n_cores = 1, verbose = TRUE) {
   set.seed(seed)
   perm_labels <- replicate(n_perm, sample(sample_region), simplify = FALSE)
   
   one_perm <- function(perm_region) {
      fit <- fit_one_config(hap, perm_region, REGIONS,
                            min_obs_per_region_snp, min_region_n,
                            sample_ids = colnames(hap))
      if (is.null(fit$Q)) return(NULL)
      fit$Q[perm_region == west_label, east_label]
   }
   
   if (n_cores > 1) {
      out <- parallel::mclapply(perm_labels, one_perm, mc.cores = n_cores)
   } else {
      out <- vector("list", n_perm)
      for (p in seq_len(n_perm)) {
         out[[p]] <- one_perm(perm_labels[[p]])
         if (verbose && p %% 25 == 0) message(sprintf("  permutation %d/%d", p, n_perm))
      }
   }
   null_vec <- unlist(out)
   list(null_distribution = null_vec,
        q95 = unname(quantile(null_vec, 0.95, na.rm = TRUE)),
        q99 = unname(quantile(null_vec, 0.99, na.rm = TRUE)),
        mean = mean(null_vec, na.rm = TRUE), sd = sd(null_vec, na.rm = TRUE),
        n_perm_effective = sum(!sapply(out, is.null)))
}

if (RUN_PERMUTATION) {
   log_msg(sprintf("--- B. Permutation null (n_perm=%d) ---", N_PERM))
   perm_res <- permutation_null_threshold(raw$hap, raw$sample_region, REGIONS,
                                          WEST_LABEL, EAST_LABEL, n_perm = N_PERM,
                                          n_cores = N_CORES)
   write_tsv(tibble(q_east_null = perm_res$null_distribution),
             file.path(output_dir, "permutation_null.tsv"))
   
   log_msg(sprintf(
      "Null Q_east (no true structure): mean=%.3f sd=%.3f | 95th pct=%.3f | 99th pct=%.3f",
      perm_res$mean, perm_res$sd, perm_res$q95, perm_res$q99))
   
   log_msg(sprintf(
      "Protocol threshold (%.2f) is %s the permutation 95th percentile (%.3f)",
      IMPORT_THRESH, if (IMPORT_THRESH > perm_res$q95) "ABOVE" else "BELOW", perm_res$q95))
   
   p_perm <- ggplot() +
      geom_histogram(data = tibble(x = perm_res$null_distribution),
                     aes(x = x, y = after_stat(density)), bins = 40,
                     fill = "grey70", colour = "white") +
      geom_density(data = west_qd, aes(x = .data[[EAST_LABEL]]),
                   colour = "#2E86AB", linewidth = 1) +
      geom_vline(xintercept = IMPORT_THRESH, linetype = "dashed", colour = "red") +
      geom_vline(xintercept = perm_res$q95, linetype = "dotted", colour = "black") +
      annotate("text", x = IMPORT_THRESH, y = 0, label = "protocol", angle = 90,
               vjust = -0.5, hjust = 0, size = 3, colour = "red") +
      annotate("text", x = perm_res$q95, y = 0, label = "null 95th pct", angle = 90,
               vjust = 1.3, hjust = 0, size = 3) +
      labs(title = "Observed Q_east (western samples) vs. permutation null",
           subtitle = "Grey = Q_east under shuffled region labels; blue = observed western distribution",
           x = "Q_east", y = "Density") +
      theme_minimal()
   
   ggsave(file.path(output_dir, "permutation_null.pdf"), p_perm, width = 6.5, height = 4.5, dpi = 600)
} else {
   log_msg("--- B. Permutation null: SKIPPED (run_permutation=FALSE) ---")
   perm_res <- NULL
}

# ============================================================================
# C. BETA-MIXTURE THRESHOLD — data-driven crossing point, if bimodal
# ============================================================================
log_msg("--- C. Beta-mixture threshold ---")

.mom_params_weighted <- function(v, w) {
   w <- w / sum(w)
   m <- sum(w * v); s2 <- sum(w * (v - m)^2); s2 <- max(s2, 1e-6)
   common <- m * (1 - m) / s2 - 1
   c(max(m * common, 0.1), max((1 - m) * common, 0.1))
}

# Two-component beta mixture via EM with method-of-moments M-step. This is
# an approximation (no closed-form MLE for beta mixtures) intended to locate
# a plausible crossing point, not for formal inference — if you need a rigorous
# fit, refit with mixtools::betamix or flexmix::flexmix(FLXMRbeta()) and use
# this as a sanity check / starting value.
fit_beta_mixture_2 <- function(x, max_iter = 500, tol = 1e-6, seed = 1) {
   set.seed(seed)
   x <- pmin(pmax(x[is.finite(x)], 1e-4), 1 - 1e-4)
   med <- median(x)
   pi_ <- c(0.7, 0.3)
   p1 <- .mom_params_weighted(x, as.numeric(x <= med) + 1e-6)
   p2 <- .mom_params_weighted(x, as.numeric(x >  med) + 1e-6)
   ll_prev <- -Inf; it <- 0L
   for (it in seq_len(max_iter)) {
      d1 <- pi_[1] * dbeta(x, p1[1], p1[2])
      d2 <- pi_[2] * dbeta(x, p2[1], p2[2])
      denom <- pmax(d1 + d2, 1e-300)
      gam1 <- d1 / denom; gam2 <- d2 / denom
      pi_ <- c(mean(gam1), mean(gam2))
      p1 <- .mom_params_weighted(x, gam1)
      p2 <- .mom_params_weighted(x, gam2)
      ll <- sum(log(denom))
      if (abs(ll - ll_prev) < tol) break
      ll_prev <- ll
   }
   list(pi = pi_, params1 = p1, params2 = p2, loglik = ll_prev, iter = it)
}

find_mixture_crossing <- function(fit, grid = seq(0.001, 0.999, by = 0.001)) {
   d1 <- fit$pi[1] * dbeta(grid, fit$params1[1], fit$params1[2])
   d2 <- fit$pi[2] * dbeta(grid, fit$params2[1], fit$params2[2])
   mean1 <- fit$params1[1] / sum(fit$params1)
   mean2 <- fit$params2[1] / sum(fit$params2)
   diff_d <- d1 - d2
   sc <- which(diff_d[-1] * diff_d[-length(diff_d)] < 0)
   if (length(sc) == 0) return(NA_real_)
   rng <- sort(c(mean1, mean2))
   cand <- grid[sc]
   cand <- cand[cand > rng[1] & cand < rng[2]]
   if (length(cand) == 0) return(NA_real_)
   cand[1]
}

mix_fit <- tryCatch(fit_beta_mixture_2(west_qd[[EAST_LABEL]]), error = function(e) NULL)
mix_threshold <- if (!is.null(mix_fit)) find_mixture_crossing(mix_fit) else NA_real_
log_msg(sprintf("Beta-mixture crossing-point threshold: %s",
                if (is.na(mix_threshold)) "not identifiable (likely unimodal distribution)"
                else sprintf("%.3f", mix_threshold)))

if (!is.null(mix_fit) && !is.na(mix_threshold)) {
   grid <- seq(0.001, 0.999, by = 0.002)
   dens_df <- tibble(
      x = grid,
      component1 = mix_fit$pi[1] * dbeta(grid, mix_fit$params1[1], mix_fit$params1[2]),
      component2 = mix_fit$pi[2] * dbeta(grid, mix_fit$params2[1], mix_fit$params2[2])) %>%
      pivot_longer(c(component1, component2), names_to = "component", values_to = "density")
   
   p_mix <- ggplot() +
      geom_histogram(data = west_qd, aes(x = .data[[EAST_LABEL]], y = after_stat(density)),
                     bins = 25, fill = "grey85", colour = "white") +
      geom_line(data = dens_df, aes(x = x, y = density, colour = component), linewidth = 1) +
      geom_vline(xintercept = mix_threshold, linetype = "dashed", colour = "red") +
      geom_vline(xintercept = IMPORT_THRESH, linetype = "dotted", colour = "black") +
      annotate("text", x = mix_threshold, y = 0, label = "mixture crossing", angle = 90,
               vjust = -0.4, hjust = 0, size = 3, colour = "red") +
      labs(title = "Two-component beta mixture fit to western Q_east",
           subtitle = "Dashed red = data-driven crossing point; dotted black = protocol threshold (0.20)",
           x = "Q_east", y = "Density", colour = NULL) +
      theme_minimal()
   ggsave(file.path(output_dir, "mixture_threshold.pdf"), p_mix, width = 6.5, height = 4.5, dpi = 600)
}

write_tsv(tibble(method = "beta_mixture_crossing", threshold = mix_threshold),
          file.path(output_dir, "mixture_threshold.tsv"))

# ============================================================================
# D. BOOTSTRAP-OVER-SNPs — per-sample Q_east uncertainty from a finite panel
# ============================================================================
bootstrap_qeast_ci <- function(hap, sample_region, REGIONS, sample_ids,
                               east_label, B = 200,
                               min_obs_per_region_snp = 2, min_region_n = 3,
                               seed = 1, n_cores = 1, verbose = TRUE) {
   set.seed(seed)
   L <- nrow(hap)
   boot_idx <- replicate(B, sample.int(L, L, replace = TRUE), simplify = FALSE)
   
   one_boot <- function(idx) {
      fit <- fit_one_config(hap[idx, , drop = FALSE], sample_region, REGIONS,
                            min_obs_per_region_snp, min_region_n,
                            sample_ids = sample_ids)
      fit$Q
   }
   
   if (n_cores > 1) {
      boots <- parallel::mclapply(boot_idx, one_boot, mc.cores = n_cores)
   } else {
      boots <- vector("list", B)
      for (b in seq_len(B)) {
         boots[[b]] <- one_boot(boot_idx[[b]])
         if (verbose && b %% 25 == 0) message(sprintf("  bootstrap %d/%d", b, B))
      }
   }
   boots <- boots[!sapply(boots, is.null)]
   east_mat <- sapply(boots, function(Q) Q[, east_label])
   rownames(east_mat) <- sample_ids
   ci <- t(apply(east_mat, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
   colnames(ci) <- c("q_east_lo95", "q_east_median", "q_east_hi95")
   as_tibble(ci, rownames = "SampleID")
}

if (RUN_BOOTSTRAP) {
   log_msg(sprintf("--- D. Bootstrap-over-SNPs (B=%d) ---", N_BOOT))
   boot_ci <- bootstrap_qeast_ci(raw$hap, raw$sample_region, REGIONS, raw$sample_ids,
                                 EAST_LABEL, B = N_BOOT, n_cores = N_CORES)
   boot_ci <- boot_ci %>%
      left_join(Q_df %>% select(SampleID, region, point_est = all_of(EAST_LABEL)),
                by = "SampleID") %>%
      mutate(ci_excludes_thresh = q_east_lo95 > IMPORT_THRESH |
                q_east_hi95 < IMPORT_THRESH,
             flag_uncertain = !ci_excludes_thresh & region == WEST_LABEL)
   
   write_tsv(boot_ci, file.path(output_dir, "bootstrap_qeast_ci.tsv"))
   
   n_unstable <- sum(boot_ci$flag_uncertain, na.rm = TRUE)
   log_msg(sprintf(
      "Bootstrap: %d western samples have a 95%% CI on Q_east straddling the %.2f threshold (i.e. flag status is not robust to SNP resampling)",
      n_unstable, IMPORT_THRESH))
   
   west_boot <- boot_ci %>% filter(region == WEST_LABEL) %>% arrange(point_est)
   west_boot <- west_boot %>% mutate(SampleID = factor(SampleID, levels = SampleID))
   p_boot <- ggplot(west_boot, aes(x = SampleID, y = point_est)) +
      geom_errorbar(aes(ymin = q_east_lo95, ymax = q_east_hi95), width = 0,
                    colour = "grey50") +
      geom_point(aes(colour = flag_uncertain), size = 1.6) +
      geom_hline(yintercept = IMPORT_THRESH, linetype = "dashed", colour = "red") +
      scale_colour_manual(values = c(`TRUE` = "#D7263D", `FALSE` = "#2E86AB"),
                          labels = c(`TRUE` = "CI straddles threshold",
                                     `FALSE` = "stable call"),
                          name = NULL) +
      labs(title = "Bootstrap 95% CI on Q_east, western samples",
           subtitle = sprintf("%d/%d samples ordered by point estimate; uncertainty from SNP-panel resampling",
                              nrow(west_boot), nrow(west_boot)),
           x = NULL, y = "Q_east") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
   
   ggsave(file.path(output_dir, "bootstrap_qeast_ci.pdf"), p_boot, width = 8, height = 4.5, dpi = 600)
} else {
   log_msg("--- D. Bootstrap-over-SNPs: SKIPPED (run_bootstrap=FALSE) ---")
   boot_ci <- NULL
}

# ============================================================================
# E. PARAMETER SWEEP — MIN_REGION_N x min_obs_per_region_snp
# ============================================================================
log_msg("--- E. Parameter sweep: MIN_REGION_N x min_obs_per_region_snp ---")

parameter_sweep <- function(hap, sample_region, REGIONS, sample_ids,
                            west_label, east_label,
                            min_region_n_grid = c(3, 5, 8, 10, 15, 20),
                            min_obs_grid = c(2, 3, 5, 10),
                            threshold = 0.20) {
   configs <- expand.grid(min_region_n = min_region_n_grid,
                          min_obs_per_region_snp = min_obs_grid)
   pmap_dfr(configs, function(min_region_n, min_obs_per_region_snp) {
      fit <- fit_one_config(hap, sample_region, REGIONS,
                            min_obs_per_region_snp, min_region_n,
                            sample_ids = sample_ids)
      if (is.null(fit$Q)) {
         return(tibble(min_region_n = min_region_n,
                       min_obs_per_region_snp = min_obs_per_region_snp,
                       n_snps_used = fit$n_snps_used, n_flagged = NA_integer_,
                       n_west = NA_integer_, error = fit$error))
      }
      west_mask <- sample_region == west_label
      tibble(min_region_n = min_region_n,
             min_obs_per_region_snp = min_obs_per_region_snp,
             n_snps_used = fit$n_snps_used,
             n_flagged = sum(fit$Q[west_mask, east_label] > threshold),
             n_west = sum(west_mask), error = NA_character_)
   })
}

sweep_param_df <- parameter_sweep(raw$hap, raw$sample_region, REGIONS, raw$sample_ids,
                                  WEST_LABEL, EAST_LABEL, threshold = IMPORT_THRESH)

write_tsv(sweep_param_df, file.path(output_dir, "parameter_sweep.tsv"))
log_msg(sprintf("Parameter sweep: %d/%d configs produced a usable fit",
                sum(!is.na(sweep_param_df$n_flagged)), nrow(sweep_param_df)))

p_param <- ggplot(sweep_param_df %>% filter(!is.na(n_flagged)),
                  aes(x = factor(min_obs_per_region_snp), y = factor(min_region_n),
                      fill = n_flagged)) +
   geom_tile() +
   geom_text(aes(label = n_flagged), colour = "white", size = 3) +
   scale_fill_viridis_c(name = "N flagged\n(western)") +
   labs(title = "Importation calls across MIN_REGION_N x SNP-coverage filter",
        subtitle = sprintf("Fixed at Q_east > %.2f; cells with NA = config infeasible", IMPORT_THRESH),
        x = "min_obs_per_region_snp", y = "min_region_n (LOO floor)") +
   theme_minimal()

ggsave(file.path(output_dir, "parameter_sweep.pdf"), p_param, width = 6, height = 4.5, dpi = 600)

# ============================================================================
# Combined reviewer-facing summary
# ============================================================================
log_msg("--- Summary: candidate thresholds across methods ---")

thresh_summary <- tibble(
   method = c("Protocol (a priori)",
              if (RUN_PERMUTATION) "Permutation null, 95th pct" else NA,
              if (RUN_PERMUTATION) "Permutation null, 99th pct" else NA,
              "Beta-mixture crossing point"),
   threshold = c(IMPORT_THRESH,
                 if (RUN_PERMUTATION) perm_res$q95 else NA,
                 if (RUN_PERMUTATION) perm_res$q99 else NA,
                 mix_threshold)) %>%
   filter(!is.na(method))
log_msg(paste(capture.output(print(thresh_summary)), collapse = "\n"))

write_tsv(thresh_summary, file.path(output_dir, "threshold_comparison_summary.tsv"))

saveRDS(list(Q_df_baseline = Q_df, threshold_sweep = sweep_df,
             permutation = perm_res, mixture_fit = mix_fit,
             mixture_threshold = mix_threshold, bootstrap_ci = boot_ci,
             parameter_sweep = sweep_param_df, thresh_summary = thresh_summary),
        file.path(output_dir, "phase3_01b_sensitivity_results.rds"))

writeLines(log_lines, file.path(output_dir, "phase3_01b_summary.txt"))
message("\nDone. Outputs in ", output_dir)

