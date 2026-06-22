# =============================================================================
# Phase 3 — Step 3.1b: Reusable core functions for sensitivity analysis
# =============================================================================
# Companion LIBRARY to phase3_01_supervised_admixture.R. Source this file
# from phase3_01b_threshold_sensitivity.R (or any other script).
#
# Why this exists:
#   The original script's data-loading, frequency-computation, and
#   per-sample LOO-EM steps are written as one linear block with hardcoded
#   parameters (IMPORT_THRESH, MIN_REGION_N, the implicit min-observations-
#   per-region-SNP filter, FREQ_EPS, EM tolerances). A sensitivity analysis
#   needs to re-run that same pipeline dozens to thousands of times under
#   different parameter values, different permuted region labels (for a
#   null distribution), and different resampled SNP sets (for bootstrap
#   uncertainty) — without copy-pasting the pipeline body each time.
#
# NAMING CONVENTION (see note in the calling script): all functions below
# take REGIONS as an explicit argument and use it directly for indexing —
# they never hardcode "west"/"east" in lowercase, so they are robust
# regardless of how you resolve the case-mismatch in the original script's
# sections 4-5.
# =============================================================================

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

suppressPackageStartupMessages({
   library(tidyverse)
})

# ----- Paths ----------------------------------------------------------------
geno_in    <- geno_arg
phase22_in <- phase22_arg
meta_in    <- meta_arg
cent_in    <- cent_arg

# -----------------------------------------------------------------------------
# 1. Raw data loading (mirrors phase3_01 section 1) — NO filtering applied.
#    Filtering (min_obs_per_region_snp, min_region_n) is deferred to
#    fit_one_config() so that the same loaded matrix can be re-filtered
#    under many different parameter configurations without re-reading disk.
# -----------------------------------------------------------------------------
load_phase3_01_raw_inputs <- function(geno_in, phase22_in, meta_in, REGIONS) {
   
   phase22 <- readRDS(phase22_in)
   sample_ids <- phase22$eigenvec$SampleID
   n <- length(sample_ids)
   
   geno_lines <- readLines(geno_in)
   L <- length(geno_lines)
   geno_pseudo <- matrix(0L, nrow = L, ncol = n)
   for (i in seq_len(L)) {
      chars <- strsplit(geno_lines[i], "", fixed = TRUE)[[1]]
      stopifnot(length(chars) == n)
      geno_pseudo[i, ] <- as.integer(chars)
   }
   
   hap <- geno_pseudo
   hap[hap == 0L] <- 0L
   hap[hap == 2L] <- 1L
   hap[geno_pseudo == 9L | geno_pseudo == 1L] <- NA_integer_
   storage.mode(hap) <- "integer"
   colnames(hap) <- sample_ids
   
   meta <- readRDS(meta_in)
   meta_aln <- meta %>% filter(SampleID %in% sample_ids) %>%
      arrange(match(SampleID, sample_ids))
   stopifnot(all(meta_aln$SampleID == sample_ids))
   sample_region <- as.character(meta_aln$region)
   stopifnot(all(sample_region %in% REGIONS))
   
   region_n <- table(factor(sample_region, levels = REGIONS))
   message("Region sample sizes (full data): ",
           paste(sprintf("%s=%d", names(region_n), region_n), collapse = ", "))
   message(sprintf("Loaded .geno: %d SNPs x %d samples (no filtering applied yet)",
                   L, n))
   
   list(hap = hap, sample_ids = sample_ids, sample_region = sample_region,
        region_n = region_n, meta_aln = meta_aln, L_total = L)
}

# -----------------------------------------------------------------------------
# 2. Per-region allele counts/frequencies + SNP coverage filter
#    (mirrors phase3_01 section 2; min_obs_per_region_snp is now a parameter,
#    not a hardcoded "2")
# -----------------------------------------------------------------------------
compute_region_freqs <- function(hap, sample_region, REGIONS,
                                 min_obs_per_region_snp = 2) {
   L <- nrow(hap)
   sum_rl <- matrix(0L, nrow = L, ncol = length(REGIONS))
   n_rl   <- matrix(0L, nrow = L, ncol = length(REGIONS))
   for (r_idx in seq_along(REGIONS)) {
      sub <- hap[, sample_region == REGIONS[r_idx], drop = FALSE]
      sum_rl[, r_idx] <- rowSums(sub, na.rm = TRUE)
      n_rl[, r_idx]   <- rowSums(!is.na(sub))
   }
   F_full <- sum_rl / pmax(n_rl, 1L)
   F_full[n_rl == 0L] <- NA_real_
   ok_snp <- rowSums(n_rl >= min_obs_per_region_snp) == length(REGIONS)
   list(sum_rl = sum_rl, n_rl = n_rl, F_full = F_full, ok_snp = ok_snp)
}

# -----------------------------------------------------------------------------
# 3. EM core — identical maths to the original em_fit(); FREQ_EPS and the
#    convergence tolerances are now explicit arguments instead of globals.
# -----------------------------------------------------------------------------
clip_freqs <- function(F, freq_eps = 1e-6) pmax(pmin(F, 1 - freq_eps), freq_eps)

em_fit <- function(a, F, em_max_iter = 5000, em_tol_q = 1e-8, em_tol_ll = 1e-8) {
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
   it <- 0L
   for (it in seq_len(em_max_iter)) {
      pa <- as.vector(F_alt %*% Q)
      pr <- as.vector(F_ref %*% Q)
      ll <- sum(log(pa)) + sum(log(pr))
      gam_alt <- sweep(F_alt, 2, Q, `*`) / pa
      gam_ref <- sweep(F_ref, 2, Q, `*`) / pr
      Q_new <- (colSums(gam_alt) + colSums(gam_ref)) / nL
      Q_new <- Q_new / sum(Q_new)
      step <- max(abs(Q_new - Q))
      ll_rel <- if (it == 1) Inf else abs(ll - ll_prev) / max(1, abs(ll_prev))
      Q <- Q_new
      ll_prev <- ll
      if (ll_rel < em_tol_ll || step < em_tol_q) { converged <- TRUE; break }
   }
   list(Q = Q, iter = it, ll = ll_prev, step = step, converged = converged)
}

# -----------------------------------------------------------------------------
# 4. Per-sample LOO + EM across all samples (mirrors phase3_01 section 3)
# -----------------------------------------------------------------------------
fit_loo_admixture <- function(hap, sample_region, REGIONS, sum_rl, n_rl, F_full,
                              freq_eps = 1e-6, em_max_iter = 5000,
                              em_tol_q = 1e-8, em_tol_ll = 1e-8,
                              sample_ids = NULL) {
   n <- ncol(hap)
   if (is.null(sample_ids)) sample_ids <- colnames(hap)
   if (is.null(sample_ids)) sample_ids <- as.character(seq_len(n))
   
   Q_all <- matrix(NA_real_, nrow = n, ncol = length(REGIONS),
                   dimnames = list(sample_ids, REGIONS))
   iters <- integer(n); final_step <- numeric(n); converged <- logical(n)
   
   for (i in seq_len(n)) {
      r_i <- match(sample_region[i], REGIONS)
      F_loo <- F_full
      mask <- !is.na(hap[, i])
      if (any(mask)) {
         new_sum <- sum_rl[mask, r_i] - hap[mask, i]
         new_n   <- n_rl[mask, r_i] - 1L
         safe <- new_n > 0L
         F_loo[mask, r_i][safe] <- new_sum[safe] / new_n[safe]
      }
      F_loo_clip <- clip_freqs(F_loo, freq_eps)
      res <- em_fit(hap[, i], F_loo_clip, em_max_iter, em_tol_q, em_tol_ll)
      Q_all[i, ]    <- res$Q
      iters[i]      <- res$iter
      final_step[i] <- res$step
      converged[i]  <- res$converged
   }
   list(Q = Q_all, iters = iters, final_step = final_step, converged = converged)
}

# -----------------------------------------------------------------------------
# 5. ONE-CALL WRAPPER — the function every sensitivity analysis below calls
#    repeatedly with different arguments. Given an already-loaded hap matrix
#    and region labels, runs freq computation -> SNP filter -> LOO-EM under
#    any parameter configuration, and fails soft (returns an error string,
#    not a stop()) if a config is statistically unsafe — important because
#    a parameter sweep WILL hit some invalid configs by design.
# -----------------------------------------------------------------------------
fit_one_config <- function(hap, sample_region, REGIONS,
                           min_obs_per_region_snp = 2,
                           min_region_n = 3,
                           freq_eps = 1e-6,
                           em_max_iter = 5000,
                           em_tol_q = 1e-8,
                           em_tol_ll = 1e-8,
                           sample_ids = NULL) {
   region_n <- table(factor(sample_region, levels = REGIONS))
   if (any(region_n < min_region_n)) {
      return(list(Q = NULL, n_snps_used = NA_integer_, error = sprintf(
         "region_n below min_region_n=%d: %s", min_region_n,
         paste(sprintf("%s=%d", names(region_n), region_n), collapse = ", "))))
   }
   freqs <- compute_region_freqs(hap, sample_region, REGIONS,
                                 min_obs_per_region_snp)
   if (sum(freqs$ok_snp) < 10) {
      return(list(Q = NULL, n_snps_used = sum(freqs$ok_snp), error = sprintf(
         "only %d SNPs pass min_obs_per_region_snp=%d; too few to fit",
         sum(freqs$ok_snp), min_obs_per_region_snp)))
   }
   hap_sub <- hap[freqs$ok_snp, , drop = FALSE]
   fit <- fit_loo_admixture(hap_sub, sample_region, REGIONS,
                            freqs$sum_rl[freqs$ok_snp, , drop = FALSE],
                            freqs$n_rl[freqs$ok_snp, , drop = FALSE],
                            freqs$F_full[freqs$ok_snp, , drop = FALSE],
                            freq_eps, em_max_iter, em_tol_q, em_tol_ll,
                            sample_ids)
   c(fit, list(n_snps_used = sum(freqs$ok_snp), error = NA_character_))
}

# -----------------------------------------------------------------------------
# 6. Attach metadata to a fitted Q matrix (mirrors phase3_01 lines 242-249)
# -----------------------------------------------------------------------------
build_Q_df <- function(Q, meta_aln, REGIONS,
                       season_levels = c("early", "peak", "late")) {
   as_tibble(Q, rownames = "SampleID") %>%
      left_join(meta_aln %>% select(SampleID, region, Location, VillageCode,
                                    transmission_season, collection_month,
                                    latitude, longitude),
                by = "SampleID") %>%
      mutate(region = factor(region, levels = REGIONS),
             transmission_season = factor(transmission_season,
                                          levels = season_levels))
}

