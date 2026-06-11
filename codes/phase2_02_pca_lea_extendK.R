# =============================================================================
# Phase 2 — Step 2.2 (LEA variant) — cross-K aggregation
# =============================================================================
# Aggregates per-K sNMF cross-entropy values produced by phase2_02_pca_lea.R
# (one job per K) and produces the CE-vs-K curve + K summary table. This is
# Snakemake plan Task 5 Pattern A — no snmfProject I/O here.
#
# Input  : --ce_files=path1,path2,...   per-K ce.tsv files
#          --output_dir=...             where to write aggregated outputs
#
# Output : <output_dir>/ce_vs_K.pdf    cross-entropy curve
#          <output_dir>/K_summary.tsv  best run + min CE per K
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}

output_dir   <- .parse_arg("output_dir", "results/phase2_02_lea")
ce_files_str <- .parse_arg("ce_files",   "")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
})

ce_files <- if (nzchar(ce_files_str)) {
   strsplit(ce_files_str, ",", fixed = TRUE)[[1]]
} else {
   # interactive default: glob per-K dirs
   Sys.glob(file.path(output_dir, "K*", "ce.tsv"))
}

stopifnot(length(ce_files) > 0)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read each per-K ce.tsv (K column present, written by phase2_02_pca_lea.R)
ce_df <- map_dfr(ce_files, ~ read_tsv(.x, show_col_types = FALSE))
stopifnot(all(c("K", "run", "cross_entropy") %in% names(ce_df)))

best_runs <- ce_df %>%
   group_by(K) %>%
   slice_min(cross_entropy, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   arrange(K)

optimal_K <- best_runs$K[which.min(best_runs$cross_entropy)]

write_tsv(best_runs, file.path(output_dir, "K_summary.tsv"))

message(sprintf("Aggregated K = %s (optimal = %d, CE = %.4f, run %d)",
                paste(range(best_runs$K), collapse = ".."),
                optimal_K,
                best_runs$cross_entropy[best_runs$K == optimal_K],
                best_runs$run[best_runs$K == optimal_K]))

p_ce <- ggplot(best_runs, aes(x = K, y = cross_entropy)) +
   geom_line() + geom_point(size = 3) +
   geom_point(data = best_runs %>% filter(K == optimal_K),
              colour = "red", size = 4) +
   geom_jitter(data = ce_df, width = 0.12, alpha = 0.4, 
               size = 1.6, colour = "grey40") +
   scale_x_continuous(breaks = sort(unique(best_runs$K))) +
   labs(title = sprintf("sNMF cross-entropy, K = %d..%d (optimal K = %d)",
                        min(best_runs$K), max(best_runs$K), optimal_K),
        subtitle = "Line = best run per K; points = all runs",
        x = "K", y = "Cross-entropy") +
   theme_minimal()

ggsave(file.path(output_dir, "ce_vs_K.pdf"), p_ce,
       width = 7, height = 4, dpi = 600)

message("Done. Aggregated outputs in ", output_dir)
