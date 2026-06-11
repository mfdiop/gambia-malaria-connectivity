# =============================================================================
# Phase 2 — Step 2.3: Baseline relatedness landscape (hmmIBD)
# =============================================================================
# Per notes/objective3_proposal.md §6 Phase 2 Step 2.3:
#   - Run hmmIBD on all pairwise sample combinations to produce a genome-wide
#     pairwise IBD matrix.
#   - Classify pairs:
#       Highly related   (fract_sites_IBD > 0.5)    candidate recent transmission
#       Moderately rel.  (0.125 <= IBD <= 0.5)      2nd-4th degree relatives
#       Unrelated        (IBD < 0.125)              background
#   - Visualise the spatial distribution of highly related pairs:
#       chord diagram (region-level)  and  geographic network (village-level).
#   - Question: are high-IBD pairs spatially clustered within regions, or do
#     cross-regional high-IBD pairs form a directional pattern?
#
# Sample set : monoclonal subset (results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz)
#              hmmIBD assumes haploid genotypes; polyclonal samples re-enter
#              downstream via dcifer (Phase 3.2 / 4.1).
#
# hmmIBD encoding (haploid):
#   "0/0"            -> 0  (ref)
#   "1/1"            -> 1  (alt)
#   "0/1" / "1/0"    -> -1 (het in monoclonal: treated as missing)
#   "./." / NA / "." -> -1
#
# Chromosomes: nuclear only, Pf3D7_01_v3..Pf3D7_14_v3, renumbered 1..14 for
# hmmIBD (which expects integer chromosomes).
#
# Output : results/phase2_03/
#            hmmibd_input.txt                 hmmIBD genotype input
#            hmmibd.hmm.txt                   per-segment IBD calls (hmmIBD)
#            hmmibd.hmm_fract.txt             per-pair summary (hmmIBD)
#            pair_ibd_classified.tsv          all pairs + class + region info
#            ibd_matrix.tsv                   144 x 144 fract_sites_IBD matrix
#            ibd_distribution.pdf             histogram of pairwise IBD
#            ibd_within_vs_cross_region.pdf   boxplot
#            ibd_chord_region.pdf             chord diagram of high-IBD pairs
#            ibd_geographic_network.pdf       village-level high-IBD network
#            phase2_03_summary.txt            human-readable log
#            phase2_03_results.rds            persisted state
#
# Tools  : hmmIBD (built into tools/hmmIBD/hmmIBD), R packages vcfR, tidyverse,
#          circlize, ggrepel.
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir   <- .parse_arg("output_dir",   "results/phase2_03")
mono_vcf     <- .parse_arg("mono_vcf",     "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
meta_arg     <- .parse_arg("meta",         "results/phase1_03/sample_master_metadata.rds")
centroid_arg <- .parse_arg("centroid",     "results/phase1_03/village_centroids.tsv")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(vcfR)
   library(circlize)
   library(ggrepel)
})

# ----- Paths ----------------------------------------------------------------
vcf_in       <- mono_vcf
meta_in      <- meta_arg
centroid_in  <- centroid_arg

out_dir      <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hmmibd_bin   <- Sys.getenv("HMMIBD", "tools/hmmIBD/hmmIBD")
stopifnot(file.exists(hmmibd_bin))

input_path   <- file.path(out_dir, "hmmibd_input.txt")
out_prefix   <- file.path(out_dir, "hmmibd")

# ----- Parameters -----------------------------------------------------------
IBD_HIGH    <- 0.5     # proposal
IBD_MOD_LO  <- 0.125   # proposal
MAX_ITER    <- 5       # hmmIBD default
PF_CHRS     <- sprintf("Pf3D7_%02d_v3", 1:14)
SEED        <- 42
set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load VCF + metadata; restrict to nuclear chromosomes
# ============================================================================
log_msg("=== Phase 2 Step 2.3: hmmIBD baseline relatedness ===")
vcf  <- read.vcfR(vcf_in, verbose = FALSE)
meta <- readRDS(meta_in)
vcf_samples <- colnames(vcf@gt)[-1]
log_msg(sprintf("Monoclonal VCF samples: %d", length(vcf_samples)))

nuc <- vcf@fix[, "CHROM"] %in% PF_CHRS
vcf <- vcf[nuc, ]
log_msg(sprintf("Nuclear SNPs (chrs 1..14): %d", nrow(vcf@fix)))

# ============================================================================
# 2. Encode haploid genotypes for hmmIBD
# ============================================================================
gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)   # SNPs x samples
encode_hmm <- function(x) {
   out <- rep(-1L, length(x))
   out[x %in% c("0", "0/0", "0|0")] <- 0L
   out[x %in% c("1", "1/1", "1|1")] <- 1L
   out
}
geno <- matrix(encode_hmm(as.vector(gt_chr)),
               nrow = nrow(gt_chr), ncol = ncol(gt_chr),
               dimnames = dimnames(gt_chr))
storage.mode(geno) <- "integer"

miss <- mean(geno == -1L)
log_msg(sprintf("Genotype matrix: %d SNPs x %d samples (missing/het %.2f%%)",
                nrow(geno), ncol(geno), 100 * miss))
stopifnot(miss < 0.5)   # encoding sanity gate

# ============================================================================
# 3. Build sorted hmmIBD input (chrom int, pos int, then samples)
# ============================================================================
chrom_int <- as.integer(factor(vcf@fix[, "CHROM"], levels = PF_CHRS))
pos_int   <- as.integer(vcf@fix[, "POS"])
ord <- order(chrom_int, pos_int)

hmm_df <- as_tibble(geno[ord, , drop = FALSE]) %>%
   mutate(chrom = chrom_int[ord], pos = pos_int[ord], .before = 1)

write_tsv(hmm_df, input_path)
log_msg(sprintf("Wrote hmmIBD input: %s (%d SNPs)", input_path, nrow(hmm_df)))

# ============================================================================
# 4. Run hmmIBD (single population, all pairs)
# ============================================================================
log_msg(sprintf("Running %s ...", hmmibd_bin))
run_log_path <- file.path(out_dir, "hmmibd_run.log")
res <- system2(hmmibd_bin,
               args = c("-i", input_path,
                        "-o", out_prefix,
                        "-m", MAX_ITER),
               stdout = TRUE, stderr = TRUE)
writeLines(res, run_log_path)
status <- attr(res, "status")
if (!is.null(status) && status != 0) {
   writeLines(res); stop(sprintf("hmmIBD failed (exit %d)", status))
}
log_msg(sprintf("hmmIBD finished. Log: %s", run_log_path))

# ============================================================================
# 5. Parse per-pair output + classify
# ============================================================================
fract_path <- paste0(out_prefix, ".hmm_fract.txt")
stopifnot(file.exists(fract_path))
pair_df <- read_tsv(fract_path, show_col_types = FALSE)
log_msg(sprintf("Pairs analysed: %d (expected n*(n-1)/2 = %d)",
                nrow(pair_df), choose(length(vcf_samples), 2)))

meta_min <- meta %>% select(SampleID, region, Location, VillageCode,
                            latitude, longitude)

pair_df <- pair_df %>%
   rename(SampleID1 = sample1, SampleID2 = sample2) %>%
   left_join(meta_min %>% rename_with(~ paste0(., "1")),
             by = "SampleID1") %>%
   left_join(meta_min %>% rename_with(~ paste0(., "2")),
             by = "SampleID2") %>%
   mutate(
      region1_chr = as.character(region1),
      region2_chr = as.character(region2),
      ibd_class = case_when(
         fract_sites_IBD >  IBD_HIGH                                  ~ "high",
         fract_sites_IBD >= IBD_MOD_LO & fract_sites_IBD <= IBD_HIGH ~ "moderate",
         TRUE                                                         ~ "unrelated"),
      pair_type = ifelse(region1_chr == region2_chr, "within-region", "cross-region"),
      region_pair = ifelse(region1_chr < region2_chr,
                           paste(region1_chr, region2_chr, sep = " <-> "),
                           paste(region2_chr, region1_chr, sep = " <-> ")))

write_tsv(pair_df, file.path(out_dir, "pair_ibd_classified.tsv"))

n_high <- sum(pair_df$ibd_class == "high")
n_mod  <- sum(pair_df$ibd_class == "moderate")
log_msg(sprintf("Pair classification: high = %d, moderate = %d, unrelated = %d",
                n_high, n_mod, sum(pair_df$ibd_class == "unrelated")))
log_msg(sprintf("Median fract_sites_IBD: %.4f",
                median(pair_df$fract_sites_IBD, na.rm = TRUE)))

# ============================================================================
# 6. Build symmetric IBD matrix
# ============================================================================
samp <- sort(unique(c(pair_df$SampleID1, pair_df$SampleID2)))
ibd_mat <- matrix(NA_real_, nrow = length(samp), ncol = length(samp),
                  dimnames = list(samp, samp))
diag(ibd_mat) <- 1
ibd_mat[cbind(match(pair_df$SampleID1, samp),
              match(pair_df$SampleID2, samp))] <- pair_df$fract_sites_IBD
ibd_mat[cbind(match(pair_df$SampleID2, samp),
              match(pair_df$SampleID1, samp))] <- pair_df$fract_sites_IBD
write_tsv(as_tibble(ibd_mat, rownames = "SampleID"),
          file.path(out_dir, "ibd_matrix.tsv"))
# Standalone RDS so downstream rules can depend on a single declared file
saveRDS(ibd_mat, file.path(out_dir, "ibd_matrix.rds"))

# ============================================================================
# 7. Distribution + within vs cross-region IBD
# ============================================================================
p_hist <- ggplot(pair_df, aes(x = fract_sites_IBD)) +
   geom_histogram(bins = 80, fill = "steelblue", colour = "white") +
   geom_vline(xintercept = c(IBD_MOD_LO, IBD_HIGH),
              linetype = "dashed", colour = "red") +
   annotate("text", x = IBD_MOD_LO + 0.005, y = Inf,
            label = "moderate (>= 0.125)", hjust = 0, vjust = 1.4,
            size = 3, colour = "red") +
   annotate("text", x = IBD_HIGH + 0.005, y = Inf,
            label = "high (> 0.5)", hjust = 0, vjust = 1.4,
            size = 3, colour = "red") +
   labs(title = "Pairwise IBD distribution (hmmIBD, monoclonal)",
        x = "fract_sites_IBD", y = "Number of pairs") +
   theme_minimal()
ggsave(file.path(out_dir, "ibd_distribution.pdf"), p_hist,
       width = 7.5, height = 4.5, dpi = 600)

p_box <- ggplot(pair_df, aes(x = pair_type, y = fract_sites_IBD,
                              fill = pair_type)) +
   geom_violin(alpha = 0.4) +
   geom_boxplot(width = 0.18, outlier.alpha = 0.3) +
   scale_fill_manual(values = c("within-region" = "steelblue",
                                "cross-region"  = "tomato"),
                     guide = "none") +
   geom_hline(yintercept = c(IBD_MOD_LO, IBD_HIGH),
              linetype = "dashed", colour = "grey40") +
   labs(title = "Within-region vs cross-region pairwise IBD",
        x = NULL, y = "fract_sites_IBD") +
   theme_minimal()
ggsave(file.path(out_dir, "ibd_within_vs_cross_region.pdf"), p_box,
       width = 6, height = 4.5, dpi = 600)

# Wilcoxon test for the cross-vs-within shift
wt <- wilcox.test(fract_sites_IBD ~ pair_type, data = pair_df)
log_msg(sprintf("Wilcoxon within vs cross-region: W = %.0f, p = %.3g",
                wt$statistic, wt$p.value))

# ============================================================================
# 8. Chord diagram: region-level adjacency of high-IBD pairs
# ============================================================================
high_pairs <- pair_df %>% filter(ibd_class == "high")
log_msg(sprintf("High-IBD pairs (n = %d): within = %d, cross-region = %d",
                nrow(high_pairs),
                sum(high_pairs$pair_type == "within-region"),
                sum(high_pairs$pair_type == "cross-region")))

if (nrow(high_pairs) > 0) {
   region_adj <- high_pairs %>%
      mutate(r1 = pmin(region1_chr, region2_chr),
             r2 = pmax(region1_chr, region2_chr)) %>%
      count(r1, r2, name = "n_pairs")

   pdf(file.path(out_dir, "ibd_chord_region.pdf"), width = 6, height = 6)
   circos.clear()
   chordDiagram(
      region_adj %>% select(r1, r2, n_pairs),
      grid.col = c(west = "#D7263D", central = "#F46036",
                   east = "#2E86AB"),
      annotationTrack = c("name", "grid"),
      preAllocateTracks = list(track.height = 0.04))
   title("High-IBD pairs (fract_sites_IBD > 0.5) — region chord")
   dev.off()
} else {
   log_msg("No high-IBD pairs — chord diagram skipped.")
   file.create(file.path(out_dir, "ibd_chord_region.pdf"))
}

# ============================================================================
# 9. Geographic network of high-IBD pairs (village centroids)
# ============================================================================
cent_raw <- read_tsv(centroid_in, show_col_types = FALSE)
# Dedupe to one centroid per Location (the centroid file has one row per
# VillageCode, and several Locations aggregate multiple VillageCodes — see
# phase1_03_sample_annotation.R). Within-Location centroids are averaged
# weighted by sample count.
cent <- cent_raw %>%
   group_by(Location, region) %>%
   summarise(lon = weighted.mean(lon, n),
             lat = weighted.mean(lat, n),
             n   = sum(n), .groups = "drop")

if (nrow(high_pairs) > 0) {
   # Within-village high-IBD pairs (loc_a == loc_b): annotate node, no edge
   within_village_counts <- high_pairs %>%
      filter(Location1 == Location2) %>%
      count(Location1, name = "n_within_village") %>%
      rename(Location = Location1)

   edges <- high_pairs %>%
      filter(Location1 != Location2) %>%
      mutate(loc_a = pmin(Location1, Location2),
             loc_b = pmax(Location1, Location2)) %>%
      count(loc_a, loc_b, name = "n_high_pairs") %>%
      left_join(cent %>% select(Location, lon_a = lon, lat_a = lat),
                by = c("loc_a" = "Location")) %>%
      left_join(cent %>% select(Location, lon_b = lon, lat_b = lat),
                by = c("loc_b" = "Location"))

   nodes <- cent %>%
      left_join(within_village_counts, by = "Location") %>%
      mutate(n_within_village = replace_na(n_within_village, 0L))

   p_geo <- ggplot() +
      geom_segment(data = edges,
                   aes(x = lon_a, y = lat_a, xend = lon_b, yend = lat_b,
                       linewidth = n_high_pairs),
                   colour = "tomato", alpha = 0.7, lineend = "round") +
      geom_point(data = nodes, aes(x = lon, y = lat, size = n,
                                   fill = region),
                 shape = 21, colour = "black") +
      geom_text_repel(data = nodes,
                      aes(x = lon, y = lat,
                          label = sprintf("%s\n(n=%d, %d within-village)",
                                          Location, n, n_within_village)),
                      size = 3, max.overlaps = 20) +
      scale_size_continuous(name = "n samples", range = c(3, 10)) +
      scale_linewidth_continuous(name = "high-IBD pairs", range = c(0.3, 3)) +
      scale_fill_manual(values = c(west = "#D7263D", central = "#F46036",
                                   east = "#2E86AB")) +
      coord_quickmap() +
      labs(title = "Geographic network of high-IBD pairs",
           subtitle = sprintf("Edges link village pairs with >=1 hmmIBD pair > %.2f",
                              IBD_HIGH),
           x = "Longitude", y = "Latitude") +
      theme_minimal()
   ggsave(file.path(out_dir, "ibd_geographic_network.pdf"), p_geo,
          width = 10, height = 7, dpi = 600)
} else {
   log_msg("No high-IBD pairs — geographic network skipped.")
   file.create(file.path(out_dir, "ibd_geographic_network.pdf"))
}

# ============================================================================
# 10. Persist + log
# ============================================================================
saveRDS(list(pair_df    = pair_df,
             ibd_matrix = ibd_mat,
             n_high     = n_high,
             n_moderate = n_mod,
             thresholds = c(high = IBD_HIGH, mod_lo = IBD_MOD_LO),
             wilcoxon   = list(W = unname(wt$statistic), p = wt$p.value)),
        file.path(out_dir, "phase2_03_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase2_03_summary.txt"))
message("\nDone. Outputs in ", out_dir)
