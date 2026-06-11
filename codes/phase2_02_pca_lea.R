# =============================================================================
# Phase 2 — Step 2.2 (LEA variant): Population structure (PCA + sNMF)
# =============================================================================
# In-R alternative to the PLINK2 + ADMIXTURE pipeline in
# phase2_02_pca_admixture.R. Same scientific question (proposal §6 Step 2.2):
#   - PCA on LD-pruned SNPs (r2 < 0.1 in 50 kb windows)
#   - Ancestry estimation for K = 2..8, optimal K from cross-entropy
#   - Question: discrete west/east clusters or continuous gradient?
#
# Why an R-native variant:
#   - Avoids the PLINK2 / ADMIXTURE binary dependency.
#   - LEA::snmf handles missing data natively, no imputation step.
#   - SNPRelate is already in proposal §7 (PCA / FST) and pulls the LD-prune
#     and PCA steps into the same session that runs the ancestry inference.
#
# Sample set : monoclonal subset (results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz)
#              sNMF, like ADMIXTURE, assumes unambiguous genotypes. Polyclonal
#              samples re-enter downstream via dcifer.
#
# Encoding for haploid P. falciparum:
#   The .geno format LEA expects is dosage in {0, 1, 2, 9}. Haploid calls are
#   mapped 0 -> 0 (ref) and 1 -> 2 (alt), i.e. pseudo-homozygous-diploid. This
#   is the standard trick for running diploid-model ancestry tools on Pf and
#   mirrors what PLINK does when it reads haploid genotypes for ADMIXTURE.
#   sNMF is then called with ploidy = 2; the loss surface is identical up to a
#   factor of 2 in allele-count units, so the inferred Q matrices are directly
#   comparable to the ADMIXTURE output.
#
# Chromosomes: restricted to Pf3D7_01_v3..Pf3D7_14_v3 (nuclear). Apicoplast
# and mitochondrial scaffolds dropped before genotype-matrix construction.
#
# Output : results/phase2_02_lea/
#            pca_eigenvec.tsv                  sample x PC1..PC10
#            pca_eigenval.tsv                  eigenvalues + variance explained
#            pca_pc1_pc2.pdf                   PC1 vs PC2 by region
#            pca_pc1_pc2_by_village.pdf        coloured by village
#            pca_pc1_pc2_by_season.pdf         coloured by transmission season
#            pca_scree.pdf                     scree
#            snmf/                             LEA snmfProject + .geno + .snmf
#            snmf_cross_entropy.tsv            cross-entropy per (K, repetition)
#            snmf_cross_entropy.pdf            CE vs K curve (best run per K)
#            snmf_Q_K{2..8}.tsv                ancestry proportions (best run)
#            snmf_barplot_K{2..8}.pdf          stacked-bar ancestry per K
#            snmf_barplot_optimalK.pdf         bar plot at min-CE K
#            phase2_02_lea_summary.txt         human-readable log
#
# Tools  : R only — vcfR, SNPRelate, LEA, tidyverse, scales, RColorBrewer
# =============================================================================

# --- snakemake interface (no-op when sourced interactively) -----------------
.args <- commandArgs(trailingOnly = TRUE)
.parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}
output_dir  <- .parse_arg("output_dir", "results/phase2_02_lea/K3")
k_single    <- as.integer(.parse_arg("k", 3))
# k_multiple  <- as.integer(.parse_arg("k", 2:12))
mono_vcf    <- .parse_arg("mono_vcf", "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz")
meta_arg    <- .parse_arg("meta",     "results/phase1_03/sample_master_metadata.rds")
# --- end snakemake interface ------------------------------------------------

suppressPackageStartupMessages({
   library(tidyverse)
   library(vcfR)
   library(SNPRelate)
   library(LEA)
   library(scales)
   library(RColorBrewer)
})

theme_fig <- function() {
   theme_minimal() +
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

# ----- Paths ----------------------------------------------------------------
vcf_in   <- mono_vcf
meta_in  <- meta_arg

out_dir   <- output_dir
snmf_dir  <- file.path(out_dir, "snmf")
dir.create(snmf_dir, recursive = TRUE, showWarnings = FALSE)

gds_path  <- file.path(out_dir, "gambia_mono.gds")
geno_path <- file.path(snmf_dir, "gambia_mono_pruned.geno")

# ----- Parameters -----------------------------------------------------------
LD_WINDOW_BP <- 50000           # 50 kb window (proposal)
LD_R2        <- 0.1             # r^2 threshold (proposal)
LD_R_THRESH  <- sqrt(LD_R2)     # SNPRelate ld.threshold is |r|, not r^2
K_RANGE      <- k_single        # single K per Snakemake job
N_REPS       <- 10              # snmf repetitions per K
N_PC         <- 10
ALPHA        <- 10              # sNMF regularisation (LEA default)
SEED         <- 42
PLOIDY       <- 2               # pseudo-diploid encoding

PF_CHRS <- sprintf("Pf3D7_%02d_v3", 1:14)

set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# ============================================================================
# 1. Load VCF + metadata, restrict to nuclear chromosomes
# ============================================================================
log_msg("=== Phase 2 Step 2.2 (LEA variant): PCA + sNMF ===")
vcf  <- read.vcfR(vcf_in, verbose = FALSE)
meta <- readRDS(meta_in)

vcf_samples <- colnames(vcf@gt)[-1]
log_msg(sprintf("Monoclonal VCF samples: %d", length(vcf_samples)))
log_msg(sprintf("Sites (all chromosomes): %d", nrow(vcf@fix)))

chrom <- vcf@fix[, "CHROM"]
nuc <- chrom %in% PF_CHRS
vcf <- vcf[nuc, ]
log_msg(sprintf("Sites after restricting to %d nuclear chromosomes: %d",
                length(PF_CHRS), nrow(vcf@fix)))

# ============================================================================
# 2. Build numeric genotype matrix (samples x SNPs) with pseudo-haploid encoding
# ============================================================================
gt_chr <- extract.gt(vcf, element = "GT", as.numeric = FALSE)   # SNPs x samples

# The phase1_02 monoclonal VCF keeps the original diploid GT encoding (Pf
# monoclonal calls present as 0/0 or 1/1; a small fraction of sites are 0/1
# within an overall-monoclonal sample, treated here as missing for the
# pseudo-haploid encoding required by sNMF / ADMIXTURE).
#
# Mapping:
#   "0", "0/0", "0|0"   -> 0 (homo ref)
#   "1", "1/1", "1|1"   -> 2 (homo alt, pseudo-homozygous-diploid)
#   "0/1", "1/0", ...   -> 9 (het in a monoclonal sample: drop)
#   ".", "./.", NA      -> 9
encode_pseudo_hap <- function(x) {
   out <- rep(9L, length(x))
   out[x %in% c("0", "0/0", "0|0")] <- 0L
   out[x %in% c("1", "1/1", "1|1")] <- 2L
   out
}
n_het_string <- sum(gt_chr %in% c("0/1", "1/0", "0|1", "1|0"), na.rm = TRUE)
geno_mat <- matrix(encode_pseudo_hap(as.vector(gt_chr)),
                   nrow = nrow(gt_chr), ncol = ncol(gt_chr),
                   dimnames = dimnames(gt_chr))
storage.mode(geno_mat) <- "integer"

miss_frac <- mean(geno_mat == 9L)
log_msg(sprintf("Genotype matrix: %d SNPs x %d samples (missing %.2f%%, of which %d are het calls in monoclonal samples)",
                nrow(geno_mat), ncol(geno_mat),
                100 * miss_frac, n_het_string))

# Defensive: implausibly high missingness almost certainly means encoding mismatch.
if (miss_frac > 0.5) {
   sample_vals <- head(sort(unique(as.vector(gt_chr))), 20)
   stop(sprintf(
      "Aborting: genotype matrix is %.1f%% missing — likely a GT-encoding mismatch.\n  Unique GT values in VCF: %s",
      100 * miss_frac, paste(sample_vals, collapse = ", ")))
}

# Position info aligned with rows of geno_mat
snp_chrom <- vcf@fix[, "CHROM"]
snp_pos   <- as.integer(vcf@fix[, "POS"])
snp_id    <- seq_len(nrow(geno_mat))   # SNPRelate uses integer snp.id

# ============================================================================
# 3. Build GDS and LD-prune via SNPRelate
# ============================================================================
# Re-encode for SNPRelate which uses 0/1/2/3 with 3 = missing on the file side
# but the in-memory function snpgdsCreateGeno expects 0/1/2/NA.
geno_for_gds <- geno_mat
geno_for_gds[geno_for_gds == 9L] <- NA
# SNPRelate expects samples x SNPs (snpfirstdim = FALSE) or SNPs x samples
# (snpfirstdim = TRUE). We use samples x SNPs.
geno_for_gds <- t(geno_for_gds)
storage.mode(geno_for_gds) <- "integer"

# Map Pf chromosomes to integers 1..14 (SNPRelate is content-agnostic but
# downstream Manhattan-style plots prefer numeric chromosomes).
chr_int <- as.integer(factor(snp_chrom, levels = PF_CHRS))

if (file.exists(gds_path)) file.remove(gds_path)
snpgdsCreateGeno(gds_path,
                 genmat        = geno_for_gds,
                 sample.id     = vcf_samples,
                 snp.id        = snp_id,
                 snp.rs.id     = paste(snp_chrom, snp_pos, sep = ":"),
                 snp.chromosome = chr_int,
                 snp.position  = snp_pos,
                 snpfirstdim   = FALSE)
gds <- snpgdsOpen(gds_path)
on.exit(snpgdsClose(gds), add = TRUE)

log_msg(sprintf("LD pruning: r2 < %.2f, %d bp window", LD_R2, LD_WINDOW_BP))
prune_res <- snpgdsLDpruning(gds,
                             method        = "corr",
                             ld.threshold  = LD_R_THRESH,
                             slide.max.bp  = LD_WINDOW_BP,
                             autosome.only = FALSE,
                             verbose       = FALSE)
pruned_snp_id <- unlist(prune_res, use.names = FALSE)
log_msg(sprintf("LD-pruned SNPs retained: %d / %d",
                length(pruned_snp_id), nrow(geno_mat)))

# ============================================================================
# 4. PCA via SNPRelate on the pruned set
# ============================================================================
pca <- snpgdsPCA(gds, snp.id = pruned_snp_id,
                 autosome.only = FALSE, num.thread = 1, verbose = FALSE)

eigenval <- tibble(
   PC  = seq_len(N_PC),
   eig = pca$eigenval[seq_len(N_PC)],
   var_explained = pca$varprop[seq_len(N_PC)])

write_tsv(eigenval, file.path(out_dir, "pca_eigenval.tsv"))

eigenvec <- as_tibble(pca$eigenvect[, seq_len(N_PC)],
                      .name_repair = ~ paste0("PC", seq_along(.))) %>%
   mutate(SampleID = pca$sample.id, .before = 1)

write_tsv(eigenvec, file.path(out_dir, "pca_eigenvec.tsv"))

pca_df <- eigenvec %>% left_join(meta, by = "SampleID")

plot_pca <- function(df, colour_var, title) {
   ggplot(df, aes(x = PC1, y = PC2, colour = .data[[colour_var]])) +
      geom_point(size = 2.4, alpha = 0.85) +
      labs(
         # title = title,
           x = sprintf("PC1 (%.1f%%)", 100 * eigenval$var_explained[1]),
           y = sprintf("PC2 (%.1f%%)", 100 * eigenval$var_explained[2]),
           colour = colour_var) +
      theme_fig()
}

ggsave(file.path(out_dir, "pca_pc1_pc2.pdf"),
       plot_pca(pca_df, "region", "PCA (sNMF / SNPRelate) — by region"),
       width = 7, height = 5.5, dpi = 600)
ggsave(file.path(out_dir, "pca_pc1_pc2_by_village.pdf"),
       plot_pca(pca_df, "Location", "PCA — by village"),
       width = 7.5, height = 5.5, dpi = 600)
ggsave(file.path(out_dir, "pca_pc1_pc2_by_season.pdf"),
       plot_pca(pca_df, "transmission_season",
                "PCA — by transmission season"),
       width = 7, height = 5.5, dpi = 600)

p_scree <- ggplot(eigenval, aes(x = PC, y = var_explained)) +
   geom_col(fill = "steelblue") +
   geom_text(aes(label = percent(var_explained, accuracy = 0.1)),
             vjust = -0.4, size = 3) +
   scale_y_continuous(labels = percent) +
   labs(title = "PCA scree plot (SNPRelate, LD-pruned)",
        x = "PC", y = "Variance explained") +
   theme_minimal()
ggsave(file.path(out_dir, "pca_scree.pdf"), p_scree,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 5. Write .geno file (rows = SNPs, cols = samples) for the pruned set
# ============================================================================
# LEA::write.geno expects a matrix with rows = individuals, columns = SNPs.
# Reconstruct it from the pruned-SNP rows of geno_mat (which is SNPs x samples).
geno_pruned <- t(geno_mat[pruned_snp_id, , drop = FALSE])   # samples x SNPs
storage.mode(geno_pruned) <- "integer"
if (file.exists(geno_path)) file.remove(geno_path)
LEA::write.geno(geno_pruned, geno_path)
log_msg(sprintf("Wrote .geno: %s (%d samples x %d SNPs)",
                geno_path, nrow(geno_pruned), ncol(geno_pruned)))

# ============================================================================
# 6. sNMF across K = 2..8, with masked-genotype cross-entropy
# ============================================================================
# Remove any prior snmfProject for these inputs (LEA stores state on disk).
snmf_project_file <- file.path(snmf_dir, "gambia_mono_pruned.snmfProject")
if (file.exists(snmf_project_file)) {
   LEA::remove.snmfProject(snmf_project_file)
}

log_msg(sprintf("Running sNMF: K = %s, reps = %d, alpha = %d, entropy = TRUE",
                paste(range(K_RANGE), collapse = ".."), N_REPS, ALPHA))

project <- LEA::snmf(
   geno_path,
   K           = K_RANGE,
   repetitions = N_REPS,
   alpha       = ALPHA,
   entropy     = TRUE,
   ploidy      = PLOIDY,
   seed        = SEED,
   project     = "new",
   CPU         = 1)

# Cross-entropy per (K, run)
ce_rows <- lapply(K_RANGE, function(k) {
   ce <- LEA::cross.entropy(project, K = k)
   tibble(K = k, run = seq_along(ce), cross_entropy = as.numeric(ce))
})
ce_df <- bind_rows(ce_rows)
write_tsv(ce_df, file.path(out_dir, "snmf_cross_entropy.tsv"))

best_runs <- ce_df %>% group_by(K) %>%
   slice_min(cross_entropy, n = 1, with_ties = FALSE) %>% ungroup()

optimal_K <- best_runs$K[which.min(best_runs$cross_entropy)]
log_msg(sprintf("Optimal K (min cross-entropy) = %d (CE = %.4f, run = %d)",
                optimal_K,
                best_runs$cross_entropy[best_runs$K == optimal_K],
                best_runs$run[best_runs$K == optimal_K]))

p_ce <- ggplot(best_runs, aes(x = K, y = cross_entropy)) +
   geom_line() + geom_point(size = 3) +
   geom_point(data = best_runs %>% filter(K == optimal_K),
              colour = "red", size = 4) +
   geom_jitter(data = ce_df, width = 0.12, alpha = 0.6, size = 1.8,
               colour = "grey40") +
   scale_x_continuous(breaks = K_RANGE) +
   labs(title = sprintf("sNMF cross-entropy (optimal K = %d)", optimal_K),
        subtitle = sprintf("%d repetitions per K; line = best run per K",
                           N_REPS),
        x = "K", y = "Cross-entropy") +
   theme_minimal() +
   theme(
      legend.position = "none",
      plot.title        = element_text(size = 15, color = "black", face = "bold"),
      axis.title        = element_text(size = 14, color = "black", face = "bold"),
      axis.text         = element_text(size = 12, color = "black"),
      axis.line         = element_line(linewidth = 1, colour = "black", lineend = "square"),
      axis.ticks        = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.22, "cm")
   )

ggsave(file.path(out_dir, "snmf_cross_entropy.pdf"), plot = p_ce,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 7. Q matrices + ancestry bar plots
# ============================================================================
plot_order_df <- meta %>%
   filter(SampleID %in% vcf_samples) %>%
   mutate(region = factor(region, levels = c("Western", "Central", "Eastern"))) %>%
   arrange(region, Location, collection_month, SampleID)

stacked_barplot <- function(q_df, K, title) {
   q_long <- q_df %>%
      pivot_longer(starts_with("Q"), names_to = "cluster", values_to = "prop") %>%
      mutate(cluster = factor(cluster, levels = paste0("Q", seq_len(K))),
             SampleID = factor(SampleID, levels = plot_order_df$SampleID)) %>%
      left_join(meta %>% select(SampleID, region) %>%
                   mutate(region = factor(region,
                                          levels = c("Western", "Central", "Eastern"))),
                by = "SampleID")
   pal <- if (K <= 8) brewer.pal(max(3, K), "Set2")[seq_len(K)]
          else colorRampPalette(brewer.pal(8, "Set2"))(K)
   ggplot(q_long, aes(x = SampleID, y = prop, fill = cluster)) +
      geom_col(width = 1) +
      facet_grid(~ region, scales = "free_x", space = "free_x", switch = "x") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = pal) +
      labs(title = title, x = NULL, y = "Ancestry proportion", fill = "Clusters") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(size = 14, color = "black", face = "bold"),
            panel.spacing.x = unit(0.2, "lines"),
            strip.placement = "outside")
}

for (k in K_RANGE) {
   best_run_k <- best_runs$run[best_runs$K == k]
   q_mat <- LEA::Q(project, K = k, run = best_run_k)
   q_df <- as_tibble(q_mat, .name_repair = ~ paste0("Q", seq_along(.))) %>%
      mutate(SampleID = vcf_samples, .before = 1)
   
   write_tsv(q_df, file.path(out_dir, sprintf("snmf_Q_K%d.tsv", k)))

   p <- stacked_barplot(q_df, k,
                        sprintf("sNMF ancestry, K = %d (best of %d runs, CE = %.4f)",
                                k, N_REPS,
                                best_runs$cross_entropy[best_runs$K == k]))
   ggsave(file.path(out_dir, sprintf("snmf_barplot_K%d.pdf", k)),
          p, width = 12, height = 3.5, dpi = 600)
}

# ============================================================================
# 7b. Canonical filenames for downstream Snakemake aggregation
# ============================================================================
# Each per-K Snakemake job writes a stable `qmatrix.tsv` (best run's Q) and
# `ce.tsv` (cross-entropy across this K's runs). The aggregator (extendK)
# reads ce.tsv across K dirs to build the CE-vs-K curve without touching the
# snmfProject file — see plan Task 5 Pattern A.
q_best <- read_tsv(file.path(out_dir, sprintf("snmf_Q_K%d.tsv", k_single)),
                   show_col_types = FALSE)
write_tsv(q_best, file.path(out_dir, "qmatrix.tsv"))
write_tsv(ce_df, file.path(out_dir, "ce.tsv"))

# ============================================================================
# 8. Persist + log
# ============================================================================
saveRDS(list(eigenvec    = eigenvec,
             eigenval    = eigenval,
             cross_entropy = ce_df,
             best_runs   = best_runs,
             optimal_K   = optimal_K,
             pruned_snps = length(pruned_snp_id),
             total_snps  = nrow(geno_mat)),
        file.path(out_dir, "phase2_02_lea_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase2_02_lea_summary.txt"))
message("\nDone. Outputs in ", out_dir)
