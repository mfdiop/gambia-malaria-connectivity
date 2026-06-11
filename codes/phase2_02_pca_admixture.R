# =============================================================================
# Phase 2 — Step 2.2: Population structure (PCA + ADMIXTURE)
# =============================================================================
# Per notes/objective3_proposal.md §6 Phase 2 Step 2.2:
#   - PCA on LD-pruned SNPs (PLINK2, r2 < 0.1 in 50 kb windows)
#   - ADMIXTURE K = 2..8, 10-fold cross-validation, pick optimal K (min CV error)
#   - Question: do western and eastern Gambia form discrete clusters, or a
#     continuous gradient consistent with ongoing gene flow?
#
# Sample set : monoclonal subset (results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz)
#              ADMIXTURE assumes unrelated, unambiguous genotype calls.
#              Polyclonal samples re-enter downstream via dcifer.
#
# P. falciparum quirks handled here:
#   - Non-standard chromosome names (Pf3D7_01_v3 .. Pf3D7_14_v3) — PLINK2 reads
#     them with --allow-extra-chr; .bim chromosomes renamed to 1..14 for
#     ADMIXTURE (requires integer chromosomes).
#   - Haploid calls — PLINK2 keeps them as such; ADMIXTURE treats homozygous
#     diploid genotypes equivalently for allele-frequency estimation, which is
#     the standard pseudo-haploid encoding used in Pf population genetics.
#   - Mitochondria / apicoplast / hypervariable scaffolds dropped (kept only
#     Pf3D7_01_v3..Pf3D7_14_v3).
#
# Output : results/phase2_02/
#            plink/                              PLINK2 intermediate files
#            admixture/                          ADMIXTURE outputs (Q, P, log)
#            pca_eigenvec.tsv                    sample x PC1..PC10
#            pca_eigenval.tsv                    eigenvalues
#            pca_pc1_pc2.pdf                     PC1 vs PC2 by region
#            pca_pc1_pc2_by_village.pdf          coloured by village
#            pca_pc1_pc2_by_season.pdf           coloured by transmission season
#            pca_scree.pdf                       eigenvalue scree
#            admixture_cv_error.tsv              CV error per K
#            admixture_cv_error.pdf              CV vs K curve
#            admixture_Q_K{2..8}.tsv             per-sample ancestry proportions
#            admixture_barplot_K{2..8}.pdf       stacked-bar ancestry per K
#            admixture_barplot_optimalK.pdf      bar plot at min-CV K
#            phase2_02_summary.txt               human-readable log
#
# Tools  : PLINK2 (>= 2.0), ADMIXTURE v1.3.
#          Override binary paths with env vars PLINK2 and ADMIXTURE if not on PATH.
# =============================================================================

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

suppressPackageStartupMessages({
   library(tidyverse)
   library(scales)
})

# ----- Paths ----------------------------------------------------------------
vcf_in   <- "results/phase1_02/gambia_2014_qc_monoclonal.vcf.gz"
meta_in  <- "results/phase1_03/sample_master_metadata.rds"

out_dir   <- "results/phase2_02"
plink_dir <- file.path(out_dir, "plink")
adm_dir   <- file.path(out_dir, "admixture")
dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(adm_dir,   recursive = TRUE, showWarnings = FALSE)

# ----- External binaries ----------------------------------------------------
plink2_bin    <- Sys.getenv("PLINK2",    "plink2")
# plink2_bin    <- Sys.getenv("PLINK2",    "plink")
admixture_bin <- Sys.getenv("ADMIXTURE", "admixture")

check_bin <- function(bin, label) {
   ok <- nzchar(Sys.which(bin)) || file.exists(bin)
   if (!ok) stop(sprintf(
      "%s binary not found ('%s'). Install it or set the %s env var to its full path.",
      label, bin, toupper(label)))
}

check_bin(plink2_bin,    "plink2")
check_bin(admixture_bin, "admixture")

# ----- Parameters -----------------------------------------------------------
LD_WINDOW   <- "50kb"     # LD-pruning window (proposal)
LD_STEP     <- 1          # variant step
LD_R2       <- 0.2        # r^2 threshold (proposal)
K_RANGE     <- 2:8        # ADMIXTURE K range (proposal)
CV_FOLDS    <- 10         # ADMIXTURE cross-validation folds (proposal)
N_PC        <- 10         # PCs to extract
SEED        <- 42

PF_CHRS <- sprintf("Pf3D7_%02d_v3", 1:14)   # nuclear chromosomes only

set.seed(SEED)

log_lines <- c()
log_msg <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

run <- function(cmd, args, log_file = NULL) {
   log_msg(sprintf("$ %s %s", cmd, paste(args, collapse = " ")))
   res <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
   if (!is.null(log_file)) writeLines(res, log_file)
   status <- attr(res, "status")
   if (!is.null(status) && status != 0) {
      writeLines(res)
      stop(sprintf("%s failed (exit %d)", cmd, status))
   }
   invisible(res)
}

# ============================================================================
# 1. Load metadata
# ============================================================================
log_msg("=== Phase 2 Step 2.2: PCA + ADMIXTURE ===")
meta <- readRDS(meta_in)
log_msg(sprintf("Metadata rows: %d", nrow(meta)))

# ============================================================================
# 2. VCF -> PLINK pgen, restrict to nuclear chromosomes, biallelic SNPs
# ============================================================================
pfile_all <- file.path(plink_dir, "gambia_mono")
run(plink2_bin, c(
   "--vcf", vcf_in,
   "--allow-extra-chr",
   "--chr",  paste(PF_CHRS, collapse = ","),
   "--max-alleles", "2",
   "--min-alleles", "2",
   "--snps-only",
   "--set-all-var-ids", "@:#",
   "--rm-dup", "exclude-all",
   "--make-pgen",
   "--out", pfile_all),
   log_file = file.path(plink_dir, "01_import.log"))

# ============================================================================
# 3. LD pruning (r2 < 0.1, 50 kb window)
# ============================================================================
prune_prefix <- file.path(plink_dir, "ld_prune")
run(plink2_bin, c(
   "--pfile", pfile_all,
   "--allow-extra-chr",
   "--indep-pairwise", LD_WINDOW, LD_STEP, LD_R2,
   "--out", prune_prefix),
   log_file = file.path(plink_dir, "02_indep_pairwise.log"))

pfile_pruned <- file.path(plink_dir, "gambia_mono_pruned")
run(plink2_bin, c(
   "--pfile", pfile_all,
   "--allow-extra-chr",
   "--extract", paste0(prune_prefix, ".prune.in"),
   "--make-pgen",
   "--out", pfile_pruned),
   log_file = file.path(plink_dir, "03_extract_pruned.log"))

n_pruned <- length(readLines(paste0(prune_prefix, ".prune.in")))
log_msg(sprintf("LD-pruned SNPs retained: %d", n_pruned))

# ============================================================================
# 4. PCA
# ============================================================================
pca_prefix <- file.path(plink_dir, "pca")
run(plink2_bin, c(
   "--pfile", pfile_pruned,
   "--allow-extra-chr",
   "--pca", N_PC,
   "--out", pca_prefix),
   log_file = file.path(plink_dir, "04_pca.log"))

eigenvec <- read.table(paste0(pca_prefix, ".eigenvec"),
                       header = TRUE, comment.char = "", check.names = FALSE)

# PLINK2 writes "#IID" or "#FID IID"; normalise to SampleID
names(eigenvec)[1] <- sub("^#", "", names(eigenvec)[1])

if ("IID" %in% names(eigenvec)) {
   eigenvec <- eigenvec %>% rename(SampleID = IID) %>% select(-any_of("FID"))
} else {
   names(eigenvec)[1] <- "SampleID"
}

eigenval <- tibble(
   PC  = seq_len(N_PC),
   eig = scan(paste0(pca_prefix, ".eigenval"), quiet = TRUE)) %>%
   mutate(var_explained = eig / sum(eig))

write_tsv(eigenvec, file.path(out_dir, "pca_eigenvec.tsv"))
write_tsv(eigenval, file.path(out_dir, "pca_eigenval.tsv"))

pca <- eigenvec %>% left_join(meta, by = "SampleID")
log_msg(sprintf("PCA samples with metadata: %d / %d",
                sum(!is.na(pca$region)), nrow(pca)))

pca <- pca  %>% 
   mutate(region = recode(region, "west" = "Western", 
                          "central" = "Central", "east" = "Eastern"
   )) %>% 
   rename(., Regions = region)

# PC1 vs PC2 by region
plot_pca <- function(df, colour_var, title, file) {
   ggplot(df, aes(x = PC1, y = PC2, colour = .data[[colour_var]])) +
      geom_point(size = 2.7, alpha = 1) +
      labs(
         # title = title,
         x = sprintf("PC1 (%.1f%%)", 100 * eigenval$var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", 100 * eigenval$var_explained[2]),
         colour = colour_var) +
      theme_fig()
}

ggsave(file.path(out_dir, "pca_pc1_pc2.pdf"),
       plot_pca(pca, "Regions", "PCA — coloured by region", NULL),
       width = 7, height = 5.5, dpi = 600)
ggsave(file.path(out_dir, "pca_pc1_pc2_by_village.pdf"),
       plot_pca(pca, "Location", "PCA — coloured by village", NULL),
       width = 7.5, height = 5.5, dpi = 600)
ggsave(file.path(out_dir, "pca_pc1_pc2_by_season.pdf"),
       plot_pca(pca, "transmission_season",
                "PCA — coloured by transmission season", NULL),
       width = 7, height = 5.5, dpi = 600)

# Scree
p_scree <- ggplot(eigenval, aes(x = PC, y = var_explained)) +
   geom_col(fill = "steelblue") +
   geom_text(aes(label = percent(var_explained, accuracy = 0.1)),
             vjust = -0.4, size = 3) +
   scale_y_continuous(labels = percent) +
   labs(title = "PCA scree plot", x = "PC", y = "Variance explained") +
   theme_minimal()
ggsave(file.path(out_dir, "pca_scree.pdf"), p_scree,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 5. Export bed for ADMIXTURE (integer chromosomes 1..14)
# ============================================================================
# ADMIXTURE requires .bed/.bim/.fam with numeric chromosomes. Use PLINK2 to
# emit bed, then rewrite the .bim chrom column.
bed_prefix <- file.path(adm_dir, "gambia_mono_pruned")
run(plink2_bin, c(
   "--pfile", pfile_pruned,
   "--allow-extra-chr",
   "--make-bed",
   "--out", bed_prefix),
   log_file = file.path(adm_dir, "05_make_bed.log"))

bim <- read.table(paste0(bed_prefix, ".bim"),
                  sep = "\t", stringsAsFactors = FALSE)
chr_map <- setNames(as.character(1:14), PF_CHRS)
bim$V1 <- chr_map[bim$V1]
stopifnot(!any(is.na(bim$V1)))
write.table(bim, paste0(bed_prefix, ".bim"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
log_msg("Renamed .bim chromosomes Pf3D7_NN_v3 -> 1..14 for ADMIXTURE")

# ============================================================================
# 6. ADMIXTURE K = 2..8 with 10-fold CV
# ============================================================================
# ADMIXTURE writes outputs to the current working directory; run from adm_dir.
bed_basename <- basename(paste0(bed_prefix, ".bed"))
cv_lines <- list()
old_wd <- getwd()
setwd(adm_dir)
tryCatch({
   for (K in K_RANGE) {
      log_msg(sprintf("ADMIXTURE K = %d (cv = %d folds, seed = %d)",
                      K, CV_FOLDS, SEED))
      log_file <- sprintf("admixture_K%d.log", K)
      run(admixture_bin,
          c(sprintf("--cv=%d", CV_FOLDS),
            "-s", SEED,
            bed_basename, K),
          log_file = log_file)
      cv_lines[[as.character(K)]] <- grep("^CV error", readLines(log_file),
                                          value = TRUE)
   }
}, finally = setwd(old_wd))

# Parse CV error
parse_cv <- function(line, K) {
   v <- as.numeric(sub(".*:\\s*", "", line))
   tibble(K = K, cv_error = v)
}
cv_df <- bind_rows(mapply(parse_cv, cv_lines, as.integer(names(cv_lines)),
                          SIMPLIFY = FALSE))
write_tsv(cv_df, file.path(out_dir, "admixture_cv_error.tsv"))

optimal_K <- cv_df$K[which.min(cv_df$cv_error)]
log_msg(sprintf("Optimal K (min CV error) = %d (CV = %.4f)",
                optimal_K, min(cv_df$cv_error)))

p_cv <- ggplot(cv_df, aes(x = K, y = cv_error)) +
   geom_line() + geom_point(size = 3) +
   geom_point(data = cv_df %>% filter(K == optimal_K),
              colour = "red", size = 4) +
   scale_x_continuous(breaks = K_RANGE) +
   labs(title = sprintf("ADMIXTURE cross-validation error (optimal K = %d)",
                        optimal_K),
        x = "K", y = sprintf("%d-fold CV error", CV_FOLDS)) +
   theme_minimal()

ggsave(file.path(out_dir, "admixture_cv_error.pdf"), p_cv,
       width = 6, height = 4, dpi = 600)

# ============================================================================
# 7. Build Q tables and bar plots per K
# ============================================================================
fam <- read.table(paste0(bed_prefix, ".fam"), stringsAsFactors = FALSE)
sample_order <- fam$V2   # IID column

# Sample ordering for plotting: region (west/central/east) -> village -> month
plot_order_df <- meta %>%
   filter(SampleID %in% sample_order) %>%
   mutate(region = factor(region, levels = c("Western", "Central", "Eastern"))) %>%
   arrange(region, Location, collection_month, SampleID)

stacked_barplot <- function(q_df, K, title) {
   q_long <- q_df %>%
      pivot_longer(starts_with("Q"), names_to = "cluster", values_to = "prop") %>%
      mutate(cluster = factor(cluster, levels = paste0("Q", seq_len(K))),
             SampleID = factor(SampleID, levels = plot_order_df$SampleID))
   
   region_levels <- meta %>%
      filter(SampleID %in% q_df$SampleID) %>%
      mutate(region = factor(region, levels = c("Western", "Central", "Eastern"))) %>%
      arrange(match(SampleID, plot_order_df$SampleID))
   
   ggplot(q_long, aes(x = SampleID, y = prop, fill = cluster)) +
      geom_col(width = 1) +
      facet_grid(~ region_levels$region[match(SampleID, region_levels$SampleID)],
                 scales = "free_x", space = "free_x", switch = "x") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = title, x = NULL, y = "Ancestry proportion", fill = "Clusters") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.title = element_text(size = 14, color = "black", face = "bold"),
            panel.spacing.x = unit(0.2, "lines"),
            strip.placement = "outside")
}

for (K in K_RANGE) {
   q_file <- file.path(adm_dir, sprintf("%s.%d.Q",
                                        tools::file_path_sans_ext(bed_basename), K))
   q <- read.table(q_file)
   names(q) <- paste0("Q", seq_len(K))
   q <- bind_cols(SampleID = sample_order, q)
   write_tsv(q, file.path(out_dir, sprintf("admixture_Q_K%d.tsv", K)))

   p <- stacked_barplot(q, K,
                        sprintf("ADMIXTURE ancestry, K = %d (samples ordered by region/village)",
                                K))
   ggsave(file.path(out_dir, sprintf("admixture_barplot_K%d.pdf", K)),
          p, width = 12, height = 3.5, dpi = 600)
}

# Highlighted optimal-K bar plot
q_opt <- read_tsv(file.path(out_dir, sprintf("admixture_Q_K%d.tsv", optimal_K)),
                  show_col_types = FALSE)
p_opt <- stacked_barplot(q_opt, optimal_K,
                         sprintf("ADMIXTURE ancestry at optimal K = %d",
                                 optimal_K))
ggsave(file.path(out_dir, "admixture_barplot_optimalK.pdf"),
       p_opt, width = 12, height = 4, dpi = 600)

# ============================================================================
# 8. Persist + log
# ============================================================================
saveRDS(list(eigenvec = eigenvec, eigenval = eigenval,
             cv = cv_df, optimal_K = optimal_K,
             pruned_snps = n_pruned),
        file.path(out_dir, "phase2_02_results.rds"))

writeLines(log_lines, file.path(out_dir, "phase2_02_summary.txt"))
message("\nDone. Outputs in ", out_dir)
