#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(bsseq)
  library(methylSig)
  library(ggplot2)
})

# ============================================================
# User settings
# ============================================================

# ---- Input files ----
cpg_bed_path   <- "~/software/wgbs_tools/wgbs_tools/references/hg38/CpG.bed.gz"
#coseg_bed_path <- "coseg_1.bed"

# Directory containing paired bigWig files:
#   *_FractionalMethylation.bigwig
#   *_ReadCoverage.bigwig
bigwig_dir <- "roadmap/sample_1"

# Sample metadata CSV:
# must contain columns:
#   sample_id, group
# where sample_id matches the GSM/sample prefix in bigWig filenames
#
# Example:
# sample_id,group
# GSM5652181,case
# GSM5652182,case
# GSM5652183,control
# GSM5652184,control
sample_sheet_path <- "roadmap/sample_1/sample_1_sheet.csv"

# Chromosomes to include. Set to NULL to use all chromosomes in CpG BED.
chroms_to_use <- c("chr2")

# Segment counts to test
segment_counts <- c(500, 1000, 2000, 5000, 10000, 20000)

# Number of random repetitions per segment count
n_random_reps <- 3

# methylSig settings
fdr_threshold     <- 0.05
min_cov_per_cpg   <- 1
n_cores           <- 4

# Output files
results_csv <- "dm_fraction_sig_results_adi-cd8.csv"
plot_pdf    <- "dm_fraction_sig_plot_adi-cd8.pdf"

# Random seed for reproducibility
set.seed(42)

# ============================================================
# Helper functions
# ============================================================

read_coseg_bed <- function(path) {
  # Skip first track header line if present
  first_line <- readLines(path, n = 1)
  skip_n <- 0
  #skip_n <- if (length(first_line) > 0 && grepl("^track", first_line[1])) 1 else 0
  
  dt <- fread(
    path,
    skip = skip_n,
    header = FALSE,
    col.names = c("chrom", "start", "end", "n_tracks", "co_count")
  )
  
  dt[, `:=`(
    start    = as.integer(start),
    end      = as.integer(end),
    n_tracks = as.integer(n_tracks),
    co_count = as.integer(co_count)
  )]
  
  dt
}

make_sample_pairs <- function(bigwig_dir, sample_sheet) {
  frac_files <- list.files(
    bigwig_dir,
    pattern = "_FractionalMethylation\\.bigwig$",
    full.names = TRUE
  )
  cov_files <- list.files(
    bigwig_dir,
    pattern = "_ReadCoverage\\.bigwig$",
    full.names = TRUE
  )
  
  if (length(frac_files) == 0) stop("No fractional methylation bigWig files found.")
  if (length(cov_files) == 0) stop("No read coverage bigWig files found.")
  
  frac_ids <- sub("_WGBS_FractionalMethylation\\.bigwig$", "", basename(frac_files))
  cov_ids  <- sub("_WGBS_ReadCoverage\\.bigwig$", "", basename(cov_files))
  
  dt_frac <- data.table(sample_id = frac_ids, frac_bw = frac_files)
  dt_cov  <- data.table(sample_id = cov_ids,  cov_bw  = cov_files)
  
  dt <- merge(dt_frac, dt_cov, by = "sample_id", all = FALSE)
  dt <- merge(dt, sample_sheet, by = "sample_id", all = FALSE)
  
  if (nrow(dt) == 0) {
    stop("No matched sample_id across sample_sheet and bigWig file pairs.")
  }
  
  dt
}

import_bigwig_at_cpgs <- function(bw_path, cpg_gr) {
  bw <- import(bw_path, which = reduce(cpg_gr))
  hits <- findOverlaps(cpg_gr, bw)
  
  vals <- rep(NA_real_, length(cpg_gr))
  vals[queryHits(hits)] <- bw$score[subjectHits(hits)]
  vals
}

build_bsseq_from_bigwigs <- function(cpg_dt, sample_pairs, chroms_to_use = NULL,
                                     min_cov_per_cpg = 1) {
  if (!is.null(chroms_to_use)) {
    cpg_dt <- cpg_dt[chrom %in% chroms_to_use]
  }
  
  cpg_dt <- copy(cpg_dt)
  setorder(cpg_dt, chrom, start)
  
  cpg_pos <- as.integer(cpg_dt$start)
  
  cpg_gr <- GRanges(
    seqnames = cpg_dt$chrom,
    ranges   = IRanges(start = cpg_pos, end = cpg_pos + 1L)
  )
  
  n_cpg <- nrow(cpg_dt)
  n_samples <- nrow(sample_pairs)
  
  M_mat   <- matrix(NA_integer_, nrow = n_cpg, ncol = n_samples)
  Cov_mat <- matrix(NA_integer_, nrow = n_cpg, ncol = n_samples)
  
  colnames(M_mat)   <- sample_pairs$sample_id
  colnames(Cov_mat) <- sample_pairs$sample_id
  
  chrom_index <- split(seq_len(n_cpg), cpg_dt$chrom)
  
  for (i in seq_len(n_samples)) {
    sid <- sample_pairs$sample_id[i]
    message("Importing sample: ", sid)
    
    frac_vals_all <- rep(NA_real_, n_cpg)
    cov_vals_all  <- rep(NA_real_, n_cpg)
    
    for (chrom_i in names(chrom_index)) {
      idx <- chrom_index[[chrom_i]]
      cpg_gr_chr <- cpg_gr[idx]
      
      message("  chromosome: ", chrom_i)
      
      frac_vals_all[idx] <- import_bigwig_at_cpgs(sample_pairs$frac_bw[i], cpg_gr_chr)
      cov_vals_all[idx]  <- import_bigwig_at_cpgs(sample_pairs$cov_bw[i],  cpg_gr_chr)
    }
    
    finite_frac <- frac_vals_all[is.finite(frac_vals_all)]
    if (length(finite_frac) > 0 && max(finite_frac, na.rm = TRUE) > 1.5) {
      message("Detected fraction values > 1 for sample ", sid, "; dividing by 100.")
      frac_vals_all <- frac_vals_all / 100
    }
    
    cov_vals_all[!is.finite(cov_vals_all)] <- NA_real_
    frac_vals_all[!is.finite(frac_vals_all)] <- NA_real_
    
    cov_int <- as.integer(round(cov_vals_all))
    cov_int[!is.finite(cov_int)] <- 0L
    cov_int[cov_int < min_cov_per_cpg] <- 0L
    
    frac_vals_all[!is.finite(frac_vals_all)] <- 0
    
    M_int <- as.integer(round(frac_vals_all * cov_int))
    
    # clamp into valid count range
    M_int[!is.finite(M_int)] <- 0L
    M_int[M_int < 0] <- 0L
    
    over_idx <- which(M_int > cov_int)
    if (length(over_idx) > 0) M_int[over_idx] <- cov_int[over_idx]
    
    M_mat[, i]   <- M_int
    Cov_mat[, i] <- cov_int
  }
  
  keep <- rowSums(Cov_mat > 0) > 0
  
  M_keep   <- M_mat[keep, , drop = FALSE]
  Cov_keep <- Cov_mat[keep, , drop = FALSE]
  
  # make absolutely sure assay column names are present
  colnames(M_keep)   <- sample_pairs$sample_id
  colnames(Cov_keep) <- sample_pairs$sample_id
  
  bs <- BSseq(
    chr = as.character(seqnames(cpg_gr))[keep],
    pos = cpg_pos[keep],
    M   = M_keep,
    Cov = Cov_keep,
    sampleNames = sample_pairs$sample_id
  )
  
  # add sample metadata after construction
  colData(bs)$group <- sample_pairs$group
  
  bs
}

make_coseg_segments_one_chrom <- function(cpg_pos, coseg_chr_dt, K, chrom) {
  n_cpg <- length(cpg_pos)
  if (n_cpg < 2) stop("Need at least 2 CpGs.")
  if (K < 1) stop("K must be >= 1.")
  if (K > n_cpg) stop("K cannot exceed number of CpGs.")
  
  if (K == 1) {
    return(GRanges(
      seqnames = chrom,
      ranges = IRanges(start = min(cpg_pos), end = max(cpg_pos) + 1)
    ))
  }
  
  # coseg file should have n_cpg - 1 adjacent intervals
  if (nrow(coseg_chr_dt) != (n_cpg - 1)) {
    stop(sprintf(
      "Chrom %s: expected %d coseg rows, found %d. CpG positions and coseg file are inconsistent.",
      chrom, n_cpg - 1, nrow(coseg_chr_dt)
    ))
  }
  
  n_tracks <- unique(coseg_chr_dt$n_tracks)
  if (length(n_tracks) != 1) {
    stop("Expected one unique n_tracks per chromosome in coseg file.")
  }
  
  sim <- coseg_chr_dt$co_count / n_tracks
  boundary_strength <- 1 - sim  # larger = stronger cut
  
  cut_idx <- sort(order(boundary_strength, decreasing = TRUE)[seq_len(K - 1)])
  
  starts_idx <- c(1, cut_idx + 1)
  ends_idx   <- c(cut_idx, n_cpg)
  
  GRanges(
    seqnames = chrom,
    ranges = IRanges(
      start = cpg_pos[starts_idx],
      end   = cpg_pos[ends_idx] + 1
    )
  )
}

make_random_segments_one_chrom <- function(cpg_pos, K, chrom) {
  n_cpg <- length(cpg_pos)
  if (n_cpg < 2) stop("Need at least 2 CpGs.")
  if (K < 1) stop("K must be >= 1.")
  if (K > n_cpg) stop("K cannot exceed number of CpGs.")
  
  if (K == 1) {
    return(GRanges(
      seqnames = chrom,
      ranges = IRanges(start = min(cpg_pos), end = max(cpg_pos) + 1)
    ))
  }
  
  cut_idx <- sort(sample(seq_len(n_cpg - 1), K - 1, replace = FALSE))
  
  starts_idx <- c(1, cut_idx + 1)
  ends_idx   <- c(cut_idx, n_cpg)
  
  GRanges(
    seqnames = chrom,
    ranges = IRanges(
      start = cpg_pos[starts_idx],
      end   = cpg_pos[ends_idx] + 1
    )
  )
}

allocate_segments_across_chroms <- function(cpg_dt, K_total, chroms_to_use = NULL) {
  if (!is.null(chroms_to_use)) {
    cpg_dt <- cpg_dt[chrom %in% chroms_to_use]
  }
  
  counts <- cpg_dt[, .(n_cpg = .N), by = chrom]
  counts <- counts[n_cpg >= 1]
  if (nrow(counts) == 0) stop("No CpGs available after filtering chromosomes.")
  
  # Allocate segments roughly proportional to CpG counts
  counts[, raw_k := K_total * n_cpg / sum(n_cpg)]
  counts[, K := pmax(1L, as.integer(floor(raw_k)))]
  
  # Adjust to match exact total K_total
  diff_k <- K_total - sum(counts$K)
  
  if (diff_k > 0) {
    counts[, frac_part := raw_k - floor(raw_k)]
    add_idx <- order(counts$frac_part, decreasing = TRUE)[seq_len(diff_k)]
    counts$K[add_idx] <- counts$K[add_idx] + 1L
  } else if (diff_k < 0) {
    counts[, frac_part := raw_k - floor(raw_k)]
    subtract_needed <- abs(diff_k)
    can_subtract <- which(counts$K > 1L)
    if (length(can_subtract) < subtract_needed) {
      stop("Cannot allocate requested total segment count across chromosomes.")
    }
    sub_idx <- can_subtract[order(counts$frac_part[can_subtract], decreasing = FALSE)[seq_len(subtract_needed)]]
    counts$K[sub_idx] <- counts$K[sub_idx] - 1L
  }
  
  counts[, .(chrom, K)]
}

build_segments_coseg <- function(cpg_dt, coseg_dt, K_total, chroms_to_use = NULL) {
  alloc <- allocate_segments_across_chroms(cpg_dt, K_total, chroms_to_use)
  
  seg_list <- vector("list", nrow(alloc))
  
  for (i in seq_len(nrow(alloc))) {
    chrom_i <- alloc$chrom[i]
    K_i     <- alloc$K[i]
    
    cpg_chr   <- cpg_dt[chrom == chrom_i]
    coseg_chr <- coseg_dt[chrom == chrom_i]
    
    seg_list[[i]] <- make_coseg_segments_one_chrom(
      cpg_pos = cpg_chr$start,
      coseg_chr_dt = coseg_chr,
      K = K_i,
      chrom = chrom_i
    )
  }
  
  do.call(c, seg_list)
}

build_segments_random <- function(cpg_dt, K_total, chroms_to_use = NULL) {
  alloc <- allocate_segments_across_chroms(cpg_dt, K_total, chroms_to_use)
  
  seg_list <- vector("list", nrow(alloc))
  
  for (i in seq_len(nrow(alloc))) {
    chrom_i <- alloc$chrom[i]
    K_i     <- alloc$K[i]
    
    cpg_chr <- cpg_dt[chrom == chrom_i]
    
    seg_list[[i]] <- make_random_segments_one_chrom(
      cpg_pos = cpg_chr$start,
      K = K_i,
      chrom = chrom_i
    )
  }
  
  do.call(c, seg_list)
}

run_diff_methylsig_tile_table <- function(bs, seg_gr, fdr_threshold = 0.05, n_cores = 1) {
  tiled <- tile_by_regions(bs = bs, gr = seg_gr)
  
  diff_gr <- diff_methylsig(
    bs = tiled,
    group_column = "group",
    comparison_groups = c(case = "case", control = "control"),
    disp_groups = c(case = TRUE, control = TRUE),
    local_window_size = 0,
    t_approx = TRUE,
    n_cores = n_cores
  )
  
  meta <- as.data.table(as.data.frame(mcols(diff_gr)))
  meta_names <- names(meta)
  
  fdr_col <- NULL
  for (nm in c("fdr", "FDR", "padj", "qvalue", "adj_pvalue")) {
    if (nm %in% meta_names) {
      fdr_col <- nm
      break
    }
  }
  
  if (is.null(fdr_col)) {
    stop(
      "Could not find an FDR-like column in diff_methylsig output. Available columns are: ",
      paste(meta_names, collapse = ", ")
    )
  }
  
  p_col <- NULL
  for (nm in c("pvalue", "p_value", "p.val", "pval", "PValue")) {
    if (nm %in% meta_names) {
      p_col <- nm
      break
    }
  }
  
  tile_dt <- data.table(
    chrom = as.character(seqnames(diff_gr)),
    start = start(diff_gr),
    end   = end(diff_gr),
    width = width(diff_gr)
  )
  
  tile_dt[, fdr := meta[[fdr_col]]]
  
  if (!is.null(p_col)) {
    tile_dt[, pvalue := meta[[p_col]]]
  } else {
    tile_dt[, pvalue := NA_real_]
  }
  
  tile_dt[, significant := !is.na(fdr) & fdr < fdr_threshold]
  tile_dt[, signal_score := -log10(pmax(fdr, .Machine$double.xmin))]
  tile_dt[, contribution_to_fraction := as.numeric(significant) / .N]
  
  tile_dt
}

# ============================================================
# Main
# ============================================================

message("Reading sample sheet...")
sample_sheet <- fread(sample_sheet_path)
required_cols <- c("sample_id", "group")
if (!all(required_cols %in% names(sample_sheet))) {
  stop("sample_sheet.csv must contain columns: sample_id, group")
}
sample_sheet[, group := as.character(group)]

if (!all(sort(unique(sample_sheet$group)) == c("case", "control"))) {
  stop("sample_sheet$group must contain exactly two labels: 'case' and 'control'")
}

message("Reading CpG BED...")
cpg_dt <- fread(cpg_bed_path, col.names = c("chrom", "start", "end"))
cpg_dt[, `:=`(start = as.integer(start), end = as.integer(end))]

if (!is.null(chroms_to_use)) {
  cpg_dt <- cpg_dt[chrom %in% chroms_to_use]
}
setorder(cpg_dt, chrom, start)

# temporary speed test
#cpg_dt <- cpg_dt[1:1000]

message("Matching fraction/coverage bigWig pairs...")
sample_pairs <- make_sample_pairs(bigwig_dir, sample_sheet)

message("Building BSseq object from bigWigs...")
bs <- build_bsseq_from_bigwigs(
  cpg_dt = cpg_dt,
  sample_pairs = sample_pairs,
  chroms_to_use = chroms_to_use,
  min_cov_per_cpg = min_cov_per_cpg
)

bs_dt <- data.table(
  chrom = as.character(seqnames(bs)),
  start = start(granges(bs))
)
setorder(bs_dt, chrom, start)

message("Running segmentation benchmark...")
results <- list()
tile_results <- list()
idx <- 1
tile_idx <- 1

for (K in segment_counts) {
  message("K = ", K, " | random baseline")
  for (rep_i in seq_len(n_random_reps)) {
    message("  random replicate ", rep_i, "/", n_random_reps)
    
    seg_rand <- build_segments_random(
      cpg_dt = bs_dt,
      K_total = K,
      chroms_to_use = chroms_to_use
    )
    
    tile_dt <- run_diff_methylsig_tile_table(
      bs = bs,
      seg_gr = seg_rand,
      fdr_threshold = fdr_threshold,
      n_cores = n_cores
    )
    
    frac_sig_rand <- mean(tile_dt$significant, na.rm = TRUE)
    
    results[[idx]] <- data.table(
      method = "random",
      segment_count = K,
      replicate = rep_i,
      frac_sig = frac_sig_rand
    )
    idx <- idx + 1
    
    tile_dt[, `:=`(
      method = "random",
      segment_count = K,
      replicate = rep_i,
      tile_id = seq_len(.N)
    )]
    tile_results[[tile_idx]] <- tile_dt
    tile_idx <- tile_idx + 1
  }
}

# results_dt = one row per K/replicate, used for original fraction plot
# tile_results_dt  = one row per tile, used for the new tile-level plots
results_dt <- rbindlist(results)
tile_results_dt <- rbindlist(tile_results)

message("Saving raw results...")
fwrite(results_dt, results_csv)
tile_results_csv <- sub("\\.csv$", "_tile_level.csv", results_csv)
fwrite(tile_results_dt, tile_results_csv)

summary_dt <- results_dt[
  , .(
    mean_frac_sig = mean(frac_sig, na.rm = TRUE),
    sd_frac_sig   = sd(frac_sig, na.rm = TRUE)
  ),
  by = .(method, segment_count)
]

summary_dt[, ymin := pmax(0, mean_frac_sig - sd_frac_sig)]
summary_dt[, ymax := pmin(1, mean_frac_sig + sd_frac_sig)]

message("Making plot...")

# plot 1: fraction of segments with FDR>0.05 vs. segment number
p <- ggplot(summary_dt, aes(x = segment_count, y = mean_frac_sig)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of segments",
    y = "Fraction of segments with FDR < 0.05",
    title = "Differential methylation with random tiles"
  ) +
  theme_bw()

ggsave(plot_pdf, p, width = 8, height = 5)

# plot 2: number of segments vs. -log10(FDR)
tile_signal_hist_pdf <- sub("\\.pdf$", "_tile_signal_distribution.pdf", plot_pdf)

p_tile_hist <- ggplot(tile_results_dt, aes(x = signal_score)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = -log10(fdr_threshold), linetype = "dashed") +
  facet_wrap(~ segment_count, scales = "free_y") +
  labs(
    x = "-log10(FDR)",
    y = "Number of random tiles",
    title = "Distribution of tile-level statistical signal",
    subtitle = paste0("Dashed line = FDR < ", fdr_threshold)
  ) +
  theme_bw()

ggsave(tile_signal_hist_pdf, p_tile_hist, width = 10, height = 6)

# plot 3: segment signal vs. segment size
tile_signal_width_pdf <- sub("\\.pdf$", "_tile_signal_vs_width.pdf", plot_pdf)

p_tile_width <- ggplot(tile_results_dt, aes(x = width, y = signal_score)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
  facet_wrap(~ segment_count, scales = "free") +
  labs(
    x = "Tile genomic width",
    y = "-log10(FDR)",
    title = "Tile-level statistical signal versus tile size",
    subtitle = paste0("Dashed line = FDR < ", fdr_threshold)
  ) +
  theme_bw()

ggsave(tile_signal_width_pdf, p_tile_width, width = 10, height = 6)

message("Done.")
message("Wrote: ", results_csv)
message("Wrote: ", tile_results_csv)
message("Wrote: ", plot_pdf)
message("Wrote: ", tile_signal_hist_pdf)
message("Wrote: ", tile_signal_width_pdf)