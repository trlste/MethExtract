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
cpg_bed_path   <- "CpG.bed.gz"

# Directory containing paired bigWig files:
#   *_FractionalMethylation.bigwig
#   *_ReadCoverage.bigwig
bigwig_dir <- "sample_1"

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
sample_sheet_path <- "sample_1/sample_1_sheet.csv"

# Chromosomes to include. Set to NULL to use all chromosomes in CpG BED.
chroms_to_use <- c("chr2")

coseg_bed_paths <- c(
  "output_20_per.bedgraph",
  "output_40_per.bedgraph"
)

# methylSig settings
fdr_threshold     <- 0.05
min_cov_per_cpg   <- 1
n_cores           <- 8

# Output files
results_csv <- "dm_results_cd4-cd8_0714_20_40.csv"

# Random seed for reproducibility
set.seed(42)

# ============================================================
# Helper functions
# ============================================================

read_coseg_regions <- function(path, chroms_to_use = NULL) {
  first_line <- readLines(path, n = 1)
  
  skip_n <- if (length(first_line) > 0 && grepl("^track", first_line[1])) {
    1
  } else {
    0
  }
  
  dt <- fread(
    path,
    skip = skip_n,
    header = FALSE,
    col.names = c("chrom", "start", "end", "score")
  )
  
  dt[, `:=`(
    start = as.integer(start),
    end   = as.integer(end),
    score = as.numeric(score)
  )]
  
  if (!is.null(chroms_to_use)) {
    dt <- dt[chrom %in% chroms_to_use]
  }
  
  dt <- dt[end > start]
  setorder(dt, chrom, start, end)
  
  gr <- GRanges(
    seqnames = dt$chrom,
    ranges = IRanges(
      start = dt$start + 1L, # not sure whether it's 0- or 1-based
      end = dt$end
    )
  )
  
  mcols(gr)$coseg_score <- dt$score
  
  gr
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
  
  # These values are methylation percentages, not 0-1 fractions.
  required_meth_cols <- c("meth_case", "meth_control")
  if (!all(required_meth_cols %in% meta_names)) {
    stop(
      "Could not find expected methylSig methylation columns: ",
      paste(required_meth_cols, collapse = ", "),
      ". Available columns are: ",
      paste(meta_names, collapse = ", ")
    )
  }
  
  tile_dt[, mean_methylation_case := meta[["meth_case"]]]
  tile_dt[, mean_methylation_control := meta[["meth_control"]]]
  
  if ("meth_diff" %in% meta_names) {
    tile_dt[, methylation_difference := meta[["meth_diff"]]]
  } else {
    tile_dt[, methylation_difference := mean_methylation_case - mean_methylation_control]
  }
  
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

message("Running differential methylation on pre-defined co-segmentation regions...")

results <- list()
tile_results <- list()
idx <- 1
tile_idx <- 1

for (coseg_path in coseg_bed_paths) {
  
  coseg_name <- tools::file_path_sans_ext(basename(coseg_path))
  
  message("Processing region file: ", coseg_path)
  
  seg_gr <- read_coseg_regions(
    path = coseg_path,
    chroms_to_use = chroms_to_use
  )
  
  message("  Number of regions: ", length(seg_gr))
  
  tile_dt <- run_diff_methylsig_tile_table(
    bs = bs,
    seg_gr = seg_gr,
    fdr_threshold = fdr_threshold,
    n_cores = n_cores
  )
  
  frac_sig <- mean(tile_dt$significant, na.rm = TRUE)
  
  results[[idx]] <- data.table(
    method = "coseg_predefined",
    coseg_file = coseg_name,
    n_regions = nrow(tile_dt),
    frac_sig = frac_sig,
    n_significant = sum(tile_dt$significant, na.rm = TRUE),
    mean_width = mean(tile_dt$width, na.rm = TRUE),
    median_width = median(tile_dt$width, na.rm = TRUE),
    mean_signal_score = mean(tile_dt$signal_score, na.rm = TRUE),
    median_signal_score = median(tile_dt$signal_score, na.rm = TRUE),
    mean_abs_methylation_difference = mean(abs(tile_dt$methylation_difference), na.rm = TRUE),
    median_abs_methylation_difference = median(abs(tile_dt$methylation_difference), na.rm = TRUE)
  )
  
  idx <- idx + 1
  
  tile_dt[, `:=`(
    method = "coseg_predefined",
    coseg_file = coseg_name,
    tile_id = seq_len(.N)
  )]
  
  tile_results[[tile_idx]] <- tile_dt
  tile_idx <- tile_idx + 1
}

results_dt <- rbindlist(results)
tile_results_dt <- rbindlist(tile_results)

message("Saving results...")

metrics_csv <- results_csv
tile_results_csv <- sub("\\.csv$", "_tile_level.csv", results_csv)

fwrite(results_dt, metrics_csv)
fwrite(tile_results_dt, tile_results_csv)

message("Done.")
message("Wrote metrics: ", metrics_csv)
message("Wrote tile-level results: ", tile_results_csv)
