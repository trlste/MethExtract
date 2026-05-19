#!/usr/bin/env Rscript
# Apply a segmentation to a list of .beta tracks.
#
# For each input track, writes:
#   <name>_meth.bedGraph  - mean methylation per interval
#   <name>_cov.bedGraph   - total coverage (sum of total counts) per interval
#
# Also writes a combined TSV with all samples side by side.

suppressPackageStartupMessages({
  library(data.table)
})

# ---- configuration ----
seg_file   <- "output.bedgraph"
cpg_index  <- "wgbs_tools/references/hg38/CpG.bed.gz"
chrom_of_interest <- "chr2"
out_dir    <- "segmentation_out"
combined_tsv <- file.path(out_dir, "segmentation_methylation.tsv")

# Named list of .beta files. Names become the output prefixes.
beta_files <- list(
  TCD4 = "/net/dali/home/mscbio/yur28/directed_study/betadir/GSM5652296_Blood-T-Naive-CD4-Z0000041E.hg38.beta",
  TCD8 = "/net/dali/home/mscbio/yur28/directed_study/betadir/GSM5652298_Blood-T-Naive-CD8-Z0000041H.hg38.beta"
  # add more here: SampleX = "path/to/SampleX.beta", ...
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helpers ----
read_beta <- function(path) {
  fsize <- file.info(path)$size
  if (is.na(fsize)) stop("Cannot stat .beta file: ", path)
  if (fsize %% 2L != 0L) stop("Odd .beta file size, not a valid beta file: ", path)
  raw <- readBin(path, what = "raw", n = fsize)
  list(
    meth  = as.integer(raw[seq(1L, fsize, by = 2L)]),
    total = as.integer(raw[seq(2L, fsize, by = 2L)]),
    n     = fsize %/% 2L
  )
}

aggregate_sample <- function(beta, seg) {
  m_cum <- c(0L, cumsum(beta$meth))
  t_cum <- c(0L, cumsum(beta$total))
  has   <- !is.na(seg$first_cpg)

  sum_m <- rep(NA_integer_, nrow(seg))
  sum_t <- rep(NA_integer_, nrow(seg))
  sum_m[has] <- m_cum[seg$last_cpg[has] + 1L] - m_cum[seg$first_cpg[has]]
  sum_t[has] <- t_cum[seg$last_cpg[has] + 1L] - t_cum[seg$first_cpg[has]]

  mean_meth <- ifelse(is.na(sum_t) | sum_t == 0, NA_real_, sum_m / sum_t)
  list(mean = mean_meth, cov = sum_t)
}

write_bedgraph <- function(seg, values, path, track_name, description,
                           view_limits = NULL) {
  keep <- !is.na(values)
  vals <- values[keep]
  fmt  <- if (is.integer(vals) || all(vals == floor(vals))) "%d" else "%.6f"
  bg <- data.frame(
    chrom = seg$chrom[keep],
    start = seg$start[keep],
    end   = seg$end[keep],
    value = sprintf(fmt, vals)
  )
  header <- sprintf(
    'track type=bedGraph name="%s" description="%s" visibility=full %s',
    track_name,
    description,
    if (!is.null(view_limits))
      sprintf('autoScale=off viewLimits=%g:%g', view_limits[1], view_limits[2])
    else
      'autoScale=on'
  )
  writeLines(header, path)
  write.table(
    bg,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

# ---- 1. Load segmentation ----
seg <- fread(
  seg_file,
  sep = "\t",
  skip = 1,
  header = FALSE,
  col.names = c("chrom", "start", "end", "score")
)
seg <- seg[chrom == chrom_of_interest]
seg[, `:=`(start = as.integer(start), end = as.integer(end))]
setorder(seg, start)

# ---- 2. Load CpG index, restrict to chromosome ----
cpg <- fread(
  cmd = sprintf("zcat %s | awk '$1==\"%s\"'", cpg_index, chrom_of_interest),
  sep = "\t", header = FALSE,
  col.names = c("chrom", "start", "cpg_idx")
)
cpg[, start := start - 1L]
setorder(cpg, start)
cpg_pos <- cpg$start
cpg_idx <- cpg$cpg_idx

# print(file.info(beta_files[[1]])$size / 2)
# print(nrow(cpg))
# print(nrow(cpg) == file.info(beta_files[[1]])$size / 2)

# ---- 3. Map each interval to a CpG index range ----
lo <- findInterval(seg$start - 1L, cpg_pos, left.open = FALSE) + 1L
hi <- findInterval(seg$end   - 1L, cpg_pos, left.open = FALSE)

seg[, first_cpg := ifelse(lo > hi, NA_integer_, cpg_idx[lo])]
seg[, last_cpg  := ifelse(lo > hi, NA_integer_, cpg_idx[hi])]
seg[, n_cpg     := ifelse(lo > hi, 0L, hi - lo + 1L)]

# ---- 4. Process each track ----
out <- seg[, .(chrom, start, end, n_cpg)]

for (nm in names(beta_files)) {
  path <- beta_files[[nm]]
  message("Processing ", nm, " from ", path)

  # Derive output prefix from the .beta filename (strip directory and .beta extension)
  base_name <- sub("\\.beta$", "", basename(path))

  beta <- read_beta(path)
  agg  <- aggregate_sample(beta, seg)

  # Add to combined table (still keyed by short list name)
  out[, paste0(nm, "_mean") := round(agg$mean, 6)]
  out[, paste0(nm, "_cov")  := agg$cov]

  # Write per-track bedGraphs using full .beta basename
  meth_path <- file.path(out_dir, paste0(base_name, "_meth.bedGraph"))
  cov_path  <- file.path(out_dir, paste0(base_name, "_cov.bedGraph"))

  write_bedgraph(
    seg, agg$mean, meth_path,
    track_name  = paste0(nm, " mean methylation"),
    description = paste0("Mean methylation over segmentation intervals (", base_name, ")"),
    view_limits = c(0, 1)
  )
  write_bedgraph(
    seg, agg$cov, cov_path,
    track_name  = paste0(nm, " coverage"),
    description = paste0("Total coverage over segmentation intervals (", base_name, ")")
  )

  message("  wrote ", meth_path)
  message("  wrote ", cov_path)
}

# ---- 5. Write combined TSV ----
fwrite(out, combined_tsv, sep = "\t")
message("Wrote combined table: ", combined_tsv, " (", nrow(out), " intervals)")