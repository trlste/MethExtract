#!/usr/bin/env Rscript
# Build a global methylation segmentation directly from one or more wgbs_tools
# .beta files. This replaces MethExtract README steps 1--3.
#
# Example:
#   Rscript run_pipeline.R \
#     --beta /data/A.hg38.beta /data/B.hg38.beta \
#     --cpg-bed wgbs_tools/references/hg38/CpG.bed.gz \
#     --chrom chr2 --lambda 0.5 --keep-frac 0.10 \
#     --output output.bedgraph
#
# Add --debug to save the per-sample L0 segments, raw co-segmentation scores,
# and the chosen boundaries. The final --output is always written.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(l01segmentation)
})

usage <- function(status = 0L) {
  message(
    "Usage:\n",
    "  Rscript run_pipeline.R --beta FILE [FILE ...] --cpg-bed CpG.bed.gz [--chrom chr2] [--lambda 0.5] [--keep-frac 0.10 | --target-islands N] [--all-covered-cpgs] [--output output.bedgraph] [--debug] [--debug-dir pipeline_debug]\n\n",
    "--beta accepts one or more wgbs_tools .beta files.\n",
    "--keep-frac is the fraction of original CpG islands retained after co-segmentation.\n",
    "It defaults to 0.10, matching cofrequency.R's nrow(bed) %/% 10 choice.\n",
    "By default, segmentation uses every second covered CpG, matching meth.R; use --all-covered-cpgs to opt out."
  )
  quit(status = status)
}

parse_args <- function(argv) {
  out <- list(
    beta = character(), cpg_bed = NULL, chrom = "chr2", lambda = 0.5,
    keep_frac = 0.10, target_islands = NULL, output = "output.bedgraph",
    debug = FALSE, debug_dir = "pipeline_debug", every_other_covered_cpg = TRUE
  )
  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("-h", "--help")) usage()
    if (arg == "--debug") {
      out$debug <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--all-covered-cpgs") {
      out$every_other_covered_cpg <- FALSE
      i <- i + 1L
      next
    }
    if (arg == "--beta") {
      i <- i + 1L
      start <- i
      while (i <= length(argv) && !startsWith(argv[[i]], "--")) i <- i + 1L
      if (i == start) stop("--beta requires at least one file.")
      out$beta <- c(out$beta, argv[start:(i - 1L)])
      next
    }
    if (!arg %in% c("--cpg-bed", "--chrom", "--lambda", "--keep-frac",
                    "--target-islands", "--output", "--debug-dir")) {
      stop("Unknown argument: ", arg)
    }
    if (i == length(argv)) stop(arg, " requires a value.")
    value <- argv[[i + 1L]]
    switch(arg,
      "--cpg-bed" = out$cpg_bed <- value,
      "--chrom" = out$chrom <- value,
      "--lambda" = out$lambda <- as.numeric(value),
      "--keep-frac" = out$keep_frac <- as.numeric(value),
      "--target-islands" = out$target_islands <- as.integer(value),
      "--output" = out$output <- value,
      "--debug-dir" = out$debug_dir <- value
    )
    i <- i + 2L
  }

  if (!length(out$beta)) stop("At least one --beta file is required.")
  if (is.null(out$cpg_bed)) stop("--cpg-bed is required.")
  if (!all(file.exists(out$beta))) stop("Missing beta file(s): ", paste(out$beta[!file.exists(out$beta)], collapse = ", "))
  if (!file.exists(out$cpg_bed)) stop("CpG BED does not exist: ", out$cpg_bed)
  if (is.na(out$lambda) || out$lambda <= 0) stop("--lambda must be positive.")
  if (is.na(out$keep_frac) || out$keep_frac <= 0 || out$keep_frac > 1) stop("--keep-frac must be in (0, 1].")
  if (!is.null(out$target_islands) && (is.na(out$target_islands) || out$target_islands < 1L)) {
    stop("--target-islands must be a positive integer.")
  }
  out
}

read_beta <- function(path) {
  n_bytes <- file.info(path)$size
  if (is.na(n_bytes) || n_bytes %% 2L != 0L) stop("Invalid .beta file: ", path)
  x <- readBin(path, what = "raw", n = n_bytes)
  list(
    meth = as.integer(x[seq.int(1L, n_bytes, by = 2L)]),
    total = as.integer(x[seq.int(2L, n_bytes, by = 2L)])
  )
}

load_cpg_chrom <- function(path, chrom) {
  # Reading only the requested chromosome avoids loading the full hg38 CpG index.
  cmd <- sprintf("zcat %s | awk '$1 == \"%s\"'", shQuote(normalizePath(path)), chrom)
  x <- fread(cmd = cmd, header = FALSE)
  if (ncol(x) < 3L || !nrow(x)) stop("No CpGs found for ", chrom, " in ", path)
  setnames(x, names(x)[1:3], c("chrom", "start", "cpg_idx"))
  x <- x[, .(chrom = as.character(chrom), start = as.integer(start) - 1L,
             cpg_idx = as.integer(cpg_idx))]
  if (anyNA(x) || any(x$cpg_idx < 1L)) stop("CpG BED must have a 1-based CpG index in column 3.")
  setorder(x, start)
  x
}

aggregate_granges <- function(gr, start_index, end_index) {
  GRanges(
    seqnames = seqnames(gr)[1L],
    ranges = IRanges(start = start(gr)[start_index], end = end(gr)[end_index])
  )
}

segment_sample <- function(beta, cpg, lambda, every_other_covered_cpg) {
  if (max(cpg$cpg_idx) > length(beta$total)) {
    stop("The beta file is shorter than the CpG index for the requested chromosome.")
  }
  meth <- beta$meth[cpg$cpg_idx]
  cov <- beta$total[cpg$cpg_idx]
  kept <- which(cov > 0L)
  # meth.R explicitly retains its 2nd, 4th, ... imported covered CpG entries.
  if (every_other_covered_cpg) kept <- kept[seq.int(2L, length(kept), by = 2L)]
  if (length(kept) < 2L) stop("Fewer than two covered CpGs on ", cpg$chrom[1L], ".")

  cpg_gr <- GRanges(cpg$chrom[1L], IRanges(start = cpg$start[kept] + 1L, width = 1L))
  frac <- meth[kept] / cov[kept]
  spacing <- diff(start(cpg_gr))
  # Same inverse-distance weighting used by meth.R; guard against duplicate sites.
  weights <- c(1, 1 / pmax(1, as.numeric(spacing)))
  fit <- fusedsegmentation(
    meth[kept], lambda2 = lambda, C = cov[kept], weight = weights,
    objective = "binomial"
  )
  seg <- aggregate_granges(cpg_gr, fit$start, fit$end)
  mcols(seg)$score <- as.numeric(fit$value)
  mcols(seg)$segment_id <- seq_along(seg)
  list(segments = seg, covered = kept, fraction = frac, coverage = cov[kept])
}

co_segment_edges <- function(cpg, segmentations) {
  n_cpg <- nrow(cpg)
  cpg_gr <- GRanges(cpg$chrom[1L], IRanges(start = cpg$start + 1L, width = 1L))
  co_count <- integer(n_cpg - 1L)

  for (sample_name in names(segmentations)) {
    seg <- segmentations[[sample_name]]$segments
    hits <- findOverlaps(cpg_gr, seg, type = "within")
    # cosegmentation.R compares the scores exported in each L0 BigWig, rather
    # than an internal segment ID. Preserve that comparison exactly.
    values <- rep(NA_real_, n_cpg)
    values[queryHits(hits)] <- mcols(seg)$score[subjectHits(hits)]
    left <- values[-n_cpg]
    right <- values[-1L]
    co_count <- co_count + as.integer(!is.na(left) & !is.na(right) & left == right)
  }

  data.table(
    chrom = cpg$chrom[-n_cpg], start = cpg$start[-n_cpg], end = cpg$start[-1L],
    co_count = co_count, n_tracks = length(segmentations)
  )
}

choose_merged_edges <- function(edges, n_islands) {
  # cofrequency.R's string labels are only bookkeeping for connected
  # components. Scores never change during merging, so direct selection gives
  # the same merge decisions without serializing intermediate labels.
  n_edges <- nrow(edges)
  if (n_islands < 1L || n_islands > n_edges + 1L) {
    stop("Invalid target number of islands.")
  }
  n_to_merge <- n_edges + 1L - n_islands
  merged <- rep(FALSE, n_edges)
  if (n_to_merge) {
    merged[head(order(-edges$co_count, seq_len(n_edges)), n_to_merge)] <- TRUE
  }
  merged
}

build_legacy_output_intervals <- function(edges, merged, n_tracks) {
  # Each non-merged edge is preserved by cofrequency.R. Its left and right
  # labels resolve to the starts of adjacent merged components. Construct those
  # same intervals directly, without materializing the labels themselves.
  n_edges <- nrow(edges)
  cpg_positions <- c(edges$start[1L], edges$end)
  boundaries <- which(!merged)
  if (!length(boundaries)) {
    return(data.table(chrom = character(), start = integer(), end = integer(),
                      score = numeric()))
  }
  left_component_start <- c(1L, head(boundaries, -1L) + 1L)
  data.table(
    chrom = edges$chrom[1L],
    start = cpg_positions[left_component_start],
    end = cpg_positions[boundaries + 1L],
    score = edges$co_count[boundaries] / n_tracks
  )
}

write_final_bedgraph <- function(blocks, path) {
  final <- blocks[start < end, .(chrom, start, end, score)]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    'track type=bedGraph name="CpG co-segmentation (score/total)" description="Fraction of samples where adjacent CpGs were co-segmented" visibility=full autoScale=off viewLimits=0:1 color=31,119,180',
    path
  )
  final_to_write <- copy(final)
  final_to_write[, score := sprintf("%.6f", score)]
  fwrite(final_to_write, path, sep = "\t", col.names = FALSE, append = TRUE)
  final
}

write_debug <- function(debug_dir, cpg, segmentations, edges, merged, blocks) {
  dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
  for (sample_name in names(segmentations)) {
    x <- segmentations[[sample_name]]
    frac <- data.table(
      chrom = cpg$chrom[x$covered], start = cpg$start[x$covered],
      end = cpg$start[x$covered] + 1L, fractional_methylation = x$fraction,
      coverage = x$coverage
    )
    seg <- x$segments
    l0 <- data.table(
      chrom = as.character(seqnames(seg)), start = start(seg) - 1L, end = end(seg),
      segment_id = mcols(seg)$segment_id, methylation = mcols(seg)$score
    )
    fwrite(frac, file.path(debug_dir, paste0(sample_name, "_covered_cpgs.tsv")), sep = "\t")
    fwrite(l0, file.path(debug_dir, paste0(sample_name, "_L0_segments.bed")), sep = "\t")
  }
  fwrite(edges, file.path(debug_dir, "cosegmentation_raw.tsv"), sep = "\t")
  fwrite(cbind(edges, merged = merged), file.path(debug_dir, "merge_decisions.tsv"), sep = "\t")
  fwrite(blocks, file.path(debug_dir, "global_blocks.tsv"), sep = "\t")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  cpg <- load_cpg_chrom(args$cpg_bed, args$chrom)
  sample_names <- make.unique(sub("\\.beta$", "", basename(args$beta)))
  names(args$beta) <- sample_names

  message("Loaded ", nrow(cpg), " CpGs on ", args$chrom, ".")
  segmentations <- vector("list", length(args$beta))
  names(segmentations) <- sample_names
  for (sample_name in sample_names) {
    message("Segmenting ", sample_name, " ...")
    segmentations[[sample_name]] <- segment_sample(
      read_beta(args$beta[[sample_name]]), cpg, args$lambda,
      args$every_other_covered_cpg
    )
    message("  produced ", length(segmentations[[sample_name]]$segments), " L0 segments.")
  }

  message("Computing adjacent-CpG co-segmentation ...")
  edges <- co_segment_edges(cpg, segmentations)
  n_islands <- args$target_islands %||% max(1L, floor(nrow(edges) * args$keep_frac))
  n_islands <- min(n_islands, nrow(edges) + 1L)
  message("Squashing to ", n_islands, " islands ...")
  merged <- choose_merged_edges(edges, n_islands)
  blocks <- build_legacy_output_intervals(edges, merged, length(segmentations))
  final <- write_final_bedgraph(blocks, args$output)
  if (args$debug) write_debug(args$debug_dir, cpg, segmentations, edges, merged, blocks)
  message("Wrote ", args$output, " (", nrow(final), " intervals).")
}

`%||%` <- function(x, y) if (is.null(x)) y else x
main()
