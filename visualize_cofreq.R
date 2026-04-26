#!/usr/bin/env Rscript
# Convert quasi-BED co-segmentation file to IGV-compatible bedGraph.
#
# Usage:
#   Rscript quasi_bed_to_bedgraph.R input.bed [output.bedgraph]
#
#   If output path is omitted, writes to stdout.

parse_pos <- function(field) {
  as.integer(sub("\\+.*", "", field))
}

convert <- function(in_path, out_path = NULL) {
  lines <- readLines(in_path)

  out_lines <- character(0)

  for (line in lines) {
    # Pass through blank lines
    if (nchar(trimws(line)) == 0) {
      out_lines <- c(out_lines, line)
      next
    }

    # Rewrite track/browser header lines
    if (grepl("^track", line)) {
      out_lines <- c(out_lines, paste0(
        'track type=bedGraph name="CpG co-segmentation (score/total)"',
        ' description="Fraction of tracks where CpG_i and CpG_i+1 are co-segmented"',
        ' visibility=full autoScale=off viewLimits=0:1 color=31,119,180'
      ))
      next
    }
    if (grepl("^browser", line)) {
      out_lines <- c(out_lines, line)
      next
    }

    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 5) {
      message("WARNING: skipping malformed line: ", line)
      next
    }

    chrom <- parts[1]
    start <- parse_pos(parts[2])
    end   <- parse_pos(parts[3])
    score <- as.integer(parts[4])
    total <- as.integer(parts[5])

    if (is.na(total) || total == 0) {
      message("WARNING: total=0 or NA, skipping: ", line)
      next
    }
    if (start >= end) {
      message("WARNING: start >= end after parsing (", start, " >= ", end, "), skipping: ", line)
      next
    }

    value <- score / total
    out_lines <- c(out_lines, sprintf("%s\t%d\t%d\t%.6f", chrom, start, end, value))
  }

  if (is.null(out_path)) {
    cat(paste(out_lines, collapse = "\n"), "\n", sep = "")
  } else {
    writeLines(out_lines, out_path)
    message("Written to: ", out_path)
  }
}

in_path  <- "cofrequency.bed"
out_path <- "output.bedgraph"

if (!file.exists(in_path)) {
  stop("Input file not found: ", in_path)
}

convert(in_path, out_path)
