library(rtracklayer)
library(GenomicRanges)
library(data.table)

# --- Arguments ---
bigwig_dir <- "roadmap/segments"
cpg_bed_path <- "wgbs_tools/references/hg38/CpG.bed.gz"
outfile <- "cosegmentation.bed"
chroms_to_use <- c("chr2")  # or NULL to use all

# --- Load CpG positions ---
cat("Loading CpG BED...\n")
cpg_bed <- fread(cpg_bed_path, col.names = c("chrom", "start", "end"))

bigwig_files <- list.files(bigwig_dir, pattern = "\\.bigwig$", full.names = TRUE)
cat(sprintf("Found %d bigWig files\n", length(bigwig_files)))

chroms <- if (!is.null(chroms_to_use)) chroms_to_use else unique(cpg_bed$chrom)

results <- list()

for (chrom_name in chroms) {
  cat(sprintf("Processing %s...\n", chrom_name))
  
  positions <- cpg_bed[chrom == chrom_name, start]

  #positions <- positions[1:10]

  n_cpgs <- length(positions)
  
  if (n_cpgs < 2) next
  
  co_count <- integer(n_cpgs - 1)

  #print(co_count)
  
  for (bw_path in bigwig_files) {
    bw <- import(bw_path, which = GRanges(chrom_name, IRanges(1, max(positions) + 1)))
    
    cpg_gr <- GRanges(chrom_name, IRanges(positions, positions + 1))
    hits <- findOverlaps(cpg_gr, bw)
    
    values <- rep(NA_real_, n_cpgs)
    values[queryHits(hits)] <- bw$score[subjectHits(hits)]

    #print(values)
    
    left  <- values[-n_cpgs]
    right <- values[-1]
    
    both_covered <- !is.na(left) & !is.na(right)
    same_segment <- both_covered & (left == right)
    
    co_count <- co_count + as.integer(same_segment)
  }
  
  results[[chrom_name]] <- data.table(
    chrom       = chrom_name,
    start       = positions[-n_cpgs],
    end         = positions[-1],
    co_count    = co_count
  )
}

out_dt <- rbindlist(results)
out_dt[, n_tracks := length(bigwig_files)]

# --- Write BED file ---
bed_out <- out_dt[, .(chrom, start, end, name = n_tracks, score = co_count)]

header <- paste0(
  'track name="CpG co-segmentation" ',
  'description="Score = number of counts where CpG_i and CpG_{i+1} ',
  'are in the same segment. Name = how many tracks used. Start/end spans the interval between two neighboring CpGs." ',
  'type=bedGraph'
)

writeLines(header, outfile)
fwrite(bed_out, outfile, sep = "\t", col.names = FALSE, append = TRUE)
cat(sprintf("Written to %s\n", outfile))
