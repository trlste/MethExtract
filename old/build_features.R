library(rtracklayer)
library(BiocParallel)

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
beta_folder  <- "betas"
cpg_bed_path <- "CpG.bed"
lbeta        <- FALSE   # set TRUE if using .lbeta (uint16) format

# ----------------------------------------------------------------------
# Load CpG reference positions (one row per CpG site, genome-ordered)
# ----------------------------------------------------------------------
cpg_bed <- import(cpg_bed_path)

print("done importing cpg bed file")

# ----------------------------------------------------------------------
# Load BED segments (5-column format: chrom, start, end, + 2 unused columns)
# ----------------------------------------------------------------------
seg_df   <- read.table("kNN_Components_clusters.bed", header = FALSE, sep = "\t",
                        col.names = c("chrom", "start", "end", "V4", "V5"))
segments <- GRanges(seqnames = seg_df$chrom,
                    ranges   = IRanges(start = seg_df$start, end = seg_df$end))
seg_names <- paste0(seqnames(segments), ":", start(segments), "-", end(segments))

print("done loading segments")
# ----------------------------------------------------------------------
# Read one .beta (or .lbeta) file and return meth + total vectors
# ----------------------------------------------------------------------
read_beta_file <- function(path, lbeta = FALSE) {
    dtype <- if (lbeta) "integer" else "raw"  # raw = uint8, integer = uint16

    con <- file(path, "rb")
    raw_data <- readBin(con, what = dtype, n = 2 * length(cpg_bed), size = if (lbeta) 2L else 1L, signed = FALSE)
    close(con)

    mat   <- matrix(raw_data, ncol = 2, byrow = TRUE)
    meth  <- as.numeric(mat[, 1])
    total <- as.numeric(mat[, 2])
    list(meth = meth, total = total)
}

# ----------------------------------------------------------------------
# Compute mean meth, total, and beta over segments for one .beta file
# ----------------------------------------------------------------------
get_segment_values <- function(beta_path, lbeta = FALSE) {
    vals  <- read_beta_file(beta_path, lbeta)
    meth  <- vals$meth
    total <- vals$total

    # CpGs with total == 0 are missing (NaN), consistent with wgbstools
    covered <- total > 0
    beta_vals <- ifelse(covered, meth / total, NA_real_)

    # Find overlaps between CpG sites and segments
    hits <- findOverlaps(segments, cpg_bed)

    seg_meth  <- rep(NA_real_, length(segments))
    seg_total <- rep(NA_real_, length(segments))
    seg_beta  <- rep(NA_real_, length(segments))

    for (i in seq_along(segments)) {
        cpg_idx <- subjectHits(hits[queryHits(hits) == i])
        cpg_idx <- cpg_idx[covered[cpg_idx]]   # restrict to covered CpGs
        if (length(cpg_idx) == 0) next
        seg_meth[i]  <- mean(meth[cpg_idx],      na.rm = TRUE)
        seg_total[i] <- mean(total[cpg_idx],     na.rm = TRUE)
        seg_beta[i]  <- mean(beta_vals[cpg_idx], na.rm = TRUE)
    }

    list(meth = seg_meth, total = seg_total, beta = seg_beta)
}

# ----------------------------------------------------------------------
# Build matrices across all .beta files
# ----------------------------------------------------------------------
ext         <- if (lbeta) "\\.lbeta$" else "\\.beta$"
beta_files  <- list.files(beta_folder, pattern = ext, full.names = TRUE)
beta_files
sample_names <- sub("_.*", "", basename(beta_files))   # mirrors Python prefix logic

bpparam <- MulticoreParam(workers = parallel::detectCores() - 1)
results <- bplapply(beta_files, get_segment_values, lbeta = lbeta, BPPARAM = bpparam)

build_matrix <- function(slot) {
    mat <- do.call(cbind, lapply(results, `[[`, slot))
    colnames(mat) <- sample_names
    rownames(mat) <- seg_names
    mat
}

assays <- list(
    beta  = build_matrix("beta"),
    M     = build_matrix("meth"),
    total = build_matrix("total")
)

assays

# Save to disk — reload later with: assays <- readRDS("assays.rds")
saveRDS(assays, file = "assays.rds")