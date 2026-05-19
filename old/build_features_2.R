############################################################
# Segment methylation values using tile_by_regions
############################################################

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install(c("methylSig","bsseq","GenomicRanges"), ask = FALSE)

suppressPackageStartupMessages({
    library(methylSig)
    library(bsseq)
    library(GenomicRanges)
    library(rtracklayer)
    library(BiocParallel)
    library(data.table)
})

############################################################
# Function to read binary .beta files
############################################################

read_beta_file <- function(path, n_cpg, lbeta = FALSE) {
    fname <- path
    N <- file.info(fname)$size
    mat <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)

    #dtype <- "integer"

    #con <- file(path, "rb")

    #raw_data <- readBin(
    #    con,
    #    what = dtype,
    #    n = 2 * n_cpg,
    #    size = if (lbeta) 2 else 1,
    #    signed = FALSE
    #)

    #close(con)

    #mat <- matrix(raw_data, ncol = 2, byrow = TRUE)

    meth  <- as.numeric(mat[,1])
    total <- as.numeric(mat[,2])

    list(meth = meth, total = total)
}

############################################################
# Build BSseq object from beta files
############################################################

build_bsseq_from_beta <- function(beta_files, cpg_bed, lbeta = FALSE) {

    n_cpg <- length(cpg_bed)
    n_samples <- length(beta_files)

    M   <- matrix(0, nrow = n_cpg, ncol = n_samples)
    Cov <- matrix(0, nrow = n_cpg, ncol = n_samples)

    for (i in seq_along(beta_files)) {

        vals <- read_beta_file(beta_files[i], n_cpg, lbeta)

        M[,i]   <- vals$meth
        Cov[,i] <- vals$total
    }

    bs <- BSseq(
        chr = as.character(seqnames(cpg_bed)),
        pos = start(cpg_bed),
        M = M,
        Cov = Cov,
        sampleNames = basename(beta_files)
    )

    return(bs)
}

############################################################
# Compute segment methylation values
############################################################

get_segment_values_tile <- function(beta_files, cpg_bed, segments, lbeta = FALSE) {

    message("Building BSseq object...")
    bs <- build_bsseq_from_beta(beta_files, cpg_bed, lbeta)

    message("Aggregating CpGs into segments...")
    tiled <- tile_by_regions(
        bs = bs,
        gr = segments
    )

    # Count NAs
    sum(is.na(tiled))
    #sum(is.na(seg_meth))

    # See which regions have any NA (row-wise)
    #any_na_beta <- apply(seg_beta, 1, function(x) any(is.na(x)))
    #sum(any_na_beta)  # number of regions with at least one NA sample

    # If you have multiple samples, see NA count per sample (column-wise)
    #colSums(is.na(seg_beta))
    #colSums(is.na(seg_meth))

    message("Extracting region counts...")

    region_M   <- getCoverage(tiled, type = "M")
    region_Cov <- getCoverage(tiled, type = "Cov")

    region_beta <- region_M / region_Cov

    list(
        meth  = region_M,
        total = region_Cov,
        beta  = region_beta,
        bs_tiled = tiled
    )
}

############################################################
# Example usage
############################################################

beta_folder  <- "betas"
cpg_bed_path <- "./wgbs_tools/references/hg38/CpG.bed.gz"

#cpg_bed <- import(cpg_bed_path)
cpg_df  <- fread("./wgbs_tools/references/hg38/CpG.bed.gz",
                  header = FALSE, sep = "\t",
                  col.names = c("chrom", "start", "end"))

# BED is 0-based, convert start to 1-based for GRanges
cpg_bed <- GRanges(
    seqnames = cpg_df$chrom,
    ranges   = IRanges(start = cpg_df$start + 1,
                       end   = cpg_df$start + 1)  # single-base sites
)

ext         <- "\\.beta$"
beta_files  <- list.files(beta_folder, pattern = ext, full.names = TRUE)

seg_df   <- read.table("kNN_Components_clusters.bed", header = FALSE, sep = "\t",
                        col.names = c("chrom", "start", "end", "V4", "V5"))
segments <- GRanges(seqnames = seg_df$chrom,
                    ranges   = IRanges(start = seg_df$start, end = seg_df$end))

print("done loading segments")

# Compute segment methylation
seg_results <- get_segment_values_tile(
    beta_files = beta_files,
    cpg_bed = cpg_bed,
    segments = segments,
    lbeta = FALSE
)

# Access results
seg_beta  <- seg_results$beta
seg_meth  <- seg_results$meth
seg_total <- seg_results$total

# Count NAs
#sum(is.na(seg_beta))
#sum(is.na(seg_meth))

# See which regions have any NA (row-wise)
#any_na_beta <- apply(seg_beta, 1, function(x) any(is.na(x)))
#sum(any_na_beta)  # number of regions with at least one NA sample

# If you have multiple samples, see NA count per sample (column-wise)
#colSums(is.na(seg_beta))
#colSums(is.na(seg_meth))

# Find a region with NaN beta in sample 1
na_idx <- which(is.nan(seg_beta[, 1]))[1]  # first NaN region

# Check its total coverage in sample 1
seg_total[na_idx, 1]  # should be 0

# And its meth count
seg_meth[na_idx, 1]   # should also be 0

print("Indices:")
segments[na_idx]
print("For the following samples, this segment is NA:")
which(is.nan(seg_beta[na_idx, ]))

assays <- list(
    beta  = seg_beta,
    M     = seg_meth,
    total = seg_total
)

# Save to disk — reload later with: assays <- readRDS("assays.rds")
saveRDS(assays, file = "assays.rds")