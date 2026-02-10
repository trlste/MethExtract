library(rtracklayer)
library(GenomicRanges)

current_path <- getwd()

files <- list.files(
  path = paste0(current_path, "/roadmap"),
  pattern = "*_Seg_L0\\.5\\.bigwig$",
  full.names = TRUE
)
print(files)

all_bps <- unlist(lapply(files, function(f) {
  gr <- import(f)
  end(gr)[-length(gr)]
}), use.names = FALSE)

#bps <- end(gr)[-length(gr)]
D <- dist(all_bps, method = "manhattan")
hc <- hclust(D)
hotspot_clusters <- cutree(hc, h = 5000)

# Step 4: Summarize clusters as genomic ranges
hotspot_df <- data.frame(
  start = tapply(all_bps, hotspot_clusters, min),
  end   = tapply(all_bps, hotspot_clusters, max),
  count = as.numeric(table(hotspot_clusters))
)

hotspot_gr <- GRanges(
  seqnames = "chr2",  # change if multiple chromosomes
  ranges = IRanges(start = hotspot_df$start, end = hotspot_df$end),
  count = hotspot_df$count
)

# Optional: Save as BED
out_file  <- file.path(current_path, "clusters_L0.5_adipocytes_T-CD4_T-CD8.bed")
export(hotspot_gr, con=out_file, format="Bed")
