args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_file <- args[2]

#input_file <- "GSM5652176_WGBS_Seg_L0.1.csv"
dm <- as.matrix(read.csv(input_file, header=FALSE))
dm_dist <- as.dist(dm)
hc <- hclust(dm_dist, method="average")

clusters_png <- paste0(out_file, "_clusters.png")
png(clusters_png, width=800, height=600)
plot(hc, main="Breakpoint Clustering", xlab="", sub="")
dev.off()

heatmap_png <- paste0(out_file, "_heatmap.png")
png(heatmap_png, width=800, height=800)
heatmap(dm)
dev.off()

clusters_cut <- paste0(out_file, "_clusters_cut.csv")
clusters <- cutree(hc, h = 10000)
write.csv(clusters, clusters_cut)