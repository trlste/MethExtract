library(data.table)

for (gran in c("40", "20", "10", "5", "2.5", "1.25", "0.625")) {
  
  seg_name <- paste0("output_", gran, ".bedgraph")
  seg <- fread(seg_name, header = TRUE, select = 1:3)
  
  colnames(seg) <- c("chr", "start", "end")
  
  seg$width <- seg$end - seg$start
  seg$log_width <- log10(seg$width)
  
  cut1 <- 40
  cut2 <- 3700
  
  seg$length_group <- cut(
    seg$width,
    breaks = c(-Inf, cut1, cut2, Inf),
    labels = c("short", "medium", "long")
  )
  
  table(seg$length_group)
  
  tapply(seg$width, seg$length_group, summary)
  
  for (subgroup in c("short", "medium", "long")) {
    x <- seg[seg$length_group == subgroup, c("chr", "start", "end")]
    
    write.table(
      x,
      paste0("output_", gran, "_", subgroup, ".bed"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
}
