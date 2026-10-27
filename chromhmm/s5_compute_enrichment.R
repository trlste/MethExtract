library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
lambda <- args[1]

grans <- c("40", "20", "10", "5", "2.5", "1.25", "0.625")
groups <- c("short", "medium", "long")

# ChromHMM state order
state_levels <- c(
  "1_TssA",
  "2_TssAFlnk",
  "3_TxFlnk",
  "4_Tx",
  "5_TxWk",
  "6_EnhG",
  "7_Enh",
  "8_ZNF/Rpts",
  "9_Het",
  "10_TssBiv",
  "11_BivFlnk",
  "12_EnhBiv",
  "13_ReprPC",
  "14_ReprPCWk",
  "15_Quies"
)

summarize_overlap <- function(file, total_bp) {

  x <- fread(file, header = FALSE)

  # V7 = ChromHMM state
  # V8 = overlap bp
  out <- x[V8 > 0,
           .(overlap_bp = sum(V8)),
           by = .(state = V7)]

  # Add states that happened to have zero overlap
  out <- merge(
    data.table(state = state_levels),
    out,
    by = "state",
    all.x = TRUE
  )

  out[is.na(overlap_bp), overlap_bp := 0]

  out[, fraction := overlap_bp / total_bp]

  out
}

get_total_bp <- function(file) {
  x <- fread(file, header = FALSE, select = 1:3)
  sum(x$V3 - x$V2)
}



for (gran in grans) {
  
  total_bp <- sapply(
    groups,
    function(g)
      get_total_bp(
        paste0("output_", gran, "_", g, ".bed")
      )
  )
  
  total_bp
  
  observed_list <- lapply(groups, function(g) {
    
    #overlap/observed/gran_${gran}/${group}_E038.tsv
    x <- summarize_overlap(
      paste0("overlap/observed/gran_", gran, "/", g, "_E038.tsv"),
      total_bp[g]
    )
    
    x[, group := g]
    x
  })
  
  observed <- rbindlist(observed_list)
  
  random_list <- list()
  
  for (g in groups) {
    
    #overlap/random/gran_${gran}/${group}/random_${i}_E038.tsv
    files <- list.files(
      paste0("overlap/random/gran_", gran, "/", g),
      pattern = "\\.tsv$",
      full.names = TRUE
    )
    
    for (i in seq_along(files)) {
      
      x <- summarize_overlap(files[i], total_bp[g])
      
      x[, `:=`(
        group = g,
        randomization = i
      )]
      
      random_list[[length(random_list) + 1]] <- x
    }
  }
  
  random_dt <- rbindlist(random_list)
  
  table(random_dt$group, random_dt$randomization)
  
  random_summary <- random_dt[, .(
    random_mean = mean(fraction),
    random_sd   = sd(fraction),
    
    random_q025 = quantile(fraction, 0.025),
    random_q50  = quantile(fraction, 0.50),
    random_q975 = quantile(fraction, 0.975)
  ), by = .(group, state)]
  
  result <- merge(
    observed,
    random_summary,
    by = c("group", "state")
  )
  
  result[, fold_enrichment := fraction / random_mean]
  result[, log2_enrichment := log2(fraction / random_mean)]
  
  p_list <- list()
  
  for (g in groups) {
    
    for (s in state_levels) {
      
      obs <- observed[group == g & state == s, fraction]
      rnd <- random_dt[group == g & state == s, fraction]
      
      p_list[[length(p_list) + 1]] <- data.table(
        group = g,
        state = s,
        p_enrich = (1 + sum(rnd >= obs)) / (length(rnd) + 1),
        p_deplete = (1 + sum(rnd <= obs)) / (length(rnd) + 1)
      )
    }
  }
  
  empirical_p <- rbindlist(p_list)
  
  result <- merge(result, empirical_p, by = c("group", "state"))
  result[, group := factor(group, levels = c("short", "medium", "long"))]
  result[, state := factor(state, levels = rev(state_levels))]
  
  gran_title <- paste0("ChromHMM enrichment by segment length")
  gran_subtitle <- paste0(gran, "%, L = ", lambda, "E038")
  gran_save_file <- paste0("chromhmm_enrichment_L", lambda, "_", gran, "_E038.pdf")
    
  p <- ggplot(
    result,
    aes(
      x = group,
      y = state,
      fill = log2_enrichment
    )
  ) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", log2_enrichment)), size = 4.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = gran_title,
      subtitle = gran_subtitle,
      x = "Segment length group",
      y = "ChromHMM state",
      fill = "log2 enrichment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  ggsave(gran_save_file, p, width = 8, height = 10, dpi = 300)
}

