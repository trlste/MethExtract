library(collections)

squash_cpg_edges_pq <- function(bed, n_islands) {
  # bed columns required: chrom, start, end, score
  # each row is an edge between two neighboring CpG islands
  # repeatedly squash the edge with the largest score
  # until only n_islands remain
  #
  # current #islands = current #edges + 1, for a single chain

  req <- c("chrom", "start", "end", "score")
  if (!is.data.frame(bed)) stop("bed must be a data.frame")
  if (!all(req %in% names(bed))) {
    stop("bed must contain columns: chrom, start, end, score")
  }

  if (length(n_islands) != 1L || is.na(n_islands) || n_islands < 1) {
    stop("n_islands must be a single positive integer")
  }
  n_islands <- as.integer(n_islands)

  k <- nrow(bed)
  initial_islands <- k + 1L

  if (n_islands > initial_islands) {
    stop(sprintf(
      "n_islands = %d is larger than the current number of islands = %d",
      n_islands, initial_islands
    ))
  }

  if (n_islands == initial_islands) {
    return(bed[, req, drop = FALSE])
  }

  # Keep labels as character, because merged labels become strings like "10296+10360"
  chrom <- as.character(bed$chrom)
  left  <- as.character(bed$start)
  right <- as.character(bed$end)
  score <- bed$score

  n_edges <- k

  # Doubly linked list over edges
  prev_edge <- c(NA_integer_, seq_len(k - 1L))
  next_edge <- c(seq_len(k - 1L) + 1L, NA_integer_)
  alive     <- rep(TRUE, k)

  alive_edges <- k

  # Max-priority queue of edge ids by score.
  # If your installed PriorityQueue is min-first instead of max-first,
  # change priority = score[i] to priority = -score[i],
  # and similarly adjust the optional tie-break below.
  pq <- priority_queue()

  # Optional deterministic tie-break: prefer smaller row index when scores tie
  # by adding a tiny decreasing adjustment.
  push_edge <- function(i) {
    pq$push(
      i,
      priority = score[i] - 1e-12 * i
    )
  }

  for (i in seq_len(k)) {
    push_edge(i)
  }

  make_label <- function(a, b) paste0(a, "+", b)

  # Number of squashes needed:
  # initial islands = k + 1
  # each squash removes 1 edge and reduces islands by 1
  target_edges <- n_islands - 1L

  while (alive_edges > target_edges) {
    # Pop until we find a live edge
    repeat {
      if (pq$size() == 0L) {
        stop("Priority queue exhausted before reaching target n_islands.")
      }
      e <- pq$pop()
      if (alive[e]) break
    }

    merged_label <- make_label(left[e], right[e])

    p <- prev_edge[e]
    q <- next_edge[e]

    # Update left neighboring edge: (... -> left[e]) becomes (... -> merged_label)
    if (!is.na(p)) {
      right[p]    <- merged_label
      next_edge[p] <- q
    }

    # Update right neighboring edge: (right[e] -> ...) becomes (merged_label -> ...)
    if (!is.na(q)) {
      left[q]     <- merged_label
      prev_edge[q] <- p
    }

    # Remove edge e
    alive[e] <- FALSE
    alive_edges <- alive_edges - 1L
  }

  # Reconstruct surviving edges in chain order
  live_ids <- which(alive)

  if (length(live_ids) == 0L) {
    # One island remains, so there are zero edges
    return(data.frame(
      chrom = character(0),
      start = character(0),
      end   = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  head <- live_ids[is.na(prev_edge[live_ids])]
  if (length(head) != 1L) {
    stop("Internal error: expected exactly one live head edge.")
  }

  out <- vector("list", length(live_ids))
  cur <- head[1L]
  idx <- 1L

  while (!is.na(cur)) {
    out[[idx]] <- data.frame(
      chrom = chrom[cur],
      start = left[cur],
      end   = right[cur],
      score = score[cur],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    cur <- next_edge[cur]
  }

  do.call(rbind, out)
}

track_line <- readLines("cosegmentation.bed", n = 1)

bed <- read.table(
  "cosegmentation.bed",
  sep = "\t",
  header = FALSE,
  skip = 1,
  stringsAsFactors = FALSE,
  col.names = c("chrom", "start", "end", "name", "score")
)

extra_value <- bed$name[1]

# Run the function
result <- squash_cpg_edges_pq(bed, n_islands = nrow(bed) %/% 10)

result$name <- extra_value

out_file <- "cofrequency.bed"

writeLines(track_line, out_file)

# Save as BED file
write.table(
  result,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)