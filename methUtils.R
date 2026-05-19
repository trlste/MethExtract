# Step 1: Fetch chromosome sizes for hg38
fetch_chrom_sizes <- function(assembly = "hg38") {
  url <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/", assembly, "/bigZips/", assembly, ".chrom.sizes")
  chrom_sizes <- read.table(url, header = FALSE, col.names = c("chrom", "size"))
  return(chrom_sizes)
}



mid_point_ranges=function(gr, chrom_end) {
  # Compute midpoints for each range
  midpoints <- start(gr) + (end(gr) - start(gr)) %/% 2
  
  # Create new start and end positions
  new_starts <- c(1, midpoints)+1
  new_ends <- c(midpoints, chrom_end)
  
  # Construct new GRanges object
  new_gr <- GRanges(
    seqnames = seqnames(gr)[1], # Assume all ranges are on the same chromosome
    ranges = IRanges(start = new_starts, end = new_ends)
  )
  
  return(new_gr)
}


aggregate_granges <- function(gr, startIndex, endIndex) {
  # Convert the aggregation index to a GRanges object
  agg_ranges <- GRanges(
    seqnames = seqnames(gr)[1], # Assume single chromosome for simplicity
    ranges = IRanges(
      start = start(gr)[startIndex],
      end = end(gr)[endIndex]
    )
  )
  agg_ranges
}

fill_gaps <- function(gr, chrom_length, proportional=T) {
  # Ensure the ranges are sorted
  gr <- sort(gr)
  
  grEndsNew=grEnds=end(gr)
  grStartsNew=grStarts=start(gr)
  
  # Adjust the first range to start at 1
  grStartsNew[1] <- 1
  
  # Adjust the last range to end at chrom_length
  grEndsNew[length(gr)] <- chrom_length
  
  # Iterate over the ranges and adjust the ends and starts
  for (i in seq_len(length(gr) - 1)) {
    if(proportional){
      
      gap=grStarts[i+1]-grEnds[i]
      widthLeft=grEnds[i]-grStarts[i]+1
      widthRight=grEnds[i+1]-grStarts[i+1]+1
      
      proportion=widthLeft/(widthLeft+widthRight)
      offset=floor(gap*proportion)
      midpoint=
        grEnds[i]+offset
    }
    else{
      # Calculate the midpoint of the gap between current and next range
      midpoint <- round((grEnds[i] + grStarts[i + 1]) / 2)
    }
    # Adjust the end of the current range and start of the next range
    grStartsNew[i+1] <- midpoint+1
    grEndsNew[i] <- midpoint
  }
    #print((dat<-cbind(grStartsNew, grEndsNew))[1:10,])

  start(gr)=grStartsNew
  end(gr)=grEndsNew

  gr
}

