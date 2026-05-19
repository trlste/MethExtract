library(rtracklayer)
source("methUtils.R")
#devtools::install_github("https://github.com/boooooogey/l01segmentation")

library(l01segmentation)

cell_lines=sub("_.*", "", list.files(pattern="*Fraction*",path = "roadmap/"))


# Get hg38 chromosome sizes
chrom_sizes <- fetch_chrom_sizes()
par(mfrow=c(3,2))
for(cell_line in cell_lines){
  for (lambda in c(0.5)){
    
    dir <- "roadmap"
    frac_file <- file.path(dir, paste0(cell_line, "_WGBS_FractionalMethylation.bigwig"))
    cov_file  <- file.path(dir, paste0(cell_line, "_WGBS_ReadCoverage.bigwig"))
    
    
    #do for a single chromosome
    i=2
    chrom <- chrom_sizes$chrom[i]
    chrom_length <- chrom_sizes$size[i]
    
    # Define the chromosome range
    region <- GRanges(chrom, IRanges(1, chrom_length))
    data_frac <- import(frac_file, which = region)
    data_cov <- import(cov_file, which = region)
    
    # the scores appear twice so we take every other one
    ii=seq(2,length(data_frac$score),2) 
    data_frac=data_frac[ii,]
    data_cov=data_cov[ii,]
    meth=data_frac$score*data_cov$score
    
    #compute the distance between sites
    dd=diff(data_frac@ranges@start)
    
    #perform segmentation
    #lambda2 controls  the number of segments
    #weights are optional, weighing here by invese distance
    out=fusedsegmentation(meth, lambda2 = lambda, C=data_cov$score,weight = c(1,1/dd), objective = "binomial")
    out_file  <- file.path(dir, paste0(cell_line, "_WGBS_chr2_Seg", "_L", lambda,".bigwig"))
    
    #compute the numer of singletons
    sum((out$end-out$start) ==1)
    
    #collapse ranges that got the same value
    aggRanges=aggregate_granges(data_frac, out$start, out$end)
    length(aggRanges)
    
    aggRanges$score=c(out$value)

    seqlengths(aggRanges) <- chrom_sizes[match(seqlevels(aggRanges), chrom_sizes[,1]),2]
    hist(aggRanges$score, 200, main=lambda)
    export(aggRanges, con=out_file, format="BigWig")

    #fill the gaps to a continuous track
    #regionsMid=fill_gaps(aggRanges, chrom_length)
    
    #regionsMid$score=c(out$value)
    
    #regionsMid[1:10,]
    #seqlengths(regionsMid) <- chrom_sizes[match(seqlevels(regionsMid), chrom_sizes[,1]),2]
    #hist(regionsMid$score, 200,main=lambda)
    #print(regionsMid)
    #export(regionsMid, con=out_file, format="BigWig")
  }
}
