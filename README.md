# MethExtract

First extract methylation ratio and counts information from .beta files.
If genome is not initiated, intiate the genome using wgbs_tools:
```
wgbs_tools init_genome hg38
```
## Step 1
extract the information from .beta files:
```
python3 extract_files.py /path/to/beta/XXXX.hg38.beta --indir /path/to/wgbs_tools/wgbs_tools/references/hg38 
```
## Step 2
compute L0 segmentation via meth.R, we used lambda = 0.5. meth.R expects inputs to be in roadmap/:
```
Rscript meth.R
```
TODO: automatically collect the results of L0 segmentation into roadmap/segments/ .

## Step 3
Run scripts to obtain new global segmentation `output.bedgraph`
```
Rscript cosegmentation.R
Rscript cofrequency.R
Rscript visualize_cofreq.R
```
## Step 4
Apply new global segmentation to files. Results will be under segmentation_out/
```
Rscript resegment.R
```

# Evaluation of segmentation methods

## Step 1

Use R package bsseq to do differential methylation analysis, given cosegments or random tiles.
Prepare a list of samples involved `sample_sheet.csv` (please refer to the annotation).
```
Rscript diff_methyl_random.R # for random tiles
Rscript diff_methyl_coseg.R # for cosegments
```

# Step 2

Run R scripts to plot metrics of cosegments and random tiles
There are three components in the plot: cosegments, random tiles, reference `GSE186458_blocks.s207.hg38.bed.gz` (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458).
Adjust cut-offs in the scripts.
```
Rscript plot_three.R
```

# Compute distance matrix

Two ways to cluster breakpoints across samples.

To use only native R, run
```
Rscript cluster.R
```
change the path/file names in the script to run on different lambdas. This outputs a .bed file representing the clusters of breakpoints with distance cutoff at 5000.

To use a combination of pyBigWig / R and explore how different lambda influences the breakpoints, run compute_breakpoints.py and cluster_breakpoints.R sequentially for each sample:
```
# This will output the distance matrix of breakpoints
python compute_breakpoints.py -i roadmap/GSM5652176_WGBS_Seg_L0.1.bigwig -c 2 -o matrix/GSM5652176_WGBS_Seg_L0.1.csv

# This will output the cluster plot, the heatmap of distance matrix, and the .csv containing clusters cut by clusters <- cutree(hc, h = 10000)
Rscript cluster_breakpoints.R matrix/GSM5652176_WGBS_Seg_L0.1.csv matrix/GSM5652176_WGBS_Seg_L0.1
```

Here is a shell script to perform the loop:
```
# Perform matrix computing and clustering for lambda=[0.1, 0.2, 0.5, 1, 2, 5] in sample=["GSM5652176", "GSM5652296", "GSM5652297"]
# Change as you need
sh compute_matrix_loop.sh
```
