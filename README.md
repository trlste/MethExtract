# MethExtract

This script is for extracting methylation ratio and counts information from .beta files.
If the genome is not initiated, please intiate the genome using wgbs_tools:
```
wgbs_tools init_genome hg38
```
Then extract the information from .beta files:
```
python3 extract_files.py /path/to/beta/XXXX.hg38.beta --indir /path/to/wgbs_tools/wgbs_tools/references/hg38 
```

# Compute distance matrix

To explore how different lambda influences the breakpoints, run compute_breakpoints.py and cluster_breakpoints.R sequentially for each sample:
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
