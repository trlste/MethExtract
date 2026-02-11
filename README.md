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
