#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pyBigWig
import pandas as pd
import re

# assuming we have run wgbstools init_genome
def main():
    parser = argparse.ArgumentParser(
        description="Convert CpG beta file to methylation and coverage bigWig tracks"
    )

    parser.add_argument(
        "--beta",
        required=True,
        help="Input beta file (binary: methylated count, total coverage)"
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Input directory containing CpG.bed.gz, and CpG.chrome.size. If you have run init_genome, \
            probably under [INSTALLATION_PATH]/wgbs_tools/references/[GENOME]"
    )

    parser.add_argument(
        "--outdir",
        default=f'{os.getcwd()}/roadmap',
        help="Output directory"
    )

    args = parser.parse_args()

    outdir = args.outdir if args.outdir else f'{os.getcwd()}/roadmap'
    os.makedirs(args.outdir, exist_ok=True)
    indir = args.indir
    cpg_bed_path = os.path.join(indir, "CpG.bed.gz")
    chrom_sizes_path = os.path.join(indir, "CpG.chrome.size")

    # -------------------------------
    # Derive prefix from beta filename
    # -------------------------------
    prefix = re.sub(r"_.*", "", os.path.basename(args.beta))
    print(prefix)

    # ----------------------------------------------------
    # 1. Read beta file
    # ----------------------------------------------------
    content = np.fromfile(args.beta, dtype=np.uint8).reshape((-1, 2))
    print(f"Number of CpGs in {args.beta}: {content.shape[0]}")
    meth_count = content[:, 0].astype(float)
    total_cov = content[:, 1].astype(float)

    print(f"Coverage summary for {args.beta} using reference {args.indir}:")
    print(f"  Mean coverage: {total_cov.mean():.2f}")
    print(f"  Max coverage:  {total_cov.max():.0f}")
    print(f"  CpGs with coverage > 0: {(total_cov > 0).sum()}")
    print()

    # Calculate methylation fraction
    meth_frac = np.zeros_like(meth_count)
    mask = total_cov > 0
    meth_frac[mask] = meth_count[mask] / total_cov[mask]

    print(f"Methylation fraction summary (covered CpGs only):")
    if mask.any():
        print(f"  Mean: {meth_frac[mask].mean():.3f}")
        print(f"  Min:  {meth_frac[mask].min():.3f}")
        print(f"  Max:  {meth_frac[mask].max():.3f}")
    else:
        print("  WARNING: No CpGs with coverage > 0")
    print()

    # ----------------------------------------------------
    # 2. Load CpG positions
    # ----------------------------------------------------
    cpg_bed = pd.read_csv(
        cpg_bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end"]
    )

    if len(cpg_bed) != content.shape[0]:
        print("WARNING: CpG BED length does not match beta file length!")
        print("   This will cause incorrect genomic mapping.")
    else:
        print("✔ CpG BED length matches beta file")
    print()

    # ----------------------------------------------------
    # 3. Read chromosome sizes
    # ----------------------------------------------------
    chrom_sizes = {}
    with open(chrom_sizes_path) as f:
        for line in f:
            chrom, size = line.strip().split()
            chrom_sizes[chrom] = int(size)

    # -------------------------------
    # 4. Create bigWig files with new suffixes
    # -------------------------------
    meth_bw_path = os.path.join(outdir, f"{prefix}_WGBS_FractionalMethylation.bigwig")
    cov_bw_path = os.path.join(outdir, f"{prefix}_WGBS_ReadCoverage.bigwig")

    bw_meth = pyBigWig.open(meth_bw_path, "w")
    bw_cov = pyBigWig.open(cov_bw_path, "w")

    bw_meth.addHeader(list(chrom_sizes.items()))
    bw_cov.addHeader(list(chrom_sizes.items()))

    # ----------------------------------------------------
    # 5. Write data by chromosome
    # ----------------------------------------------------
    for chrom in cpg_bed["chrom"].unique():
        idx = cpg_bed["chrom"] == chrom
        cov_mask = total_cov[idx] > 0

        if cov_mask.any():
            positions = cpg_bed.loc[idx, "start"].values[cov_mask]

            bw_meth.addEntries(
                [chrom] * cov_mask.sum(),
                starts=positions.tolist(),
                ends=(positions + 1).tolist(),
                values=meth_frac[idx][cov_mask].tolist(),
            )

            bw_cov.addEntries(
                [chrom] * cov_mask.sum(),
                starts=positions.tolist(),
                ends=(positions + 1).tolist(),
                values=total_cov[idx][cov_mask].tolist(),
            )

    bw_meth.close()
    bw_cov.close()

main()