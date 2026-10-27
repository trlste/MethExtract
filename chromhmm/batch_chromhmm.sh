#!/bin/bash
#SBATCH --job-name=chrhm
#SBATCH -p dept_cpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output chrhm_enrichment_1026.out

set -euo pipefail

PROJECT="/net/dali/home/mscbio/jiw417/project/methyl_seg"
SCRATCH="/scr/${USER}/chrhm_${SLURM_JOB_ID}"

echo "Job started on $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Scratch directory: ${SCRATCH}"
date

# Load conda

eval "$(conda shell.bash hook)"
# Activate conda environment
conda activate methyl

mkdir -p "${SCRATCH}"

cleanup() {
    echo "Cleaning scratch directory: ${SCRATCH}"
    rm -rf "${SCRATCH}"
}
trap cleanup EXIT

# Copy and paste scripts and E038 annotation file

cp "${PROJECT}/chromhmm/chromhmm/E038_15_coreMarks_hg38lift_mnemonics.bed.gz" "${SCRATCH}"
cp "${PROJECT}/chromhmm/hg38_chr2.genome" "${SCRATCH}"
cp "${PROJECT}/chromhmm/batch_script_chromhmm"/* "${SCRATCH}"
cd "${SCRATCH}"

# Run each loop by lambda

for lambda in 5 1 '5e-1' '5e-2'
do
	cp "${PROJECT}/cofreq/coseg_${lambda}"/output_*.bedgraph "${SCRATCH}"

	Rscript s1_divide_segments.R
	sh s2_make_100_random.sh
	sh s3_compute_overlap_observe.sh
	sh s4_compute_overlap_random.sh
	Rscript s5_compute_enrichment.R "${lambda}"

	rm -r random
	rm -r overlap
	rm -r output_*.bedgraph

	echo "L = ${lambda} done."
	date
done

echo "All done. Copying files..."
date
cp "${SCRATCH}"/chromhmm_enrichment_L*.pdf "${PROJECT}/chromhmm/chromhmm_enrichment_pdf"

cd "${PROJECT}"

echo "Job finished."

