#!/bin/bash
#SBATCH --job-name=beta_to_bigwig_array
#SBATCH --output=logs/beta_to_bigwig_%A_%a.out
#SBATCH --error=logs/beta_to_bigwig_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=dept_cpu
#SBATCH --array=1-100%10

# Job description:
# Convert beta files to bigWig format using SLURM array jobs
# Each array task processes one beta file
# The %10 limits to 10 concurrent jobs - adjust as needed

# ============================================================================
# CONFIGURATION - MODIFY THESE VARIABLES FOR YOUR SETUP
# ============================================================================

# Path to the original conversion script
ORIGINAL_SCRIPT="/net/dali/home/mscbio/yur28/directed_study/extract_files.py"

# Directory containing beta files
BETA_DIR="/net/dali/home/mscbio/yur28/directed_study/betadir"

# Reference directory (contains CpG.bed.gz and CpG.chrome.size)
REF_DIR="/net/dali/home/mscbio/yur28/directed_study/wgbs_tools/references/hg38"

# Output directory
OUTPUT_DIR="/net/dali/home/mscbio/yur28/directed_study/roadmap"

# File pattern to match (e.g., "*.beta" or "sample_*.beta")
PATTERN="*.beta"

# ============================================================================
# SETUP
# ============================================================================

# Create logs directory if it doesn't exist
mkdir -p logs

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Get list of beta files (only on first task to avoid race conditions)
FILE_LIST="${OUTPUT_DIR}/.beta_file_list_${SLURM_ARRAY_JOB_ID}.txt"

if [ ! -f "${FILE_LIST}" ]; then
    # Create file list
    find ${BETA_DIR} -name "${PATTERN}" -type f | sort > ${FILE_LIST}
fi

# Wait for file list to be created
while [ ! -f "${FILE_LIST}" ]; do
    sleep 1
done

# Get the specific file for this array task
BETA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${FILE_LIST})

# Check if we got a valid file
if [ -z "${BETA_FILE}" ]; then
    echo "ERROR: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Print job information
echo "===================================================================="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "Started: $(date)"
echo "===================================================================="
echo ""
echo "Processing file: ${BETA_FILE}"
echo "Reference directory: ${REF_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "===================================================================="
echo ""

# ============================================================================
# LOAD REQUIRED MODULES (adjust for your cluster)
# ============================================================================

# Example for common module systems:
# module load python/3.9
# module load numpy/1.21
# module load pandas/1.3

# OR activate a conda environment:
# source /path/to/miniconda3/etc/profile.d/conda.sh
# conda activate wgbs_tools

# OR activate a virtual environment:
# source /path/to/venv/bin/activate

# ============================================================================
# RUN THE ORIGINAL SCRIPT
# ============================================================================

python3 ${ORIGINAL_SCRIPT} \
    --beta ${BETA_FILE} \
    --indir ${REF_DIR} \
    --outdir ${OUTPUT_DIR}

# Capture exit status
EXIT_STATUS=$?

# ============================================================================
# CLEANUP AND SUMMARY
# ============================================================================

echo ""
echo "===================================================================="
echo "Job finished: $(date)"
echo "Exit status: ${EXIT_STATUS}"
echo "===================================================================="

# If this is the last task, clean up the file list
TOTAL_TASKS=$(wc -l < ${FILE_LIST})
if [ ${SLURM_ARRAY_TASK_ID} -eq ${TOTAL_TASKS} ]; then
    # Give other jobs time to read the file
    sleep 10
    rm -f ${FILE_LIST}
fi

# Exit with the status from the Python script
exit ${EXIT_STATUS}
