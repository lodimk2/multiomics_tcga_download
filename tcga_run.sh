#!/bin/bash
#SBATCH --job-name=tcga_kirp
#SBATCH --output=tcga-%j.out
#SBATCH --error=tcga-%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lodimk2@vcu.edu

# Set project code - change this for different cancer types
PROJECT="TCGA-KIRP"

# Set base directory for output
BASE_DIR="/lustre/home/lodimk2/0325_biolinks/01202025_test"

# Load required modules
module load R/4.3.1

# Create logs directory if it doesn't exist
mkdir -p "$BASE_DIR/logs"

# Create a timestamp for job tracking
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$BASE_DIR/logs/slurm_${PROJECT}_${TIMESTAMP}.log"

echo "Job started at: $(date)" | tee -a "$LOG_FILE"
echo "Project: $PROJECT" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

# Execute R script (logs will be created inside logs directory by the R script itself)
Rscript TCGA_download.R -p $PROJECT -d "$BASE_DIR" | tee -a "$LOG_FILE"

# Job completion message
echo "Job completed at: $(date)"
exit 0
