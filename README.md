# TCGA Biolinks Download Code for Multiomics
This pipeline downloads miRNA, mRNA, and DNA Methylation data for TCGA Project Codes. It will also report summary statistics in the log file to determine how many patients have multiple modalities of data available. 
## Usage Instructions and Requirements 

### Required R packages (R version 4.3.1)

Install:
- TCGAbiolinks: '2.30.4'
- SummarizedExperiment: '1.32.0'

### Run Instructions 

To run the download, execute ```TCGA_download.R```. For HPRC environments, please use ```tcga_run.sh```. 

The two inputs required for ```TCGA_download.R``` are the TCGA project code and an output directory. An example command for this is below:

```Rscript TCGA_download.R -p $PROJECT -d "$BASE_DIR" | tee -a "$LOG_FILE"```. 

To run the SLURM script simply change the command arguments for the $PROJECT, $BASE_DIR and $LOG_FILE to the desired inputs. 
