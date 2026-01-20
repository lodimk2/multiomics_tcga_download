#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(DT)
  library(SummarizedExperiment)
  library(openxlsx)
  library(optparse)
})

# Main function to download and process TCGA data
download_tcga_data <- function(project_code, base_dir = getwd(), log_file = NULL) {
  # Create logs directory if it doesn't exist
  logs_dir <- file.path(base_dir, "logs")
  if (!dir.exists(logs_dir)) {
    dir.create(logs_dir, recursive = TRUE)
  }
  
  # If log file is NULL, set default in logs directory
  if (is.null(log_file)) {
    log_file <- file.path(
      logs_dir,
      paste0(project_code, "_download_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
    )
  } else if (!grepl("^/", log_file)) {
    # If log_file is not an absolute path, put it in logs directory
    log_file <- file.path(logs_dir, log_file)
  }
  
  log_con <- file(log_file, "a")
  
  # Redefine cat function to write to log file
  log_message <- function(...) {
    msg <- paste0(...)
    # Write to console
    cat(msg, "\n")
    # Write to log file
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), msg, "\n",
        file = log_con, append = TRUE)
  }
  
  log_message(sprintf("Log file created at: %s", log_file))
  
  # Register cleanup function to close log file when done
  on.exit({
    if (!is.null(log_file) && exists("log_con") && isOpen(log_con)) {
      close(log_con)
      # NOTE: can't call log_message after close(log_con) safely for file logging,
      # but console output is fine.
      cat("Log file closed\n")
    }
  }, add = TRUE)
  
  # Create project directory if it doesn't exist
  project_dir <- file.path(base_dir, paste0(project_code, "_Data"))
  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE)
    log_message(sprintf("Created directory: %s", project_dir))
  }
  
  # Set working directory to project directory
  old_dir <- getwd()
  on.exit({
    setwd(old_dir)  # Reset working directory when function exits
    log_message(sprintf("Reset working directory to: %s", old_dir))
  }, add = TRUE)
  
  setwd(project_dir)
  log_message(sprintf("Working directory set to: %s", project_dir))
  
  # Log start time
  start_time <- Sys.time()
  log_message(sprintf("Processing %s data started at %s", project_code, start_time))
  
  # Initialize data containers
  rna_tumor_data_mat  <- NULL
  rna_normal_data_mat <- NULL
  meth_combined_data  <- NULL
  mirna_tumor_data_mat <- NULL
  clinical_data_df    <- NULL
  
  # Helper function to extract patient IDs
  extract_patient_id <- function(sample_names) {
    valid_samples <- grep("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_names, value = TRUE)
    unique(substr(valid_samples, 1, 12))
  }
  
  # Try-catch wrapper for GDC operations
  safe_gdc_operation <- function(expr, operation_name) {
    result <- NULL
    start_op_time <- Sys.time()
    log_message(sprintf("Starting %s...", operation_name))
    
    tryCatch({
      result <- eval(expr)
      end_op_time <- Sys.time()
      duration <- difftime(end_op_time, start_op_time, units = "mins")
      log_message(sprintf("Completed %s in %.2f minutes", operation_name, as.numeric(duration)))
      return(result)
    }, error = function(e) {
      log_message(sprintf("ERROR in %s: %s", operation_name, e$message))
      return(NULL)
    }, warning = function(w) {
      log_message(sprintf("WARNING in %s: %s", operation_name, w$message))
      # Return whatever was computed (if any) despite warning
      return(invisible(NULL))
    })
  }
  
  # 1. Download RNA-Seq Data - Normal Samples
  log_message("=== RNA-Seq Normal Samples ===")
  normal_query <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      sample.type = c(
        "Blood Derived Normal",
        "Solid Tissue Normal",
        "Buccal Cell Normal",
        "EBV Immortalized Normal",
        "Bone Marrow Normal"
      )
    )),
    operation_name = "RNA normal query"
  )
  
  if (!is.null(normal_query)) {
    safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCdownload(query = normal_query)),
      operation_name = "RNA normal download"
    )
    
    normal_data <- safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCprepare(query = normal_query)),
      operation_name = "RNA normal prepare"
    )
    
    if (!is.null(normal_data)) {
      rna_normal_data_mat <- as.data.frame(SummarizedExperiment::assay(normal_data))
      normal_file <- paste0("RNA_", project_code, "_Normal.csv")
      write.csv(rna_normal_data_mat, file = normal_file)
      log_message(sprintf("%s NORMAL RNA PATIENTS: %d (saved to %s)",
                          project_code, ncol(rna_normal_data_mat), normal_file))
    }
  } else {
    log_message("No normal RNA-Seq data found for this project")
  }
  
  # 2. Download RNA-Seq Data - Tumor Samples
  log_message("=== RNA-Seq Tumor Samples ===")
  tumor_query <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c(
        "Primary Tumor", "Recurrent Tumor", "Metastatic",
        "Additional - New Primary", "Primary Blood Derived Cancer - Peripheral Blood"
      )
    )),
    operation_name = "RNA tumor query"
  )
  
  if (!is.null(tumor_query)) {
    safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCdownload(query = tumor_query)),
      operation_name = "RNA tumor download"
    )
    
    tumor_data <- safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCprepare(query = tumor_query)),
      operation_name = "RNA tumor prepare"
    )
    
    if (!is.null(tumor_data)) {
      rna_tumor_data_mat <- as.data.frame(SummarizedExperiment::assay(tumor_data))
      tumor_file <- paste0("RNA_", project_code, "_Tumor.csv")
      write.csv(rna_tumor_data_mat, file = tumor_file)
      log_message(sprintf("%s TUMOR RNA PATIENTS: %d (saved to %s)",
                          project_code, ncol(rna_tumor_data_mat), tumor_file))
    }
  } else {
    log_message("No tumor RNA-Seq data found for this project")
  }
  
  # 3. Download DNA Methylation Data - 27K Platform
  log_message("=== DNA Methylation 27K Platform ===")
  meth_27k_query <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "DNA Methylation",
      data.type = "Methylation Beta Value",
      platform = "Illumina Human Methylation 27",
      sample.type = c(
        "Primary Tumor", "Recurrent Tumor", "Metastatic",
        "Additional - New Primary", "Primary Blood Derived Cancer - Peripheral Blood"
      )
    )),
    operation_name = "Methylation 27K query"
  )
  
  meth_27k_data_mat <- NULL
  if (!is.null(meth_27k_query)) {
    safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCdownload(query = meth_27k_query, files.per.chunk = 10, method = "api")),
      operation_name = "Methylation 27K download"
    )
    
    meth_27k_data <- safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCprepare(query = meth_27k_query)),
      operation_name = "Methylation 27K prepare"
    )
    
    if (!is.null(meth_27k_data)) {
      meth_27k_data_mat <- as.data.frame(SummarizedExperiment::assay(meth_27k_data))
      log_message(sprintf("%s TUMOR METHYLATION 27K PATIENTS: %d",
                          project_code, ncol(meth_27k_data_mat)))
    }
  } else {
    log_message("No methylation 27K data found for this project")
  }
  
  # 4. Download DNA Methylation Data - 450K Platform
  log_message("=== DNA Methylation 450K Platform ===")
  meth_450k_query <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "DNA Methylation",
      data.type = "Methylation Beta Value",
      platform = "Illumina Human Methylation 450",
      sample.type = c(
        "Primary Tumor", "Recurrent Tumor", "Metastatic",
        "Additional - New Primary", "Primary Blood Derived Cancer - Peripheral Blood"
      )
    )),
    operation_name = "Methylation 450K query"
  )
  
  meth_450k_data_mat <- NULL
  if (!is.null(meth_450k_query)) {
    safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCdownload(query = meth_450k_query, files.per.chunk = 10, method = "api")),
      operation_name = "Methylation 450K download"
    )
    
    meth_450k_data <- safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCprepare(query = meth_450k_query)),
      operation_name = "Methylation 450K prepare"
    )
    
    if (!is.null(meth_450k_data)) {
      meth_450k_data_mat <- as.data.frame(SummarizedExperiment::assay(meth_450k_data))
      log_message(sprintf("%s TUMOR METHYLATION 450K PATIENTS: %d",
                          project_code, ncol(meth_450k_data_mat)))
    }
  } else {
    log_message("No methylation 450K data found for this project")
  }
  
  # Combine methylation data if available
  log_message("=== Combining Methylation Data ===")
  if (!is.null(meth_27k_data_mat) && !is.null(meth_450k_data_mat)) {
    meth_combined_data <- merge(meth_27k_data_mat, meth_450k_data_mat, by = "row.names", all = TRUE)
    rownames(meth_combined_data) <- meth_combined_data$Row.names
    meth_combined_data$Row.names <- NULL
    log_message("Combined both 27K and 450K methylation data")
  } else if (!is.null(meth_27k_data_mat)) {
    meth_combined_data <- meth_27k_data_mat
    log_message("Using only 27K methylation data")
  } else if (!is.null(meth_450k_data_mat)) {
    meth_combined_data <- meth_450k_data_mat
    log_message("Using only 450K methylation data")
  } else {
    log_message("No methylation data available for this project")
  }
  
  if (!is.null(meth_combined_data)) {
    meth_file <- paste0("Methylation_", project_code, "_Tumor.csv")
    write.csv(meth_combined_data, file = meth_file)
    log_message(sprintf("%s TOTAL METHYLATION PATIENTS: %d (saved to %s)",
                        project_code, ncol(meth_combined_data), meth_file))
  }
  
  # 5. Download miRNA Data
  log_message("=== miRNA Expression Data ===")
  mirna_query <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery(
      project = project_code,
      data.category = "Transcriptome Profiling",
      data.type = "miRNA Expression Quantification",
      sample.type = c(
        "Primary Tumor", "Recurrent Tumor", "Metastatic",
        "Additional - New Primary", "Primary Blood Derived Cancer - Peripheral Blood"
      )
    )),
    operation_name = "miRNA query"
  )
  
  if (!is.null(mirna_query)) {
    safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCdownload(query = mirna_query, files.per.chunk = 10, method = "api")),
      operation_name = "miRNA download"
    )
    
    mirna_data <- safe_gdc_operation(
      expr = quote(TCGAbiolinks::GDCprepare(query = mirna_query)),
      operation_name = "miRNA prepare"
    )
    
    if (!is.null(mirna_data)) {
      mirna_data_mat <- as.data.frame(mirna_data)
      
      # TCGAbiolinks miRNA format often includes columns like read_count_*.
      mirna_tumor_data_mat <- mirna_data_mat %>%
        dplyr::select(dplyr::starts_with("read_count")) %>%
        dplyr::rename_with(~ sub("read_count_", "", .))
      
      if ("miRNA_ID" %in% colnames(mirna_data_mat)) {
        log_message("Using miRNA data format with miRNA_ID column detected")
      } else {
        log_message("Using standard miRNA data format")
      }
      
      mirna_file <- paste0("miRNA_", project_code, "_Tumor.csv")
      write.csv(mirna_tumor_data_mat, file = mirna_file)
      log_message(sprintf("%s TUMOR miRNA PATIENTS: %d (saved to %s)",
                          project_code, ncol(mirna_tumor_data_mat), mirna_file))
    }
  } else {
    log_message("No miRNA data found for this project")
  }
  
  # 6. Download Clinical Data (ADDED BACK)
  log_message("=== Clinical Data ===")
  clinical_data <- safe_gdc_operation(
    expr = quote(TCGAbiolinks::GDCquery_clinic(project = project_code, type = "clinical")),
    operation_name = "Clinical query"
  )
  
  if (!is.null(clinical_data)) {
    clinical_data_df <- as.data.frame(clinical_data)
    clinical_file <- paste0("Clinical_", project_code, ".csv")
    write.csv(clinical_data_df, file = clinical_file, row.names = FALSE)
    log_message(sprintf("%s CLINICAL ROWS: %d (saved to %s)",
                        project_code, nrow(clinical_data_df), clinical_file))
  } else {
    log_message("No clinical data returned for this project")
  }
  
  # Generate patient overlap statistics
  log_message("=== Generating Patient Overlap Statistics ===")
  
  # Extract patient IDs for available data types
  rna_patients  <- if (!is.null(rna_tumor_data_mat))  extract_patient_id(colnames(rna_tumor_data_mat)) else character(0)
  meth_patients <- if (!is.null(meth_combined_data))  extract_patient_id(colnames(meth_combined_data)) else character(0)
  mirna_patients <- if (!is.null(mirna_tumor_data_mat)) extract_patient_id(colnames(mirna_tumor_data_mat)) else character(0)
  
  log_message(sprintf("Number of RNA patients: %d", length(rna_patients)))
  log_message(sprintf("Number of Methylation patients: %d", length(meth_patients)))
  log_message(sprintf("Number of miRNA patients: %d", length(mirna_patients)))
  
  # List available data types
  available_types <- c(
    ifelse(length(rna_patients) > 0, "RNA", NA),
    ifelse(length(meth_patients) > 0, "Methylation", NA),
    ifelse(length(mirna_patients) > 0, "miRNA", NA)
  )
  available_types <- available_types[!is.na(available_types)]
  log_message(sprintf("Available data types: %s", paste(available_types, collapse = ", ")))
  
  # Get union of all patients
  all_patients <- unique(c(rna_patients, meth_patients, mirna_patients))
  log_message(sprintf("Total unique patients across all modalities: %d", length(all_patients)))
  
  # Log simple overlap statistics
  log_message("=== Patient Overlap Statistics ===")
  
  # RNA-specific patients
  rna_only <- character(0)
  meth_only <- character(0)
  mirna_only <- character(0)
  rna_meth <- 0
  rna_mirna <- 0
  meth_mirna <- 0
  three_way <- 0
  
  if (length(rna_patients) > 0) {
    rna_only <- setdiff(rna_patients, union(meth_patients, mirna_patients))
    log_message(sprintf("RNA Only: %d", length(rna_only)))
  }
  
  if (length(meth_patients) > 0) {
    meth_only <- setdiff(meth_patients, union(rna_patients, mirna_patients))
    log_message(sprintf("Methylation Only: %d", length(meth_only)))
  }
  
  if (length(mirna_patients) > 0) {
    mirna_only <- setdiff(mirna_patients, union(rna_patients, meth_patients))
    log_message(sprintf("miRNA Only: %d", length(mirna_only)))
  }
  
  # Pairwise overlaps (excluding three-way)
  if (length(rna_patients) > 0 && length(meth_patients) > 0) {
    rna_meth <- length(intersect(rna_patients, meth_patients))
  }
  if (length(rna_patients) > 0 && length(mirna_patients) > 0) {
    rna_mirna <- length(intersect(rna_patients, mirna_patients))
  }
  if (length(meth_patients) > 0 && length(mirna_patients) > 0) {
    meth_mirna <- length(intersect(meth_patients, mirna_patients))
  }
  
  if (length(rna_patients) > 0 && length(meth_patients) > 0 && length(mirna_patients) > 0) {
    three_way <- length(intersect(intersect(rna_patients, meth_patients), mirna_patients))
    log_message(sprintf("RNA & Methylation & miRNA: %d", three_way))
  }
  
  # Subtract the three-way overlap from pairwise to keep categories disjoint
  if (three_way > 0) {
    if (rna_meth > 0)   rna_meth   <- rna_meth - three_way
    if (rna_mirna > 0)  rna_mirna  <- rna_mirna - three_way
    if (meth_mirna > 0) meth_mirna <- meth_mirna - three_way
  }
  
  if (length(rna_patients) > 0 && length(meth_patients) > 0) {
    log_message(sprintf("RNA & Methylation: %d", rna_meth))
  }
  if (length(rna_patients) > 0 && length(mirna_patients) > 0) {
    log_message(sprintf("RNA & miRNA: %d", rna_mirna))
  }
  if (length(meth_patients) > 0 && length(mirna_patients) > 0) {
    log_message(sprintf("Methylation & miRNA: %d", meth_mirna))
  }
  
  # Create overlap summary outputs
  overlap_data <- data.frame(
    Category = c("RNA Only", "Methylation Only", "miRNA Only",
                 "RNA & Methylation", "RNA & miRNA", "Methylation & miRNA",
                 "All Three"),
    Count = c(length(rna_only), length(meth_only), length(mirna_only),
              rna_meth, rna_mirna, meth_mirna, three_way),
    stringsAsFactors = FALSE
  )
  
  overlap_file <- paste0("Patient_Overlap_", project_code, ".csv")
  write.csv(overlap_data, file = overlap_file, row.names = FALSE)
  log_message(sprintf("Patient overlap statistics written to: %s", overlap_file))
  
  # Create a simpler all-modality overlap file that lists which modalities each patient has
  all_modality_overlap <- data.frame(
    Patient_ID = all_patients,
    RNA_Seq = all_patients %in% rna_patients,
    Methylation = all_patients %in% meth_patients,
    miRNA = all_patients %in% mirna_patients
  )
  all_modality_overlap$Modality_Count <- rowSums(all_modality_overlap[, c("RNA_Seq", "Methylation", "miRNA")])
  all_modality_overlap <- all_modality_overlap[order(all_modality_overlap$Patient_ID), ]
  
  all_modality_file <- paste0("All_Modality_Overlap_", project_code, ".csv")
  write.csv(all_modality_overlap, file = all_modality_file, row.names = FALSE)
  log_message(sprintf("All modality overlap written to: %s", all_modality_file))
  
  # Modality count summary
  modality_counts <- table(all_modality_overlap$Modality_Count)
  log_message("--- Modality Count Summary ---")
  for (i in 1:3) {
    count <- ifelse(as.character(i) %in% names(modality_counts), as.integer(modality_counts[as.character(i)]), 0)
    pct <- if (length(all_patients) > 0) (count / length(all_patients) * 100) else 0
    log_message(sprintf("Patients with %d modalities: %d (%.1f%%)", i, count, pct))
  }
  
  # Log end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(sprintf("Processing completed at %s", end_time))
  log_message(sprintf("Total processing time: %.2f minutes", as.numeric(duration)))
  
  # Return a structured object (fixes the undefined patient_overlap_stats in your version)
  return(list(
    project_code = project_code,
    files = list(
      rna_tumor = if (!is.null(rna_tumor_data_mat)) paste0("RNA_", project_code, "_Tumor.csv") else NA_character_,
      rna_normal = if (!is.null(rna_normal_data_mat)) paste0("RNA_", project_code, "_Normal.csv") else NA_character_,
      methylation_tumor = if (!is.null(meth_combined_data)) paste0("Methylation_", project_code, "_Tumor.csv") else NA_character_,
      mirna_tumor = if (!is.null(mirna_tumor_data_mat)) paste0("miRNA_", project_code, "_Tumor.csv") else NA_character_,
      clinical = if (!is.null(clinical_data_df)) paste0("Clinical_", project_code, ".csv") else NA_character_,
      overlap_summary = overlap_file,
      overlap_by_patient = all_modality_file,
      log_file = log_file
    ),
    patient_ids = list(
      rna = rna_patients,
      methylation = meth_patients,
      mirna = mirna_patients,
      all = all_patients
    ),
    overlap_summary = overlap_data,
    clinical = clinical_data_df
  ))
  
  # --- Venn diagram of patient overlap (RNA / Methylation / miRNA) ---
  venn_png <- paste0("Venn_Patient_Overlap_", project_code, ".png")
  
  # Only attempt if at least 2 modalities exist (a 1-set Venn isn't very informative)
  n_nonempty <- sum(c(length(rna_patients) > 0, length(meth_patients) > 0, length(mirna_patients) > 0))
  
  if (n_nonempty >= 2) {
    # Use VennDiagram package
    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
      log_message("VennDiagram package not installed; skipping Venn PNG export.")
    } else {
      log_message(sprintf("Exporting Venn diagram to: %s", venn_png))
      
      # Prepare sets (named list)
      venn_sets <- list(
        RNA = unique(rna_patients),
        Methylation = unique(meth_patients),
        miRNA = unique(mirna_patients)
      )
      
      # Optionally drop empty sets to avoid errors/odd plots
      venn_sets <- venn_sets[sapply(venn_sets, length) > 0]
      
      # Create PNG device
      grDevices::png(filename = venn_png, width = 1600, height = 1600, res = 200)
      
      # Draw Venn
      grid::grid.newpage()
      venn_grob <- VennDiagram::venn.diagram(
        x = venn_sets,
        filename = NULL,          # return grob instead of writing directly
        main = paste0(project_code, " Patient Overlap by Modality"),
        main.cex = 1.4,
        category.cex = 1.2,
        cex = 1.3,
        cat.pos = 0,
        margin = 0.1
      )
      grid::grid.draw(venn_grob)
      
      # Close device
      grDevices::dev.off()
      
      log_message("Venn diagram export completed.")
    }
  } else {
    log_message("Fewer than two modalities have patient IDs; skipping Venn PNG export.")
  }
  
}


# Parse command line arguments
option_list <- list(
  make_option(c("-p", "--project"), type="character", default=NULL, 
              help="TCGA project code (e.g., TCGA-BRCA, TCGA-GBM)", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=getwd(), 
              help="Base directory to save data [default= current directory]", metavar="character"),
  make_option(c("-l", "--log"), type="character", default=NULL, 
              help="Log file path (optional)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$project)) {
  cat("Error: Project code is required.\n")
  cat("Usage: Rscript tcga_download.R -p TCGA-GBM [-d /path/to/save] [-l logfile.txt]\n")
  quit(status = 1)
}

# Set default log file if not provided
if (is.null(opt$log)) {
  opt$log <- file.path(opt$dir, paste0(opt$project, "_download_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
}

# Start the log file
cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Starting TCGA data download script\n", file = opt$log)
cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Project: ", opt$project, "\n", file = opt$log, append = TRUE)
cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Base directory: ", opt$dir, "\n", file = opt$log, append = TRUE)
cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Log file: ", opt$log, "\n", file = opt$log, append = TRUE)

# Execute the download
tryCatch({
  result <- download_tcga_data(opt$project, opt$dir, opt$log)
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "Script completed successfully\n", file = opt$log, append = TRUE)
  quit(status = 0)
}, error = function(e) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), "FATAL ERROR: ", e$message, "\n", file = opt$log, append = TRUE)
  cat("FATAL ERROR: ", e$message, "\n")
  quit(status = 1)
})