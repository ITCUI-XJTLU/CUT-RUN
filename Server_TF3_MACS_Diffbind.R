# --- Setup and Configuration ---
# Clear workspace (optional, uncomment if needed)
# rm(list = ls())
# gc()

set.seed(1)

# --- Paths (MODIFY THESE AS NEEDED) ---
# Define the root directory for your project where data, results, and scripts are located.
user_project_dir <- "/projects/gqilab/DAESC_GPU/data/simulation" # Base project directory
source_file_main <- file.path(user_project_dir, "scripts", "ChIPtest_source.R") # Should contain ALL helpers: simulate_tf, simulate_hm, run_MACS3_server, find_overlap, calculate_fdr_recall, AND run_macs3_diffbind_pipeline

# --- Load Libraries and Source Files ---
if (file.exists(source_file_main)) {
  source(source_file_main)
  cat("Sourced main helper and pipeline functions from:", source_file_main, "\n")
} else {
  stop("CRITICAL: Main source file (ChIPtest_source.R) not found at: ", source_file_main)
}

# Make sure the MACS3+DiffBind analysis function is defined
if (!exists("run_macs3_diffbind_pipeline") || !is.function(run_macs3_diffbind_pipeline)) {
  stop("CRITICAL: The analysis function 'run_macs3_diffbind_pipeline' is not defined. Please ensure it's in ChIPtest_source.R.")
}
if (!exists("simulate_tf") || !is.function(simulate_tf) || !exists("simulate_hm") || !is.function(simulate_hm) ) {
  stop("CRITICAL: Simulation functions 'simulate_tf' and/or 'simulate_hm' are not defined. Please ensure they are in ChIPtest_source.R.")
}
# Ensure other necessary functions (run_MACS3_server, find_overlap, calculate_fdr_recall) are also loaded from source_file_main.

# --- Global Simulation Settings ---
num_iterations <- 100 # Total number of simulation iterations to run (e.g., 100)
IS_TF_SIMULATION <- TRUE # SET TO TRUE for TF simulation, FALSE for Histone Mark simulation

# --- MACS3+DiffBind Pipeline Specific Global Parameters ---
PIPELINE_PEAK_WIDTH_TF <- 150
PIPELINE_PEAK_WIDTH_HM <- 1000
PIPELINE_GSIZE <- 300678703
PIPELINE_FRAGLEN <- 100

# --- Define Simulation Name Tag ---
# This tag will define the subfolder names under Data/ and Results/
# Example: "HM_MyModel_MACS3_DiffBind" or "TF_SpecificParams_MACS3_DiffBind"
simulation_model_tag <- if (IS_TF_SIMULATION) "TF3_MACS_DiffBind" else "TF3_MACS_DiffBind" # <<<< MODIFY THIS for each distinct simulation model/setting

# --- Define Base Directories for this Specific Simulation Model ---
# Data will be stored in /projects/gqilab/DAESC_GPU/data/simulation/Data/[simulation_model_tag]/iter_X/
# Results will be stored in /projects/gqilab/DAESC_GPU/data/simulation/Results/[simulation_model_tag]/iter_X_analysis_files/
# Aggregated results will be in /projects/gqilab/DAESC_GPU/data/simulation/Results/[simulation_model_tag]/

model_specific_simulated_data_base_dir <- file.path(user_project_dir, "Data", simulation_model_tag)
model_specific_analysis_results_base_dir <- file.path(user_project_dir, "Results", simulation_model_tag)

# --- Create Base Model-Specific Directories if they don't exist ---
if (!dir.exists(model_specific_simulated_data_base_dir)) {
  dir.create(model_specific_simulated_data_base_dir, recursive = TRUE)
  cat("Created base directory for simulated data for this model:", model_specific_simulated_data_base_dir, "\n")
}
if (!dir.exists(model_specific_analysis_results_base_dir)) {
  dir.create(model_specific_analysis_results_base_dir, recursive = TRUE)
  cat("Created base directory for analysis results for this model:", model_specific_analysis_results_base_dir, "\n")
}

# --- Simulation Parameters (customize based on IS_TF_SIMULATION) ---
if (IS_TF_SIMULATION) {
  sim_function_to_call <- simulate_tf
  sim_params <- list(
    fraglen = 100,
    base.mu = 120,
    up.mu = 180,
    down.mu.1 = 60,
    down.mu.2 = 60,
    bins_vis = 15,
    count_length_vis = 200
  )
  analysis_peak_width_for_pipeline <- PIPELINE_PEAK_WIDTH_TF
} else { # HM Simulation
  sim_function_to_call <- simulate_hm
  sim_params <- list(
    true.width = 2000,
    base.mu = 200,
    bins = 15,
    count_length = 200
  )
  analysis_peak_width_for_pipeline <- PIPELINE_PEAK_WIDTH_HM
}

# --- Results Aggregation Setup ---
all_iterations_performance_metrics <- vector("list", num_iterations)
simulation_run_log <- data.frame(
  iteration = integer(num_iterations),
  iteration_tag_pipeline = character(num_iterations),
  status = character(num_iterations),
  error_message = character(num_iterations),
  sim_data_iter_dir = character(num_iterations),         # Path to specific iteration's simulated data
  analysis_results_iter_dir = character(num_iterations), # Path to specific iteration's analysis results
  auc_value = numeric(num_iterations),
  start_time = as.POSIXct(rep(NA, num_iterations)),
  end_time = as.POSIXct(rep(NA, num_iterations)),
  duration_seconds = numeric(num_iterations),
  stringsAsFactors = FALSE
)
simulation_run_log$error_message <- NA_character_
simulation_run_log$status <- NA_character_
simulation_run_log$auc_value <- NA_real_


# --- Main Simulation Loop ---
cat(paste0("\n--- Starting simulation batch of ", num_iterations, " iterations for model: ", simulation_model_tag, " (", Sys.time(), ") ---\n"))

for (iter_num in 1:num_iterations) {
  cat(paste0("\n--- Starting Iteration: ", iter_num, "/", num_iterations, " (", Sys.time(), ") ---\n"))
  iter_start_time <- Sys.time()
  current_status <- "Pending"
  current_error_msg <- NA_character_
  current_auc <- NA_real_
  pipeline_iteration_tag <- paste0("iter", iter_num) # Tag for pipeline function, used to name subfolder within iter_analysis_parent_dir
  
  # --- Define and Create Iteration-Specific Directories ---
  # Directory for storing simulated data for this iteration (e.g. .../Data/HM_Model1.../iter_1/)
  iter_sim_data_dir <- file.path(model_specific_simulated_data_base_dir, paste0("iter_", iter_num))
  # Parent directory for analysis results of this iteration (e.g. .../Results/HM_Model1.../iter_1_analysis_output_base/)
  # The run_macs3_diffbind_pipeline will create its own "analysis_run_iterX" subfolder inside this.
  iter_analysis_output_base_dir <- file.path(model_specific_analysis_results_base_dir, paste0("iter_", iter_num, "_analysis_output_base"))
  
  tryCatch({
    if (!dir.exists(iter_sim_data_dir)) {
      dir.create(iter_sim_data_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(iter_analysis_output_base_dir)) {
      dir.create(iter_analysis_output_base_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cat("Simulated data for iter", iter_num, "will be in:", iter_sim_data_dir, "\n")
    cat("MACS3+DiffBind analysis output base for iter", iter_num, "will be in:", iter_analysis_output_base_dir, "\n")
  }, error = function(e_dir) {
    stop(paste("CRITICAL: Failed to ensure directories for iteration", iter_num, ". Error:", conditionMessage(e_dir)))
  })
  
  analysis_results_obj <- NULL
  
  tryCatch({
    # Step 1: Simulate data
    cat("Step 1: Running data simulation for iteration", iter_num, "...\n")
    current_sim_args <- c(list(data_path = iter_sim_data_dir), sim_params)
    sim_output <- do.call(sim_function_to_call, current_sim_args)
    cat("Data simulation completed for iteration", iter_num, ". Output in:", iter_sim_data_dir, "\n")
    
    # Step 2: Run MACS3 + DiffBind Pipeline Analysis
    cat("Step 2: Running MACS3+DiffBind Pipeline for iteration", iter_num, "...\n")
    
    analysis_results_obj <- run_macs3_diffbind_pipeline(
      data_directory = iter_sim_data_dir,
      output_directory_base = iter_analysis_output_base_dir, # The pipeline will create "analysis_run_iterX" inside this
      iteration_tag = pipeline_iteration_tag, # This tag is for the pipeline's internal naming
      is_tf = IS_TF_SIMULATION,
      peak_width = analysis_peak_width_for_pipeline,
      count_bin=10,
      source_file_path = source_file_main, # Path to ChIPtest_source.R (for find_overlap etc.)
      gsize = PIPELINE_GSIZE,
      fraglen = PIPELINE_FRAGLEN
    )
    cat("MACS3+DiffBind Pipeline analysis completed for iteration", iter_num, ".\n")
    
    if (!is.null(analysis_results_obj) && !is.null(analysis_results_obj$metrics_df)) {
      current_metrics_df <- analysis_results_obj$metrics_df
      current_metrics_df$main_simulation_iter_num <- iter_num
      all_iterations_performance_metrics[[iter_num]] <- current_metrics_df
      current_status <- "Success"
      if ("auc_value" %in% colnames(current_metrics_df)) {
        current_auc <- current_metrics_df$auc_value[1]
      } else if (!is.null(analysis_results_obj$auc_value)) {
        current_auc <- analysis_results_obj$auc_value
      }
    } else {
      warning("No performance metrics dataframe returned from pipeline for iteration ", iter_num)
      all_iterations_performance_metrics[[iter_num]] <- NULL
      current_status <- "CompletedWithWarning_NoMetrics"
      current_auc <- NA_real_
    }
    
  }, error = function(e_main) {
    cat("ERROR in simulation/analysis for iteration", iter_num, ": ", conditionMessage(e_main), "\n")
    current_status <<- "Failed"
    current_error_msg <<- conditionMessage(e_main)
    all_iterations_performance_metrics[[iter_num]] <- NULL
    current_auc <<- NA_real_
  })
  
  iter_end_time <- Sys.time()
  iter_duration_secs <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
  
  simulation_run_log[iter_num, "iteration"] <- iter_num
  simulation_run_log[iter_num, "iteration_tag_pipeline"] <- pipeline_iteration_tag
  simulation_run_log[iter_num, "status"] <- current_status
  simulation_run_log[iter_num, "error_message"] <- current_error_msg
  simulation_run_log[iter_num, "sim_data_iter_dir"] <- iter_sim_data_dir
  # Store the path where the specific iteration's analysis files were saved by run_macs3_diffbind_pipeline
  simulation_run_log[iter_num, "analysis_results_iter_dir"] <- if(!is.null(analysis_results_obj) && !is.null(analysis_results_obj$output_path)) analysis_results_obj$output_path else NA_character_
  simulation_run_log[iter_num, "auc_value"] <- current_auc
  simulation_run_log[iter_num, "start_time"] <- iter_start_time
  simulation_run_log[iter_num, "end_time"] <- iter_end_time
  simulation_run_log[iter_num, "duration_seconds"] <- iter_duration_secs
  
  cat(paste0("--- Iteration ", iter_num, " finished. Status: ", current_status,
             ". Duration: ", round(iter_duration_secs, 2), " seconds. AUC: ", round(current_auc, 4) ," ---\n"))
}

# --- Post-Simulation Processing ---
cat("\n--- Aggregating results from all iterations for model: ", simulation_model_tag, " ---\n")

valid_metrics_list <- all_iterations_performance_metrics[!sapply(all_iterations_performance_metrics, is.null)]

if (length(valid_metrics_list) > 0) {
  final_performance_summary_df <- do.call(rbind, valid_metrics_list)
  # Save aggregated metrics directly into the model_specific_analysis_results_base_dir
  summary_filename <- file.path(model_specific_analysis_results_base_dir, paste0("Aggregated_Metrics_", simulation_model_tag, ".csv"))
  tryCatch({
    write.csv(final_performance_summary_df, summary_filename, row.names = FALSE)
    cat("Aggregated performance metrics saved to:", summary_filename, "\n")
  }, error = function(e_save_metrics) {
    cat("Error saving aggregated performance metrics:", conditionMessage(e_save_metrics), "\n")
  })
  
  # Optional: Boxplot of AUCs from the log, saved into model_specific_analysis_results_base_dir
  if("auc_value" %in% names(simulation_run_log) && sum(!is.na(simulation_run_log$auc_value)) > 0) {
    # auc_plot_filename <- file.path(model_specific_analysis_results_base_dir, paste0("AUC_Boxplot_", simulation_model_tag, ".png"))
    # png(auc_plot_filename, width=800, height=600)
    #   boxplot(simulation_run_log$auc_value[!is.na(simulation_run_log$auc_value)],
    #           main=paste("AUC Values for", simulation_model_tag), ylab="AUC",
    #           ylim=c(min(0.4, min(simulation_run_log$auc_value, na.rm=TRUE)), 1))
    # dev.off()
    # cat("AUC boxplot saved to:", auc_plot_filename, "\n")
  }
  
} else {
  cat("No performance metrics were successfully collected to aggregate for model: ", simulation_model_tag, ".\n")
}

# Save the detailed simulation run log into model_specific_analysis_results_base_dir
log_filename <- file.path(model_specific_analysis_results_base_dir, paste0("Run_Log_", simulation_model_tag, ".csv"))
tryCatch({
  write.csv(simulation_run_log, log_filename, row.names = FALSE)
  cat("Simulation run log saved to:", log_filename, "\n")
}, error = function(e_save_log) {
  cat("Error saving simulation log:", conditionMessage(e_save_log), "\n")
})

cat(paste0("\n--- Simulation batch for model: ", simulation_model_tag, " Finished (", Sys.time(), ") ---\n"))