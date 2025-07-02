# --- Setup and Configuration ---
# Clear workspace (optional, uncomment if needed)
# rm(list = ls())
# gc()

set.seed(1)

# --- Paths (MODIFY THESE AS NEEDED) ---
# Define the root directory for your project where data, results, and scripts are located.
# It's good practice to use absolute paths or project-relative paths via `here::here()` if using RStudio Projects.
# .libPaths("/Users/cuitengfei/Graduate/Research/Cancer/SCAW/methods") # Ensure this path is correct or manage libraries globally
user_project_dir <- "/projects/gqilab/DAESC_GPU/data/simulation" # <<<< MODIFY THIS
source_file <- file.path(user_project_dir, "scripts", "ChIPtest_source.R") # Should contain simulate_hm, simulate_tf, and functions needed by PePr pipeline

# --- Load Libraries and Source Files ---
if (file.exists(source_file)) {
  source(source_file)
  cat("Sourced file:", source_file, "\n")
} else {
  stop("CRITICAL: Source file (e.g., ChIPtest_source.R) not found at: ", source_file)
}

# Make sure the PePr analysis function is defined (e.g., sourced or defined in this script)
if (!exists("run_pepr_pipeline") || !is.function(run_pepr_pipeline)) {
  stop("CRITICAL: The PePr analysis function 'run_pepr_pipeline' is not defined. Please source or define it.")
}
if (!exists("simulate_tf") || !is.function(simulate_tf) || !exists("simulate_hm") || !is.function(simulate_hm) ) {
  stop("CRITICAL: Simulation functions 'simulate_tf' and/or 'simulate_hm' are not defined. Please ensure they are in the sourced file.")
}


# --- Global Simulation Settings ---
num_iterations <- 100 # Total number of simulation iterations to run (e.g., 100)
IS_TF_SIMULATION <- FALSE # SET TO TRUE for TF simulation, FALSE for Histone Mark simulation

# --- PePr Pipeline Specific Global Parameters (can be overridden in function call if needed) ---
# These are passed to run_pepr_pipeline_v3
PEP_PEAK_WIDTH_TF <- 1000 # peak_width for find_overlap in TF case
PEP_PEAK_WIDTH_HM <- 500  # peak_width for find_overlap in HM case
PEP_GSIZE <- 300678703
PEP_FRAGLEN <- 100 # This is used for PePr's shift size calculation (fraglen/2L)

# --- Define Base Directories for this Simulation Batch ---
simulation_name_tag <- if (IS_TF_SIMULATION) "TF_PePr_Sim" else "HM1_PePr_Sim"

master_simulated_data_dir <- file.path(user_project_dir, "Data", simulation_name_tag)
master_analysis_results_dir <- file.path(user_project_dir, "Results", simulation_name_tag)

# --- Create Base Directories if they don't exist ---
if (!dir.exists(master_simulated_data_dir)) {
  dir.create(master_simulated_data_dir, recursive = TRUE)
  cat("Created master directory for simulated data:", master_simulated_data_dir, "\n")
}
if (!dir.exists(master_analysis_results_dir)) {
  dir.create(master_analysis_results_dir, recursive = TRUE)
  cat("Created master directory for analysis results:", master_analysis_results_dir, "\n")
}

# --- Simulation Parameters (customize based on IS_TF_SIMULATION) ---
# These parameters are passed to your simulate_tf or simulate_hm functions
if (IS_TF_SIMULATION) {
  sim_function_to_call <- simulate_tf
  # Parameters for simulate_tf from your example:
  # (data_path, fraglen, base.mu, up.mu, down.mu.1, down.mu.2, bins_vis, count_length_vis)
  sim_params <- list(
    fraglen = 150,      # Fragment length for simulation (distinct from PePr's fraglen for shift)
    base.mu = 30,
    up.mu = 90,
    down.mu.1 = 10,
    down.mu.2 = 10,
    bins_vis = 15,      # For visualization part within simulate_tf
    count_length_vis = 200 # For visualization part
  )
  analysis_peak_width_for_pepr <- PEP_PEAK_WIDTH_TF
} else { # HM Simulation
  sim_function_to_call <- simulate_hm
  # Parameters for simulate_hm from your example:
  # (true.width, base.mu, data_path, bins, count_length)
  sim_params <- list(
    true.width = 500,
    base.mu = 60,
    bins = 15,          # For visualization part within simulate_hm
    count_length = 200  # For visualization part
  )
  analysis_peak_width_for_pepr <- PEP_PEAK_WIDTH_HM
}

# --- Results Aggregation Setup ---
all_iterations_performance_metrics <- vector("list", num_iterations)
simulation_run_log <- data.frame(
  iteration = integer(num_iterations),
  iteration_tag_pepr = character(num_iterations),
  status = character(num_iterations),
  error_message = character(num_iterations),
  sim_data_dir = character(num_iterations),
  analysis_output_base_dir_pepr = character(num_iterations), # Base dir passed to PePr function
  analysis_output_specific_dir_pepr = character(num_iterations), # Actual dir created by PePr function
  start_time = as.POSIXct(rep(NA, num_iterations)),
  end_time = as.POSIXct(rep(NA, num_iterations)),
  duration_seconds = numeric(num_iterations),
  stringsAsFactors = FALSE
)
simulation_run_log$error_message <- NA_character_
simulation_run_log$status <- NA_character_


# --- Main Simulation Loop ---
cat(paste0("\n--- Starting simulation batch of ", num_iterations, " iterations for ", simulation_name_tag, " (", Sys.time(), ") ---\n"))

for (iter_num in 1:num_iterations) {
  cat(paste0("\n--- Starting Iteration: ", iter_num, "/", num_iterations, " (", Sys.time(), ") ---\n"))
  iter_start_time <- Sys.time()
  current_status <- "Pending"
  current_error_msg <- NA_character_
  pepr_iteration_tag <- paste0("iter", iter_num) # Tag for PePr function
  
  # --- Define and Create Iteration-Specific Directories ---
  # Directory for storing simulated data for this iteration
  iter_sim_data_dir <- file.path(master_simulated_data_dir, paste0("iter_", iter_num))
  # Parent directory for analysis results of this iteration (PePr function will create its own subfolder within this)
  iter_analysis_parent_dir_for_pepr <- file.path(master_analysis_results_dir, paste0("iter_", iter_num, "_analysis_files"))
  
  tryCatch({
    if (!dir.exists(iter_sim_data_dir)) {
      dir.create(iter_sim_data_dir, recursive = TRUE, showWarnings = FALSE)
    }
    # PePr function expects output_dir_base and creates its own structure within it.
    # So, iter_analysis_parent_dir_for_pepr is output_dir_base for PePr function for this iteration.
    if (!dir.exists(iter_analysis_parent_dir_for_pepr)) {
      dir.create(iter_analysis_parent_dir_for_pepr, recursive = TRUE, showWarnings = FALSE)
    }
    cat("Simulated data for iter", iter_num, "will be in:", iter_sim_data_dir, "\n")
    cat("PePr analysis output base for iter", iter_num, "will be in:", iter_analysis_parent_dir_for_pepr, "\n")
  }, error = function(e_dir) {
    stop(paste("CRITICAL: Failed to ensure directories for iteration", iter_num, ". Error:", conditionMessage(e_dir)))
  })
  
  analysis_results_obj <- NULL # Initialize
  
  tryCatch({
    # Step 1: Simulate data
    cat("Step 1: Running data simulation for iteration", iter_num, "...\n")
    # Pass data_path (iter_sim_data_dir) to the simulation function
    current_sim_args <- c(list(data_path = iter_sim_data_dir), sim_params)
    
    # The simulation function should create BAMs, _log.txt, and .RData files in iter_sim_data_dir
    # We assume it doesn't return anything critical for the next step, or handles its own saving.
    do.call(sim_function_to_call, current_sim_args)
    cat("Data simulation completed for iteration", iter_num, ". Output in:", iter_sim_data_dir, "\n")
    
    # Step 2: Run PePr Pipeline Analysis
    cat("Step 2: Running PePr Pipeline (run_pepr_pipeline) for iteration", iter_num, "...\n")
    
    analysis_results_obj <- run_pepr_pipeline(
      input_data_dir = iter_sim_data_dir,
      output_dir_base = iter_analysis_parent_dir_for_pepr, # PePr function will make its own sub-folder here
      is_tf = IS_TF_SIMULATION,
      iteration_tag = pepr_iteration_tag,
      peak_width = analysis_peak_width_for_pepr,
      gsize = PEP_GSIZE,
      fraglen = PEP_FRAGLEN,
      save_initial_pepr_results = FALSE,
      save_sites_overlap_csv = FALSE,
      save_final_filtered_results_csv = TRUE, # Save main results per iteration
      save_metrics_csv = TRUE,                # Save metrics per iteration
      save_roc_data_csv = TRUE                # Save ROC data per iteration
    )
    cat("PePr Pipeline analysis completed for iteration", iter_num, ".\n")
    
    if (!is.null(analysis_results_obj) && !is.null(analysis_results_obj$metrics)) {
      current_metrics_df <- analysis_results_obj$metrics
      # The 'iteration_tag' from PePr function should match our 'pepr_iteration_tag'
      # The metrics df from PePr already contains 'iteration_tag', 'peak_width', 'is_tf_analysis', 'runtime_seconds'
      # We add the main loop's iter_num for cross-referencing if needed.
      current_metrics_df$main_simulation_iter_num <- iter_num 
      all_iterations_performance_metrics[[iter_num]] <- current_metrics_df
      current_status <- "Success"
    } else {
      warning("No performance metrics returned from PePr pipeline for iteration ", iter_num)
      all_iterations_performance_metrics[[iter_num]] <- NULL
      current_status <- "CompletedWithWarning_NoMetrics"
      if (!is.null(analysis_results_obj) && !is.null(analysis_results_obj$error_message_from_pepr)) { # Assuming PePr func might return an error field
        current_error_msg <- analysis_results_obj$error_message_from_pepr
      }
    }
    
  }, error = function(e_main) {
    cat("ERROR in simulation/analysis for iteration", iter_num, ": ", conditionMessage(e_main), "\n")
    # Use <<- to assign to variables in the loop's parent environment (i.e., the for loop's variables)
    current_status <<- "Failed" 
    current_error_msg <<- conditionMessage(e_main)
    all_iterations_performance_metrics[[iter_num]] <- NULL
  }) # End tryCatch for simulation/analysis
  
  iter_end_time <- Sys.time()
  iter_duration_secs <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
  
  # Log details for the current iteration
  simulation_run_log[iter_num, "iteration"] <- iter_num
  simulation_run_log[iter_num, "iteration_tag_pepr"] <- pepr_iteration_tag
  simulation_run_log[iter_num, "status"] <- current_status
  simulation_run_log[iter_num, "error_message"] <- current_error_msg
  simulation_run_log[iter_num, "sim_data_dir"] <- iter_sim_data_dir
  simulation_run_log[iter_num, "analysis_output_base_dir_pepr"] <- iter_analysis_parent_dir_for_pepr
  simulation_run_log[iter_num, "analysis_output_specific_dir_pepr"] <- if(!is.null(analysis_results_obj) && !is.null(analysis_results_obj$output_directory)) analysis_results_obj$output_directory else NA_character_
  simulation_run_log[iter_num, "start_time"] <- iter_start_time
  simulation_run_log[iter_num, "end_time"] <- iter_end_time
  simulation_run_log[iter_num, "duration_seconds"] <- iter_duration_secs
  
  cat(paste0("--- Iteration ", iter_num, " finished. Status: ", current_status,
             ". Duration: ", round(iter_duration_secs, 2), " seconds. ---\n"))
  
  # gc() # Optional: Call gc() to free up memory after each iteration
} # End main simulation loop

# --- Post-Simulation Processing ---
cat("\n--- Aggregating results from all iterations ---\n")

valid_metrics_list <- all_iterations_performance_metrics[!sapply(all_iterations_performance_metrics, is.null)]

if (length(valid_metrics_list) > 0) {
  # Assuming each element in valid_metrics_list is a data.frame (from PePr function's $metrics output)
  final_performance_summary_df <- do.call(rbind, valid_metrics_list)
  
  # Save the aggregated performance metrics
  summary_filename <- file.path(master_analysis_results_dir, paste0("Aggregated_Metrics_", simulation_name_tag, ".csv"))
  tryCatch({
    write.csv(final_performance_summary_df, summary_filename, row.names = FALSE)
    cat("Aggregated performance metrics saved to:", summary_filename, "\n")
  }, error = function(e_save_metrics) {
    cat("Error saving aggregated performance metrics:", conditionMessage(e_save_metrics), "\n")
  })
  
  # Optional: Generate a summary plot or further analysis of final_performance_summary_df
  # For example, boxplots of AUC values:
  # if("auc_value" %in% names(final_performance_summary_df) && !all(is.na(final_performance_summary_df$auc_value))) {
  #   boxplot_filename <- file.path(master_analysis_results_dir, paste0("AUC_Boxplot_", simulation_name_tag, ".png"))
  #   png(boxplot_filename)
  #   boxplot(final_performance_summary_df$auc_value, main=paste("AUC Values for", simulation_name_tag), ylab="AUC")
  #   dev.off()
  #   cat("AUC boxplot saved to:", boxplot_filename, "\n")
  # }
  
} else {
  cat("No performance metrics were successfully collected to aggregate.\n")
  final_performance_summary_df <- data.frame() 
}

# Save the detailed simulation run log
log_filename <- file.path(master_analysis_results_dir, paste0("Run_Log_", simulation_name_tag, ".csv"))
tryCatch({
  write.csv(simulation_run_log, log_filename, row.names = FALSE)
  cat("Simulation run log saved to:", log_filename, "\n")
}, error = function(e_save_log) {
  cat("Error saving simulation log:", conditionMessage(e_save_log), "\n")
})

cat(paste0("\n--- Simulation batch for ", simulation_name_tag, " Finished (", Sys.time(), ") ---\n"))




























