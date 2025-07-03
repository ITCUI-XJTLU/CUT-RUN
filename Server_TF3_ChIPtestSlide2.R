# --- Setup and Configuration ---
# Clear workspace (optional, uncomment if needed)
# rm(list = ls())
# gc()

# Set seed for reproducibility for the entire simulation batch
set.seed(1)

# --- Paths (MODIFY THESE AS NEEDED) ---
user_home_dir <- "/projects/gqilab/DAESC_GPU/data/simulation" # Adjust to your project root
source_file <- file.path(user_home_dir, "scripts/ChIPtest_source.R") # Contains simulate_hm, simulate_tf, run_chiptest_slide2_analysis, etc.

# --- Load Libraries and Source Files ---
# Ensure necessary packages are loaded here if not handled by the sourced file or functions
# For example: library(csaw); library(edgeR); library(data.table); etc.
# However, run_chiptest_slide2_analysis has requireNamespace checks.

if (file.exists(source_file)) {
  source(source_file)
  cat("Sourced file:", source_file, "\n")
} else {
  stop("CRITICAL: Source file (ChIPtest_source.R) not found: ", source_file)
}

# --- Global Simulation Settings ---
num_iterations <- 100 # Total number of simulation iterations to run
IS_TF_SIMULATION <- TRUE # SET TO TRUE for TF simulation, FALSE for Histone Mark simulation

# --- Define Base Directories for this Simulation Batch ---
# These will contain iteration-specific subfolders.
simulation_name_tag <- if (IS_TF_SIMULATION) "TF3_ChIPtestSlide2" else "HM4_ChIPtestSlide2"

# Base directory for storing simulated data (each iteration in a subfolder)
master_simulated_data_dir <- file.path(user_home_dir, "Data", paste0(simulation_name_tag, "_SimData"))

# Base directory for storing analysis results (each iteration's analysis output in a subfolder)
master_analysis_results_dir <- file.path(user_home_dir, "Results", paste0(simulation_name_tag, "_SimResults"))

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
if (IS_TF_SIMULATION) {
  sim_function_to_call <- simulate_tf
  # Parameters for simulate_tf: data_path, fraglen, base.mu, up.mu, down.mu.1, down.mu.2, bins_vis, count_length_vis
  sim_params <- list(
    fraglen = 100,         # Example
    base.mu = 120,          # Example
    up.mu = 180,            # Example
    down.mu.1 = 60,        # Example for group A's "down" regions compared to baseline for other group
    down.mu.2 = 60,        # Example for group B's "down" regions compared to baseline for other group
    bins_vis = 15,         # Example for visualization part within simulate_tf
    count_length_vis = 200 # Example for visualization part
  )
  analysis_peak_width <- 1000 # peak_width for run_chiptest_slide2_analysis when is.tf=TRUE
  # chiptest_params for run_chiptest_slide2_analysis can use its internal defaults for TF
  # or be specified here if you want to override:
  # specific_chiptest_params_for_run <- list(bins_tf = 25, counts_tf = 120, ...)
} else { # HM Simulation
  sim_function_to_call <- simulate_hm
  # Parameters for simulate_hm: true.width, base.mu, data_path, bins, count_length
  sim_params <- list(
    true.width = 500,
    base.mu = 30,
    bins = 15,          # For visualization part within simulate_hm
    count_length = 200  # For visualization part
  )
  analysis_peak_width <- 500 # peak_width for run_chiptest_slide2_analysis when is.tf=FALSE
  # chiptest_params for run_chiptest_slide2_analysis can use its internal defaults for HM
  # or be specified here:
  # specific_chiptest_params_for_run <- list(bins_hm = 15, counts_hm = 120, ...)
}


# --- Results Aggregation Setup ---
all_iterations_performance_metrics <- vector("list", num_iterations)
simulation_run_log <- data.frame(
  iteration = integer(num_iterations),
  status = character(num_iterations),
  error_message = character(num_iterations),
  sim_data_dir = character(num_iterations),
  analysis_output_parent_dir = character(num_iterations),
  start_time = as.POSIXct(rep(NA, num_iterations)),
  end_time = as.POSIXct(rep(NA, num_iterations)),
  duration_seconds = numeric(num_iterations),
  stringsAsFactors = FALSE
)
simulation_run_log$error_message <- NA_character_
simulation_run_log$status <- NA_character_


# --- Main Simulation Loop ---
cat(paste0("\n--- Starting simulation batch of ", num_iterations, " iterations for ", simulation_name_tag, " (", Sys.time(), ") ---\n"))

for (iter in 1:num_iterations) {
  cat(paste0("\n--- Starting Iteration: ", iter, "/", num_iterations, " (", Sys.time(), ") ---\n"))
  iter_start_time <- Sys.time()
  current_status <- "Pending"
  current_error_msg <- NA_character_
  
  # --- Define and Create Iteration-Specific Directories ---
  iter_sim_data_dir <- file.path(master_simulated_data_dir, paste0("iter_", iter))
  iter_analysis_output_parent_dir <- file.path(master_analysis_results_dir, paste0("iter_", iter)) # Analysis function will create its own subfolder inside this
  
  tryCatch({
    if (!dir.exists(iter_sim_data_dir)) {
      dir.create(iter_sim_data_dir, recursive = TRUE, showWarnings = FALSE)
      cat("Created simulation data directory for iteration", iter, ":", iter_sim_data_dir, "\n")
    } else {
      cat("Simulation data directory for iteration", iter, "already exists:", iter_sim_data_dir, ". Will use/overwrite.\n")
      # Consider cleaning if necessary: unlink(iter_sim_data_dir, recursive = TRUE); dir.create(...)
    }
    if (!dir.exists(iter_analysis_output_parent_dir)) {
      dir.create(iter_analysis_output_parent_dir, recursive = TRUE, showWarnings = FALSE)
      cat("Created analysis output parent directory for iteration", iter, ":", iter_analysis_output_parent_dir, "\n")
    } else {
      cat("Analysis output parent directory for iteration", iter, "already exists:", iter_analysis_output_parent_dir, ". Analysis will write into it.\n")
    }
  }, error = function(e_dir) {
    stop(paste("CRITICAL: Failed to create directories for iteration", iter, ". Error:", conditionMessage(e_dir)))
  })
  
  tryCatch({
    # Step 1: Simulate data
    cat("Step 1: Running data simulation for iteration", iter, "...\n")
    # Pass data_path to the simulation function
    current_sim_args <- c(list(data_path = iter_sim_data_dir), sim_params)
    sim_results_obj <- do.call(sim_function_to_call, current_sim_args)
    # sim_results_obj might contain paths to BAM files, log file etc., if needed later, though run_chiptest_slide2_analysis reconstructs them.
    cat("Data simulation completed for iteration", iter, ". Output in:", iter_sim_data_dir, "\n")
    
    # Step 2: Run ChIPtest Slide-Window 2 analysis
    cat("Step 2: Running run_chiptest_slide2_analysis for iteration", iter, "...\n")
    
    # Prepare arguments for run_chiptest_slide2_analysis
    # The function uses internal defaults for many params, we only need to provide the core ones.
    # If you defined 'specific_chiptest_params_for_run' above, pass it via 'chiptest_params' argument.
    analysis_args <- list(
      input_data_dir = iter_sim_data_dir,
      output_parent_dir = iter_analysis_output_parent_dir,
      is.tf = IS_TF_SIMULATION,
      peak_width = analysis_peak_width,
      chiptest_params = list(
        bins_tf = 8,
        counts_tf = 120,
        filter_tf = 0,
        window_width_csaw_tf = 10,
        bins_hm = 10,
        counts_hm = 120,
        filter_hm = 0,
        window_width_csaw_hm = 150
      )
    )
    
    analysis_results <- do.call(run_chiptest_slide2_analysis, analysis_args)
    cat("run_chiptest_slide2_analysis completed for iteration", iter, ".\n")
    
    if (!is.null(analysis_results) && !is.null(analysis_results$performance_metrics)) {
      # Add iteration number to the performance metrics df for easier aggregation if not already present
      current_metrics <- analysis_results$performance_metrics
      if(is.data.frame(current_metrics) && !"iteration" %in% colnames(current_metrics)){
        current_metrics$iteration <- iter
      }
      all_iterations_performance_metrics[[iter]] <- current_metrics
      current_status <- "Success"
    } else {
      warning("No performance metrics returned from run_chiptest_slide2_analysis for iteration ", iter)
      all_iterations_performance_metrics[[iter]] <- NULL # Mark as NULL if no metrics
      if (!is.null(analysis_results$error)) {
        current_status <- paste0("CompletedWithWarning_NoMetrics (Error: ", analysis_results$error,")")
        current_error_msg <- analysis_results$error
      } else {
        current_status <- "CompletedWithWarning_NoMetrics"
      }
    }
    
  }, error = function(e_main) {
    cat("ERROR in simulation/analysis for iteration", iter, ": ", conditionMessage(e_main), "\n")
    current_status <<- "Failed" # Use <<- to assign to 'current_status' in the loop's environment
    current_error_msg <<- conditionMessage(e_main)
    all_iterations_performance_metrics[[iter]] <- NULL
  }) # End tryCatch for simulation/analysis
  
  iter_end_time <- Sys.time()
  iter_duration_secs <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
  
  # Log details for the current iteration
  simulation_run_log[iter, "iteration"] <- iter
  simulation_run_log[iter, "status"] <- current_status
  simulation_run_log[iter, "error_message"] <- current_error_msg
  simulation_run_log[iter, "sim_data_dir"] <- iter_sim_data_dir
  simulation_run_log[iter, "analysis_output_parent_dir"] <- iter_analysis_output_parent_dir
  simulation_run_log[iter, "start_time"] <- iter_start_time
  simulation_run_log[iter, "end_time"] <- iter_end_time
  simulation_run_log[iter, "duration_seconds"] <- iter_duration_secs
  
  cat(paste0("--- Iteration ", iter, " finished. Status: ", current_status,
             ". Duration: ", round(iter_duration_secs, 2), " seconds. ---\n"))
  
  # gc() # Optional: Call gc() to free up memory after each iteration
} # End main simulation loop

# --- Post-Simulation Processing ---
cat("\n--- Aggregating results from all iterations ---\n")

# Combine all valid performance metrics into a single data frame
# Filter out NULL entries (failed or no-metric iterations)
valid_metrics_list <- all_iterations_performance_metrics[!sapply(all_iterations_performance_metrics, is.null)]

if (length(valid_metrics_list) > 0) {
  # If each element is already a data frame (as it should be from run_chiptest_slide2_analysis)
  final_performance_summary_df <- do.call(rbind, valid_metrics_list)
  
  # Save the aggregated performance metrics
  # Use master_analysis_results_dir for aggregated files, or a new specific summary dir
  summary_filename <- file.path(master_analysis_results_dir, paste0("Aggregated_Metrics_", simulation_name_tag, ".csv"))
  tryCatch({
    write.csv(final_performance_summary_df, summary_filename, row.names = FALSE)
    cat("Aggregated performance metrics saved to:", summary_filename, "\n")
  }, error = function(e_save_metrics) {
    cat("Error saving aggregated performance metrics:", conditionMessage(e_save_metrics), "\n")
  })
} else {
  cat("No performance metrics were successfully collected to aggregate.\n")
  final_performance_summary_df <- data.frame() # Create an empty data frame if no results
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