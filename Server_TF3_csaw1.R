# --- Setup and Configuration ---
# Clear workspace (optional, uncomment if needed)
# rm(list = ls())
# gc()

# Set seed for reproducibility for the entire simulation batch
set.seed(1)

# --- Paths (MODIFY THESE AS NEEDED) ---
user_home_dir <- "/projects/gqilab/DAESC_GPU/data/simulation" # Adjust to your project root
source_file <- file.path(user_home_dir, "scripts/ChIPtest_source.R") # Contains simulate_hm, simulate_tf, run_csaw_differential_analysis, etc.

# --- Load Libraries and Source Files ---
if (file.exists(source_file)) {
  source(source_file)
  cat("Sourced file:", source_file, "\n")
} else {
  stop("CRITICAL: Source file (ChIPtest_source.R) not found: ", source_file)
}

# --- Global Simulation Settings ---
num_iterations <- 100      # Total simulation iterations PER window_width setting
IS_TF_SIMULATION <- TRUE  # SET TO TRUE for TF simulation, FALSE for Histone Mark simulation

# --- Define window widths to test (remains the same) ---
# window_widths_to_test <- c(50, 150, 250, 350, 550, 1000)
window_widths_to_test <- c(10)

# --- Simulation Parameters (customize based on IS_TF_SIMULATION) ---
simulation_name_prefix <- if (IS_TF_SIMULATION) "TF3" else "HM"

if (IS_TF_SIMULATION) {
  sim_function_to_call <- simulate_tf
  # Parameters for simulate_tf: data_path, fraglen, base.mu, up.mu, down.mu.1, down.mu.2, bins_vis, count_length_vis
  # Ensure these parameters are appropriate for your simulate_tf function
  sim_params_all_iterations <- list(
    fraglen = 100,
    base.mu = 120,
    up.mu = 180,
    down.mu.1 = 60,
    down.mu.2 = 60,
    bins_vis = 15,        # For visualization part within simulate_tf
    count_length_vis = 200 # For visualization part
  )
  # Base Parameters for run_csaw_differential_analysis for TF
  csaw_params_base_all_iterations <- list(
    is.tf = TRUE,
    peak_width = 500, # TF peak width for CSAW analysis; check if simulate_tf internal true.width matches
    grouping_vector = c("A", "A", "B", "B"),
    csaw_count_filter_value = 20,
    csaw_extension_length = 100, # TF-specific extension?
    csaw_binned_width_bg = 8000,
    csaw_merge_tolerance = 100,
    csaw_merge_max_width = 5000,
    save_results_df = TRUE,
    save_roc_data = TRUE
  )
} else { # HM Simulation
  sim_function_to_call <- simulate_hm
  # Parameters for simulate_hm: true.width, base.mu, data_path, bins, count_length
  sim_params_all_iterations <- list(
    true.width = 500,
    base.mu = 30,
    bins = 15,          # For visualization part within simulate_hm
    count_length = 200  # For visualization part
  )
  # Base Parameters for run_csaw_differential_analysis for HM
  csaw_params_base_all_iterations <- list(
    is.tf = FALSE,
    peak_width = 500, # HM peak width for CSAW analysis; matches simulate_hm$true.width
    grouping_vector = c("A", "A", "B", "B"),
    csaw_count_filter_value = 20,
    csaw_extension_length = 100, # HM-specific extension?
    csaw_binned_width_bg = 8000,
    csaw_merge_tolerance = 100,
    csaw_merge_max_width = 5000,
    save_results_df = TRUE,
    save_roc_data = TRUE
  )
}

cat(paste0("\n=== STARTING ENTIRE CSAW SIMULATION SUITE for ", simulation_name_prefix, " (", Sys.time(), ") ===\n"))

# --- Outer Loop for Different CSAW Window Widths ---
for (current_csaw_ww in window_widths_to_test) {
  cat(paste0("\n\n--- Processing for ", simulation_name_prefix, " with CSAW Window Width: ", current_csaw_ww, " (", Sys.time(), ") ---\n"))
  
  # --- Update csaw_params for the current window width ---
  current_run_csaw_params <- csaw_params_base_all_iterations # This now correctly has is.tf set
  current_run_csaw_params$window_width <- current_csaw_ww     # Set the looping window_width
  
  # --- Define Paths based on current_csaw_ww and simulation type ---
  data_dir_suffix <- paste0(simulation_name_prefix, "_csaw_ww", current_csaw_ww, "_data")
  results_dir_suffix <- paste0(simulation_name_prefix, "_csaw_ww", current_csaw_ww, "_results")
  
  scenario_simulation_data_base_dir <- file.path(user_home_dir, "Data", data_dir_suffix)
  scenario_main_csaw_results_dir <- file.path(user_home_dir, "Results", results_dir_suffix)
  
  # --- Create Directories for the current scenario if they don't exist ---
  if (!dir.exists(scenario_simulation_data_base_dir)) {
    dir.create(scenario_simulation_data_base_dir, recursive = TRUE)
    cat("Created base simulation data directory:", scenario_simulation_data_base_dir, "\n")
  }
  if (!dir.exists(scenario_main_csaw_results_dir)) {
    dir.create(scenario_main_csaw_results_dir, recursive = TRUE)
    cat("Created main CSAW results directory:", scenario_main_csaw_results_dir, "\n")
  }
  
  # --- Results Aggregation Setup (re-initialize for each window_width scenario) ---
  all_iterations_performance_metrics <- vector("list", num_iterations)
  simulation_run_log <- data.frame(
    iteration = integer(num_iterations),
    csaw_window_width = integer(num_iterations),
    status = character(num_iterations),
    error_message = character(num_iterations),
    sim_data_dir = character(num_iterations),
    analysis_output_dir = character(num_iterations),
    start_time = as.POSIXct(rep(NA, num_iterations)),
    end_time = as.POSIXct(rep(NA, num_iterations)),
    duration_seconds = numeric(num_iterations),
    stringsAsFactors = FALSE
  )
  simulation_run_log$error_message <- NA_character_
  simulation_run_log$status <- NA_character_
  
  cat(paste0("\n--- Starting simulation batch of ", num_iterations, " iterations for ",
             simulation_name_prefix, ", ww=", current_csaw_ww, " (", Sys.time(), ") ---\n"))
  
  # --- Main Simulation Loop (for iterations) ---
  for (iter in 1:num_iterations) {
    cat(paste0("\n--- Starting Iteration: ", iter, "/", num_iterations,
               " for ", simulation_name_prefix, ", ww=", current_csaw_ww, " (", Sys.time(), ") ---\n"))
    iter_start_time <- Sys.time()
    current_status <- "Pending"
    current_error_msg <- NA_character_
    
    # --- Define and Create Iteration-Specific Data Directory ---
    iteration_specific_data_dir <- file.path(scenario_simulation_data_base_dir, paste0("iter_", iter))
    
    tryCatch({
      if (!dir.exists(iteration_specific_data_dir)) {
        dir.create(iteration_specific_data_dir, recursive = TRUE, showWarnings = FALSE)
        cat("Created data directory for iter ", iter, ": ", iteration_specific_data_dir, "\n")
      } else {
        cat("Data directory for iter ", iter, " already exists: ", iteration_specific_data_dir, ". Will use/overwrite.\n")
      }
    }, error = function(e_dir_create) {
      stop(paste("CRITICAL: Failed to create iteration-specific data directory for iter", iter,
                 ". Path:", iteration_specific_data_dir, ". Error:", conditionMessage(e_dir_create)))
    })
    
    tryCatch({
      # Step 1: Simulate data
      cat("Step 1: Running data simulation (", simulation_name_prefix, ") for iteration", iter, "...\n")
      current_sim_args_iter <- c(list(data_path = iteration_specific_data_dir), sim_params_all_iterations)
      sim_output_obj <- do.call(sim_function_to_call, current_sim_args_iter)
      cat("Data simulation completed for iteration", iter, ". Output in:", iteration_specific_data_dir, "\n")
      
      # Step 2: Run CSAW differential analysis
      cat("Step 2: Running run_csaw_differential_analysis for iteration", iter, "...\n")
      # run_csaw_differential_analysis takes 'output_dir' where it saves iteration-specific files if its internal flags are TRUE
      # and where the aggregated metrics will also be eventually stored by this script.
      analysis_args_iter <- c(
        list(
          input_data_dir = iteration_specific_data_dir,
          output_dir = scenario_main_csaw_results_dir, # CSAW iteration files (if any) go here
          iteration_number = iter # Crucial for run_csaw_differential_analysis to name its own files
        ),
        current_run_csaw_params # Contains the correct is.tf and current_csaw_ww
      )
      analysis_results <- do.call(run_csaw_differential_analysis, analysis_args_iter)
      cat("run_csaw_differential_analysis completed for iteration", iter, ".\n")
      
      if (!is.null(analysis_results) && !is.null(analysis_results$performance)) {
        metrics_df <- analysis_results$performance
        # Add csaw_window_width to the metrics for easier aggregation if not already present
        if(is.data.frame(metrics_df) && !"csaw_window_width" %in% colnames(metrics_df)){
          metrics_df$csaw_window_width <- current_csaw_ww
        }
        all_iterations_performance_metrics[[iter]] <- metrics_df
        current_status <- "Success"
      } else {
        warning("No performance metrics returned from run_csaw_differential_analysis for iter ", iter)
        all_iterations_performance_metrics[[iter]] <- NULL
        if (!is.null(analysis_results$runtime_seconds) && is.null(analysis_results$performance)) { # It ran but no metrics
          current_status <- "CompletedWithWarning_NoMetrics"
        } else { # Could be an error object from tryCatch if run_csaw itself failed badly
          current_status <- "Failed_AnalysisError"
        }
      }
      
    }, error = function(e_main_iter) {
      cat("ERROR in simulation/analysis for iter", iter, ": ", conditionMessage(e_main_iter), "\n")
      current_status <<- "Failed"
      current_error_msg <<- conditionMessage(e_main_iter)
      all_iterations_performance_metrics[[iter]] <- NULL
    })
    
    iter_end_time <- Sys.time()
    iter_duration_secs <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
    
    simulation_run_log[iter, "iteration"] <- iter
    simulation_run_log[iter, "csaw_window_width"] <- current_csaw_ww
    simulation_run_log[iter, "status"] <- current_status
    simulation_run_log[iter, "error_message"] <- current_error_msg
    simulation_run_log[iter, "sim_data_dir"] <- iteration_specific_data_dir
    simulation_run_log[iter, "analysis_output_dir"] <- scenario_main_csaw_results_dir # Main dir for this scenario
    simulation_run_log[iter, "start_time"] <- iter_start_time
    simulation_run_log[iter, "end_time"] <- iter_end_time
    simulation_run_log[iter, "duration_seconds"] <- iter_duration_secs
    
    cat(paste0("--- Iteration ", iter, " for ", simulation_name_prefix, ", ww=", current_csaw_ww, " finished. Status: ", current_status,
               ". Duration: ", round(iter_duration_secs, 2), " seconds. ---\n"))
    # gc()
  } # End main iteration loop
  
  # --- Post-Processing for the current CSAW window_width scenario ---
  cat(paste0("\n--- Post-processing results for ", simulation_name_prefix, ", ww=", current_csaw_ww, " ---\n"))
  valid_metrics_list <- all_iterations_performance_metrics[!sapply(all_iterations_performance_metrics, is.null)]
  
  if (length(valid_metrics_list) > 0) {
    final_performance_summary_df_ww <- do.call(rbind, valid_metrics_list)
    summary_filename <- file.path(scenario_main_csaw_results_dir, paste0("Aggregated_Metrics_", data_dir_suffix, ".csv")) # Filename includes ww
    tryCatch({
      write.csv(final_performance_summary_df_ww, summary_filename, row.names = FALSE)
      cat("Aggregated performance metrics for current scenario saved to:", summary_filename, "\n")
    }, error = function(e_save) {
      cat("Error saving aggregated performance metrics:", conditionMessage(e_save), "\n")
    })
  } else {
    cat("No performance metrics were successfully collected to aggregate for this scenario.\n")
  }
  
  log_filename_ww <- file.path(scenario_main_csaw_results_dir, paste0("Run_Log_", data_dir_suffix, ".csv")) # Filename includes ww
  tryCatch({
    write.csv(simulation_run_log, log_filename_ww, row.names = FALSE)
    cat("Simulation run log for current scenario saved to:", log_filename_ww, "\n")
  }, error = function(e_save_log) {
    cat("Error saving simulation log:", conditionMessage(e_save_log), "\n")
  })
  
  cat(paste0("\n--- Finished processing for ", simulation_name_prefix, ", CSAW Window Width: ", current_csaw_ww, " (", Sys.time(), ") ---\n"))
  
} # End outer loop for CSAW window_widths_to_test

cat(paste0("\n\n=== ENTIRE CSAW SIMULATION SUITE for ", simulation_name_prefix, " FINISHED (", Sys.time(), ") ===\n"))