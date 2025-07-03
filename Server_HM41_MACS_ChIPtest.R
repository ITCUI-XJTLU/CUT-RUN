# --- 主迭代运行脚本 (调用 run_macs3_chiptest_pipeline) ---

# --- 0. 设置与配置 ---
# 清理工作空间 (可选, 如果需要取消注释)
# rm(list = ls())
# gc()

set.seed(1)

# --- 路径定义 (请根据您的实际情况修改) ---
user_project_dir <- "/projects/gqilab/DAESC_GPU/data/simulation" # 项目根目录
# 主源文件，应包含 simulate_tf, simulate_hm, run_macs3_chiptest_pipeline, 及所有它们依赖的辅助函数
source_file_main_helpers_and_pipeline <- file.path(user_project_dir, "scripts", "ChIPtest_source.R")

# --- 1. 加载库与源文件 ---
if (file.exists(source_file_main_helpers_and_pipeline)) {
  source(source_file_main_helpers_and_pipeline)
  cat("已加载主源文件 (包含模拟、ChIPtest流程及辅助函数):", source_file_main_helpers_and_pipeline, "\n")
} else {
  stop("关键错误: 主源文件 (ChIPtest_source.R) 未在指定路径找到: ", source_file_main_helpers_and_pipeline)
}

# 确保 run_macs3_chiptest_pipeline 函数已定义
if (!exists("run_macs3_chiptest_pipeline") || !is.function(run_macs3_chiptest_pipeline)) {
  stop("关键错误: 'run_macs3_chiptest_pipeline' 函数未定义。请确保它在 ChIPtest_source.R 中。")
}
# 确保模拟函数已定义
if (!exists("simulate_tf") || !is.function(simulate_tf) || !exists("simulate_hm") || !is.function(simulate_hm) ) {
  stop("关键错误: 模拟函数 'simulate_tf' 和/或 'simulate_hm' 未定义。请确保它们在 ChIPtest_source.R 中。")
}

# --- 2. 全局模拟设置 ---
num_iterations <- 100 # <<<< 要运行的总迭代次数 (例如: 5, 100)
IS_TF_SIMULATION <- FALSE # <<<< TRUE 进行TF模拟, FALSE 进行组蛋白修饰(HM)模拟

# --- 3. run_macs3_chiptest_pipeline 函数的全局参数 ---
# 这些参数将传递给 run_macs3_chiptest_pipeline
# (注意: data_base_dir, output_base_dir, run_tag 会在循环中动态设置或基于 model_tag)
PIPELINE_PEAK_WIDTH_TF <- 150  # TF情况下用于 find_overlap 的 peak_width
PIPELINE_PEAK_WIDTH_HM <- 500  # HM情况下用于 find_overlap 的 peak_width (与您函数定义中的默认值一致)
PIPELINE_GSIZE <- 300678703
PIPELINE_FRAGLEN <- 100
PIPELINE_NUM_CORES_MCLAPPLY <- 8
PIPELINE_NUM_CORES_FINETUNE <- 5
# ChIPtest相关的微调参数 (可以使用函数中的默认值，或在此处覆盖)
# PIPELINE_CHIPTEST_BAND_VALUES <- c(5,10,30,60)
# PIPELINE_CHIPTEST_QUANTILE_VALUES <- c(-1, 0.0001, 0.01, 0.02, 0.1, 0.2, 0.4, 0.5, 0.99, 1.2, 1.5, 1.8, 2)
# PIPELINE_CHIPTEST_VAR_EST <- 1
# PIPELINE_CHIPTEST_VAR_THRED <- 0.01


# --- 4. 定义此模拟批次的模型标识符和基础目录 ---
# <<<< 为每个不同的模拟模型/参数集修改此标识符 >>>>
simulation_model_identifier <- if (IS_TF_SIMULATION) "TF_Model_ChIPtest_A" else "HM41_MACS_ChIPtest"

# 模拟数据的基础目录: user_project_dir/Data/[simulation_model_identifier]/iter_X/
model_specific_sim_data_parent_dir <- file.path(user_project_dir, "Data", simulation_model_identifier)
# 分析结果的基础目录: user_project_dir/Results/[simulation_model_identifier]/
# (run_macs3_chiptest_pipeline 会在此目录下根据 run_tag 创建子目录)
model_specific_results_parent_dir <- file.path(user_project_dir, "Results", simulation_model_identifier)

# 创建这些父级目录 (如果它们不存在)
if (!dir.exists(model_specific_sim_data_parent_dir)) {
  dir.create(model_specific_sim_data_parent_dir, recursive = TRUE)
  cat("已创建模型特定的模拟数据父目录:", model_specific_sim_data_parent_dir, "\n")
}
if (!dir.exists(model_specific_results_parent_dir)) {
  dir.create(model_specific_results_parent_dir, recursive = TRUE)
  cat("已创建模型特定的分析结果父目录:", model_specific_results_parent_dir, "\n")
}

# --- 5. 模拟函数参数配置 (基于 IS_TF_SIMULATION) ---
if (IS_TF_SIMULATION) {
  sim_function_to_call <- simulate_tf
  # simulate_tf 的参数 (请确保这些与您 ChIPtest_source.R 中定义的函数参数匹配)
  sim_params_list <- list(
    fraglen = 150,
    base.mu = 30,
    up.mu = 90,
    down.mu.1 = 10,
    down.mu.2 = 10,
    bins_vis = 15,
    count_length_vis = 200
  )
  pipeline_peak_width_to_use <- PIPELINE_PEAK_WIDTH_TF
} else { # HM 模拟
  sim_function_to_call <- simulate_hm
  # simulate_hm 的参数
  sim_params_list <- list(
    true.width = 500,
    base.mu = 200,
    bins = 15,
    count_length = 200
  )
  pipeline_peak_width_to_use <- PIPELINE_PEAK_WIDTH_HM
}

# --- 6. 结果聚合设置 ---
all_iterations_metrics_summary <- vector("list", num_iterations)
overall_run_log <- data.frame(
  iteration_number = integer(num_iterations),
  pipeline_run_tag = character(num_iterations), # run_tag 传递给 pipeline 函数
  iteration_status = character(num_iterations),
  error_details = character(num_iterations),
  simulated_data_path_iter = character(num_iterations),   # 存储模拟数据的迭代子目录
  analysis_results_path_iter = character(num_iterations), # pipeline 函数输出的特定迭代结果目录
  iteration_start_time = as.POSIXct(rep(NA, num_iterations)),
  iteration_end_time = as.POSIXct(rep(NA, num_iterations)),
  iteration_duration_seconds = numeric(num_iterations),
  stringsAsFactors = FALSE
)
overall_run_log$error_details <- NA_character_
overall_run_log$iteration_status <- "Not Started"

# --- 7. 主模拟与分析循环 ---
cat(paste0("\n--- 开始模拟批次: ", num_iterations, " 次迭代, 模型: ", simulation_model_identifier, " (", Sys.time(), ") ---\n"))

for (i in 1:num_iterations) {
  cat(paste0("\n--- 开始迭代: ", i, "/", num_iterations, " (", Sys.time(), ") ---\n"))
  iter_loop_start_time <- Sys.time() # 记录迭代循环的开始时间
  
  current_iter_status <- "Pending_Simulation"
  current_iter_error <- NA_character_
  
  # 为当前迭代生成唯一的 run_tag (用于 pipeline 函数内部创建子目录)
  # 并定义迭代特定的模拟数据目录和分析输出基础目录
  current_pipeline_run_tag <- paste0(simulation_model_identifier, "_iter", i)
  
  # 迭代特定的模拟数据目录 (例如: .../Data/HM_Model_ChIPtest_B/iter_1/)
  iter_specific_sim_data_dir <- file.path(model_specific_sim_data_parent_dir, paste0("iter_", i))
  
  # run_macs3_chiptest_pipeline 的 output_base_dir 将是 model_specific_results_parent_dir
  # run_macs3_chiptest_pipeline 的 run_tag 将是 current_pipeline_run_tag
  # 这意味着最终的迭代结果会存放在:
  # model_specific_results_parent_dir/macs_chiptest_analysis_[current_pipeline_run_tag]/
  
  pipeline_output_path_for_this_iter <- NA_character_ # 初始化
  
  tryCatch({
    # 确保迭代特定的模拟数据目录存在
    if (!dir.exists(iter_specific_sim_data_dir)) {
      dir.create(iter_specific_sim_data_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cat("迭代 ", i, " 的模拟数据将存放于: ", iter_specific_sim_data_dir, "\n")
    cat("迭代 ", i, " 的分析结果将基于父目录: ", model_specific_results_parent_dir, " 使用run_tag: ", current_pipeline_run_tag, "\n")
    
    # --- 步骤 1: 数据模拟 ---
    cat("步骤 1: 开始为迭代 ", i, " 进行数据模拟...\n")
    # 将 iter_specific_sim_data_dir 作为 data_path 参数传递给模拟函数
    simulation_args <- c(list(data_path = iter_specific_sim_data_dir), sim_params_list)
    
    # 执行模拟函数 (它应该在 iter_specific_sim_data_dir 中创建BAM, log, RData文件)
    simulation_output_details <- do.call(sim_function_to_call, simulation_args)
    cat("迭代 ", i, " 的数据模拟完成。输出位于: ", iter_specific_sim_data_dir, "\n")
    current_iter_status <- "Pending_Analysis"
    
    # --- 步骤 2: 运行 MACS3 + ChIPtest 分析流程 ---
    cat("步骤 2: 开始为迭代 ", i, " 运行 run_macs3_chiptest_pipeline 分析...\n")
    
    analysis_output <- run_macs3_chiptest_pipeline(
      data_base_dir = iter_specific_sim_data_dir, # 这是包含BAM/log/RData的目录
      output_base_dir = model_specific_results_parent_dir, # 分析结果的父目录
      run_tag = current_pipeline_run_tag,         # 本次运行的特定标签
      peak_width = pipeline_peak_width_to_use,
      is_tf = IS_TF_SIMULATION,
      source_file_path = source_file_main_helpers_and_pipeline,
      gsize = PIPELINE_GSIZE,
      fraglen = PIPELINE_FRAGLEN,
      num_cores_mclapply = PIPELINE_NUM_CORES_MCLAPPLY,
      num_cores_finetune = PIPELINE_NUM_CORES_FINETUNE
      # 可以按需传递其他 chiptest_* 参数
    )
    cat("迭代 ", i, " 的 run_macs3_chiptest_pipeline 分析完成。\n")
    
    if (!is.null(analysis_output) && !is.null(analysis_output$metrics_ChIPtest_df)) {
      metrics_from_iter <- analysis_output$metrics_ChIPtest_df
      metrics_from_iter$simulation_model_id <- simulation_model_identifier # 添加模型标识
      metrics_from_iter$iteration_num_overall <- i                   # 添加整体迭代号
      all_iterations_metrics_summary[[i]] <- metrics_from_iter
      current_iter_status <- "Success"
    } else {
      warning("迭代 ", i, " 的 run_macs3_chiptest_pipeline 未返回有效的 metrics_ChIPtest_df。")
      all_iterations_metrics_summary[[i]] <- NULL # 或一个有正确列名的空数据框
      current_iter_status <- "Completed_NoMetrics"
    }
    
    if(!is.null(analysis_output) && !is.null(analysis_output$output_path)){
      pipeline_output_path_for_this_iter <- analysis_output$output_path
    }
    
  }, error = function(e) {
    cat("错误发生在迭代 ", i, " 的模拟或分析阶段: ", conditionMessage(e), "\n")
    current_iter_status <<- paste0("Failed_", current_iter_status) # 记录失败发生在哪一步
    current_iter_error <<- conditionMessage(e)
    all_iterations_metrics_summary[[i]] <- NULL
  }) # 结束 tryCatch
  
  iter_loop_end_time <- Sys.time() # 记录迭代循环的结束时间
  iter_loop_duration_secs <- as.numeric(difftime(iter_loop_end_time, iter_loop_start_time, units = "secs"))
  
  # --- 记录本次迭代的日志信息 ---
  overall_run_log[i, "iteration_number"] <- i
  overall_run_log[i, "pipeline_run_tag"] <- current_pipeline_run_tag
  overall_run_log[i, "iteration_status"] <- current_iter_status
  overall_run_log[i, "error_details"] <- current_iter_error
  overall_run_log[i, "simulated_data_path_iter"] <- iter_specific_sim_data_dir
  overall_run_log[i, "analysis_results_path_iter"] <- pipeline_output_path_for_this_iter
  overall_run_log[i, "iteration_start_time"] <- iter_loop_start_time
  overall_run_log[i, "iteration_end_time"] <- iter_loop_end_time
  overall_run_log[i, "iteration_duration_seconds"] <- iter_loop_duration_secs
  
  cat(paste0("--- 迭代 ", i, " 完成。状态: ", current_iter_status,
             "。耗时: ", round(iter_loop_duration_secs, 2), " 秒. ---\n"))
  
  # gc() # 可选: 每次迭代后进行垃圾回收
} # 结束主模拟与分析循环

# --- 8. 模拟后处理与结果汇总 ---
cat("\n--- 开始汇总所有迭代的结果 (模型: ", simulation_model_identifier, ") ---\n")

# 过滤掉因迭代失败可能产生的NULL元素
valid_metrics <- all_iterations_metrics_summary[!sapply(all_iterations_metrics_summary, is.null)]

if (length(valid_metrics) > 0) {
  # 假设 all_iterations_metrics_summary 中的每个元素都是一个数据框
  aggregated_metrics_df <- do.call(rbind, valid_metrics)
  
  # 将汇总的性能指标保存到模型特定的结果目录中
  aggregated_metrics_filename <- file.path(model_specific_results_parent_dir, paste0("Aggregated_Metrics_", simulation_model_identifier, ".csv"))
  tryCatch({
    write.csv(aggregated_metrics_df, aggregated_metrics_filename, row.names = FALSE)
    cat("汇总的性能指标已保存至:", aggregated_metrics_filename, "\n")
  }, error = function(e_save) {
    cat("保存汇总性能指标时发生错误:", conditionMessage(e_save), "\n")
  })
  
} else {
  cat("未能成功收集到任何性能指标进行汇总 (模型: ", simulation_model_identifier, ").\n")
}

# 保存详细的运行日志到模型特定的结果目录中
overall_log_filename <- file.path(model_specific_results_parent_dir, paste0("Overall_Run_Log_", simulation_model_identifier, ".csv"))
tryCatch({
  write.csv(overall_run_log, overall_log_filename, row.names = FALSE)
  cat("详细的运行日志已保存至:", overall_log_filename, "\n")
}, error = function(e_save_log) {
  cat("保存运行日志时发生错误:", conditionMessage(e_save_log), "\n")
})

cat(paste0("\n--- 模拟批次 (模型: ", simulation_model_identifier, ") 已完成 (", Sys.time(), ") ---\n"))









































