library(GenomicRanges)
library(SummarizedExperiment)
library(metapod)
library(InteractionSet)
library(diffHic)
library(csaw)
library(xscss)
library(DiffBind)

require(xscss)
require(csaw)
require(edgeR)

library(dplyr)
# library(ggplot2)
library(data.table)
library(parallel)
library(ChIPtest)

library(pROC)
# library(ggplot2)
library(gridExtra)
library(ROCR)
library(parallel)
library(dplyr)
library(tidyr)
library(LaplacesDemon) # KL 
library(parallel)
library(entropy)
library(data.table)

grouping <- c("A","A","B","B")
design<-model.matrix(~factor(grouping))
# dump2file <- function(id, cutoff, result) {
#   write.table(file=result.file, data.frame(id, cutoff, 1-result$overlap/result$found, result$recall), 
#               sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
# }

my_assessChIP <- function (observed, known, tol = 200, checkfc = TRUE) 
{
  # browser()
  obs <- read.table(observed, header = TRUE)
  # browser()
  if (checkfc) {
    up.o <- obs$logFC > 0
    if (is.null(up.o)) {
      stop("need a log-FC field in the table of observed sites")
    }
  }
  else {
    up.o <- rep(TRUE, nrow(obs))
  }
  obranges <- GRanges(obs$chr, IRanges(obs$start, obs$end))
  if (!is.na(tol)) {
    out <- mergeWindows(obranges, sign = up.o, tol = tol)
    obranges <- out$region
    all.signs <- logical(length(obranges))
    all.signs[out$id] <- up.o
    up.o <- all.signs
  }
  kx <- read.table(known, header = TRUE)
  if (checkfc) {
    up.t <- kx$logFC > 0
    if (is.null(up.t)) {
      stop("need a log-FC field in the table of known sites")
    }
  }
  else {
    up.t <- rep(TRUE, nrow(kx))
  }
  if (!nrow(kx)) {
    stop("no known sites to check against")
  }
  kranges <- GRanges(kx[, 1], IRanges(kx[, 2], kx[, 3]))
  if (is.null(kx$name)) {
    kranges$name <- 1:length(kranges)
  }
  else {
    kranges$name <- kx$name
  }
  known.up <- kranges[up.t]
  known.down <- kranges[!up.t]
  u.olap <- findOverlaps(known.up, obranges[up.o])
  d.olap <- findOverlaps(known.down, obranges[!up.o])
  recall <- length(unique(known.up$name[queryHits(u.olap)])) + 
    length(unique(known.down$name[queryHits(d.olap)]))
  overlapped <- length(unique(subjectHits(u.olap))) + length(unique(subjectHits(d.olap)))
  found <- length(obranges)
  return(list(overlap = overlapped, found = found, recall = recall))
}


my_est.c <- function (data1, data4, max1 = 5, max4 = 5) 
{
  # browser()
  M = 0
  Dmean = 0
  dat1 = apply(data1, 1, max)
  dat4 = apply(data4, 1, max)
  
  # change the way to select sites, we decide to use the difference
  # loc41 = which(dat1 <= max1 & dat4 <= max4)
  loc41 = which(abs(dat1 - dat4) < max1)
  
  for (k in 1:dim(data1)[1]) {
    n = dim(data1)[2]
    L1 = data1[k, ]
    L4 = data4[k, ]
    Diff = L4 - L1
    D41_upp = c(Diff[2:n])
    D41_low = c(Diff[1:(n - 1)])
    M[k] = sum(D41_upp * D41_low)/n
    Dmean[k] = mean(Diff)
  }
  tao = mean(M[loc41])
  return(tao)
}

my_TS_twosample <- function (data1, data4, tao, band, quant, var.est = 1, 
                             var.thred = 0.01) 
{
  # browser()
  n = dim(data1)[2]
  hwidth = band/n
  x = c(1:n)/n
  Snw <- matrix(0, nrow = n, ncol = n)
  In <- diag(rep(1, n))
  for (j in 1:n) {
    y <- In[, j]
    Snw[, j] <- ksmooth(x, y, kernel = "normal", bandwidth = hwidth, 
                        x.points = x)$y
  }
  ad_df = sum(diag(Snw))
  Sev = 0
  Suv = 0
  Dsum = 0
  Deql = 0
  Dnun = 0
  sigma1 = 0
  sigma4 = 0
  M = 0
  sigma41 = 0
  sigma_1_sq = 0
  sigma_4_sq = 0
  sigma1_fg = 0
  sigma4_fg = 0
  sig4sig1 = 0
  sigma_unequal = 0
  Ts_yvec = 0
  pairT = 0
  pairT.sigma = 0
  sigma = 0
  Xg = 0
  # browser()
  for (k in 1:dim(data1)[1]) {
    n = dim(data1)[2]
    L1 = data1[k, ]
    L4 = data4[k, ]
    Diff = L4 - L1
    d41 = ksmooth(x, Diff, "normal", bandwidth = hwidth)
    Xg[k] = sqrt(sum((d41$y - Diff)^2)/(n - ad_df))
    Ts_yvec[k] = mean((d41$y)^2)
    D41_upp = c(Diff[2:n])
    D41_low = c(Diff[1:(n - 1)])
    M[k] = sum(D41_upp * D41_low)/(n - 1)
    sigma1[k] = mean((L1[2:(length(L1))] - L1[1:(length(L1) - 
                                                   1)])^2)/2
    sigma4[k] = mean((L4[2:(length(L4))] - L4[1:(length(L4) - 
                                                   1)])^2)/2
    sigma41[k] = (sigma1[k] + sigma4[k])^2 + 4 * M[k] * (sigma1[k] + 
                                                           sigma4[k])
    sigma_1_sq[k] = sum((L1[2:(length(L1) - 2)] - L1[1:(length(L1) - 
                                                          3)])^2 * (L1[4:(length(L1))] - L1[3:(length(L1) - 
                                                                                                 1)])^2)/(4 * (n - 3))
    sigma_4_sq[k] = sum((L4[2:(length(L4) - 2)] - L4[1:(length(L4) - 
                                                          3)])^2 * (L4[4:(length(L4))] - L4[3:(length(L4) - 
                                                                                                 1)])^2)/(4 * (n - 3))
    sigma1_fg[k] = sum((L4[1:(n - 2)] - L1[1:(n - 2)]) * 
                         (c(L4[1], L4[1:(n - 3)]) - c(L1[1], L1[1:(n - 3)])) * 
                         (L1[3:n] - L1[2:(n - 1)])^2)/(2 * (n - 3))
    sigma4_fg[k] = sum((L4[1:(n - 2)] - L1[1:(n - 2)]) * 
                         (c(L4[1], L4[1:(n - 3)]) - c(L1[1], L1[1:(n - 3)])) * 
                         (L4[3:n] - L4[2:(n - 1)])^2)/(2 * (n - 3))
    sig4sig1[k] = sum((L1[2:n] - L1[1:(n - 1)])^2 * (L4[2:n] - 
                                                       L4[1:(n - 1)])^2)/(4 * (n - 3))
    sigma_unequal[k] = (sigma_1_sq[k] + 4 * sigma1_fg[k]) + 
      (sigma_4_sq[k] + 4 * sigma4_fg[k]) + 2 * sig4sig1[k]
    Sev[k] = sigma41[k]
    Suv[k] = sigma_unequal[k]
  }
  # Suv[which(Suv <= 0)] = min(Suv[which(Suv > 0)])
  # Sev[which(Sev <= 0)] = min(Sev[which(Sev > 0)])
  # Xg[which(Xg <= 0)] = min(Xg[which(Xg > 0)])
  Suv[which(Suv <= 0)] <- median(Suv[which(Suv > 0)])
  Sev[which(Sev <= 0)] <- median(Sev[which(Sev > 0)])
  Xg[which(Xg <= 0)] <- median(Xg[which(Xg > 0)])
  
  
  Sev_a0 = sqrt(Sev)
  Suv_a0 = sqrt(Suv)
  CHQBC_1_adjB = Xg
  
  if(var.est==1){
    Sev <- Sev 
    Suv <- Suv
    if (quant[1] < 0){
      Sev_a0 = sqrt(Sev)
    }else if(quant[1] > 0 && quant[1] < 1){
      Sev_a0 = sqrt(Sev) + quantile(sqrt(Sev), quant[1]) 
    }else{
      Sev_a0 = sqrt(Sev) + quantile(sqrt(Sev), 0.99) * quant[1]
    }
    
    if (quant[2] < 0){
      Suv_a0 = sqrt(Suv) 
    }else if(quant[2] > 0 && quant[2] < 1){
      Suv_a0 = sqrt(Suv) + quantile(sqrt(Suv), quant[2]) 
    }else{
      Suv_a0 = sqrt(Suv) + quantile(sqrt(Suv), 0.99)*quant[2] 
    }
    
    if (quant[3] < 0){
      CHQBC_1_adjB = Xg
    }else if(quant[3] > 0 && quant[3] < 1){
      CHQBC_1_adjB = Xg + quantile(Xg, quant[3])
    }else{
      CHQBC_1_adjB = Xg + quantile(Xg, 0.99)*quant[3]
    }
    
  }else{
    Sev_a0[Sev < var.thred] <- sqrt(Sev[Sev < var.thred]) + median(sqrt(Sev))
    Suv_a0[Suv < var.thred] <- sqrt(Suv[Suv < var.thred]) + median(sqrt(Suv))
  }
  
  Dsum = M
  Deql = (Dsum - tao)/(Sev_a0/sqrt(n - 1))
  Dnun = (Dsum - tao)/(Suv_a0/sqrt(n - 1))
  
  Tsb = 1/(2 * (band * 0.37) * sqrt(pi))
  vu = sqrt(2 * pi)/(2 * pi * n * (band * 0.37))
  Amax = Snw %*% Snw
  eigenvalue = as.numeric(eigen(Amax)$values)
  d = sum(eigenvalue)^2/sum(eigenvalue^2)
  delta = sum(eigenvalue)/(n * d)
  
  Test.adj = (Ts_yvec)/(CHQBC_1_adjB^2)
  Test.adj_1 = (Test.adj/(delta * d))
  Test.adj_2 = sign(Test.adj_1) * abs(Test.adj_1)^(1/3)
  
  #original setting (I think it does not correct, so I try my own idea on 2025.4.15)
  #Test.adj_3 = Test.adj_2 - mean(Test.adj_2) + (1 - 2/(9 * d))
  #TS_kn = (( Test.adj_3 - (1 - 2/(9 * d)))/sqrt(2/(9 * d)))
  
  # my new idea
  print("new idea")
  Test.adj_3 = Test.adj_2 - (1 - 2/(9 * d)) 
  TS_kn = Test.adj_3 /sqrt(2/(9 * d))
  
  return(a = list(sigma1 = sigma1, sigma4 = sigma4, TS_kn = TS_kn, 
                  Ts_yvec = Ts_yvec, Dsum = Dsum, Deql = Deql, Dnun = Dnun, 
                  Sev = Sev, Suv = Suv, Xg = Xg))
}

###################################################################################
# 先定义两个抽取 count sequence 的函数
process_peak_conA <- function(x, df, width) {
  start_row <- x - (width / 2)
  end_row <- x + (width / 2)
  return(c(df[start_row, 1], df[end_row, 2], df[start_row:end_row, 3]))
}
process_peak_conB <- function(x, df, width) {
  start_row <- x - (width / 2)
  end_row <- x + (width / 2)
  return(c(df[start_row, 1], df[end_row, 2], df[start_row:end_row, 3]))
}
process_peak <- function(x, df, width) {
  start_row <- x - (width / 2)
  end_row <- x + (width / 2)
  return(c(df[start_row, 1], df[end_row, 2], df[start_row:end_row, 3]))
}

################################################################################
## 为找到的sites定义overlap标签的函数

# 
# find_overlap <- function(is.tf, pos.1, pos.3, peak_width, lfile, csaw_sites) {
#   print("new overlap function (modified to keep unmatched all_sites)")
#   # browser()
#   # 1. 构建 all_sites 区间
#   if (!is.tf) {
#     all_sites <- data.frame(start = pos.1 - peak_width / 2,
#                             end = pos.3 + peak_width / 2)  # for HM
#   } else {
#     all_sites <- data.frame(start = pos.1 - peak_width / 2,
#                             end = pos.1 + peak_width / 2)  # for TF
#   }
#   
#   # 2. 读取真实 peak 信息
#   hist_log <- read.table(lfile, header = TRUE, sep = "\t")
#   if (!is.tf) {
#     True_sites <- hist_log %>%
#       dplyr::group_by(name) %>%
#       dplyr::summarise(
#         chr = dplyr::first(chr),
#         start = min(start),
#         end = max(end)
#       )
#   } else {
#     True_sites <- hist_log
#   }
# 
#   # 3. 转换为 data.table 并设置 key
#   data.table::setDT(all_sites)
#   data.table::setDT(True_sites)
#   data.table::setDT(csaw_sites)
#   data.table::setkey(csaw_sites, start, end)
#   data.table::setkey(True_sites, start, end)
# 
#   # 4. 判断 all_sites 是否与真实 peak 重叠
#   result <- data.table::foverlaps(all_sites, True_sites,
#                                   by.x = c("start", "end"),
#                                   by.y = c("start", "end"),
#                                   type = "any", nomatch = 0L)
#   all_sites$overlap <- ifelse(all_sites$start %in% result$i.start, 1, 0)
# 
#   # 5. 判断 all_sites 与 csaw_sites 的重叠
#   overlap_csaw <- data.table::foverlaps(all_sites, csaw_sites,
#                                         by.x = c("start", "end"),
#                                         by.y = c("start", "end"),
#                                         type = "any", nomatch = NA)
# 
#   data.table::setnames(overlap_csaw, old = c("i.start", "i.end"), new = c("sim_start", "sim_end"))
#   
#   # 7. 对未匹配到 csaw_sites 的行（index 为 NA）填充默认值
#   overlap_csaw[is.na(index), `:=`(
#     Pvalues = 1,
#     FDR = 1
#   )]
#   return(overlap_csaw)
# }
# 
# 

find_overlap <- function(is.tf, pos.1, pos.3, peak_width, lfile, csaw_sites) {
  # 1. 构建 all_sites 区间
  if (!is.tf) {
    all_sites <- data.frame(start = pos.1 - peak_width / 2,
                            end = pos.3 + peak_width / 2) # for HM
  } else {
    all_sites <- data.frame(start = pos.1 - peak_width / 2,
                            end = pos.1 + peak_width / 2) # for TF
  }
  
  # 2. 读取真实 peak 信息
  hist_log <- read.table(lfile, header = TRUE, sep = "\t")
  if (!is.tf) {
    True_sites <- hist_log %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(
        chr = dplyr::first(chr),
        start = min(start),
        end = max(end)
      )
  } else {
    True_sites <- hist_log
  }
  
  # 3. 转换为 data.table 并设置 key
  data.table::setDT(all_sites)
  data.table::setDT(True_sites)
  data.table::setDT(csaw_sites)
  data.table::setkey(all_sites, start, end)
  data.table::setkey(True_sites, start, end)
  data.table::setkey(csaw_sites, start, end)
  
  # 4. 判断 all_sites 与 True_sites 的重叠
  overlap_result <- data.table::foverlaps(all_sites,True_sites,
                                          by.x = c("start", "end"),
                                          by.y = c("start", "end"),
                                          type = "any", nomatch = NA)
  
  data.table::setnames(overlap_result,
                       old = c("i.start", "i.end", "start", "end"),
                       new = c("sim_start", "sim_end", "true_start", "true_end"))
  overlap_result[, overlap_true := ifelse(is.na(true_start), 0, 1)]
  overlap_result <- overlap_result[,c("sim_start", "sim_end", "overlap_true")]
  
  # 将 all_sites 与 True_sites 的重叠信息合并到 csaw_sites 中
  data.table::setkey(overlap_result, sim_start, sim_end)
  data.table::setkey(csaw_sites, start, end)
  
  final_results <- foverlaps(overlap_result, csaw_sites,
                               by.x = c("sim_start", "sim_end"),
                               by.y = c("start", "end"),
                               type = "any", nomatch = 0L)
  
  unmatched_from_overlap <- overlap_result[!final_results[, .(sim_start, sim_end)], on = c("sim_start", "sim_end")]
  if (nrow(unmatched_from_overlap) > 0) {
    unmatched_from_overlap_formatted <- data.table(
      start = NA_real_,
      end = NA_real_,
      index = NA_integer_,
      sim_start = unmatched_from_overlap$sim_start,
      sim_end = unmatched_from_overlap$sim_end,
      overlap_true = unmatched_from_overlap$overlap_true
    )
  } else {
    unmatched_from_overlap_formatted <- data.table()
  }
  
  matched_csaw_in_final <- final_results[!is.na(start) & !is.na(end), .(start, end)]
  setkey(matched_csaw_in_final, start, end) # 为查找设置 key
  unmatched_from_csaw <- csaw_sites[!matched_csaw_in_final, on = c("start", "end")]
  if (nrow(unmatched_from_csaw) > 0) {
    unmatched_from_csaw_formatted <- data.table(
      start = unmatched_from_csaw$start,
      end = unmatched_from_csaw$end,
      index = unmatched_from_csaw$index,
      sim_start = NA_real_,
      sim_end = NA_real_,
      overlap_true = NA_integer_ # 可以是 NA 或你认为合适的默认值，如 0
    )
  } else {
    unmatched_from_csaw_formatted <- data.table()
  }

  final_output <- rbindlist(list(final_results,
                                 unmatched_from_overlap_formatted,
                                 unmatched_from_csaw_formatted),
                            fill = TRUE)
  return(final_output)
}








###############################################################################
## 可视化函数 
finetune <- function(df, width, potential_peaks, my_quantile, my_band,
                     my_var.est, my_var.thred,
                     plot_filename_base = NULL, # New argument for plot filenames
                     num_cores = 8) { # Added num_cores as an argument
  
  # Ensure mclapply or lapply is used based on OS and num_cores
  apply_fn <- if (.Platform$OS.type != "windows" && num_cores > 1) {
    cat("Finetune: Using mclapply with", num_cores, "cores.\n")
    function(...) parallel::mclapply(..., mc.cores = num_cores)
  } else {
    if (.Platform$OS.type == "windows" && num_cores > 1) {
      cat("Finetune: mclapply not available on Windows, using sequential lapply.\n")
    } else {
      cat("Finetune: Using sequential lapply.\n")
    }
    lapply
  }
  
  # browser() # Debugging line, can be removed
  
  # Assuming df has columns: start, end, X1, X2, X3, X4, conA, conB in that order at minimum
  # df[,c(1,2,3)] would use the 3rd column (X1 counts)
  # df[,c(1,2,4)] would use the 4th column (X2 counts) etc.
  data_list_con1 <- apply_fn(potential_peaks, process_peak_conA, df = df[,c(1,2,3)], width = width)
  data_con1 <- do.call(rbind, lapply(data_list_con1, t))
  
  data_list_con2 <- apply_fn(potential_peaks, process_peak_conB, df = df[,c(1,2,4)], width = width)
  data_con2 <- do.call(rbind, lapply(data_list_con2, t))
  
  data_list_con3 <- apply_fn(potential_peaks, process_peak_conA, df = df[,c(1,2,5)], width = width)
  data_con3 <- do.call(rbind, lapply(data_list_con3, t))
  
  data_list_con4 <- apply_fn(potential_peaks, process_peak_conB, df = df[,c(1,2,6)], width = width)
  data_con4 <- do.call(rbind, lapply(data_list_con4, t))
  
  # 数据都已准备，ChIPtest做假设检验
  Data_con1 = NormTransformation(data_con1[,-c(1,2)])
  Data_con2 = NormTransformation(data_con2[,-c(1,2)])
  Data_con3 = NormTransformation(data_con3[,-c(1,2)])
  Data_con4 = NormTransformation(data_con4[,-c(1,2)])
  
  tao_12=est.c(Data_con1, Data_con2, max1=4, max4=4) # Assuming est.c is defined
  tao_34=est.c(Data_con3, Data_con4, max1=4, max4=4) # Assuming est.c is defined
  
  band=my_band
  TS_null12=my_TS_twosample(Data_con1, Data_con2, tao_12, band,
                            quant = my_quantile, var.est = my_var.est,
                            var.thred = my_var.thred) # Assuming my_TS_twosample is defined
  TS_null34=my_TS_twosample(Data_con3, Data_con4, tao_34, band,
                            quant = my_quantile, var.est = my_var.est,
                            var.thred = my_var.thred) # Assuming my_TS_twosample is defined
  
  # --- Plotting Section ---
  
  # Plot 1: Distributions of TS_kn, Deql, Dnun
  if (!is.null(plot_filename_base)) {
    plot_file_1 <- paste0(plot_filename_base, "_distributions.png")
    png(filename = plot_file_1, width = 12, height = 8, units = "in", res = 300)
    cat("Saving distribution plots to:", plot_file_1, "\n")
  }
  par(mfrow = c(2, 3))
  hist(TS_null12[["TS_kn"]], breaks = 150, freq = FALSE, border = "black",
       main = "Histogram of TS_kn of X1 X2", xlab = "TS_kn", ylab = "Frequency")
  x_norm <- seq(min(TS_null12[["TS_kn"]], -3, na.rm=TRUE), max(TS_null12[["TS_kn"]], 3, na.rm=TRUE), length = 100)
  y_norm <- dnorm(x_norm, mean = 0, sd = 1)
  lines(x_norm, y_norm, col = "red", lwd = 2)
  
  hist(TS_null12[["Deql"]], breaks = 150, freq = FALSE,
       col = "lightblue", border = "black", main = "Histogram of Deql of X1 X2", xlim = c(-5, 5),
       xlab = "Deql", ylab = "Frequency")
  x_norm2 <- seq(-5, 5, length = 100)
  y_norm2 <- dnorm(x_norm2, mean = 0, sd = 1)
  lines(x_norm2, y_norm2, col = "red", lwd = 2)
  
  hist(TS_null12[["Dnun"]], breaks = 150, freq = FALSE,
       col = "lightgreen", border = "black", main = "Histogram of Dnun of X1 X2", xlim = c(-5, 5),
       xlab = "Dnun", ylab = "Frequency")
  lines(x_norm2, y_norm2, col = "red", lwd = 2) # Using x_norm2, y_norm2 from Deql plot
  
  hist(TS_null34[["TS_kn"]], breaks = 150, freq = FALSE, border = "black",
       main = "Histogram of TS_kn of X3 X4", xlab = "TS_kn", ylab = "Frequency")
  x_norm3 <- seq(min(TS_null34[["TS_kn"]], -3, na.rm=TRUE), max(TS_null34[["TS_kn"]], 3, na.rm=TRUE), length = 100)
  y_norm3 <- dnorm(x_norm3, mean = 0, sd = 1)
  lines(x_norm3, y_norm3, col = "red", lwd = 2)
  
  
  hist(TS_null34[["Deql"]], breaks = 150, freq = FALSE,
       col = "lightblue", border = "black", main = "Histogram of Deql of X3 X4", xlim = c(-5, 5),
       xlab = "Deql", ylab = "Frequency")
  lines(x_norm2, y_norm2, col = "red", lwd = 2)
  
  hist(TS_null34[["Dnun"]], breaks = 150, freq = FALSE,
       col = "lightgreen", border = "black", main = "Histogram of Dnun of X3 X4", xlim = c(-5, 5),
       xlab = "Dnun", ylab = "Frequency")
  lines(x_norm2, y_norm2, col = "red", lwd = 2)
  if (!is.null(plot_filename_base)) {
    dev.off()
  }
  
  # Plot 2: Distributions of Xg, Sev, Suv
  if (!is.null(plot_filename_base)) {
    plot_file_2 <- paste0(plot_filename_base, "_variances.png")
    png(filename = plot_file_2, width = 12, height = 8, units = "in", res = 300)
    cat("Saving variance component plots to:", plot_file_2, "\n")
  }
  par(mfrow = c(2, 3))
  hist(TS_null12[["Xg"]],breaks = 100,main = "Histogram of Xg of X1 X2", xlab="Xg")
  hist(TS_null12[["Sev"]],breaks = 100, main = "Histogram of Sev of X1 X2", xlab="Sev")
  hist(TS_null12[["Suv"]],breaks = 100, main = "Histogram of Suv of X1 X2", xlab="Suv")
  
  hist(TS_null34[["Xg"]],breaks = 100,main = "Histogram of Xg of X3 X4", xlab="Xg")
  hist(TS_null34[["Sev"]],breaks = 100, main = "Histogram of Sev of X3 X4", xlab="Sev")
  hist(TS_null34[["Suv"]],breaks = 100, main = "Histogram of Suv of X3 X4", xlab="Suv")
  if (!is.null(plot_filename_base)) {
    dev.off()
  }
  
  # Plot 3: P-value distributions
  if (!is.null(plot_filename_base)) {
    plot_file_3 <- paste0(plot_filename_base, "_pvalues.png")
    png(filename = plot_file_3, width = 12, height = 8, units = "in", res = 300)
    cat("Saving P-value distribution plots to:", plot_file_3, "\n")
  }
  par(mfrow = c(2, 3))
  hist(pnorm(TS_null12[["TS_kn"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of TS_kn of X1 X2", xlab="P-value")
  hist(pnorm(TS_null12[["Deql"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of Deql of X1 X2", xlab="P-value")
  hist(pnorm(TS_null12[["Dnun"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of Dnun of X1 X2", xlab="P-value")
  
  hist(pnorm(TS_null34[["TS_kn"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of TS_kn of X3 X4", xlab="P-value")
  hist(pnorm(TS_null34[["Deql"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of Deql of X3 X4", xlab="P-value")
  hist(pnorm(TS_null34[["Dnun"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "P-values of Dnun of X3 X4", xlab="P-value")
  if (!is.null(plot_filename_base)) {
    dev.off()
  }
  
  # Reset par to default single plot
  par(mfrow = c(1, 1))
  
  # The function originally didn't return anything explicitly,
  # its purpose was to generate plots and potentially modify the environment
  # if TS_null12 or TS_null34 were used by the caller (they are not in this context).
  # For now, let's keep it that way, or return the TS objects if useful.
  # invisible(list(TS_null12 = TS_null12, TS_null34 = TS_null34)) # Optional return
}
# finetune <- function(df, width, potential_peaks, my_quantile, my_band, my_var.est, my_var.thred){
#   num_cores <- 8
#   # browser()
#   data_list_con1 <- mclapply(potential_peaks, process_peak_conA, df = count_bin2[,c(1,2,3)], width = width, mc.cores = num_cores)
#   data_con1 <- do.call(rbind, lapply(data_list_con1, t))
#   
#   data_list_con2 <- mclapply(potential_peaks, process_peak_conB, df = count_bin2[,c(1,2,4)], width = width, mc.cores = num_cores)
#   data_con2 <- do.call(rbind, lapply(data_list_con2, t))
#   
#   data_list_con3 <- mclapply(potential_peaks, process_peak_conA, df = count_bin2[,c(1,2,5)], width = width, mc.cores = num_cores)
#   data_con3 <- do.call(rbind, lapply(data_list_con3, t))
#   
#   data_list_con4 <- mclapply(potential_peaks, process_peak_conB, df = count_bin2[,c(1,2,6)], width = width, mc.cores = num_cores)
#   data_con4 <- do.call(rbind, lapply(data_list_con4, t))
#   
#   # 数据都已准备，ChIPtest做假设检验
#   Data_con1 = NormTransformation(data_con1[,-c(1,2)])
#   Data_con2 = NormTransformation(data_con2[,-c(1,2)])
#   Data_con3 = NormTransformation(data_con3[,-c(1,2)])
#   Data_con4 = NormTransformation(data_con4[,-c(1,2)])
#   
#   tao_12=est.c(Data_con1, Data_con2, max1=4, max4=4)
#   tao_34=est.c(Data_con3, Data_con4, max1=4, max4=4)
#   
#   # max1 <- 5
#   # max4 <- 5
#   # dat1 = apply(Data_con1, 1, max)
#   # dat4 = apply(Data_con2, 1, max)
#   # loc41 = which(dat1 <= max1 & dat4 <= max4)
#   
#   band=my_band
#   TS_null12=my_TS_twosample(Data_con1, Data_con2, tao_12, band, 
#                             quant = my_quantile, var.est = my_var.est, 
#                             var.thred = my_var.thred)
#   TS_null34=my_TS_twosample(Data_con3, Data_con4, tao_34, band, 
#                             quant = my_quantile, var.est = my_var.est, 
#                             var.thred = my_var.thred)
#   
#   # 可视化，打印三个统计量(TS_kn, Deql, Dnun)分布
#   par(mfrow = c(2, 3))
#   hist(TS_null12[["TS_kn"]], breaks = 150, freq = FALSE, border = "black",
#        main = "Histogram of TS_kn of X1 X2", xlab = "TS_kn", ylab = "Frequency")
#   x <- seq(-3, 3, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   
#   hist(TS_null12[["Deql"]], breaks = 150, freq = FALSE,
#        col = "blue", border = "black", main = "Histogram of Deql of X1 X2", xlim = c(-5, 5),
#        xlab = "Deql", ylab = "Frequency")
#   x <- seq(-5, 5, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   
#   hist(TS_null12[["Dnun"]], breaks = 150, freq = FALSE,
#        col = "blue", border = "black", main = "Histogram of Dnun of X1 X2", xlim = c(-5, 5),
#        xlab = "Dnun", ylab = "Frequency")
#   x <- seq(-5, 5, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   
#   hist(TS_null34[["TS_kn"]], breaks = 150, freq = FALSE, border = "black",
#        main = "Histogram of TS_kn of X3 X4", xlab = "TS_kn", ylab = "Frequency")
#   x <- seq(-3, 3, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   
#   hist(TS_null34[["Deql"]], breaks = 150, freq = FALSE,
#        col = "blue", border = "black", main = "Histogram of Deql of X3 X4", xlim = c(-5, 5),
#        xlab = "Deql", ylab = "Frequency")
#   x <- seq(-5, 5, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   
#   hist(TS_null34[["Dnun"]], breaks = 150, freq = FALSE,
#        col = "blue", border = "black", main = "Histogram of Dnun of X3 X4", xlim = c(-5, 5),
#        xlab = "Dnun", ylab = "Frequency")
#   x <- seq(-5, 5, length = 100)
#   y <- dnorm(x, mean = 0, sd = 1)
#   lines(x, y, col = "red", lwd = 2)
#   par(mfrow = c(1, 1))
#   
#   # 打印variantes (Sev, Suv, Xg)
#   par(mfrow = c(2, 3))
#   hist(TS_null12[["Xg"]],breaks = 100,main = "Histogram of Xg of X1 X2")
#   hist(TS_null12[["Sev"]],breaks = 100, main = "Histogram of Sev of X1 X2")
#   hist(TS_null12[["Suv"]],breaks = 100, main = "Histogram of Suv of X1 X2")
#   
#   hist(TS_null34[["Xg"]],breaks = 100,main = "Histogram of Xg of X3 X4")
#   hist(TS_null34[["Sev"]],breaks = 100, main = "Histogram of Sev of X3 X4")
#   hist(TS_null34[["Suv"]],breaks = 100, main = "Histogram of Suv of X3 X4")
#   
#   par(mfrow = c(1, 1))
#   
#   # 展示在Null下，P value的分布
#   par(mfrow = c(2, 3))
#   hist(pnorm(TS_null12[["TS_kn"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of TS_kn of X1 X2")
#   hist(pnorm(TS_null12[["Deql"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of Deql of X1 X2")
#   hist(pnorm(TS_null12[["Dnun"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of Dnun of X1 X2")
#   
#   
#   hist(pnorm(TS_null34[["TS_kn"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of TS_kn of X3 X4")
#   hist(pnorm(TS_null34[["Deql"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of Deql of X3 X4")
#   hist(pnorm(TS_null34[["Dnun"]], mean = 0, sd = 1, lower.tail = FALSE),breaks = 100, main = "Pvalues of Dnun of X3 X4")
#   
#   par(mfrow = c(1, 1))
# }

# random select wrong and right answer, plot their curves
plot_curves <- function(ChIPtest_result, data_conA, data_conB, sort_num = 10){
  metrics <- c("pvalue_TS_high", "Dnun_adjust", "Deql_adjust")
  short_names <- c("TS_kn", "Dnun", "Deql")  
  plot_list <- list()  
  
  for (j in seq_along(metrics)) {
    metric <- metrics[j]
    short_name <- short_names[j]
    
    # 按指标排序并筛选前十个 index
    top10_overlap_1 <- ChIPtest_result[order(ChIPtest_result[[metric]]) & ChIPtest_result$overlap == 1][1:sort_num, "index"][[1]]
    top10_overlap_0 <- ChIPtest_result[order(ChIPtest_result[[metric]]) & ChIPtest_result$overlap == 0][1:sort_num, "index"][[1]]
    
    # 随机选择 3 个 index
    index_overlap_1 <- sample(top10_overlap_1, 3)
    index_overlap_0 <- sample(top10_overlap_0, 3)
    
    # 合并 index 和对应的标签
    selected_indices <- c(index_overlap_1, index_overlap_0)
    # overlap_labels <- c(rep(paste0(short_name, "_R"), 3), rep(paste0(short_name, "_W"), 3))
    overlap_labels <- c(
      paste0(short_name, "_R (Index:", index_overlap_1, ")"),
      paste0(short_name, "_W (Index:", index_overlap_0, ")")
    )
    
    # 针对每个 index，生成折线图
    for (i in seq_along(selected_indices)) {
      index <- selected_indices[i]
      overlap_label <- overlap_labels[i]
      
      # 获取 data_conA 和 data_conB 中对应行的数据
      data_conA_row <- as.numeric(data_conA[index, -c(1,2)])
      data_conB_row <- as.numeric(data_conB[index, -c(1,2)])
      
      # 创建一个数据框来绘图
      df_plot <- data.frame(
        Position = 1:length(data_conA_row),
        conA = data_conA_row,
        conB = data_conB_row
      )
      
      # 创建折线图并添加到 plot_list
      p <- ggplot(df_plot, aes(x = Position)) +
        geom_line(aes(y = conA, color = "conA")) +
        geom_line(aes(y = conB, color = "conB")) +
        labs(title = overlap_label, y = "Value", x = "Position") +
        theme_minimal() +
        theme(legend.position = "none") +  # 去除图例
        scale_color_manual(values = c("conA" = "blue", "conB" = "red"))
      
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
  grid.arrange(grobs = plot_list, ncol = 6, nrow = 3)
}


########################################################################################
## 一个用于 MACS2 的函数
my_runMACS2 <- function(file, outprefix, macs.path = "macs2", threshold = NULL, 
                        gsize = "mm", fraglen = NULL, cmd.only = FALSE, extra = NULL, 
                        format = "BED", conda.env = "macs2_env", shell = "zsh") 
{
  # 更新zsh配置并激活conda环境
  conda_activate <- paste("source ~/.zshrc && conda activate", conda.env)
  
  # 生成MACS2命令
  cmds <- c(macs.path, "callpeak -t", file, "--gsize", gsize, 
            "--keep-dup=all", "-f", format, "--outdir", dirname(outprefix), 
            "-n", basename(outprefix), extra)
  
  if (!is.null(fraglen)) {
    cmds <- c(cmds, "--nomodel --extsize", fraglen)
  }
  if (!is.null(threshold)) {
    cmds <- c(cmds, "-p", threshold)
  }
  
  # 连接命令字符串
  full_cmd <- paste(conda_activate, "&&", paste(cmds, collapse = " "))
  print(full_cmd)
  if (cmd.only) {
    return(full_cmd)
  }
  
  # 使用子shell运行完整的命令字符串
  if (system(paste(shell, "-c", shQuote(full_cmd)))) {
    stop("running MACS2 failed")
  }
  return(invisible(NULL))
}

## function for HAOMER to do the peak calling 
run_HOMER <- function(file, outprefix, is.tf = TRUE,
                      conda.env = "homer_env", shell = "zsh",
                      fraglen = 100, cmd.only = FALSE) {
  # 激活 Conda 环境
  conda_activate <- paste("source ~/.zshrc && conda activate", conda.env)
  
  # 生成 tag 目录和 peak 文件路径
  tag_dir <- paste0(outprefix, "_tagdir")
  peak_file <- paste0(outprefix, "_peaks.txt")
  
  # 根据 TF/Histone 选择不同的 findPeaks 参数
  peak_style <- ifelse(is.tf, "factor", "histone")
  tbp_option <- ifelse(is.tf, "-tbp 0", "")  # 仅 TF 需要 tbp 选项
  
  # 构造命令
  make_tag_cmd <- sprintf("%s && makeTagDirectory %s %s -format sam -keepAll", 
                          conda_activate, tag_dir, file)
  find_peaks_cmd <- sprintf("%s && findPeaks %s -style %s -o %s -fragLength %d %s", 
                            conda_activate, tag_dir, peak_style, peak_file, fraglen, tbp_option)
  
  # 仅返回命令不执行
  if (cmd.only) {
    return(list(make_tag_cmd, find_peaks_cmd))
  }
  
  # 执行命令
  system(make_tag_cmd)
  system(find_peaks_cmd)
  
  # 读取 peak 文件并转换成 BED 格式
  if (file.exists(peak_file)) {
    peak_data <- read.table(peak_file, header = FALSE)
    processed_peak_file <- paste0(outprefix, "_processed_peaks.bed")
    write.table(file = processed_peak_file, peak_data[, c(2, 3, 4, 1, 8)],
                row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
    return(processed_peak_file)
  } else {
    stop("Peak calling failed: No output file generated.")
  }
}

run_MACS3 <- function(file, outprefix, is.tf = TRUE,
                      macs.path = "macs3", gsize = "hs",
                      broad_cutoff = 0.1, fraglen = 100, threshold = 0.01,
                      format = "BAM", conda.env = "macs3_env", shell = "zsh",
                      cmd.only = FALSE) {
  # 激活 Conda 环境
  conda_activate <- paste("source ~/.zshrc && conda activate", conda.env)
  
  # 处理 Histone markers（宽峰）
  if (!is.tf) {  
    cmds <- c(macs.path, "callpeak -t", file, "-g", gsize,
              "--outdir", dirname(outprefix), "-n", basename(outprefix),
              "--nomodel", "--extsize", fraglen,
              "--broad", "--broad-cutoff", broad_cutoff)
  } else {  # 处理 TF（窄峰）
    cmds <- c(macs.path, "callpeak -t", file, "-f", format, "-g", gsize,
              "--outdir", dirname(outprefix), "-n", basename(outprefix),
              "--nomodel", "--extsize", fraglen, "-B", "-q", threshold)
  }
  
  # 组装完整的 shell 命令
  full_cmd <- paste(conda_activate, "&&", paste(cmds, collapse = " "))
  print(full_cmd)
  
  # 如果 cmd.only = TRUE，则返回命令而不执行
  if (cmd.only) {
    return(full_cmd)
  }
  
  # 执行 MACS3 命令
  if (system(paste(shell, "-c", shQuote(full_cmd)))) {
    stop("Running MACS3 failed.")
  }
  
  return(invisible(NULL))
}

run_MACS3_server <- function(file, outprefix, is.tf = TRUE,
                      macs.path = "macs3", gsize = "hs",
                      broad_cutoff = 0.1, fraglen = 100, threshold = 0.01,
                      format = "BAM", shell = "bash", cmd.only = FALSE) {
  
  # 组装 MACS3 命令部分
  if (!is.tf) {  # 宽峰 for Histone markers
    cmds <- c(macs.path, "callpeak -t", shQuote(file), "-g", gsize,
              "--outdir", shQuote(dirname(outprefix)), "-n", shQuote(basename(outprefix)),
              "--nomodel", "--extsize", fraglen,
              "--broad", "--broad-cutoff", broad_cutoff)
  } else {  # 窄峰 for TF
    cmds <- c(macs.path, "callpeak -t", shQuote(file), "-f", format, "-g", gsize,
              "--outdir", shQuote(dirname(outprefix)), "-n", shQuote(basename(outprefix)),
              "--nomodel", "--extsize", fraglen, "-B", "-q", threshold)
  }
  
  # 最终命令
  full_cmd <- paste(cmds, collapse = " ")
  cat("Running command:\n", full_cmd, "\n")
  
  # 仅返回命令
  if (cmd.only) {
    return(full_cmd)
  }
  
  # 实际执行
  status <- system(paste(shell, "-c", shQuote(full_cmd)))
  if (status != 0) {
    stop("Running MACS3 failed.")
  }
  
  return(invisible(NULL))
}

########################################################################################
## new fine-tune function 
Finetune_ChIPtest <- function(potential_peaks, count_bin2, width, band_values, quantile_values,
                              my_var.est=1, my_var.thred=0.01, num_cores=5){
  # browser()
  # 数据处理
  data_list_con1 <- mclapply(potential_peaks, process_peak_conA, df = count_bin2[,c(1,2,3)], width = width, mc.cores = num_cores)
  data_con1 <- do.call(rbind, lapply(data_list_con1, t))
  
  data_list_con2 <- mclapply(potential_peaks, process_peak_conB, df = count_bin2[,c(1,2,4)], width = width, mc.cores = num_cores)
  data_con2 <- do.call(rbind, lapply(data_list_con2, t))
  
  data_list_con3 <- mclapply(potential_peaks, process_peak_conA, df = count_bin2[,c(1,2,5)], width = width, mc.cores = num_cores)
  data_con3 <- do.call(rbind, lapply(data_list_con3, t))
  
  data_list_con4 <- mclapply(potential_peaks, process_peak_conB, df = count_bin2[,c(1,2,6)], width = width, mc.cores = num_cores)
  data_con4 <- do.call(rbind, lapply(data_list_con4, t))
  
  # 标准化
  Data_con1 <- NormTransformation(data_con1[,-c(1,2)])
  Data_con2 <- NormTransformation(data_con2[,-c(1,2)])
  Data_con3 <- NormTransformation(data_con3[,-c(1,2)])
  Data_con4 <- NormTransformation(data_con4[,-c(1,2)])
  
  # 估计tao参数
  tao_12 <- my_est.c(Data_con1, Data_con2, max1=4, max4=4)
  tao_34 <- my_est.c(Data_con3, Data_con4, max1=4, max4=4)
  
  # KL散度计算函数 (与 N(0,1) 标准正态分布比较)
  calculate_kl <- function(stat) {
    x_seq <- seq(min(stat), max(stat), length.out = 1000)
    data_density <- density(stat, n = 1000)
    data_prob <- data_density$y / sum(data_density$y)
    normal_prob <- dnorm(x_seq, mean = 0, sd = 1) / sum(dnorm(x_seq, mean = 0, sd = 1))
    entropy::KL.plugin(data_prob, normal_prob)
  }
  
  # 计算统计量和KL散度
  parallel_TS <- function(q, band, Data_con1, Data_con2, Data_con3, Data_con4, tao_12, tao_34) {
    res12 <- my_TS_twosample(Data_con1, Data_con2, tao_12, band = band, 
                             quant = c(q, q, q), var.est = 1, var.thred = 0.01)
    res34 <- my_TS_twosample(Data_con3, Data_con4, tao_34, band = band, 
                             quant = c(q, q, q), var.est = 1, var.thred = 0.01)
    
    c(
      Band = band,
      Quantile = q,
      TS_kn_KL_12 = calculate_kl(res12[["TS_kn"]]),
      Dnun_KL_12 = calculate_kl(res12[["Dnun"]]),
      Deql_KL_12 = calculate_kl(res12[["Deql"]]),
      TS_kn_KL_34 = calculate_kl(res34[["TS_kn"]]),
      Dnun_KL_34 = calculate_kl(res34[["Dnun"]]),
      Deql_KL_34 = calculate_kl(res34[["Deql"]])
    )
  }
  
  # 创建参数组合
  params_list <- expand.grid(Quantile = quantile_values, Band = band_values)
  
  # 并行计算
  results <- mclapply(1:nrow(params_list), function(i) 
    parallel_TS(params_list$Quantile[i], 
                params_list$Band[i],
                Data_con1, Data_con2, Data_con3, Data_con4, tao_12, tao_34), 
    mc.cores = num_cores)
  
  # 汇总结果
  results_df <- as.data.frame(do.call(rbind, results))
  colnames(results_df) <- c("Band", "Quantile", 
                            "TS_kn_KL_12", "Dnun_KL_12", "Deql_KL_12",
                            "TS_kn_KL_34", "Dnun_KL_34", "Deql_KL_34")
  return(results_df)
}

###################################################################################################
# 计算 FDR and Recall 
calculate_fdr_recall <- function(predicted_probs, actual_labels, thresholds, actual_positives = 1000) {
  results <- vector("list", length(thresholds))
  # browser()
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    predicted <- ifelse(predicted_probs <= threshold, 1, 0)
    predicted_positives <- sum(predicted == 1)
    true_positives <- sum(predicted == 1 & actual_labels == 1)
    false_positives <- sum(predicted == 1 & actual_labels == 0)
    
    fdr <- ifelse(predicted_positives == 0, 0, false_positives / predicted_positives)
    recall <- true_positives / actual_positives
    
    results[[i]] <- data.frame(Threshold = threshold, FDR = fdr, Recall = recall)
  }
  return(do.call(rbind, results))
}

#################################################################################

simulate_tf <- function(data_path,
                        fraglen,
                        base.mu,
                        up.mu,
                        down.mu.1,
                        down.mu.2,
                        bins_vis,
                        count_length_vis) {
  # --- Set parameters ---
  # Parameters passed as arguments: data_path, fraglen, up.mu, down.mu.1, down.mu.2
  
  # Internal default parameters (from original script)
  is.tf <- TRUE
  npeaks <- 20000
  nde <- 500
  
  prior.df_val <- 20 # Renamed to avoid conflict if prior.df is used later differently
  dispersion_val <- 1 / prior.df_val
  grouping <- c("A", "A", "B", "B")
  true.width <- 500
  base.mu <- base.mu
  
  radius <- fraglen # Derived from function argument
  prior.df_sim <- 1e8 # This was the prior.df used for the chi-squared simulation
  all.fix <- "tfx" # File prefix, could be an argument too
  
  # --- Data Path and Working Directory ---
  if (!dir.exists(data_path)) {
    dir.create(data_path, recursive = TRUE)
    cat("Directory created:", data_path, "\n")
  } else {
    cat("Directory already exists:", data_path, "\n")
  }
  
  # Save current working directory to restore it later, if desired
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE) # Restore WD when function exits
  
  setwd(data_path)
  cat("Current working directory set to:", getwd(), "\n")
  
  # --- Design Matrix and Parameters ---
  design <- model.matrix(~ factor(grouping))
  # Assuming readParam is from the sourced ChIPtest_source.R
  # It might require specific columns in 'param_template.txt' or similar
  # For simplicity, ensure dedup is available or handle its absence.
  if (exists("readParam", mode = "function")) {
    xparam <- readParam(dedup = FALSE) 
  } else {
    stop("Function 'readParam' not found. Please ensure ChIPtest_source.R is sourced correctly and defines it.")
  }
  
  
  # --- Defining Peak Sets and Dispersion ---
  up.pk <- 1:nde
  down.pk <- (1:nde) + nde
  # Using prior.df_sim for dispersion calculation as in the original script
  disp <- prior.df_sim * (1/20) / rchisq(npeaks, df = prior.df_sim) # original used dispersion (1/prior.df_val)
  # which means prior.df_sim * (1/20) / rchisq(...)
  # this might need clarification if prior.df_val (20) or prior.df_sim (1e8)
  # was intended for the initial 'dispersion' factor.
  # Original: prior.df*dispersion/rchisq(...)
  # with prior.df = 1e8 and dispersion = 1/20.
  # So, disp <- prior.df_sim * dispersion_val / rchisq(npeaks, df=prior.df_sim)
  disp <- prior.df_sim * dispersion_val / rchisq(npeaks, df = prior.df_sim)
  
  
  if (!is.tf) {
    # This block is skipped as is.tf is TRUE, but kept for completeness if is.tf becomes a parameter
    type.A.1 <- type.A.2 <- type.A.3 <- !logical(npeaks)
    chosen.drop.A <- sample(1:6, nde, replace = TRUE)
    type.A.1[up.pk] <- bitwAnd(chosen.drop.A, 0x1) > 0L
    type.A.2[up.pk] <- bitwAnd(chosen.drop.A, 0x2) > 0L
    type.A.3[up.pk] <- bitwAnd(chosen.drop.A, 0x4) > 0L
    
    type.B.1 <- type.B.2 <- type.B.3 <- !logical(npeaks)
    chosen.drop.B <- sample(1:6, nde, replace = TRUE)
    type.B.1[down.pk] <- bitwAnd(chosen.drop.B, 0x1) > 0L
    type.B.2[down.pk] <- bitwAnd(chosen.drop.B, 0x2) > 0L
    type.B.3[down.pk] <- bitwAnd(chosen.drop.B, 0x4) > 0L
  }
  
  # --- Peak Positioning ---
  distances <- round(runif(npeaks, 10000, 20000))
  pos.1 <- cumsum(distances)
  if (!is.tf) {
    pos.2 <- pos.1 + true.width / 2
    pos.3 <- pos.1 + true.width
    sizes <- c(chrA = max(pos.3) + 10000)
  } else {
    sizes <- c(chrA = max(pos.1) + 10000)
  }
  chrs <- rep("chrA", npeaks)
  
  # --- Simulating Reads and Generating SAM Files ---
  fnames <- list()
  for (lib in 1:length(grouping)) {
    fname <- paste0(all.fix, "_out_", lib, ".sam")
    if (!is.tf) {
      # Histone mark simulation (skipped as is.tf is TRUE)
      if (grouping[lib] == "A") {
        drop.1 <- type.A.1
        drop.2 <- type.A.2
        drop.3 <- type.A.3
      } else {
        drop.1 <- type.B.1
        drop.2 <- type.B.2
        drop.3 <- type.B.3
      }
      peakFile(fname, chrs = chrs[drop.1], pos = pos.1[drop.1], mu = base.mu, disp = disp[drop.1],
               sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE)
      peakFile(fname, chrs = chrs[drop.2], pos = pos.2[drop.2], mu = base.mu, disp = disp[drop.2],
               sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE, append = TRUE)
      peakFile(fname, chrs = chrs[drop.3], pos = pos.3[drop.3], mu = base.mu, disp = disp[drop.3],
               sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE, append = TRUE)
    } else {
      # TF simulation (using function arguments for mu values)
      cur.mu.counts <- rep(base.mu, npeaks) # Renamed to avoid conflict with arg up.mu
      if (grouping[lib] == "A") {
        cur.mu.counts[down.pk] <- down.mu.1 # Function argument
        cur.mu.counts[up.pk] <- up.mu     # Function argument
      } else {
        cur.mu.counts[up.pk] <- down.mu.2   # Function argument
        cur.mu.counts[down.pk] <- up.mu     # Function argument
      }
      # Assuming peakFile is from the sourced ChIPtest_source.R
      if (!exists("peakFile", mode = "function")) {
        stop("Function 'peakFile' not found. Please ensure ChIPtest_source.R is sourced correctly.")
      }
      peakFile(fname, chrs = chrs, pos = pos.1, mu = cur.mu.counts, disp = disp,
               sizes = sizes, fraglen = fraglen, width = true.width, tf = TRUE)
    }
    fnames[[lib]] <- fname
  }
  
  fnames <- unlist(fnames)
  # Assuming addBackground and crunch2BAM are from the sourced ChIPtest_source.R
  if (!exists("addBackground", mode = "function")) {
    stop("Function 'addBackground' not found. Please ensure ChIPtest_source.R is sourced correctly.")
  }
  if (!exists("crunch2BAM", mode = "function")) {
    stop("Function 'crunch2BAM' not found. Please ensure ChIPtest_source.R is sourced correctly.")
  }
  
  # Using prior.df_val for addBackground as in original dispersion parameter context
  addBackground(fnames, sizes = sizes, width = 2000, rlen = 10,
                dispersion = dispersion_val, prior.df = prior.df_val, append = TRUE)
  bam.files <- crunch2BAM(fnames)
  unlink(fnames) # Remove intermediate SAM files
  
  # --- Log File Generation ---
  lfile <- paste0(all.fix, "_log.txt")
  write.table(file = lfile, data.frame(chr = chrs[up.pk], start = pos.1[up.pk] - radius, end = pos.1[up.pk] + radius, logFC = 1),
              row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(file = lfile, data.frame(chrs[down.pk], pos.1[down.pk] - radius, pos.1[down.pk] + radius, logFC = -1),
              row.names = FALSE, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
  
  pos1_rdata_file <- "pos1_vector.RData"
  save(pos.1, file = pos1_rdata_file)
  cat("Peak positions saved to:", file.path(getwd(), pos1_rdata_file), "\n")
  
  
  # --- Visualization Section ---
  cat("Starting visualization...\n")
  bins_vis <- bins_vis # Renamed from 'bins' to avoid potential conflicts
  count_length_vis <- count_length_vis
  
  # Assuming windowCounts is from sourced ChIPtest_source.R or csaw package
  if (!exists("windowCounts", mode = "function")) {
    stop("Function 'windowCounts' not found. Please ensure ChIPtest_source.R (or csaw package) is available.")
  }
  count_seq <- windowCounts(bam.files, width = bins_vis, spacing = bins_vis,
                            ext = bins_vis, filter = 0, bin = TRUE, param = xparam)
  
  counts_matrix <- assay(count_seq)
  start_positions <- start(rowRanges(count_seq))
  end_positions <- end(rowRanges(count_seq))
  count_bin2 <- data.frame(start = start_positions, end = end_positions, counts_matrix)
  
  # Dynamically get column names for replicates, assuming bam.files provides base names like X_out_1, X_out_2 etc.
  # If crunch2BAM returns names that correspond to 1,2,3,4 directly in assay, then X1,X2,X3,X4 is fine.
  # For robustness, it's better to explicitly use generated column names if they are predictable.
  # The original script uses X1, X2, X3, X4, so we'll stick to that assumption.
  # Make sure the assay(count_seq) columns are named "X1", "X2", "X3", "X4" or adjust accordingly.
  # If bam.files are e.g. "tfx_out_1.bam", "tfx_out_2.bam", colnames might be derived from these.
  # Let's assume standard csaw behavior gives colnames based on bam.files input.
  # For simplicity, if colnames are like "tfx_out_1.bam", etc., this will fail.
  # The original used X1, X2, X3, X4. This happens if no colnames are set on SummarizedExperiment.
  # Let's add a check or make it more robust if possible.
  if(!all(c("X1", "X2", "X3", "X4") %in% colnames(count_bin2))) {
    warning("Default column names X1, X2, X3, X4 not found in count_bin2. Averaging might fail or be incorrect. Actual colnames: ", paste(colnames(counts_matrix), collapse=", "))
    # Fallback: assuming the first four columns after start/end are the counts in order
    count_bin2$conA <- (count_bin2[,3] + count_bin2[,4]) / 2 
    count_bin2$conB <- (count_bin2[,5] + count_bin2[,6]) / 2
  } else {
    count_bin2$conA <- (count_bin2$X1 + count_bin2$X2) / 2
    count_bin2$conB <- (count_bin2$X3 + count_bin2$X4) / 2
  }
  
  
  hist_log <- read.table(lfile, header = TRUE, sep = "\t")
  
  # set.seed(123) # Added for reproducibility of sampling
  indices_right <- sample(1:nde, 5) # Sample from up-regulated or down-regulated peaks
  my_location_right <- pos.1[indices_right]
  
  # Sample from non-DE peaks (e.g., after 2*nde)
  # Ensure npeaks is large enough for this sampling range
  if (npeaks > (2 * nde + 5)) {
    indices_wrong <- sample((2 * nde + 1):min(npeaks, (2*nde + 1000)), 5) # sample from a pool of 1000 non-DE sites
  } else { # Fallback if not many non-DE peaks
    indices_wrong <- sample( (nde+1) : (2*nde), 5) # sample from other DE type as a placeholder for 'wrong'
    warning("Not enough non-DE peaks for 'wrong' samples based on original logic; sampling from other DE types.")
  }
  my_location_wrong <- pos.1[indices_wrong]
  
  
  my_peaks_bins_right <- round(my_location_right / bins_vis)
  my_peaks_bins_wrong <- round(my_location_wrong / bins_vis)
  
  # Assuming process_peak is from the sourced ChIPtest_source.R
  if (!exists("process_peak", mode = "function")) {
    stop("Function 'process_peak' not found. Please ensure ChIPtest_source.R is sourced correctly.")
  }
  
  # Check for mclapply and set mc.cores accordingly
  # mclapply is not available on Windows by default.
  num_cores <- 1
  if (requireNamespace("parallel", quietly = TRUE) && .Platform$OS.type != "windows") {
    num_cores <- parallel::detectCores()
    if (num_cores > 5) num_cores <- 5 # Cap at 5 as in original script for some calls
  } else if (.Platform$OS.type == "windows") {
    cat("mclapply not available on Windows, using lapply (sequential processing).\n")
  }
  
  # Using available cores or lapply for Windows
  apply_processing <- if (num_cores > 1 && .Platform$OS.type != "windows") {
    function(...) parallel::mclapply(..., mc.cores = num_cores)
  } else {
    lapply
  }
  
  # Selecting columns for conA (e.g., "start", "end", "conA")
  # The original script used indices like c(1,2,7). This depends on the structure of count_bin2.
  # If conA is column 7, conB is column 8 (after start, end, X1, X2, X3, X4)
  col_idx_conA <- which(colnames(count_bin2) == "conA")
  col_idx_conB <- which(colnames(count_bin2) == "conB")
  
  data_list_conA_right <- apply_processing(my_peaks_bins_right, process_peak, df = count_bin2[, c(1, 2, col_idx_conA)], width = count_length_vis)
  data_conA_right <- do.call(rbind, lapply(data_list_conA_right, t))
  data_list_conB_right <- apply_processing(my_peaks_bins_right, process_peak, df = count_bin2[, c(1, 2, col_idx_conB)], width = count_length_vis)
  data_conB_right <- do.call(rbind, lapply(data_list_conB_right, t))
  
  data_list_conA_wrong <- apply_processing(my_peaks_bins_wrong, process_peak, df = count_bin2[, c(1, 2, col_idx_conA)], width = count_length_vis)
  data_conA_wrong <- do.call(rbind, lapply(data_list_conA_wrong, t))
  data_list_conB_wrong <- apply_processing(my_peaks_bins_wrong, process_peak, df = count_bin2[, c(1, 2, col_idx_conB)], width = count_length_vis)
  data_conB_wrong <- do.call(rbind, lapply(data_list_conB_wrong, t))
  
  
  # Check if data frames are not empty before proceeding
  if (nrow(data_conA_right) == 0 || ncol(data_conA_right) <= 2) {
    stop("Data processing for 'right' peaks resulted in empty or insufficient data for conA.")
  }
  if (nrow(data_conB_right) == 0 || ncol(data_conB_right) <= 2) {
    stop("Data processing for 'right' peaks resulted in empty or insufficient data for conB.")
  }
  if (nrow(data_conA_wrong) == 0 || ncol(data_conA_wrong) <= 2) {
    stop("Data processing for 'wrong' peaks resulted in empty or insufficient data for conA.")
  }
  if (nrow(data_conB_wrong) == 0 || ncol(data_conB_wrong) <= 2) {
    stop("Data processing for 'wrong' peaks resulted in empty or insufficient data for conB.")
  }
  
  
  dfA_squence_right <- data_conA_right[, -c(1, 2), drop = FALSE]
  dfB_squence_right <- data_conB_right[, -c(1, 2), drop = FALSE]
  
  # Ensure there are columns to work with after removing first two
  if (ncol(dfA_squence_right) != count_length_vis +1 && ncol(dfA_squence_right) != count_length_vis){ # process_peak might return count_length or count_length+1
    warning(paste("Unexpected number of columns in dfA_squence_right:", ncol(dfA_squence_right), "Expected around:", count_length_vis))
    # Adjust Index length for data.frame creation if necessary, or error out
    # For now, assume it's count_length_vis or count_length_vis+1
    # Original code assumes 121 columns (Index 1:121 means 121 values)
    # This means process_peak likely returns count_length_vis items (e.g. 120), and the data.frame adds an Index column.
    # The data.frame creation `Index = 1:121` and `t(dfA_sequence_right)` means dfA_sequence_right should have 121 columns.
    # If count_length_vis = 120, then process_peak should return 121 data points.
    # Let's make Index dynamic to the actual number of columns from process_peak.
  }
  
  actual_seq_length_A_right = ncol(dfA_squence_right)
  actual_seq_length_B_right = ncol(dfB_squence_right)
  actual_seq_length_A_wrong = ncol(data_conA_wrong[,-c(1,2), drop=FALSE])
  actual_seq_length_B_wrong = ncol(data_conB_wrong[,-c(1,2), drop=FALSE])
  
  dfA_long_right <- data.frame(Index = 1:actual_seq_length_A_right, t(dfA_squence_right))
  colnames(dfA_long_right) <- c("Index", paste("TF_Right", seq_len(nrow(dfA_squence_right)), sep = "_"))
  dfB_long_right <- data.frame(Index = 1:actual_seq_length_B_right, t(dfB_squence_right))
  colnames(dfB_long_right) <- c("Index", paste("TF_Right", seq_len(nrow(dfB_squence_right)), sep = "_")) # Should be same series names for faceting
  
  dfA_long_right <- data.table::melt(data.table::as.data.table(dfA_long_right), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfB_long_right <- data.table::melt(data.table::as.data.table(dfB_long_right), id.vars = "Index", variable.name = "Series", value.name = "Value")
  
  dfA_long_right$Group <- "condition1"
  dfB_long_right$Group <- "condition2"
  
  dfA_squence_wrong <- data_conA_wrong[, -c(1, 2), drop = FALSE]
  dfB_squence_wrong <- data_conB_wrong[, -c(1, 2), drop = FALSE]
  
  dfA_long_wrong <- data.frame(Index = 1:actual_seq_length_A_wrong, t(dfA_squence_wrong))
  colnames(dfA_long_wrong) <- c("Index", paste("TF_Wrong", seq_len(nrow(dfA_squence_wrong)), sep = "_"))
  dfB_long_wrong <- data.frame(Index = 1:actual_seq_length_B_wrong, t(dfB_squence_wrong))
  colnames(dfB_long_wrong) <- c("Index", paste("TF_Wrong", seq_len(nrow(dfB_squence_wrong)), sep = "_"))
  
  dfA_long_wrong <- data.table::melt(data.table::as.data.table(dfA_long_wrong), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfB_long_wrong <- data.table::melt(data.table::as.data.table(dfB_long_wrong), id.vars = "Index", variable.name = "Series", value.name = "Value")
  
  dfA_long_wrong$Group <- "condition1"
  dfB_long_wrong$Group <- "condition2"
  
  df_combined <- rbind(dfA_long_right, dfB_long_right, dfA_long_wrong, dfB_long_wrong)
  
  # Ensure Series levels are correctly ordered for plotting
  series_levels <- c(
    grep("TF_Right", unique(as.character(df_combined$Series)), value = TRUE),
    grep("TF_Wrong", unique(as.character(df_combined$Series)), value = TRUE)
  )
  df_combined$Series <- factor(df_combined$Series, levels = series_levels)
  
  
  plot_object <- ggplot2::ggplot(df_combined, ggplot2::aes(x = Index, y = Value, color = Group)) +
    ggplot2::geom_line(alpha = 0.99, linewidth = 0.5) + # size deprecated, use linewidth
    ggplot2::facet_wrap(~ Series, nrow = 2, scales = "free_y") + # Added scales = "free_y" for better individual plot visibility
    ggplot2::labs(title = "Line Plot of Simulated TF ChIP-seq Counts",
                  x = "Relative Bin Index from Peak Center", y = "Average Binned Counts", color = "Condition") +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      strip.text = ggplot2::element_text(size = 10, face = "bold"), # Reduced size a bit
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1) # size deprecated
    )
  
  plot_filename <- file.path(getwd(), paste0(all.fix, "_visualization_plot.png"))
  ggplot2::ggsave(plot_filename, plot = plot_object, width = 12, height = 6, dpi = 300,bg = "white")
  cat("Visualization plot saved to:", plot_filename, "\n")
  
  cat("Simulation and visualization complete.\n")
  
  # Restore original working directory (done by on.exit)
  # setwd(original_wd) 
  
  return(list(
    log_file = file.path(getwd(), lfile),
    bam_files = file.path(getwd(), bam.files), # getwd() here is data_path
    plot_file = plot_filename,
    data_directory = getwd() 
  ))
}

# simulate_tf(data_path="/Users/cuitengfei/Graduate/Research/Cancer/SCAW_3rd/2_Server_script/Data/practice_tf_1",
#             fraglen=1000, up.mu=270, down.mu.1=90,down.mu.2=90,
#             bins_vis=25, count_length_vis=120)

simulate_hm <- function(true.width,
                        base.mu,
                        data_path,
                        bins,
                        count_length) {
  
  # --- Helper function for visualization (defined inside as per original script) ---
  process_peak_local <- function(x, df, width) { # Renamed to avoid potential global conflicts
    # Ensure x is a valid index for df
    if (x < 1 || x > nrow(df)) {
      stop(paste("Invalid peak center bin index:", x, "for data frame with", nrow(df), "rows."))
    }
    start_row <- max(1, x - (width / 2)) # Ensure start_row is not less than 1
    end_row <- min(nrow(df), x + (width / 2)) # Ensure end_row does not exceed df rows
    
    # Ensure start_row <= end_row
    if (start_row > end_row) {
      warning(paste("Calculated start_row > end_row for peak bin index:", x, ". Adjusting to a single row."))
      end_row <- start_row # Or handle as an error / return NA sequence
    }
    
    # Check if columns exist
    if (ncol(df) < 3) stop("DataFrame for process_peak_local needs at least 3 columns.")
    
    # Extract sequence; handle cases where the window is truncated by df boundaries
    sequence_values <- df[start_row:end_row, 3]
    
    # Pad if the window is smaller than width + 1
    # process_peak from original returns width+1 elements, so padding might be needed
    # The original script uses Index = 1:(count_length+1), implying (count_length+1) values.
    # If width = count_length, then we need count_length + 1 values.
    
    # Original: return(c(df[start_row, 1], df[end_row, 2], df[start_row:end_row, 3]))
    # This structure returns start_coord, end_coord, and then the sequence.
    # The subsequent code data_conA_right[,-c(1,2)] assumes the first two are not part of the sequence.
    # The length of df[start_row:end_row, 3] needs to be count_length + 1.
    # The current definition of start_row/end_row creates a window of (width+1) bins if x is central.
    
    actual_sequence <- df[start_row:end_row, 3]
    
    # The original `process_peak` returns `width+1` elements for the sequence if width is `count_length`.
    # The logic here implies the sequence extracted `df[start_row:end_row, 3]` should be of length `count_length + 1`.
    # This means `end_row - start_row + 1` should be `count_length + 1`.
    # The window is centered at `x`, with `width/2` on each side, so total `width+1` elements.
    # If `width` argument to `process_peak_local` is `count_length`, then it's correct.
    
    return(c(df[start_row, 1], df[end_row, 2], actual_sequence))
  }
  
  # --- Set parameters ---
  # Parameters passed as arguments: true.width, base.mu, data_path, bins, count_length
  
  # Internal default parameters
  is.tf <- FALSE # Fixed for histone mark simulation
  
  fraglen <- 100
  npeaks <- 20000
  nde <- 500
  rlen <- 10 # Read length, used in peakFile and addBackground
  
  prior.df_val <- 20
  dispersion_val <- 1 / prior.df_val
  grouping <- c("A", "A", "B", "B")
  # base.mu is now an argument
  # true.width is now an argument
  
  # Parameters derived from arguments or fixed for !is.tf
  all.fix <- "hist"
  radius <- true.width / 2L # true.width is an argument
  
  # --- Data Path and Working Directory ---
  if (!dir.exists(data_path)) {
    dir.create(data_path, recursive = TRUE)
    cat("Directory created:", data_path, "\n")
  } else {
    cat("Directory already exists:", data_path, "\n")
  }
  
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  
  setwd(data_path)
  cat("Current working directory set to:", getwd(), "\n")
  
  # --- csaw Parameters ---
  # Assuming readParam is from csaw package and loaded externally
  if (!requireNamespace("csaw", quietly = TRUE)) {
    stop("Package 'csaw' is required but not installed/loaded.")
  }
  xparam <- csaw::readParam(dedup = FALSE)
  
  
  # --- Generating simulated data for histone mark data ---
  up.pk <- 1:nde
  down.pk <- (1:nde) + nde
  # For dispersion, the original script used prior.df=20.
  # If a different prior.df was intended for rchisq (like the 1e8 in TF script), it needs clarification.
  # Assuming prior.df_val (20) is the one for rchisq here as well.
  disp <- prior.df_val * dispersion_val / rchisq(npeaks, df = prior.df_val)
  
  
  # This block is now the primary logic path as is.tf is FALSE
  type.A.1 <- type.A.2 <- type.A.3 <- !logical(npeaks)
  # set.seed(123) # For reproducibility of sampling differential types
  chosen.drop.A <- sample(1:6, nde, replace = TRUE) # 1,2,4 are bitmasks for dropping peak1, peak2, peak3
  type.A.1[up.pk] <- bitwAnd(chosen.drop.A, 0x1) > 0L # If bit 1 is set, peak component 1 is present for group A 'up'
  type.A.2[up.pk] <- bitwAnd(chosen.drop.A, 0x2) > 0L # If bit 2 is set, peak component 2 is present
  type.A.3[up.pk] <- bitwAnd(chosen.drop.A, 0x4) > 0L # If bit 4 is set, peak component 3 is present
  # So, up.pk for A means some components are PRESENT. If all are present, it's "fully up".
  # The log file logic seems to define DE where components are MISSING. This needs careful alignment.
  # Let's assume type.A.x = TRUE means component x IS PRESENT for up.pk in group A.
  # For group A, up.pk are those that should be "stronger" or more complete.
  # By default, all type.A.x are TRUE (all components present).
  # For up.pk, we make some components present based on chosen.drop.A.
  # This means for group A, at up.pk regions, the state is determined by chosen.drop.A.
  # If chosen.drop.A makes a component FALSE, it means it's absent.
  
  type.B.1 <- type.B.2 <- type.B.3 <- !logical(npeaks)
  chosen.drop.B <- sample(1:6, nde, replace = TRUE)
  type.B.1[down.pk] <- bitwAnd(chosen.drop.B, 0x1) > 0L
  type.B.2[down.pk] <- bitwAnd(chosen.drop.B, 0x2) > 0L
  type.B.3[down.pk] <- bitwAnd(chosen.drop.B, 0x4) > 0L
  # For group B, down.pk are those that should be "stronger" or more complete.
  # Non up.pk/down.pk regions will have all components present in both groups.
  
  # Clarification on differential logic for !is.tf:
  # - up.pk: nde regions where Group A's profile is modified by chosen.drop.A. Group B has all components.
  # - down.pk: nde regions where Group B's profile is modified by chosen.drop.B. Group A has all components.
  # The log file writes !type.A.x as logFC=-1 (lost in A) and !type.B.x as logFC=1 (gained in A, i.e. lost in B).
  # This implies:
  # For up.pk: these are regions where A is "up". We want to simulate A having more components than B, or B losing components.
  # The code sets type.A.x[up.pk] according to chosen.drop.A.
  # And type.B.x[down.pk] according to chosen.drop.B.
  # Let's re-evaluate the simulation loop's use of drop.1/2/3.
  
  # Original logic for histone:
  # if (grouping[lib]=="A") { drop.1 <- type.A.1 ...} else {drop.1 <- type.B.1 ...}
  # This means `type.A.1` defines presence/absence for group A for *all* peaks.
  # And `type.A.1[up.pk]` modifies this for the `up.pk` regions in group A.
  # So, for group A:
  #   - at up.pk indices: type.A.x[up.pk] determines component presence.
  #   - at other indices: type.A.x is !logical(npeaks) which is TRUE (all components present).
  # For group B:
  #   - at down.pk indices: type.B.x[down.pk] determines component presence.
  #   - at other indices: type.B.x is TRUE (all components present).
  # This seems to mean that `up.pk` are regions specifically altered in A, and `down.pk` are regions specifically altered in B.
  # And "alteration" means some components (determined by chosen.drop) might be absent.
  # If `chosen.drop.A` results in all components being present (e.g. chosen.drop.A = 7), then A has full signal at up.pk.
  # The log file seems to mark where components are *missing*.
  # A logFC = -1 means "down in A compared to B". This corresponds to !type.A.x.
  # A logFC = 1 means "up in A compared to B". This corresponds to !type.B.x (i.e., B is missing a component that A has).
  
  
  # --- Peak Positioning ---
  distances <- round(runif(npeaks, 10000, 20000))
  pos.1 <- cumsum(distances)
  pos.2 <- pos.1 + true.width / 2 # true.width is an argument
  pos.3 <- pos.1 + true.width   # true.width is an argument
  sizes <- c(chrA = max(pos.3) + 10000)
  chrs <- rep("chrA", npeaks)
  
  # --- Simulating Reads and Generating SAM Files ---
  fnames <- list()
  for (lib in 1:length(grouping)) {
    fname <- paste0(all.fix, "_out_", lib, ".sam")
    # This is the histone mark simulation path
    current_type.1 <- !logical(npeaks)
    current_type.2 <- !logical(npeaks)
    current_type.3 <- !logical(npeaks)
    
    if (grouping[lib] == "A") {
      current_type.1[up.pk] <- type.A.1[up.pk] # Apply specific component presence/absence for group A at up.pk
      current_type.2[up.pk] <- type.A.2[up.pk]
      current_type.3[up.pk] <- type.A.3[up.pk]
    } else { # Group B
      current_type.1[down.pk] <- type.B.1[down.pk] # Apply specific component presence/absence for group B at down.pk
      current_type.2[down.pk] <- type.B.2[down.pk]
      current_type.3[down.pk] <- type.B.3[down.pk]
    }
    
    # peakFile function needs to be available (e.g. from externally sourced script)
    if (!exists("peakFile", mode = "function")) {
      stop("Function 'peakFile' not found. Please ensure it's defined or sourced.")
    }
    
    # Call peakFile only for PRESENT components
    peakFile(fname, chrs = chrs[current_type.1], pos = pos.1[current_type.1], mu = base.mu, disp = disp[current_type.1], rlen = rlen,
             sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE)
    peakFile(fname, chrs = chrs[current_type.2], pos = pos.2[current_type.2], mu = base.mu, disp = disp[current_type.2], rlen = rlen,
             sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE, append = TRUE)
    peakFile(fname, chrs = chrs[current_type.3], pos = pos.3[current_type.3], mu = base.mu, disp = disp[current_type.3], rlen = rlen,
             sizes = sizes, fraglen = fraglen, width = true.width, tf = FALSE, append = TRUE)
    fnames[[lib]] <- fname
  }
  
  fnames <- unlist(fnames)
  
  if (!exists("addBackground", mode = "function") || !exists("crunch2BAM", mode = "function")) {
    stop("Functions 'addBackground' and/or 'crunch2BAM' not found. Please ensure they are defined or sourced.")
  }
  addBackground(fnames, sizes = sizes, width = 2000, rlen = rlen,
                dispersion = dispersion_val, prior.df = prior.df_val, append = TRUE)
  bam.files <- crunch2BAM(fnames)
  unlink(fnames)
  
  # --- Log File Generation ---
  # Log file indicates where components are ABSENT.
  # !type.A.x means component x is absent in group A at up.pk regions. LogFC = -1 (A is down).
  # !type.B.x means component x is absent in group B at down.pk regions. LogFC = 1 (A is up relative to B).
  lfile <- paste0(all.fix, "_log.txt")
  first_write <- TRUE
  
  # Regions where A is "down" (missing components at up.pk regions)
  # up.pk are the INDICES of peaks that are designated as "up-regulated in A" sites.
  # We need to check components specifically within these up.pk.
  # type.A.1[up.pk] gives the status of component 1 for these up.pk regions.
  # So, up.pk[!type.A.1[up.pk]] would be the indices within up.pk where component 1 is missing for group A.
  # The original log: which(!type.A.1) - this is over ALL peaks.
  # This implies type.A.1, type.A.2, type.A.3 are defined globally for A, and type.B for B.
  # Let's re-check the definition of type.A.x and type.B.x from original script:
  # type.A.1 <- type.A.2 <- type.A.3 <- !logical(npeaks)  (all TRUE initially)
  # type.A.1[up.pk] <- bitwAnd(chosen.drop.A, 0x1) > 0L (modified for up.pk based on chosen.drop.A)
  # This means type.A.x variables store the final state of components for group A across all peaks.
  # Similarly for type.B.x for group B.
  
  # The simulation loop used a different interpretation:
  # It created current_type.1/2/3 based on group and modified them only at up.pk for A, or down.pk for B.
  # This means outside up.pk/down.pk, current_type.x was TRUE for both groups.
  # This implies that non-DE regions have all 3 components for both groups.
  # And DE happens at up.pk (for A) and down.pk (for B).
  
  # Let's define the ground truth based on what was simulated.
  # Group A components: A1, A2, A3. Group B components: B1, B2, B3.
  # A1_final = !logical(npeaks); A1_final[up.pk] = type.A.1[up.pk] from chosen.drop.A
  # B1_final = !logical(npeaks); B1_final[down.pk] = type.B.1[down.pk] from chosen.drop.B
  # Similar for A2,A3 and B2,B3.
  
  # True states after specific assignments:
  true_A1 <- !logical(npeaks); true_A1[up.pk] <- type.A.1[up.pk]
  true_A2 <- !logical(npeaks); true_A2[up.pk] <- type.A.2[up.pk]
  true_A3 <- !logical(npeaks); true_A3[up.pk] <- type.A.3[up.pk]
  
  true_B1 <- !logical(npeaks); true_B1[down.pk] <- type.B.1[down.pk]
  true_B2 <- !logical(npeaks); true_B2[down.pk] <- type.B.2[down.pk]
  true_B3 <- !logical(npeaks); true_B3[down.pk] <- type.B.3[down.pk]
  
  # LogFC = -1 means A is down compared to B (A missing, B has)
  # LogFC = +1 means A is up compared to B (A has, B missing)
  
  # Component 1
  de_c1_Adown <- up.pk[!true_A1[up.pk] & true_B1[up.pk]] # A missing c1 at up.pk, B has c1 at up.pk
  de_c1_Aup   <- up.pk[true_A1[up.pk] & !true_B1[up.pk]] # A has c1 at up.pk, B missing c1 at up.pk
  # Add down.pk regions too for comprehensive DE
  de_c1_Adown <- union(de_c1_Adown, down.pk[!true_A1[down.pk] & true_B1[down.pk]])
  de_c1_Aup   <- union(de_c1_Aup,   down.pk[true_A1[down.pk] & !true_B1[down.pk]])
  
  # Original log logic was simpler: write.table(!type.A.1 ... logFC=-1), write.table(!type.B.1 ... logFC=1)
  # This implies log entries are for components missing in A (relative to full) or missing in B (relative to full).
  # Let's stick to the original log writing logic first, assuming 'type.A.x' and 'type.B.x' correctly define the *overall* component presence for each group.
  # The original `type.A.1[up.pk]` etc. defines the state of components for group A at `up.pk` regions. Outside `up.pk`, components are all present for A.
  # Similar for B at `down.pk` regions.
  # So, !type.A.1 means component 1 is missing for A. This happens only at up.pk if chosen.drop.A makes it so.
  # And !type.B.1 means component 1 is missing for B. This happens only at down.pk if chosen.drop.B makes it so.
  
  # Write based on actual simulated states (current_type.X used in loop)
  # The log should reflect regions that are differentially bound.
  # Example: A peak at pos.1[i]
  #   - if component 1 is present for Group A but absent for Group B --> logFC = 1
  #   - if component 1 is absent for Group A but present for Group B --> logFC = -1
  
  # Let's use the peak indices `1:npeaks`
  all_peak_indices <- 1:npeaks
  
  # DE for component 1 at pos.1
  c1_A_present <- true_A1
  c1_B_present <- true_B1
  log_c1_A_up   <- all_peak_indices[c1_A_present & !c1_B_present]
  log_c1_A_down <- all_peak_indices[!c1_A_present & c1_B_present]
  
  if (length(log_c1_A_down) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c1_A_down], start=pos.1[log_c1_A_down]-radius, end=pos.1[log_c1_A_down]+radius, logFC=-1, name=log_c1_A_down), 
                row.names=FALSE, sep="\t", quote=FALSE, col.names = first_write)
    if(first_write) first_write <- FALSE
  }
  if (length(log_c1_A_up) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c1_A_up], start=pos.1[log_c1_A_up]-radius, end=pos.1[log_c1_A_up]+radius, logFC=1, name=log_c1_A_up), 
                row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=first_write)
    if(first_write && length(log_c1_A_up)>0) first_write <- FALSE # Only set to FALSE if something was written
  }
  
  # DE for component 2 at pos.2
  c2_A_present <- true_A2
  c2_B_present <- true_B2
  log_c2_A_up   <- all_peak_indices[c2_A_present & !c2_B_present]
  log_c2_A_down <- all_peak_indices[!c2_A_present & c2_B_present]
  
  if (length(log_c2_A_down) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c2_A_down], start=pos.2[log_c2_A_down]-radius, end=pos.2[log_c2_A_down]+radius, logFC=-1, name=log_c2_A_down), 
                row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=first_write)
    if(first_write && length(log_c2_A_down)>0) first_write <- FALSE
  }
  if (length(log_c2_A_up) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c2_A_up], start=pos.2[log_c2_A_up]-radius, end=pos.2[log_c2_A_up]+radius, logFC=1, name=log_c2_A_up), 
                row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=first_write)
    if(first_write && length(log_c2_A_up)>0) first_write <- FALSE
  }
  
  # DE for component 3 at pos.3
  c3_A_present <- true_A3
  c3_B_present <- true_B3
  log_c3_A_up   <- all_peak_indices[c3_A_present & !c3_B_present]
  log_c3_A_down <- all_peak_indices[!c3_A_present & c3_B_present]
  
  if (length(log_c3_A_down) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c3_A_down], start=pos.3[log_c3_A_down]-radius, end=pos.3[log_c3_A_down]+radius, logFC=-1, name=log_c3_A_down), 
                row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=first_write)
    if(first_write && length(log_c3_A_down)>0) first_write <- FALSE
  }
  if (length(log_c3_A_up) > 0) {
    write.table(file=lfile, data.frame(chr=chrs[log_c3_A_up], start=pos.3[log_c3_A_up]-radius, end=pos.3[log_c3_A_up]+radius, logFC=1, name=log_c3_A_up), 
                row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=first_write)
    # if(first_write && length(log_c3_A_up)>0) first_write <- FALSE # No more writes after this
  }
  
  save(pos.1, file = "pos1_vector.RData")
  save(pos.2, file = "pos2_vector.RData")
  save(pos.3, file = "pos3_vector.RData")
  cat("Peak component positions saved.\n")
  
  # --- Visualization Section ---
  cat("Starting visualization...\n")
  # bins and count_length are now function arguments.
  
  count_seq <- csaw::windowCounts(bam.files, width = bins, spacing = bins,
                                  ext = bins, filter = 0, bin = TRUE, param = xparam)
  
  counts_matrix <- SummarizedExperiment::assay(count_seq)
  row_ranges <- SummarizedExperiment::rowRanges(count_seq)
  start_positions <- GenomicRanges::start(row_ranges)
  end_positions <- GenomicRanges::end(row_ranges)
  count_bin2 <- data.frame(start = start_positions, end = end_positions, counts_matrix)
  
  if (!all(c("X1", "X2", "X3", "X4") %in% colnames(count_bin2))) {
    warning("Default column names X1, X2, X3, X4 not found in count_bin2. Averaging might fail or be incorrect. Actual colnames: ", paste(colnames(counts_matrix), collapse=", "))
    count_bin2$conA <- (count_bin2[, 3] + count_bin2[, 4]) / 2
    count_bin2$conB <- (count_bin2[, 5] + count_bin2[, 6]) / 2
  } else {
    count_bin2$conA <- (count_bin2$X1 + count_bin2$X2) / 2
    count_bin2$conB <- (count_bin2$X3 + count_bin2$X4) / 2
  }
  
  if (file.exists(lfile) && file.info(lfile)$size > 0) {
    hist_log <- try(read.table(lfile, header = TRUE, sep = "\t"), silent = TRUE)
    if (inherits(hist_log, "try-error") || nrow(hist_log) == 0) {
      warning("Log file is empty or unreadable. Skipping visualization of 'right' peaks from log.")
      hist_log <- data.frame(name=integer(0)) # Empty data frame
    }
  } else {
    warning("Log file does not exist or is empty. Skipping visualization of 'right' peaks from log.")
    hist_log <- data.frame(name=integer(0)) # Empty data frame
  }
  
  # set.seed(456) # For reproducibility
  if (nrow(hist_log) > 0 && "name" %in% colnames(hist_log) && length(unique(hist_log$name)) >=5 ) {
    indices_right_log <- sample(unique(hist_log$name), 5)
    # Ensure indices are valid for pos.1 and pos.2
    indices_right_log <- indices_right_log[indices_right_log <= npeaks & indices_right_log > 0]
    if(length(indices_right_log) < 5 && length(unique(hist_log$name)) >=1){ # if not enough unique names after filtering, sample with replacement
      indices_right_log <- sample(unique(hist_log$name)[unique(hist_log$name) <= npeaks & unique(hist_log$name) > 0], 5, replace=TRUE)
    }
    
  } else {
    warning("Not enough valid entries in hist_log$name to sample 5 'right' peaks. Using fallback sampling.")
    indices_right_log <- sample(1:nde, 5, replace = TRUE) # Fallback
  }
  
  # Check if indices_right_log is not empty
  if(length(indices_right_log) == 0) { # if still empty after fallback (e.g. nde < 1)
    indices_right_log <- sample(1:min(npeaks,5), min(npeaks,5)) # sample from first few peaks
  }
  
  
  my_location_right <- (pos.1[indices_right_log] + pos.2[indices_right_log]) / 2
  
  indices_wrong <- sample(max( (2*nde + 1), (nde + 1) ):npeaks, 5, replace = TRUE) # Sample from non-DE peaks
  # Ensure indices_wrong are valid
  indices_wrong <- indices_wrong[indices_wrong <= npeaks & indices_wrong > 0]
  if(length(indices_wrong) < 5) indices_wrong <- sample(1:npeaks, 5, replace=TRUE) # Ultimate fallback
  
  my_location_wrong <- (pos.1[indices_wrong] + pos.2[indices_wrong]) / 2
  
  
  my_peaks_bins_right <- round(my_location_right / bins)
  my_peaks_bins_wrong <- round(my_location_wrong / bins)
  
  # Ensure bin indices are within bounds of count_bin2
  my_peaks_bins_right <- pmin(pmax(1, my_peaks_bins_right), nrow(count_bin2))
  my_peaks_bins_wrong <- pmin(pmax(1, my_peaks_bins_wrong), nrow(count_bin2))
  
  
  num_cores <- 1
  if (requireNamespace("parallel", quietly = TRUE) && .Platform$OS.type != "windows") {
    num_cores_detected <- parallel::detectCores()
    num_cores <- if (!is.na(num_cores_detected) && num_cores_detected > 0) min(num_cores_detected, 5) else 1
  } else if (.Platform$OS.type == "windows") {
    cat("mclapply not available on Windows, using lapply (sequential processing).\n")
  }
  
  apply_processing <- if (num_cores > 1 && .Platform$OS.type != "windows") {
    function(...) parallel::mclapply(..., mc.cores = num_cores)
  } else {
    lapply
  }
  
  col_idx_conA <- which(colnames(count_bin2) == "conA")
  col_idx_conB <- which(colnames(count_bin2) == "conB")
  if(length(col_idx_conA) == 0 || length(col_idx_conB) == 0) stop("conA or conB columns not found in count_bin2")
  
  # Note: process_peak_local width argument is `count_length`
  data_list_conA_right <- apply_processing(my_peaks_bins_right, process_peak_local, df = count_bin2[, c(1, 2, col_idx_conA)], width = count_length)
  data_conA_right <- do.call(rbind, lapply(data_list_conA_right, t))
  data_list_conB_right <- apply_processing(my_peaks_bins_right, process_peak_local, df = count_bin2[, c(1, 2, col_idx_conB)], width = count_length)
  data_conB_right <- do.call(rbind, lapply(data_list_conB_right, t))
  
  data_list_conA_wrong <- apply_processing(my_peaks_bins_wrong, process_peak_local, df = count_bin2[, c(1, 2, col_idx_conA)], width = count_length)
  data_conA_wrong <- do.call(rbind, lapply(data_list_conA_wrong, t))
  data_list_conB_wrong <- apply_processing(my_peaks_bins_wrong, process_peak_local, df = count_bin2[, c(1, 2, col_idx_conB)], width = count_length)
  data_conB_wrong <- do.call(rbind, lapply(data_list_conB_wrong, t))
  
  # Check data frames
  if (any(sapply(list(data_conA_right, data_conB_right, data_conA_wrong, data_conB_wrong), function(df) nrow(df) == 0 || ncol(df) <= 2))) {
    stop("Data processing for visualization resulted in one or more empty/insufficient data frames.")
  }
  
  dfA_sequence_right <- data_conA_right[, -c(1, 2), drop = FALSE]
  dfB_sequence_right <- data_conB_right[, -c(1, 2), drop = FALSE]
  dfA_sequence_wrong <- data_conA_wrong[, -c(1, 2), drop = FALSE]
  dfB_sequence_wrong <- data_conB_wrong[, -c(1, 2), drop = FALSE]
  
  # Original script expected count_length + 1 columns from process_peak's sequence part.
  # The Index = 1:(count_length+1) implies this.
  # The process_peak_local extracts width+1 elements (if width=count_length, then count_length+1).
  actual_seq_len <- count_length + 1 
  if(ncol(dfA_sequence_right) != actual_seq_len) warning(paste("Seq length for dfA_sequence_right is", ncol(dfA_sequence_right), "but expected", actual_seq_len))
  
  
  dfA_long_right <- data.frame(Index = 1:ncol(dfA_sequence_right), t(dfA_sequence_right))
  colnames(dfA_long_right) <- c("Index", paste("HM_Right", seq_len(nrow(dfA_sequence_right)), sep = "_"))
  dfB_long_right <- data.frame(Index = 1:ncol(dfB_sequence_right), t(dfB_sequence_right))
  colnames(dfB_long_right) <- c("Index", paste("HM_Right", seq_len(nrow(dfB_sequence_right)), sep = "_"))
  
  dfA_long_right <- data.table::melt(data.table::as.data.table(dfA_long_right), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfB_long_right <- data.table::melt(data.table::as.data.table(dfB_long_right), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfA_long_right$Group <- "condition1"
  dfB_long_right$Group <- "condition2"
  
  dfA_long_wrong <- data.frame(Index = 1:ncol(dfA_sequence_wrong), t(dfA_sequence_wrong))
  colnames(dfA_long_wrong) <- c("Index", paste("HM_Wrong", seq_len(nrow(dfA_sequence_wrong)), sep = "_"))
  dfB_long_wrong <- data.frame(Index = 1:ncol(dfB_sequence_wrong), t(dfB_sequence_wrong))
  colnames(dfB_long_wrong) <- c("Index", paste("HM_Wrong", seq_len(nrow(dfB_sequence_wrong)), sep = "_"))
  
  dfA_long_wrong <- data.table::melt(data.table::as.data.table(dfA_long_wrong), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfB_long_wrong <- data.table::melt(data.table::as.data.table(dfB_long_wrong), id.vars = "Index", variable.name = "Series", value.name = "Value")
  dfA_long_wrong$Group <- "condition1"
  dfB_long_wrong$Group <- "condition2"
  
  df_combined <- rbind(dfA_long_right, dfB_long_right, dfA_long_wrong, dfB_long_wrong)
  
  series_levels <- c(
    grep("HM_Right", unique(as.character(df_combined$Series)), value = TRUE),
    grep("HM_Wrong", unique(as.character(df_combined$Series)), value = TRUE)
  )
  df_combined$Series <- factor(df_combined$Series, levels = series_levels)
  
  plot_object <- ggplot2::ggplot(df_combined, ggplot2::aes(x = Index, y = Value, color = Group)) +
    ggplot2::geom_line(alpha = 0.99, linewidth = 0.5) +
    ggplot2::facet_wrap(~ Series, nrow = 2, scales = "free_y") +
    ggplot2::labs(title = "Line Plot of Simulated Histone Mark ChIP-seq Counts",
                  x = "Relative Bin Index from Peak Center", y = "Average Binned Counts", color = "Condition") +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  plot_filename <- file.path(getwd(), paste0(all.fix, "_visualization_plot.png"))
  ggplot2::ggsave(plot_filename, plot = plot_object, width = 12, height = 7, dpi = 300, bg = "white") # Added bg="white"
  cat("Visualization plot saved to:", plot_filename, "\n")
  
  cat("Histone mark simulation and visualization complete.\n")
  
  return(list(
    log_file = file.path(getwd(), lfile),
    bam_files = file.path(getwd(), bam.files),
    plot_file = plot_filename,
    data_directory = getwd(),
    pos1_file = file.path(getwd(), "pos1_vector.RData"),
    pos2_file = file.path(getwd(), "pos2_vector.RData"),
    pos3_file = file.path(getwd(), "pos3_vector.RData")
  ))
}

# results_hm <- simulate_hm(
#   data_path = "/Users/cuitengfei/Graduate/Research/Cancer/SCAW_3rd/2_Server_script/Data/practice_hm_1",
#   true.width = 500,
#   base.mu = 90,
#   bins = 15,
#   count_length = 200
# )

####################################################################################
run_csaw_differential_analysis <- function(
    input_data_dir,
    output_dir,
    is.tf,
    window_width,
    peak_width,
    grouping_vector = c("A", "A", "B", "B"), # Default grouping
    csaw_count_filter_value = 20,
    csaw_extension_length = 100,
    csaw_binned_width_bg = 8000,
    csaw_merge_tolerance = 100,
    csaw_merge_max_width = 5000,
    save_results_df = TRUE,        # New parameter: whether to save the detailed results CSV
    save_roc_data = TRUE,          # New parameter: whether to save the ROC data CSV
    iteration_number = 1           # New parameter: iteration identifier
) {
  
  # --- Start Timer ---
  start_time <- Sys.time()
  
  # Ensure ROCR is loaded if not handled externally, for AUC calculation
  # if (!requireNamespace("ROCR", quietly = TRUE)) {
  #   stop("Package 'ROCR' is required for AUC calculation. Please install and load it.")
  # }
  
  # --- Determine file prefixes based on is.tf ---
  all.fix <- if (is.tf) "tfx" else "hist"
  
  # --- Setup Working Directories ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Output directory created:", output_dir, "\n")
  } else {
    cat("Output directory already exists:", output_dir, "\n")
  }
  
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE) # Ensure WD is reset on exit
  
  # --- Construct Input File Paths ---
  bam_file_names <- paste0(all.fix, "_out_", 1:length(grouping_vector), ".bam")
  bam.files <- file.path(input_data_dir, bam_file_names)
  
  if (!all(sapply(bam.files, file.exists))) {
    stop("One or more BAM files not found in input_data_dir: ", input_data_dir,
         "\nExpected names like: ", paste(bam_file_names, collapse=", "))
  }
  cat("Using BAM files:\n")
  print(bam.files)
  
  lfile <- file.path(input_data_dir, paste0(all.fix, "_log.txt"))
  if (!file.exists(lfile)) {
    stop("Log file not found: ", lfile)
  }
  
  pos1_path <- file.path(input_data_dir, "pos1_vector.RData")
  if (!file.exists(pos1_path)) stop("pos1_vector.RData not found in ", input_data_dir)
  load(pos1_path, envir = environment()) 
  
  if (!is.tf) {
    pos2_path <- file.path(input_data_dir, "pos2_vector.RData")
    pos3_path <- file.path(input_data_dir, "pos3_vector.RData")
    if (!file.exists(pos2_path)) stop("pos2_vector.RData not found in ", input_data_dir)
    if (!file.exists(pos3_path)) stop("pos3_vector.RData not found in ", input_data_dir)
    load(pos2_path, envir = environment()) 
    load(pos3_path, envir = environment()) 
  }
  
  # --- CSAW Parameters and Design Matrix ---
  if (!exists("readParam", mode = "function")) {
    if (requireNamespace("csaw", quietly = TRUE)) {
      cat("Using csaw::readParam\n")
      xparam <- csaw::readParam(dedup = FALSE)
    } else {
      stop("Function 'readParam' not found and csaw package not available. Please source ChIPtest_source.R or install/load csaw.")
    }
  } else {
    cat("Using custom readParam\n")
    xparam <- readParam(dedup = FALSE) 
  }
  
  design <- model.matrix(~ factor(grouping_vector))
  
  # --- CSAW Analysis Pipeline ---
  cat("Starting CSAW analysis with window width:", window_width, "\n")
  
  setwd(output_dir)
  cat("Current working directory set to output directory:", getwd(), "\n")
  
  if (!requireNamespace("csaw", quietly = TRUE)) stop("csaw package is required.")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR package is required.")
  
  data <- csaw::windowCounts(bam.files, width = window_width, ext = csaw_extension_length,
                             param = xparam, filter = csaw_count_filter_value)
  cat("Initial windowCounts complete. Dimensions of 'data':", dim(SummarizedExperiment::assay(data)), "\n")
  
  binned <- csaw::windowCounts(bam.files, width = csaw_binned_width_bg, bin = TRUE, param = xparam)
  cat("Background estimation windowCounts complete.\n")
  
  bin.ab <- csaw::scaledAverage(binned, scale = median(csaw::getWidths(binned)) / median(csaw::getWidths(data)))
  threshold <- median(bin.ab) + log2(2) 
  cat("Calculated abundance threshold:", threshold, "\n")
  
  keep <- edgeR::aveLogCPM(asDGEList(data)) > threshold
  data <- data[keep, ]
  cat("Filtered windows based on abundance. Dimensions of 'data':", dim(SummarizedExperiment::assay(data)), "\n")
  if (nrow(data) == 0) {
    stop("No windows remained after abundance filtering. Check filter thresholds or data quality.")
  }
  
  if (!exists("analyzeQLCounts", mode = "function")) {
    stop("Function 'analyzeQLCounts' not found. Please source ChIPtest_source.R.")
  }
  tabres <- analyzeQLCounts(SummarizedExperiment::assay(data), design, totals = data$totals)
  cat("Differential binding analysis with analyzeQLCounts complete.\n")
  
  merged <- csaw::mergeWindows(SummarizedExperiment::rowRanges(data), tol = csaw_merge_tolerance, max.width = csaw_merge_max_width)
  cat("Windows merged. Number of merged regions:", length(merged$region), "\n")
  if (length(merged$region) == 0) {
    stop("No regions remained after merging windows.")
  }
  
  tabneg <- csaw::combineTests(merged$id, tabres)
  cat("Tests combined for merged regions.\n")
  
  # --- Process Results and Find Overlaps ---
  results_csaw_df <- data.frame(
    index = 1:length(merged$region),
    start = GenomicRanges::start(merged$region),
    end = GenomicRanges::end(merged$region),
    Pvalues = tabneg[["PValue"]], 
    FDR = tabneg[["FDR"]]
  )
  csaw_sites <- results_csaw_df
  
  if (!exists("find_overlap", mode = "function")) {
    stop("Function 'find_overlap' not found. Please source ChIPtest_source.R.")
  }
  
  cat("Finding overlaps with true sites...\n")
  if (is.tf) {
    csaw_sites_new <- find_overlap(is.tf, pos.1, pos.1, peak_width, lfile, csaw_sites)
  } else {
    csaw_sites_new <- find_overlap(is.tf, pos.1, pos.3, peak_width, lfile, csaw_sites)
  }
  cat("Overlap finding complete.\n")
  csaw_sites_new[is.na(overlap_true), overlap_true := 0]
  
  results_csaw_df <- csaw_sites_new
  
  # --- Calculate Performance Metrics ---
  cat("Calculating performance metrics...\n")
  if (!exists("calculate_fdr_recall", mode = "function")) {
    stop("Function 'calculate_fdr_recall' not found. Please source ChIPtest_source.R.")
  }
  fdr_recall_results <- calculate_fdr_recall(results_csaw_df$FDR, results_csaw_df$overlap,
                                             thresholds = c(0.05, 0.1, 0.2))
  
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    warning("ROCR package not available. Skipping AUC-ROC calculation.")
    auc_value <- NA
    roc_data <- data.frame(TPR=NA, FPR=NA, Threshold=NA)
  } else {
    valid_for_roc <- !is.na(results_csaw_df$FDR) & !is.na(results_csaw_df$overlap)
    if(sum(valid_for_roc) < 2 || length(unique(results_csaw_df$overlap[valid_for_roc])) < 2) {
      warning("Not enough distinct data points or classes for ROC calculation. AUC will be NA.")
      auc_value <- NA
      roc_data <- data.frame(TPR=NA, FPR=NA, Threshold=NA)
    } else {
      pred_obj <- ROCR::prediction(1 - results_csaw_df$FDR[valid_for_roc], results_csaw_df$overlap[valid_for_roc])
      perf_obj <- ROCR::performance(pred_obj, "tpr", "fpr")
      roc_data <- data.frame(
        TPR = perf_obj@y.values[[1]],
        FPR = perf_obj@x.values[[1]],
        Threshold = perf_obj@alpha.values[[1]]
      )
      auc_obj <- ROCR::performance(pred_obj, "auc")
      auc_value <- auc_obj@y.values[[1]]
    }
  }
  cat("AUC-ROC value:", auc_value, "\n")
  
  csaw_performance_metrics <- as.data.frame(fdr_recall_results) # Ensure it's a data frame
  csaw_performance_metrics$AUCROC <- auc_value
  
  # --- Calculate Runtime ---
  end_time <- Sys.time()
  runtime_diff <- difftime(end_time, start_time, units = "secs")
  runtime_seconds <- as.numeric(runtime_diff)
  
  # --- Add new columns to csaw_performance_metrics and reorder ---
  csaw_performance_metrics$iteration <- iteration_number
  csaw_performance_metrics$runtime_seconds <- runtime_seconds
  
  # Define the desired column order
  # Get original columns from fdr_recall_results (assuming it's a data frame with named columns)
  # and AUCROC
  original_metric_cols <- names(fdr_recall_results) 
  if (!("AUCROC" %in% original_metric_cols)) { # Should be added if not already by some magic
    original_metric_cols <- c(original_metric_cols, "AUCROC")
  }
  # Remove any potential duplicates if names accidentally match
  original_metric_cols <- original_metric_cols[!original_metric_cols %in% c("iteration", "runtime_seconds")]
  
  # Reorder columns
  csaw_performance_metrics <- csaw_performance_metrics[, c("iteration", "runtime_seconds", original_metric_cols)]
  
  # --- Save Results ---
  results_filename <- paste0("results_csaw_ww", window_width, "_iter", iteration_number, ".csv")
  metrics_filename <- paste0("metrics_csaw_ww", window_width, "_iter", iteration_number, ".csv")
  roc_filename <- paste0("roc_data_csaw_ww", window_width, "_iter", iteration_number, ".csv")
  
  if (save_results_df) {
    write.csv(results_csaw_df, file = results_filename, row.names = FALSE)
    cat("CSAW detailed results saved to:", file.path(getwd(), results_filename), "\n")
  } else {
    cat("Skipping saving of detailed results_csaw_df as per save_results_df=FALSE.\n")
  }
  
  # Performance metrics are always saved
  write.csv(csaw_performance_metrics, file = metrics_filename, row.names = FALSE)
  cat("Performance metrics saved to:", file.path(getwd(), metrics_filename), "\n")
  
  if (save_roc_data) {
    if(exists("roc_data") && is.data.frame(roc_data) && !all(is.na(roc_data$TPR))) {
      write.csv(roc_data, file = roc_filename, row.names = FALSE)
      cat("ROC curve data saved to:", file.path(getwd(), roc_filename), "\n")
    } else {
      cat("ROC data not available or not valid for saving, even though save_roc_data=TRUE.\n")
    }
  } else {
    cat("Skipping saving of ROC curve data as per save_roc_data=FALSE.\n")
  }
  
  cat("CSAW analysis complete. Total runtime:", format(runtime_diff), "\n")
  
  return(list(
    performance = csaw_performance_metrics,
    results_file = if(save_results_df) file.path(getwd(), results_filename) else NA,
    metrics_file = file.path(getwd(), metrics_filename),
    roc_data_file = if(save_roc_data && exists("roc_data") && is.data.frame(roc_data) && !all(is.na(roc_data$TPR))) file.path(getwd(), roc_filename) else NA,
    output_directory_used = getwd(),
    runtime_seconds = runtime_seconds,
    runtime_formatted = format(runtime_diff)
  ))
}

# analysis_results_hm_iter1 <- run_csaw_differential_analysis(
#   input_data_dir = "/Users/cuitengfei/Graduate/Research/Cancer/SCAW_3rd/2_Server_script/Data/practice_hm_1",
#   output_dir = "/Users/cuitengfei/Graduate/Research/Cancer/SCAW_3rd/2_Server_script/Results/HM_CSAW",
#   is.tf = FALSE,
#   window_width = 1000,
#   peak_width = 500,
#   grouping_vector = c("A", "A", "B", "B"),
#   csaw_count_filter_value = 20,
#   csaw_extension_length = 100,
#   save_results_df = FALSE,    # Explicitly TRUE
#   save_roc_data = FALSE,      # Explicitly TRUE
#   iteration_number = 1       # Iteration 1
# )

###################################################################################
run_chiptest_slide2_analysis <- function(
    input_data_dir,          # Full path to the specific setting directory
    output_parent_dir,       # Directory where the results folder will be created
    is.tf,                   # Boolean: TRUE for TF, FALSE for HM
    peak_width,              # Numeric: peak width used in find_overlap
    chiptest_params = list(  # List of parameters specific to this ChIPtest variant
      bins_tf = 25,
      counts_tf = 120,
      filter_tf = 0,
      window_width_csaw_tf = 150,
      bins_hm = 15,
      counts_hm = 120,
      filter_hm = 0,
      window_width_csaw_hm = 150
    )
) {
  
  # --- Start Timer ---
  overall_start_time <- Sys.time()
  
  # --- Internal Default Parameters (previously arguments) ---
  gsize_internal <- 300678703
  csaw_main_params_internal <- list(
    ext = 100,
    filter_initial_data = 20,
    bg_width = 2000,
    merge_tol = 100,
    merge_max_width = 5000
  )
  finetune_params_internal <- list(
    band_values = c(10, 30, 60),
    quantile_values = c(-1, 0.0001, 0.01, 0.02, 0.1, 0.2, 0.5, 0.99, 1.2, 1.5, 1.8),
    my_var.est = 1,
    my_var.thred = 0.01
  )
  chiptest_hyp_params_internal <- list(
    my_est_c_max1 = 4,
    my_est_c_max4 = 4,
    my_TS_var.est = 1,
    my_TS_var.thred = 0.1
  )
  num_cores_internal <- 5 # Default number of cores
  
  # --- Argument Validation (Basic) ---
  if (!dir.exists(input_data_dir)) {
    stop("Input data directory not found: ", input_data_dir)
  }
  if (!is.logical(is.tf)) stop("'is.tf' must be TRUE or FALSE.")
  if (!is.numeric(peak_width)) stop("'peak_width' must be numeric.")
  
  # --- Load required libraries ---
  required_packages <- c("csaw", "edgeR", "data.table", "ROCR", "ggplot2", "parallel")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package '", pkg, "' is required but not installed/available.", sep = ""))
    }
  }
  
  # --- Determine analysis parameters based on is.tf ---
  if (is.tf) {
    all.fix <- "tfx"
    current_chiptest_bin <- chiptest_params$bins_tf
    current_chiptest_count <- chiptest_params$counts_tf
    current_chiptest_filter <- chiptest_params$filter_tf
    current_csaw_window_width <- chiptest_params$window_width_csaw_tf
  } else {
    all.fix <- "hist"
    current_chiptest_bin <- chiptest_params$bins_hm
    current_chiptest_count <- chiptest_params$counts_hm
    current_chiptest_filter <- chiptest_params$filter_hm
    current_csaw_window_width <- chiptest_params$window_width_csaw_hm
  }
  
  # --- Manage Working Directory ---
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  
  # --- Create and Set Output Directory ---
  output_dir_name <- paste0("results_slidewin2_ChIPtest_bin", current_chiptest_bin)
  output_dir_full_path <- file.path(output_parent_dir, output_dir_name)
  
  if (!dir.exists(output_dir_full_path)) {
    dir.create(output_dir_full_path, recursive = TRUE)
    cat("Output directory created:", output_dir_full_path, "\n")
  }
  setwd(output_dir_full_path)
  cat("Current working directory set to:", getwd(), "\n")
  
  # --- Construct Input File Paths from input_data_dir ---
  bam_file_basenames <- paste0(all.fix, "_out_", 1:4, ".bam") # Assuming 4 BAM files
  bam.files <- file.path(input_data_dir, bam_file_basenames)
  
  if (!all(sapply(bam.files, file.exists))) {
    missing_bams <- bam.files[!sapply(bam.files, file.exists)]
    stop("One or more BAM files not found in '", input_data_dir, "':\n", paste(missing_bams, collapse = "\n"))
  }
  cat("Using BAM files:\n")
  print(bam.files)
  
  lfile <- file.path(input_data_dir, paste0(all.fix, "_log.txt"))
  if (!file.exists(lfile)) {
    stop("Log file not found: ", lfile)
  }
  
  pos_env <- new.env()
  load(file.path(input_data_dir, "pos1_vector.RData"), envir = pos_env)
  if (!is.tf) {
    load(file.path(input_data_dir, "pos2_vector.RData"), envir = pos_env)
    load(file.path(input_data_dir, "pos3_vector.RData"), envir = pos_env)
  }
  
  if (!exists("readParam", mode = "function")) {
    stop("Function 'readParam' not found. Ensure ChIPtest_source.R is sourced.")
  }
  # Note: gsize_internal is available here if your readParam needs it:
  # xparam <- readParam(dedup = FALSE, gsize = gsize_internal)
  xparam <- readParam(dedup = FALSE)
  
  
  # --- CSAW: Initial window counting to find regions of interest ---
  cat("Starting initial CSAW processing...\n")
  data <- csaw::windowCounts(bam.files, width = current_csaw_window_width,
                             ext = csaw_main_params_internal$ext, param = xparam, # Using internal default
                             filter = csaw_main_params_internal$filter_initial_data) # Using internal default
  binned <- csaw::windowCounts(bam.files, width = csaw_main_params_internal$bg_width, # Using internal default
                               bin = TRUE, param = xparam)
  
  bin.ab <- csaw::scaledAverage(binned, scale = median(csaw::getWidths(binned)) / median(csaw::getWidths(data)))
  threshold <- median(bin.ab) + log2(2)
  keep <- edgeR::aveLogCPM(asDGEList(data)) > threshold
  # keep <- edgeR::aveLogCPM(asDGEList(data)) > threshold
  data <- data[keep, ]
  
  if (nrow(data) == 0) {
    warning("No windows remained after initial abundance filtering in CSAW. Returning NULL.")
    return(list(performance = NULL, results_file = NA, metrics_file = NA, output_directory_used = output_dir_full_path, error = "No windows after CSAW abundance filtering"))
  }
  
  merged <- csaw::mergeWindows(SummarizedExperiment::rowRanges(data),
                               tol = csaw_main_params_internal$merge_tol, # Using internal default
                               max.width = csaw_main_params_internal$merge_max_width) # Using internal default
  regions <- merged$region
  if (length(regions) == 0) {
    warning("No regions remained after CSAW mergeWindows. Returning NULL.")
    return(list(performance = NULL, results_file = NA, metrics_file = NA, output_directory_used = output_dir_full_path, error = "No regions after CSAW mergeWindows"))
  }
  
  csaw_sites <- data.frame(start = GenomicRanges::start(regions), end = GenomicRanges::end(regions))
  csaw_sites$index <- 1:nrow(csaw_sites)
  cat("Identified", nrow(csaw_sites), "regions of interest from initial CSAW processing.\n")
  
  # --- Re-construct count sequences for ChIPtest using specified bins ---
  cat("Re-constructing count sequences for ChIPtest...\n")
  count_seq <- csaw::windowCounts(bam.files, width = current_chiptest_bin, spacing = current_chiptest_bin,
                                  ext = current_chiptest_bin, filter = current_chiptest_filter,
                                  bin = TRUE, param = xparam)
  counts_matrix <- SummarizedExperiment::assay(count_seq)
  count_bin2 <- data.frame(start = GenomicRanges::start(SummarizedExperiment::rowRanges(count_seq)),
                           end = GenomicRanges::end(SummarizedExperiment::rowRanges(count_seq)),
                           counts_matrix)
  if (ncol(counts_matrix) < 4) stop("Less than 4 count columns found after windowCounts.")
  count_bin2$conA <- (counts_matrix[,1] + counts_matrix[,2]) / 2
  count_bin2$conB <- (counts_matrix[,3] + counts_matrix[,4]) / 2
  
  
  # --- Process peaks for ChIPtest (data_conA, data_conB) ---
  cat("Processing peaks for ChIPtest input...\n")
  potential_peaks_indices <- round((csaw_sites$start + csaw_sites$end) / (2 * current_chiptest_bin))
  potential_peaks_indices <- pmin(pmax(1, potential_peaks_indices), nrow(count_bin2))
  
  if (!exists("process_peak_conA", mode = "function") || !exists("process_peak_conB", mode = "function")) {
    stop("Functions 'process_peak_conA' or 'process_peak_conB' not found.")
  }
  
  apply_fn <- if (.Platform$OS.type != "windows" && num_cores_internal > 1) { # Using internal default
    function(...) parallel::mclapply(..., mc.cores = num_cores_internal)
  } else {
    if (.Platform$OS.type == "windows" && num_cores_internal > 1) {
      cat("mclapply not available on Windows, using sequential lapply.\n")
    }
    lapply
  }
  
  data_list_conA <- apply_fn(potential_peaks_indices, process_peak_conA,
                             df = count_bin2[, c("start", "end", "conA")],
                             width = current_chiptest_count)
  data_conA <- do.call(rbind, lapply(data_list_conA, t))
  
  data_list_conB <- apply_fn(potential_peaks_indices, process_peak_conB,
                             df = count_bin2[, c("start", "end", "conB")],
                             width = current_chiptest_count)
  data_conB <- do.call(rbind, lapply(data_list_conB, t))
  
  if (nrow(data_conA) == 0 || nrow(data_conB) == 0) {
    warning("No data generated for data_conA or data_conB. ChIPtest cannot proceed. Returning NULL.")
    return(list(performance = NULL, results_file = NA, metrics_file = NA, output_directory_used = output_dir_full_path, error = "data_conA or data_conB is empty"))
  }
  
  # --- Finetune ChIPtest Null Hypothesis ---
  cat("Finetuning ChIPtest null hypothesis...\n")
  if (!exists("Finetune_ChIPtest", mode = "function") || !exists("finetune", mode = "function")) {
    stop("Functions 'Finetune_ChIPtest' or 'finetune' not found.")
  }
  finetune_results <- Finetune_ChIPtest(
    potential_peaks = potential_peaks_indices,
    count_bin2 = count_bin2,
    width = current_chiptest_count,
    band_values = finetune_params_internal$band_values, # Using internal default
    quantile_values = finetune_params_internal$quantile_values, # Using internal default
    my_var.est = finetune_params_internal$my_var.est, # Using internal default
    my_var.thred = finetune_params_internal$my_var.thred, # Using internal default
    num_cores = num_cores_internal # Using internal default
  )
  finetune_results$TS_kn_KL_mean <- rowMeans(finetune_results[, c("TS_kn_KL_12", "TS_kn_KL_34")], na.rm = TRUE)
  finetune_results$Deql_KL_mean <- rowMeans(finetune_results[, c("Deql_KL_12", "Deql_KL_34")], na.rm = TRUE)
  finetune_results$Dnun_KL_mean <- rowMeans(finetune_results[, c("Dnun_KL_12", "Dnun_KL_34")], na.rm = TRUE)
  
  quantile_TSkn_min <- finetune_results$Quantile[which.min(finetune_results$TS_kn_KL_mean)]
  band_TSkn_min <- finetune_results$Band[which.min(finetune_results$TS_kn_KL_mean)]
  quantile_Deql_min <- finetune_results$Quantile[which.min(finetune_results$Deql_KL_mean)]
  quantile_Dnun_min <- finetune_results$Quantile[which.min(finetune_results$Dnun_KL_mean)]
  
  finetune(df = count_bin2, width = current_chiptest_count,
           potential_peaks = potential_peaks_indices,
           my_quantile = c(quantile_Deql_min, quantile_Dnun_min, quantile_TSkn_min),
           my_band = band_TSkn_min, my_var.est = 1, my_var.thred = 0.1,
           plot_filename_base = "finetune_plot")
  cat("Finetuning complete. Optimal parameters selected.\n")
  
  
  # --- ChIPtest Hypothesis Testing ---
  cat("Performing ChIPtest hypothesis testing...\n")
  result_temp <- as.data.frame(data_conA[, c(1, 2)])
  colnames(result_temp) <- c("start", "end")
  result_temp$index <- 1:nrow(result_temp)
  
  if (!exists("find_overlap", mode = "function") || !exists("NormTransformation", mode = "function") ||
      !exists("my_est.c", mode = "function") || !exists("my_TS_twosample", mode = "function")) {
    stop("One or more ChIPtest core functions not found.")
  }
  
  if (is.tf) {
    sites_overlap <- find_overlap(is.tf, pos_env$pos.1, pos_env$pos.1, peak_width, lfile, result_temp)
  } else {
    sites_overlap <- find_overlap(is.tf, pos_env$pos.1, pos_env$pos.3, peak_width, lfile, result_temp)
  }
  
  Data_conA_norm <- NormTransformation(data_conA[, -c(1, 2)])
  Data_conB_norm <- NormTransformation(data_conB[, -c(1, 2)])
  tao <- my_est.c(Data_conA_norm, Data_conB_norm,
                  max1 = chiptest_hyp_params_internal$my_est_c_max1, # Using internal default
                  max4 = chiptest_hyp_params_internal$my_est_c_max4) # Using internal default
  
  TS <- my_TS_twosample(Data_conA_norm, Data_conB_norm, tao, band_TSkn_min,
                        quant = c(quantile_Deql_min, quantile_Dnun_min, quantile_TSkn_min),
                        var.est = chiptest_hyp_params_internal$my_TS_var.est, # Using internal default
                        var.thred = chiptest_hyp_params_internal$my_TS_var.thred) # Using internal default
  
  setnames(sites_overlap, old = "overlap_true", new = "overlap")
  sites_overlap[is.na(overlap), overlap := 0]
  ChIPtest_result <- copy(sites_overlap)
  ChIPtest_result_non_na_idx <- ChIPtest_result[!is.na(index)]
  ChIPtest_result_na_idx <- ChIPtest_result[is.na(index)]
  
  
  # 处理 index 非 NA 的部分
  if (nrow(ChIPtest_result_non_na_idx) > 0) {
    # 按照 index 从小到大排序 (确保正确索引 TS 列表)
    setorderv(ChIPtest_result_non_na_idx, "index")
    
    # 获取 TS 中需要引用的索引值
    ts_indices <- ChIPtest_result_non_na_idx$index
    
    # 安全性检查：确保 ts_indices 不会超出 TS 向量的范围
    max_ts_len <- length(TS[["TS_kn"]])
    # 任何超出范围的索引都会导致取值变成 NA
    ts_indices_safe <- ts_indices
    ts_indices_safe[ts_indices_safe < 1 | ts_indices_safe > max_ts_len] <- NA_integer_
    
    # 计算 P 值
    ChIPtest_result_non_na_idx[, pvalue_TS_high := pnorm(TS[["TS_kn"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
    ChIPtest_result_non_na_idx[, pvalue_Deql := pnorm(TS[["Deql"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
    ChIPtest_result_non_na_idx[, pvalue_Dnun := pnorm(TS[["Dnun"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
  } else {
    # 如果没有非 NA 的 index 行，则创建空的 data.table，并包含相应的列
    ChIPtest_result_non_na_idx <- ChIPtest_result_non_na_idx[, `:=`(pvalue_TS_high = numeric(), pvalue_Deql = numeric(), pvalue_Dnun = numeric())]
  }
  
  # 3. 处理 index 为 NA 的部分
  if (nrow(ChIPtest_result_na_idx) > 0) {
    # 将所有 P 值赋为 1
    ChIPtest_result_na_idx[, pvalue_TS_high := 1]
    ChIPtest_result_na_idx[, pvalue_Deql := 1]
    ChIPtest_result_na_idx[, pvalue_Dnun := 1]
  } else {
    # 如果没有 NA 的 index 行，则创建空的 data.table
    ChIPtest_result_na_idx <- ChIPtest_result_na_idx[, `:=`(pvalue_TS_high = numeric(), pvalue_Deql = numeric(), pvalue_Dnun = numeric())]
  }
  # 合并两部分结果
  ChIPtest_result <- rbindlist(list(ChIPtest_result_non_na_idx, ChIPtest_result_na_idx), fill = TRUE)
  ChIPtest_result$chromosome <- 'chr1'
  ChIPtest_result$TSkn_BH <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "BH")
  ChIPtest_result$Dnun_BH <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "BH")
  ChIPtest_result$Deql_BH <- p.adjust(ChIPtest_result$pvalue_Deql, method = "BH")
  ChIPtest_result$TSkn_bon <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "bonferroni")
  ChIPtest_result$Dnun_bon <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "bonferroni")
  ChIPtest_result$Deql_bon <- p.adjust(ChIPtest_result$pvalue_Deql, method = "bonferroni")
  
  # --- Select most significant p-value for overlapping true sites ---
  dt_ChIPtest_result <- data.table::setDT(as.data.frame(ChIPtest_result))
  
  ChIPtest_result_0 <- dt_ChIPtest_result[is.na(sim_start)]
  ChIPtest_result_1 <- dt_ChIPtest_result[!is.na(sim_start), .(
    index = index[1], start = start[1], end = end[1], sim_end = sim_end[1],
    overlap = overlap[1], chromosome = chromosome[1],
    pvalue_TS_high = min(pvalue_TS_high, na.rm = TRUE),
    pvalue_Deql = min(pvalue_Deql, na.rm = TRUE),
    pvalue_Dnun = min(pvalue_Dnun, na.rm = TRUE),
    TSkn_BH = min(TSkn_BH, na.rm = TRUE), Dnun_BH = min(Dnun_BH, na.rm = TRUE),
    Deql_BH = min(Deql_BH, na.rm = TRUE), TSkn_bon = min(TSkn_bon, na.rm = TRUE),
    Dnun_bon = min(Dnun_bon, na.rm = TRUE), Deql_bon = min(Deql_bon, na.rm = TRUE)
  ), by = sim_start]
  
  if(nrow(ChIPtest_result_0) > 0 && nrow(ChIPtest_result_1) > 0) {
    data.table::setcolorder(ChIPtest_result_1, names(ChIPtest_result_0))
    ChIPtest_selected <- data.table::rbindlist(list(ChIPtest_result_1, ChIPtest_result_0), use.names = TRUE, fill = TRUE)
  } else if (nrow(ChIPtest_result_1) > 0) {
    ChIPtest_selected <- ChIPtest_result_1
  } else {
    ChIPtest_selected <- ChIPtest_result_0
  }
  ChIPtest_selected <- as.data.frame(ChIPtest_selected)
  
  
  # --- Calculate Performance Metrics ---
  cat("Calculating final performance metrics...\n")
  if (!exists("calculate_fdr_recall", mode = "function")) {
    stop("Function 'calculate_fdr_recall' not found.")
  }
  
  metrics_list <- list()
  models <- c("TSkn", "Deql", "Dnun")
  controls <- c("BH", "bon")
  
  if(nrow(ChIPtest_selected) == 0 || !"overlap" %in% names(ChIPtest_selected)) {
    warning("ChIPtest_selected is empty or missing 'overlap' column. Cannot calculate performance metrics.")
    fdr_recall_all_df <- data.frame(Threshold=NA, FDR=NA, Recall=NA, model=NA, control=NA, auc_value=NA)[-1,]
  } else {
    for (model in models) {
      for (control in controls) {
        p_col_name <- paste0(model, "_", control)
        if (!p_col_name %in% names(ChIPtest_selected)) {
          warning(paste("P-value column", p_col_name, "not found in results. Skipping."))
          next
        }
        
        valid_for_roc <- !is.na(ChIPtest_selected[[p_col_name]]) & !is.na(ChIPtest_selected$overlap)
        if(sum(valid_for_roc) < 2 || length(unique(ChIPtest_selected$overlap[valid_for_roc])) < 2) {
          auc_val <- NA
          warning(paste("Not enough distinct data points for ROC calculation for", p_col_name, ". AUC will be NA."))
        } else {
          pred_obj <- ROCR::prediction(1 - ChIPtest_selected[[p_col_name]][valid_for_roc], ChIPtest_selected$overlap[valid_for_roc])
          auc_val <- ROCR::performance(pred_obj, "auc")@y.values[[1]]
        }
        
        temp_metrics <- calculate_fdr_recall(ChIPtest_selected[[p_col_name]], ChIPtest_selected$overlap, thresholds = c(0.05))
        temp_metrics_df <- as.data.frame(temp_metrics)
        temp_metrics_df$model <- model
        temp_metrics_df$control <- control
        temp_metrics_df$auc_value <- auc_val
        metrics_list[[paste0(model, "_", control)]] <- temp_metrics_df
      }
    }
    fdr_recall_all_df <- do.call(rbind, metrics_list)
  }
  
  # --- Save Final Results ---
  results_final_filename <- "results_slidewin2_ChIPtest.csv"
  metrics_final_filename <- "matrics_slidewin2_ChIPtest.csv"
  
  write.csv(ChIPtest_selected, file = results_final_filename, row.names = FALSE)
  cat("Final selected ChIPtest results saved to:", file.path(getwd(), results_final_filename), "\n")
  write.csv(fdr_recall_all_df, file = metrics_final_filename, row.names = FALSE)
  cat("Final performance metrics saved to:", file.path(getwd(), metrics_final_filename), "\n")
  
  # --- End Timer & Return ---
  overall_end_time <- Sys.time()
  overall_runtime_diff <- difftime(overall_end_time, overall_start_time, units = "secs")
  cat("ChIPtest Slide-Window 2 analysis complete. Total runtime:", format(overall_runtime_diff, digits=3), "secs\n")
  
  return(list(
    performance_metrics = fdr_recall_all_df,
    detailed_results_file = file.path(getwd(), results_final_filename),
    metrics_summary_file = file.path(getwd(), metrics_final_filename),
    output_directory = getwd(),
    runtime_seconds = as.numeric(overall_runtime_diff)
  ))
}


###################################################################################
run_pepr_pipeline <- function(
    input_data_dir,         # 直接的数据集路径, 例如: "/Users/.../Datasets/tf_setting_2"
    output_dir_base,        # 例如: "/Users/.../Results/MyPaper_PePr" (统一的父输出目录)
    is_tf,                  # TRUE 或 FALSE
    iteration_tag = "run1", # 用于输出文件名的标签/ID (替代了原 iteration_id 的部分功能)
    peak_width = 1000,      # 传递给 find_overlap
    gsize = 300678703,      # PePr 参数 (脚本中未使用，但通常是需要的)
    fraglen = 100,          # PePr -s 参数相关 (fraglen/2L)
    grouping_vector = c("A","A","B","B"),
    save_initial_pepr_results = FALSE, 
    save_sites_overlap_csv = FALSE,    
    save_final_filtered_results_csv = TRUE,
    save_metrics_csv = TRUE,
    save_roc_data_csv = TRUE
) {
  
  # --- 0. 开始计时和参数初步处理 ---
  pipeline_start_time <- Sys.time()
  cat(paste("Starting PePr Pipeline for input:", input_data_dir, ", is_tf:", is_tf, ", iteration_tag:", iteration_tag, "at", Sys.time(), "\n"))
  
  if (!is.logical(is_tf)) stop("'is_tf' must be TRUE or FALSE.")
  if (!dir.exists(input_data_dir)) {
    stop(paste("Input data directory not found:", input_data_dir))
  }
  current_dataset_path <- input_data_dir # 直接使用提供的路径
  
  # 根据 is_tf 确定文件名前缀和输出标签
  all.fix <- if (is_tf) "tfx" else "hist"
  analysis_type_label <- if (is_tf) "TF" else "HM" 
  
  # --- 1. 构建输出路径 ---
  # 输出目录 specific to this run (使用 analysis_type_label, iteration_tag, peak_width)
  # 从 input_data_dir 获取数据集名称作为输出子目录的一部分，使其更具描述性
  dataset_name_for_output <- basename(current_dataset_path) #例如 "tf_setting_2"
  
  current_run_output_dir_name <- paste0("PePr_", dataset_name_for_output, "_", analysis_type_label, "_tag", iteration_tag, "_pw", peak_width)
  current_output_path <- file.path(output_dir_base, current_run_output_dir_name)
  
  if (!dir.exists(current_output_path)) {
    dir.create(current_output_path, recursive = TRUE)
    cat("Created output directory:", current_output_path, "\n")
  } else {
    cat("Output directory already exists:", current_output_path, "\n")
  }
  
  original_wd <- getwd()
  setwd(current_output_path) 
  on.exit(setwd(original_wd), add = TRUE)
  cat("Working directory set to:", getwd(), "\n")
  
  # --- 2. 加载数据和定义文件 ---
  # 文件名现在基于 all.fix 和 current_dataset_path
  cat("Loading .RData files and defining BAM files from:", current_dataset_path, "\n")
  bam.files <- file.path(current_dataset_path, paste0(all.fix, "_out_", seq_along(grouping_vector), ".bam"))
  if(!all(sapply(bam.files, file.exists))) {
    stop("One or more BAM files not found. Checked paths under ", current_dataset_path, ": ", paste(bam.files[!sapply(bam.files, file.exists)], collapse=", "))
  }
  
  lfile <- file.path(current_dataset_path, paste0(all.fix, "_log.txt"))
  if(!file.exists(lfile)) stop(paste("Log file not found:", lfile))
  
  load_env <- new.env()
  pos1_file <- file.path(current_dataset_path, "pos1_vector.RData")
  if(!file.exists(pos1_file)) stop(paste("pos1_vector.RData not found in:", current_dataset_path))
  load(pos1_file, envir = load_env)
  
  if(!is_tf){
    pos2_file <- file.path(current_dataset_path, "pos2_vector.RData")
    pos3_file <- file.path(current_dataset_path, "pos3_vector.RData")
    if(!file.exists(pos2_file)) stop(paste("pos2_vector.RData not found in:", current_dataset_path))
    if(!file.exists(pos3_file)) stop(paste("pos3_vector.RData not found in:", current_dataset_path))
    load(pos2_file, envir = load_env)
    load(pos3_file, envir = load_env)
  }
  
  if (!exists("pos.1", envir = load_env)) stop("Variable 'pos.1' not loaded from RData files.")
  pos.1_data <- get("pos.1", envir = load_env)
  
  pos.3_data_for_findoverlap <- NULL
  if (!is_tf) {
    if (!exists("pos.3", envir = load_env)) stop("Variable 'pos.3' not loaded for non-TF analysis.")
    pos.3_data_for_findoverlap <- get("pos.3", envir = load_env)
  } else {
    pos.3_data_for_findoverlap <- pos.1_data 
  }
  
  if (!exists("readParam") || !is.function(readParam)) {
    stop("Function 'readParam' is not defined. Please source ChIPtest_source.R or define it.")
  }
  # xparam <- readParam(dedup=FALSE)
  
  # --- 3. 运行 PePr ---
  cat("Starting PePr analysis...\n")
  pepr_run_start_time <- Sys.time()
  
  pepr_native_output_subdir_name <- paste0(all.fix, "_PePr_run_files") 
  pepr_native_output_path <- file.path(getwd(), pepr_native_output_subdir_name) 
  
  if (!dir.exists(pepr_native_output_path)) {
    dir.create(pepr_native_output_path, recursive = TRUE)
    cat("Created subdirectory for PePr's native output files:", pepr_native_output_path, "\n")
  }
  
  outname_for_pepr <- file.path(pepr_native_output_subdir_name, "pepr_analysis") 
  
  first.set <- paste(shQuote(bam.files[grouping_vector=="A"]), collapse=",") 
  second.set <- paste(shQuote(bam.files[grouping_vector=="B"]), collapse=",")
  
  # pepr_command <- sprintf(
  #   "bash -c 'source %s && conda activate %s && PePr -c %s --chip2 %s -n %s -f BAM --peaktype=%s --diff -s %i && conda deactivate'",
  #   shQuote(conda_sh_path), shQuote(conda_env_name),
  #   first.set, second.set, shQuote(outname_for_pepr),
  #   ifelse(is_tf, "sharp", "broad"), as.integer(fraglen/2L)
  # )
  pepr_command <- sprintf(
    "PePr -c %s --chip2 %s -n %s -f BAM --peaktype=%s --diff -s %i",
    first.set, second.set, shQuote(outname_for_pepr),
    ifelse(is_tf, "sharp", "broad"), as.integer(fraglen / 2L)
  )
  
  cat("Executing PePr command:\n", pepr_command, "\n")
  system_status <- system(pepr_command)
  if(system_status != 0) {
    stop(paste("PePr command failed with status:", system_status, ". Check PePr logs in", pepr_native_output_path))
  }
  cat("PePr analysis completed. Time taken:", format(Sys.time() - pepr_run_start_time), "\n")
  
  # --- 4. 处理 PePr 输出 ---
  cat("Processing PePr output files...\n")
  up_peaks_file <- paste0(outname_for_pepr, "__PePr_chip1_peaks.bed")
  down_peaks_file <- paste0(outname_for_pepr, "__PePr_chip2_peaks.bed")
  
  if(!file.exists(up_peaks_file)) stop("PePr chip1 peaks file not found: ", up_peaks_file)
  if(!file.exists(down_peaks_file)) stop("PePr chip2 peaks file not found: ", down_peaks_file)
  
  collected.up <- read.table(up_peaks_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  collected.down <- read.table(down_peaks_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  if (ncol(collected.up) < 9 || ncol(collected.down) < 9) {
    stop("PePr output .bed files do not have enough columns (expected at least 9). Check V8 for p-value, V9 for FDR.")
  }
  
  all.collected <- c(
    GenomicRanges::GRanges(collected.up[,1], IRanges::IRanges(collected.up[,2] + 1, collected.up[,3])), 
    GenomicRanges::GRanges(collected.down[,1], IRanges::IRanges(collected.down[,2] + 1, collected.down[,3]))
  )
  
  curtab_pval_fdr <- rbind(collected.up[,c(8, 9)], collected.down[,c(8, 9)])
  curtab_logfc <- c(rep(1, nrow(collected.up)), rep(-1, nrow(collected.down)))
  curtab <- cbind(logFC = curtab_logfc, PValue = curtab_pval_fdr[,1], FDR = curtab_pval_fdr[,2])
  colnames(curtab) <- c("logFC", "PValue", "FDR") 
  
  results_pepr_initial <- data.frame(
    index = seq_along(all.collected),
    chr = as.character(GenomicRanges::seqnames(all.collected)),
    start = GenomicRanges::start(all.collected),
    end = GenomicRanges::end(all.collected),
    Pvalues = curtab[,"PValue"], 
    FDR = curtab[,"FDR"],
    logFC = curtab[,"logFC"]
  )
  
  # 文件名中的标识部分
  file_id_string <- paste0(analysis_type_label, "_tag", iteration_tag, "_pw", peak_width)
  
  if (save_initial_pepr_results) {
    fn <- paste0("initial_pepr_results_", file_id_string, ".csv")
    write.csv(results_pepr_initial, file.path(getwd(), fn), row.names = FALSE)
    cat("Saved initial PePr results to:", file.path(getwd(), fn), "\n")
  }
  
  # --- 5. 调用 find_overlap ---
  cat("Calling find_overlap function...\n")
  pepr_sites_for_findoverlap <- results_pepr_initial 
  
  if (!exists("find_overlap") || !is.function(find_overlap)) {
    stop("Function 'find_overlap' is not defined. Please source ChIPtest_source.R or define it.")
  }
  pos_arg2_for_findoverlap <- if(is_tf) pos.1_data else pos.3_data_for_findoverlap
  
  sites_overlap_output <- find_overlap(
    is.tf = is_tf, 
    pos.1 = pos.1_data, 
    pos.3 = pos_arg2_for_findoverlap, 
    peak_width = peak_width, 
    lfile = lfile, 
    csaw_sites = pepr_sites_for_findoverlap 
  )
  
  if (save_sites_overlap_csv) {
    fn <- paste0("sites_overlap_results_before_filter_", file_id_string, ".csv")
    write.csv(sites_overlap_output, file.path(getwd(), fn), row.names = FALSE)
    cat("Saved sites_overlap results (before dplyr filter) to:", file.path(getwd(), fn), "\n")
  }
  
  # --- 6. dplyr 筛选结果 ---
  cat("Filtering results using dplyr...\n")
  required_cols_for_dplyr <- c("FDR", "overlap", "sim_start") # 根据您的脚本，find_overlap的输出应包含这些
  if (!all(required_cols_for_dplyr %in% names(sites_overlap_output))) {
    stop(paste("Output from find_overlap is missing one or more required columns for dplyr filtering:", 
               paste(required_cols_for_dplyr[!required_cols_for_dplyr %in% names(sites_overlap_output)], collapse=", ")))
  }
  if(!is.numeric(sites_overlap_output$FDR)) {
    cat("Warning: FDR column from find_overlap is not numeric. Attempting conversion.\n")
    sites_overlap_output$FDR <- as.numeric(sites_overlap_output$FDR)
  }
  if(!is.numeric(sites_overlap_output$sim_start)) sites_overlap_output$sim_start <- as.numeric(sites_overlap_output$sim_start)
  if(!is.numeric(sites_overlap_output$overlap)) sites_overlap_output$overlap <- as.numeric(sites_overlap_output$overlap)
  
  sites_overlap_output$sim_start[is.na(sites_overlap_output$sim_start)] <- 0
  sites_overlap_output$overlap[is.na(sites_overlap_output$overlap)] <- 0
  
  filtered_results_pepr <- dplyr::bind_rows(
    sites_overlap_output %>% 
      dplyr::filter(overlap == 1 & sim_start != 0) %>%
      dplyr::group_by(sim_start) %>%
      dplyr::slice_min(FDR, with_ties = FALSE, n=1) %>% 
      dplyr::ungroup(),
    sites_overlap_output %>% 
      dplyr::filter(overlap == 0 | sim_start == 0)
  )
  
  # --- 7. 计算评估指标 ---
  cat("Calculating FDR, Recall, and AUC/ROC metrics...\n")
  if (!is.numeric(filtered_results_pepr$overlap)) {
    filtered_results_pepr$overlap <- as.numeric(filtered_results_pepr$overlap)
  }
  valid_for_roc <- !is.na(filtered_results_pepr$FDR) & !is.na(filtered_results_pepr$overlap)
  
  if (!exists("calculate_fdr_recall") || !is.function(calculate_fdr_recall)) {
    stop("Function 'calculate_fdr_recall' is not defined. Please source ChIPtest_source.R or define it.")
  }
  fdr_recall_metrics <- calculate_fdr_recall(filtered_results_pepr$FDR, filtered_results_pepr$overlap, thresholds=c(0.01, 0.05, 0.1))
  
  auc_value <- NA
  roc_curve_data <- data.frame(TPR=NA, FPR=NA, Threshold=NA)
  
  if(sum(valid_for_roc) < 2 || length(unique(filtered_results_pepr$overlap[valid_for_roc])) < 2) {
    message(paste0("Not enough data points or classes to calculate ROC/AUC (Valid N=", sum(valid_for_roc), 
                   ", Unique overlaps=", length(unique(filtered_results_pepr$overlap[valid_for_roc])), "). Skipping ROC calculation."))
  } else {
    pred <- ROCR::prediction(1 - filtered_results_pepr$FDR[valid_for_roc], filtered_results_pepr$overlap[valid_for_roc])
    perf <- ROCR::performance(pred, "tpr", "fpr")
    roc_curve_data <- data.frame(
      TPR = perf@y.values[[1]],
      FPR = perf@x.values[[1]],
      Threshold = perf@alpha.values[[1]]
    )
    auc_value <- ROCR::performance(pred, "auc")@y.values[[1]]
  }
  fdr_recall_metrics$auc_value <- auc_value
  
  # --- 8. 保存结果 ---
  cat("Saving final results and metrics...\n")
  if (save_final_filtered_results_csv) {
    fn <- paste0("final_filtered_results_PePr_", file_id_string, ".csv")
    write.csv(filtered_results_pepr, file.path(getwd(), fn), row.names = FALSE)
    cat("Saved final filtered PePr results to:", file.path(getwd(), fn), "\n")
  }
  
  pipeline_run_time_secs <- as.numeric(difftime(Sys.time(), pipeline_start_time, units = "secs"))
  # iteration_tag 替换了原来 metrics 表中的 iteration 列，但含义略有不同
  fdr_recall_metrics$iteration_tag <- iteration_tag 
  fdr_recall_metrics$input_dataset <- dataset_name_for_output # 添加数据集信息
  fdr_recall_metrics$peak_width <- peak_width
  fdr_recall_metrics$is_tf_analysis <- is_tf
  fdr_recall_metrics$runtime_seconds <- round(pipeline_run_time_secs, 2)
  
  fdr_recall_metrics <- fdr_recall_metrics[, c("iteration_tag", "input_dataset", "peak_width", "is_tf_analysis", "runtime_seconds", 
                                               names(fdr_recall_metrics)[!names(fdr_recall_metrics) %in% 
                                                                           c("iteration_tag", "input_dataset", "peak_width", "is_tf_analysis", "runtime_seconds")])]
  
  if (save_metrics_csv) {
    fn <- paste0("metrics_PePr_", file_id_string, ".csv")
    write.csv(fdr_recall_metrics, file.path(getwd(), fn), row.names = FALSE)
    cat("Saved metrics to:", file.path(getwd(), fn), "\n")
  }
  
  if (save_roc_data_csv && nrow(roc_curve_data) > 0 && !all(is.na(roc_curve_data$TPR))) {
    fn <- paste0("roc_data_PePr_", file_id_string, ".csv")
    write.csv(roc_curve_data, file.path(getwd(), fn), row.names = FALSE)
    cat("Saved ROC curve data to:", file.path(getwd(), fn), "\n")
  }
  
  # --- 9. 返回结果 ---
  cat(paste("PePr Pipeline for input:", input_data_dir, ", iteration_tag:", iteration_tag, "completed at", Sys.time(), "\n"))
  cat("Total pipeline runtime:", format(Sys.time() - pipeline_start_time), "\n\n")
  
  return(list(
    final_results = filtered_results_pepr,
    metrics = fdr_recall_metrics,
    roc_data = roc_curve_data,
    output_directory = current_output_path, # 返回实际使用的输出目录
    path_final_results_csv = if(save_final_filtered_results_csv) file.path(current_output_path, paste0("final_filtered_results_PePr_", file_id_string, ".csv")) else NA,
    path_metrics_csv = if(save_metrics_csv) file.path(current_output_path, paste0("metrics_PePr_", file_id_string, ".csv")) else NA,
    path_roc_data_csv = if(save_roc_data_csv && nrow(roc_curve_data) > 0 && !all(is.na(roc_curve_data$TPR))) file.path(current_output_path, paste0("roc_data_PePr_", file_id_string, ".csv")) else NA,
    pepr_native_output_dir = pepr_native_output_path, 
    pipeline_runtime_formatted = format(Sys.time() - pipeline_start_time),
    pipeline_runtime_seconds = pipeline_run_time_secs
  ))
}

#######################################################################################
run_macs3_diffbind_pipeline <- function(
    data_directory,
    output_directory_base,
    iteration_tag,
    count_bin=15,
    is_tf = FALSE,
    peak_width,
    source_file_path = "/projects/gqilab/DAESC_GPU/data/simulation/scripts/ChIPtest_source.R",
    gsize = 300678703,
    fraglen = 100,
    num_samples = 4, # Assuming 4 samples as per your script
    grouping_vector = c("A", "A", "B", "B"), # Default for 2 vs 2 comparison
    min_overlap_dba = 2,
    fdr_thresholds = c(0.05)
) {
  
  # --- 0. Input Validation and Setup ---
  if (!dir.exists(data_directory)) {
    stop("Data directory does not exist: ", data_directory)
  }
  if (!file.exists(source_file_path)) {
    stop("Source file does not exist: ", source_file_path)
  }
  if (length(grouping_vector) != num_samples) {
    stop("Length of grouping_vector must match num_samples.")
  }
  
  # Load required packages silently
  if (!requireNamespace("DiffBind", quietly = TRUE)) {
    stop("Package 'DiffBind' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("Package 'ROCR' is required but not installed.")
  }
  
  # Clear environment (optional, consider if this function is part of a larger script)
  # current_env <- ls(envir = .GlobalEnv) # Get current global variables
  #rm(list = ls(envir = parent.frame())) # Clear calling environment if desired, otherwise remove this line
  #gc()
  
  # Source helper functions
  message("Sourcing helper functions from: ", source_file_path)
  source(source_file_path)
  
  # --- 1. Define Paths and Variables based on inputs ---
  message("Setting up parameters for iteration_tag: ", iteration_tag)
  
  if (is_tf) {
    all.fix <- "tfx" # tfx was in your original example for tf
    # path.fix is not strictly needed if data_directory is the full specific path
  } else {
    all.fix <- "hist"
  }
  
  # Construct the specific output directory for this run
  run_output_dir <- file.path(output_directory_base, paste0("analysis_run_", iteration_tag))
  if (!dir.exists(run_output_dir)) {
    dir.create(run_output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  message("Results will be saved in: ", run_output_dir)
  
  original_wd <- getwd()
  setwd(run_output_dir)
  on.exit(setwd(original_wd), add = TRUE) # Ensure WD is reset even on error
  
  # Define BAM file paths
  bam.files <- character(num_samples)
  for (i in 1:num_samples) {
    bam.files[i] <- file.path(data_directory, paste0(all.fix, "_out_", i, ".bam"))
    if (!file.exists(bam.files[i])) {
      stop("BAM file not found: ", bam.files[i])
    }
  }
  message("Using BAM files: ", paste(basename(bam.files), collapse = ", "))
  
  # Directory for MACS3 peaks (relative to run_output_dir)
  macs3_peak_subdir <- paste0(all.fix, "_peaks_iteration_", iteration_tag)
  peakdir_full_path <- file.path(run_output_dir, macs3_peak_subdir)
  if (!dir.exists(peakdir_full_path)) {
    dir.create(peakdir_full_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Prefixes for MACS3 output files
  prefixes <- character(num_samples)
  for (i in 1:num_samples) {
    prefixes[i] <- paste0(all.fix, "_out_", i)
  }
  
  # Read MACS3 parameters (assuming readParam is from your source file)
  # If readParam is not defined in your source file, you might need to define xparam directly
  # xparam <- readParam(dedup = FALSE) # Example from your script
  # For now, assuming it's not strictly needed for run_MACS3_server if not used by it
  # If your run_MACS3_server or other functions require 'xparam', ensure it's loaded or passed.
  
  # Log file and position information
  lfile <- file.path(data_directory, paste0(all.fix, "_log.txt"))
  if (!file.exists(lfile)) {
    warning("Log file not found: ", lfile) # Changed to warning, pipeline might still run
  }
  
  message("Loading RData position vectors...")
  pos1_path <- file.path(data_directory, "pos1_vector.RData")
  if (file.exists(pos1_path)) load(pos1_path, envir = .GlobalEnv) else stop("File not found: ", pos1_path)
  
  if (!is_tf) {
    pos2_path <- file.path(data_directory, "pos2_vector.RData")
    pos3_path <- file.path(data_directory, "pos3_vector.RData")
    if (file.exists(pos2_path)) load(pos2_path, envir = .GlobalEnv) else stop("File not found: ", pos2_path)
    if (file.exists(pos3_path)) load(pos3_path, envir = .GlobalEnv) else stop("File not found: ", pos3_path)
  }
  message("Position vectors loaded.")
  
  
  # --- 2. Running MACS3 ---
  message("Starting MACS3 peak calling...")
  all.peakfiles <- vector("list", length(bam.files))
  for (x in 1:length(bam.files)) {
    # oprefix for MACS3 output: peakdir_full_path/prefix[x]
    # e.g. /path/to/output_base/analysis_run_X/hist_peaks_iteration_X/hist_out_1
    oprefix <- file.path(peakdir_full_path, prefixes[x])
    message("Running MACS3 for: ", basename(bam.files[x]), " output prefix: ", oprefix)
    
    # Assuming run_MACS3_server is defined in your source file
    # and handles output to the directory specified in oprefix.
    run_MACS3_server(
      file = bam.files[x], # Ensure your function takes full path
      outprefix = oprefix,     # Ensure your function takes full path prefix
      fraglen = fraglen,
      gsize = gsize,
      cmd.only = FALSE,
      is.tf = is_tf
    )
    
    # MACS3 output file name: oprefix_peaks.xls or oprefix_peaks.narrowPeak/broadPeak
    # Your script used _peaks.xls. Adjust if run_MACS3_server produces different common extensions.
    # Common MACS3 outputs are .narrowPeak for sharp peaks and .broadPeak for broad.
    # If run_MACS3_server generates .xls, this is fine.
    # Check what run_MACS3_server actually produces. For DiffBind, BED or specific MACS formats are common.
    # Let's assume it's _peaks.xls as per your script.
    expected_peak_file <- paste0(oprefix, "_peaks.xls")
    if (!file.exists(expected_peak_file)) {
      # Fallback to common MACS3 extensions if .xls is not found
      if (is_tf) {
        expected_peak_file <- paste0(oprefix, "_peaks.narrowPeak")
      } else {
        expected_peak_file <- paste0(oprefix, "_peaks.broadPeak")
      }
      if (!file.exists(expected_peak_file)) {
        # Try the other one just in case
        if (is_tf) {
          expected_peak_file_alt <- paste0(oprefix, "_peaks.broadPeak")
        } else {
          expected_peak_file_alt <- paste0(oprefix, "_peaks.narrowPeak")
        }
        if (file.exists(expected_peak_file_alt)) {
          expected_peak_file <- expected_peak_file_alt
        } else {
          stop("MACS3 output peak file not found after run: ", expected_peak_file,
               " nor ", if(exists("expected_peak_file_alt")) expected_peak_file_alt else "")
        }
      }
    }
    all.peakfiles[[x]] <- expected_peak_file
    message("MACS3 peak file: ", all.peakfiles[[x]])
  }
  macs.peakfiles <- unlist(all.peakfiles)
  
  # --- 3. DiffBind Analysis ---
  message("Starting DiffBind analysis...")
  sample_sheet_df <- data.frame(
    SampleID = prefixes, # Should match sample names from MACS3 runs
    Condition = grouping_vector,
    bamReads = bam.files, # Full paths to BAM files
    Peaks = macs.peakfiles, # Full paths to MACS3 peak files
    PeakCaller = "macs" # Or "narrow" / "broad" depending on MACS3 output used
  )
  message("Sample sheet for DiffBind:")
  print(sample_sheet_df)
  
  current_dba <- DiffBind::dba(sampleSheet = sample_sheet_df, minOverlap = min_overlap_dba)
  message("Counting reads...")
  current_dba <- DiffBind::dba.count(current_dba, fragmentSize = fraglen, bRemoveDuplicates = FALSE)
  
  # Define contrast based on unique values in grouping_vector
  # Assuming a two-group comparison. If more complex, this needs adjustment.
  unique_conditions <- unique(grouping_vector)
  if (length(unique_conditions) != 2) {
    stop("DiffBind contrast setup currently supports only two unique conditions in grouping_vector.")
  }
  group1_mask <- current_dba$masks[[unique_conditions[1]]]
  group2_mask <- current_dba$masks[[unique_conditions[2]]]
  
  message("Setting up contrast: ", unique_conditions[1], " vs ", unique_conditions[2])
  current_dba <- DiffBind::dba.contrast(
    current_dba,
    group1 = group1_mask, # Use the actual mask for the condition
    group2 = group2_mask, # Use the actual mask for the condition
    name1 = unique_conditions[1],
    name2 = unique_conditions[2],
    minMembers = 2 # Default, ensure at least 2 replicates in smallest group for contrast
  )
  
  message("Running differential analysis (all methods)...")
  # DBA_ALL_METHODS can be time-consuming. DBA_EDGER is often sufficient.
  # Consider making method a parameter if needed.
  current_dba <- DiffBind::dba.analyze(current_dba, method = DiffBind::DBA_ALL_METHODS)
  # current_dba <- DiffBind::dba.analyze(current_dba, method = DiffBind::DBA_EDGER) # Faster alternative
  
  message("Generating report for edgeR results...")
  # Ensure method matches one analyzed, e.g., DBA_EDGER
  # dba.show(current_dba, current_dba$masks$Consensus) will show available methods.
  # We'll use edgeR as in your script.
  if (! DiffBind::DBA_EDGER %in% current_dba$config$AnalysisMethod) {
    warning("DBA_EDGER was not part of the analysis methods run. Attempting to get report anyway.")
  }
  out_edger <- DiffBind::dba.report(current_dba, method = DiffBind::DBA_EDGER, th = 1, DataType = DBA_DATA_FRAME)
  
  if (is.null(out_edger) || nrow(out_edger) == 0) {
    stop("DiffBind dba.report for edgeR returned no results. Check contrast and analysis steps.")
  }
  
  # Convert to data.frame, GRanges to data.frame for p-value and FDR
  # dba.report with DataType=DBA_DATA_FRAME returns a data.frame directly
  results_macs3_diffbind <- data.frame(
    start = as.integer(out_edger$Start),
    end = as.integer(out_edger$End),
    p_value = as.numeric(out_edger$`p-value`),
    FDR = as.numeric(out_edger$FDR)
  )
  results_macs3_diffbind$index <- 1:nrow(results_macs3_diffbind)
  
  message("DiffBind analysis complete. Number of sites reported: ", nrow(results_macs3_diffbind))
  
  # --- 4. Overlap Analysis and Performance Metrics ---
  message("Starting overlap analysis with simulated sites...")
  macs3_sites_for_overlap <- results_macs3_diffbind[, c("index", "start", "end", "p_value", "FDR")]
  colnames(macs3_sites_for_overlap) <- c("index", "start", "end", "Pvalues", "FDR") # Match expected by find_overlap
  
  # Ensure pos.1, pos.3 (or pos.1, pos.1 for TF) are available in the global environment
  # from the RData files loaded earlier.
  if (is_tf) {
    if (!exists("pos.1")) stop("pos.1 not found in global environment. Check RData loading.")
    sites_overlap <- find_overlap(
      is.tf = TRUE, # Pass the parameter
      pos.1 = pos.1, # Argument name as per typical function definition
      pos.3 = pos.1, # Argument name as per typical function definition
      peak_width = peak_width,
      lfile = lfile,
      csaw_sites = macs3_sites_for_overlap
    )
  } else {
    if (!exists("pos.1") || !exists("pos.3")) {
      stop("pos.1 or pos.3 not found in global environment. Check RData loading.")
    }
    sites_overlap <- find_overlap(
      is.tf = FALSE, # Pass the parameter
      pos.1 = pos.1, # Argument name as per typical function definition
      pos.3 = pos.3, # Argument name as per typical function definition
      peak_width = peak_width,
      lfile = lfile,
      csaw_sites = macs3_sites_for_overlap
    )
  }
  message("Overlap analysis complete.")
  
  # Process results as per your script
  sites_overlap[is.na(sim_start), `:=` (
    Pvalues = 1,
    FDR = 1,
    overlap_true = 0
  )]
  final_results_df <- sites_overlap
  
  if (nrow(final_results_df) > 0 && "sim_start" %in% colnames(final_results_df)) {
    final_results_df <- dplyr::bind_rows(
      final_results_df %>%
        dplyr::filter(!is.na(sim_start)) %>%   # 有匹配
        dplyr::group_by(sim_start) %>%
        dplyr::slice_min(FDR, with_ties = FALSE, n = 1) %>%
        dplyr::ungroup(),
      
      final_results_df %>%
        dplyr::filter(is.na(sim_start))        # 无匹配
    )
  } else {
    message("Skipping dplyr processing for final_results_df as it's empty or missing column sim_start.")
  }
  
  
  
  message("Calculating FDR/Recall and AUC...")
  # Ensure 'overlap' column exists for these calculations
  if (nrow(final_results_df) > 0) {
    fdr_recall_metrics <- calculate_fdr_recall(
      final_results_df$FDR,
      final_results_df$overlap_true,
      thresholds = fdr_thresholds
    )
    
    # ROCR predictions: Use (1 - FDR) for scores where higher is better
    # Ensure no NA values in FDR for prediction, replace NAs if necessary, e.g. with 1
    fdr_for_pred <- final_results_df$FDR
    fdr_for_pred[is.na(fdr_for_pred)] <- 1.0 
    
    pred_obj <- ROCR::prediction(1 - fdr_for_pred, final_results_df$overlap_true)
    perf_obj <- ROCR::performance(pred_obj, "tpr", "fpr")
    auc_value <- ROCR::performance(pred_obj, "auc")@y.values[[1]]
    fdr_recall_metrics$auc_value <- auc_value
    
    roc_data <- data.frame(
      TPR = perf_obj@y.values[[1]],
      FPR = perf_obj@x.values[[1]],
      Threshold = perf_obj@alpha.values[[1]]
    )
  } else {
    message("Skipping FDR/Recall and AUC calculation as 'overlap' or 'FDR' column is missing or data frame is empty.")
    fdr_recall_metrics <- data.frame(threshold = fdr_thresholds, FDR = NA, Recall = NA, TP = NA, FP = NA, auc_value = NA)
    roc_data <- data.frame(TPR = NA, FPR = NA, Threshold = NA)
    auc_value <- NA
  }
  message("Performance metrics calculated.")
  
  # --- 5. Save Results ---
  results_filename <- paste0(all.fix, "_results_DiffBind_Homer_", iteration_tag, ".csv")
  metrics_filename <- paste0(all.fix, "_metrics_DiffBind_Homer_", iteration_tag, ".csv")
  roc_filename <- paste0(all.fix, "_roc_data_DiffBind_Homer_", iteration_tag, ".csv")
  dba_object_filename <- paste0(all.fix, "_dba_object_", iteration_tag, ".RData")
  
  
  write.csv(final_results_df, file = file.path(run_output_dir, results_filename), row.names = FALSE)
  message("Saved detailed results to: ", file.path(run_output_dir, results_filename))
  
  write.csv(fdr_recall_metrics, file = file.path(run_output_dir, metrics_filename), row.names = FALSE)
  message("Saved performance metrics to: ", file.path(run_output_dir, metrics_filename))
  
  write.csv(roc_data, file = file.path(run_output_dir, roc_filename), row.names = FALSE)
  message("Saved ROC data to: ", file.path(run_output_dir, roc_filename))
  
  #save(current_dba, file = file.path(run_output_dir, dba_object_filename))
  #message("Saved DiffBind DBA object to: ", file.path(run_output_dir, dba_object_filename))
  
  
  # --- 6. Return ---
  message("Pipeline run complete for iteration_tag: ", iteration_tag)
  return(
    list(
      results_df = final_results_df,
      metrics_df = fdr_recall_metrics,
      roc_data = roc_data,
      auc_value = auc_value,
      output_path = run_output_dir,
      dba_object_path = file.path(run_output_dir, dba_object_filename)
    )
  )
}



########################################################################################
# MACS 3 + ChIPtest 
run_macs3_chiptest_pipeline <- function(
    data_base_dir,
    output_base_dir,
    run_tag,              # 例如 "tf_setting4_chiptest_run1"
    count_bin=15,
    peak_width=500,
    is_tf = FALSE,
    source_file_path = "/projects/gqilab/DAESC_GPU/data/simulation/scripts/ChIPtest_source.R",
    gsize = 300678703,
    fraglen = 100,
    num_samples = 4,
    grouping_vector = c("A", "A", "B", "B"), # 用于 windowCounts 后续的 conA/conB 计算
    num_cores_mclapply = 8,
    num_cores_finetune = 5,
    chiptest_band_values = c(5,10,30,60),
    chiptest_quantile_values = c(-1, 0.0001, 0.01, 0.02, 0.1, 0.2, 0.4, 0.5, 0.99, 1.2, 1.5, 1.8, 2),
    chiptest_var_est = 1,
    chiptest_var_thred = 0.01
) {
  
  # --- 0. 输入校验与环境设置 ---
  if (!dir.exists(data_base_dir)) stop("Data base directory does not exist: ", data_base_dir)
  if (!dir.exists(output_base_dir)) stop("Output base directory does not exist: ", output_base_dir)
  if (!file.exists(source_file_path)) stop("Source file does not exist: ", source_file_path)
  if (length(grouping_vector) != num_samples) stop("Length of grouping_vector must match num_samples.")
  
  # 加载必要的包
  packages_to_load <- c("dplyr", "ROCR") # ggplot2, data.table, parallel 可能在source_file_path中加载或被其函数使用
  for (pkg in packages_to_load) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package '", pkg, "' is required but not installed.", sep=""))
    }
  }
  # 确保csaw或其核心函数被加载，因为windowCounts来自那里
  if (!exists("windowCounts", mode="function") && !requireNamespace("csaw", quietly = TRUE)) {
    warning("Function 'windowCounts' not found and 'csaw' package not explicitly loaded. Assuming it's in source_file_path.")
  }
  
  
  # 清理环境 (可选)
  # rm(list = ls(envir = parent.frame()))
  # gc()
  
  message("Sourcing helper functions from: ", source_file_path)
  source(source_file_path) # 预期包含 run_MACS3_server, readParam, windowCounts, process_peak_*, Finetune_ChIPtest, etc.
  
  # --- 1. 定义路径和变量 ---
  message("Setting up parameters for run_tag: ", run_tag)
  
  path.fix <- if (is_tf) "tf" else "hm"
  all.fix <- if (is_tf) "tfx" else "hist"
  
  # 输入数据目录
  # current_data_dir <- file.path(data_base_dir, paste0(path.fix, "_setting_", iteration_identifier))
  current_data_dir <- data_base_dir
  if (!dir.exists(current_data_dir)) {
    stop("Constructed data directory does not exist: ", current_data_dir)
  }
  message("Input data directory: ", current_data_dir)
  
  # 输出结果目录
  current_analysis_output_dir <- file.path(output_base_dir, paste0("macs_chiptest_analysis_", run_tag))
  if (!dir.exists(current_analysis_output_dir)) {
    dir.create(current_analysis_output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  message("Analysis results will be saved in: ", current_analysis_output_dir)
  
  original_wd <- getwd()
  setwd(current_analysis_output_dir) # 后续所有输出都基于此目录
  on.exit(setwd(original_wd), add = TRUE)
  
  # 定义BAM文件、日志文件、位置信息文件的绝对路径
  bam.files <- character(num_samples)
  for (i in 1:num_samples) {
    bam.files[i] <- file.path(current_data_dir, paste0(all.fix, "_out_", i, ".bam"))
    if (!file.exists(bam.files[i])) stop("BAM file not found: ", bam.files[i])
  }
  message("Using BAM files: ", paste(basename(bam.files), collapse = ", "))
  
  lfile <- file.path(current_data_dir, paste0(all.fix, "_log.txt"))
  if (!file.exists(lfile)) warning("Log file not found: ", lfile)
  
  message("Loading RData position vectors...")
  pos1_path <- file.path(current_data_dir, "pos1_vector.RData")
  if (file.exists(pos1_path)) load(pos1_path, envir = .GlobalEnv) else stop("File not found: ", pos1_path)
  
  if (!is_tf) {
    pos2_path <- file.path(current_data_dir, "pos2_vector.RData")
    pos3_path <- file.path(current_data_dir, "pos3_vector.RData")
    if (file.exists(pos2_path)) load(pos2_path, envir = .GlobalEnv) else stop("File not found: ", pos2_path)
    if (file.exists(pos3_path)) load(pos3_path, envir = .GlobalEnv) else stop("File not found: ", pos3_path)
  }
  message("Position vectors loaded.")
  
  # MACS3 参数 (xparam 来自 ChIPtest_source.R 的 readParam)
  if (!exists("readParam", mode="function")) stop("Function 'readParam' must be defined in source_file_path.")
  xparam <- readParam(dedup = FALSE)
  
  
  # --- 2. MACS3 Peak Calling ---
  message("Starting MACS3 peak calling...")
  # MACS3的peak输出到 current_analysis_output_dir 下的子目录
  macs3_peak_subdir_name <- paste0(all.fix, "_macs3_peaks") #相对路径名
  macs3_peak_output_path <- file.path(current_analysis_output_dir, macs3_peak_subdir_name) #绝对路径
  if (!dir.exists(macs3_peak_output_path)) {
    dir.create(macs3_peak_output_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  all.peakfiles_macs <- character(length(bam.files)) # Storing absolute paths
  macs_prefixes <- paste0(all.fix, "_out_", 1:length(bam.files))
  
  for (x in 1:length(bam.files)) {
    # oprefix for MACS3 output will be an absolute path
    oprefix_macs <- file.path(macs3_peak_output_path, macs_prefixes[x])
    message("Running MACS3 for: ", basename(bam.files[x]), ", output prefix: ", oprefix_macs)
    
    run_MACS3_server( # 确保 run_MACS3_server 接受这些参数名
      file = bam.files[x],       # 绝对路径的BAM
      outprefix = oprefix_macs,  # 绝对路径的输出前缀
      fraglen = fraglen,
      gsize = gsize,
      cmd.only = FALSE,
      is.tf = is_tf
    )
    # 检查MACS3输出文件，使用与之前 DiffBind pipeline 中类似的逻辑
    # 您的原脚本直接使用 _peaks.xls
    expected_peak_file <- paste0(oprefix_macs, "_peaks.xls")
    if (!file.exists(expected_peak_file)) {
      alt_suffix <- if (is_tf) "_peaks.narrowPeak" else "_peaks.broadPeak"
      expected_peak_file_alt_1 <- paste0(oprefix_macs, alt_suffix)
      alt_suffix_2 <- if (is_tf) "_peaks.broadPeak" else "_peaks.narrowPeak" # a less common alternative
      expected_peak_file_alt_2 <- paste0(oprefix_macs, alt_suffix_2)
      
      if(file.exists(expected_peak_file_alt_1)){
        expected_peak_file <- expected_peak_file_alt_1
      } else if (file.exists(expected_peak_file_alt_2)){
        expected_peak_file <- expected_peak_file_alt_2
      } else {
        stop(paste("MACS3 output peak file not found after run. Looked for .xls,", alt_suffix, ",", alt_suffix_2, "with prefix:", oprefix_macs))
      }
    }
    all.peakfiles_macs[x] <- expected_peak_file
    message("MACS3 peak file: ", all.peakfiles_macs[x])
  }
  
  
  # --- 3. ChIPtest 分析准备 ---
  message("Preparing data for ChIPtest analysis...")
  bed_data_list <- lapply(all.peakfiles_macs, function(file_path) {
    # 确保MACS3的 .xls 文件有正确的列数和类型，或者调整colClasses
    tryCatch(
      read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#",
                 colClasses = c("character", "integer", "integer", "integer", "numeric", "numeric", "numeric", "numeric", "character")),
      error = function(e) {
        message("Error reading MACS3 XLS file: ", file_path, ". Error: ", e$message)
        message("Attempting to read without strict colClasses, skipping header if it looks problematic.")
        # Fallback: try to read first few lines to guess header, then read data
        header_line <- try(readLines(file_path, n=20), silent=TRUE)
        is_header_present <- FALSE
        if(!inherits(header_line, "try-error")){
          # MACS2/3 XLS often has comments then header like 'chr   start   end ...'
          header_row_idx <- grep("^chr\tstart\tend", header_line, ignore.case = TRUE)
          if(length(header_row_idx) > 0){
            is_header_present <- TRUE
            # Read data skipping lines before the identified header
            return(read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", skip = header_row_idx[1]-1))
          }
        }
        # If still failing, read without header and assume standard BED-like columns + MACS specifics
        return(read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#"))
      }
    )
  })
  
  # 确保所有读入的data.frame都有统一的列名结构，特别是前几列
  standard_macs_xls_names <- c("chr", "start", "end", "length", "pileup", "log10_pval", "fold_enrichment", "log10_qval", "name")
  
  bed_data_list <- lapply(bed_data_list, function(df) {
    if (ncol(df) >= length(standard_macs_xls_names)) {
      colnames(df)[1:length(standard_macs_xls_names)] <- standard_macs_xls_names
    } else if (ncol(df) >= 3) { # Minimal BED-like
      colnames(df)[1:3] <- c("chr", "start", "end")
    } else {
      stop("MACS3 peak file does not have enough columns.")
    }
    return(df)
  })
  
  
  all_sites <- do.call(rbind, lapply(bed_data_list, function(df) { df[, c("start", "end")] }))
  all_sites$length <- all_sites$end - all_sites$start
  all_sites$medium <- (all_sites$end + all_sites$start) / 2
  
  all_sites <- all_sites %>%
    dplyr::arrange(medium) %>%
    dplyr::mutate(group = cumsum(c(0, diff(medium) > 200))) %>% # 200bp是合并邻近peaks的阈值
    dplyr::group_by(group) %>%
    dplyr::mutate(medium_avg = mean(medium)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-group)
  
  count_bin <- count_bin
  window_width <- 120 # 这个在原脚本中对于TF和HM都是120
  
  message("Performing window counts...")
  # 确保 windowCounts, start, rowRanges, assay 函数可用 (通常来自 csaw/SummarizedExperiment)
  count_seq <- windowCounts(bam.files, width=count_bin, spacing=count_bin, ext=fraglen, filter=0, bin=TRUE, param=xparam) # ext应为fraglen
  counts_matrix <- SummarizedExperiment::assay(count_seq)
  start_positions <- GenomicRanges::start(SummarizedExperiment::rowRanges(count_seq))
  end_positions <- GenomicRanges::end(SummarizedExperiment::rowRanges(count_seq))
  count_bin2 <- data.frame(start = start_positions, end = end_positions, counts_matrix)
  
  # 检查counts_matrix的列名。原脚本用X1,X2...
  # csaw::windowCounts 通常会使用bam.files的basename作为列名
  sample_names_from_counts <- colnames(counts_matrix)
  if(ncol(counts_matrix) != num_samples) stop("Mismatch between num_samples and columns in counts_matrix from windowCounts")
  
  # 假设 grouping_vector 的顺序与 bam.files 和 counts_matrix 列的顺序一致
  counts_condA <- counts_matrix[, grouping_vector == unique(grouping_vector)[1], drop=FALSE]
  counts_condB <- counts_matrix[, grouping_vector == unique(grouping_vector)[2], drop=FALSE]
  
  count_bin2$conA <- rowMeans(counts_condA)
  count_bin2$conB <- rowMeans(counts_condB)
  
  message("Extracting count sequences for ChIPtest...")
  potential_peaks_indices <- unique(all_sites$medium_avg) / count_bin # 转换为bin的索引
  
  # 确保 process_peak_conA/B 在 ChIPtest_source.R 中定义
  # 找到 conA 和 conB 在 count_bin2 中的列索引
  col_idx_conA <- which(colnames(count_bin2) == "conA")
  col_idx_conB <- which(colnames(count_bin2) == "conB")
  if(length(col_idx_conA)==0 || length(col_idx_conB)==0) stop("conA or conB columns not found in count_bin2")
  
  # 使用 parallel::mclapply (如果 ChIPtest_source.R 未加载 parallel, 需要在这里加载)
  apply_fun <- if (.Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)) {
    function(...) parallel::mclapply(..., mc.cores = num_cores_mclapply)
  } else {
    message("Using sequential lapply as mclapply is not available or num_cores_mclapply <= 1.")
    lapply
  }
  
  data_list_conA <- apply_fun(potential_peaks_indices, process_peak_conA, df = count_bin2[,c(1,2,col_idx_conA)], width = window_width)
  data_conA <- do.call(rbind, lapply(data_list_conA, t))
  
  data_list_conB <- apply_fun(potential_peaks_indices, process_peak_conB, df = count_bin2[,c(1,2,col_idx_conB)], width = window_width)
  data_conB <- do.call(rbind, lapply(data_list_conB, t))
  
  
  # --- 4. ChIPtest - Null Hypothesis Tuning ---
  message("Running Finetune_ChIPtest for null hypothesis...")
  finetune_results <- Finetune_ChIPtest(
    potential_peaks = potential_peaks_indices,
    count_bin2 = count_bin2, # 确保 Finetune_ChIPtest 能处理这个格式
    width = window_width,
    band_values = chiptest_band_values,
    quantile_values = chiptest_quantile_values,
    my_var.est = chiptest_var_est,
    my_var.thred = chiptest_var_thred,
    num_cores = num_cores_finetune
  )
  
  finetune_results$TS_kn_KL_mean <- rowMeans(finetune_results[, c("TS_kn_KL_12", "TS_kn_KL_34")], na.rm = TRUE)
  finetune_results$Deql_KL_mean <- rowMeans(finetune_results[, c("Deql_KL_12", "Deql_KL_34")], na.rm = TRUE)
  finetune_results$Dnun_KL_mean <- rowMeans(finetune_results[, c("Dnun_KL_12", "Dnun_KL_34")], na.rm = TRUE)
  
  quantile_TSkn_min <- finetune_results$Quantile[which.min(finetune_results$TS_kn_KL_mean)]
  band_TSkn_min <- finetune_results$Band[which.min(finetune_results$TS_kn_KL_mean)]
  quantile_Deql_min <- finetune_results$Quantile[which.min(finetune_results$Deql_KL_mean)]
  quantile_Dnun_min <- finetune_results$Quantile[which.min(finetune_results$Dnun_KL_mean)]
  
  message("Plotting finetuned null hypothesis (plot will be saved if finetune function does so)...")
  # 确保 finetune 函数在 source_file_path 中定义，并且能够将图保存到当前工作目录 (current_analysis_output_dir)
  if (exists("finetune", mode="function")){
    finetune(df=count_bin2, width=window_width, potential_peaks=potential_peaks_indices,
             my_quantile=c(quantile_Deql_min, quantile_Dnun_min, quantile_TSkn_min),
             my_band=band_TSkn_min, my_var.est = chiptest_var_est, my_var.thred = chiptest_var_thred)
  } else {
    warning("Function 'finetune' for plotting not found in source file.")
  }
  
  
  # --- 5. ChIPtest - Hypothesis Testing & Performance ---
  message("Preparing data for final ChIPtest hypothesis testing...")
  result_temp <- as.data.frame(data_conA[,c(1,2)]) # data_conA 行数应与 potential_peaks_indices 长度一致
  colnames(result_temp) <- c("start", "end")
  result_temp$index <- 1:nrow(result_temp)
  
  sites_overlap <- find_overlap(is_tf, pos.1, if(is_tf) pos.1 else pos.3, peak_width, lfile, result_temp)
  sites_overlap[is.na(overlap_true), overlap_true := 0]
  data.table::setnames(sites_overlap, old = "overlap_true", new = "overlap")
  
  
  ChIPtest_result <- sites_overlap[,c("index", "start", "end", "sim_start", "sim_end", "overlap")]
  
  message("Running ChIPtest hypothesis testing...")
  # 确保 NormTransformation, my_est.c, my_TS_twosample 在 source_file_path 中定义
  # ChIPtest_result$index 应该对应于 data_conA/B 的行
  valid_indices <- ChIPtest_result$index[ChIPtest_result$index <= nrow(data_conA) & ChIPtest_result$index <= nrow(data_conB)]
  if(length(valid_indices) != nrow(ChIPtest_result)){
    warning("Some indices from find_overlap results are out of bounds for data_conA/B. Filtering them.")
    ChIPtest_result <- ChIPtest_result[ChIPtest_result$index %in% valid_indices, ]
    if(nrow(ChIPtest_result) == 0) stop("No valid sites remaining after index validation for ChIPtest.")
  }
  
  Data_conA_norm <- NormTransformation(data_conA[, -c(1,2), drop=FALSE])
  Data_conB_norm <- NormTransformation(data_conB[, -c(1,2), drop=FALSE])
  
  tao <- my_est.c(Data_conA_norm, Data_conB_norm, max1=4, max4=4) # max1, max4 可作为参数
  TS <- my_TS_twosample(Data_conA_norm, Data_conB_norm, tao, band_TSkn_min,
                        quant=c(quantile_Deql_min, quantile_Dnun_min, quantile_TSkn_min),
                        var.est = chiptest_var_est, var.thred = chiptest_var_thred)
  
  # ChIPtest_result$pvalue_TS_high <- pnorm(TS[["TS_kn"]], mean = 0, sd = 1, lower.tail = FALSE)
  # ChIPtest_result$pvalue_Deql <- pnorm(TS[["Deql"]], mean = 0, sd = 1, lower.tail = FALSE)
  # ChIPtest_result$pvalue_Dnun <- pnorm(TS[["Dnun"]], mean = 0, sd = 1, lower.tail = FALSE)
  # ChIPtest_result$chromosome <- 'chrA' # 假设模拟数据都在 chrA
  # 
  # ChIPtest_result$TSkn_BH <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "BH")
  # ChIPtest_result$Dnun_BH <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "BH")
  # ChIPtest_result$Deql_BH <- p.adjust(ChIPtest_result$pvalue_Deql, method = "BH")
  # 
  # ChIPtest_result$TSkn_bon <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "bonferroni")
  # ChIPtest_result$Dnun_bon <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "bonferroni")
  # ChIPtest_result$Deql_bon <- p.adjust(ChIPtest_result$pvalue_Deql, method = "bonferroni")
  sites_overlap[is.na(overlap), overlap := 0]
  # ChIPtest_result <- copy(sites_overlap)
  ChIPtest_result_non_na_idx <- ChIPtest_result[!is.na(index)]
  ChIPtest_result_na_idx <- ChIPtest_result[is.na(index)]
  
  
  # 处理 index 非 NA 的部分
  if (nrow(ChIPtest_result_non_na_idx) > 0) {
    # 按照 index 从小到大排序 (确保正确索引 TS 列表)
    setorderv(ChIPtest_result_non_na_idx, "index")
    
    # 获取 TS 中需要引用的索引值
    ts_indices <- ChIPtest_result_non_na_idx$index
    
    # 安全性检查：确保 ts_indices 不会超出 TS 向量的范围
    max_ts_len <- length(TS[["TS_kn"]])
    # 任何超出范围的索引都会导致取值变成 NA
    ts_indices_safe <- ts_indices
    ts_indices_safe[ts_indices_safe < 1 | ts_indices_safe > max_ts_len] <- NA_integer_
    
    # 计算 P 值
    ChIPtest_result_non_na_idx[, pvalue_TS_high := pnorm(TS[["TS_kn"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
    ChIPtest_result_non_na_idx[, pvalue_Deql := pnorm(TS[["Deql"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
    ChIPtest_result_non_na_idx[, pvalue_Dnun := pnorm(TS[["Dnun"]][ts_indices_safe], mean = 0, sd = 1, lower.tail = FALSE)]
  } else {
    # 如果没有非 NA 的 index 行，则创建空的 data.table，并包含相应的列
    ChIPtest_result_non_na_idx <- ChIPtest_result_non_na_idx[, `:=`(pvalue_TS_high = numeric(), pvalue_Deql = numeric(), pvalue_Dnun = numeric())]
  }
  
  # 3. 处理 index 为 NA 的部分
  if (nrow(ChIPtest_result_na_idx) > 0) {
    # 将所有 P 值赋为 1
    ChIPtest_result_na_idx[, pvalue_TS_high := 1]
    ChIPtest_result_na_idx[, pvalue_Deql := 1]
    ChIPtest_result_na_idx[, pvalue_Dnun := 1]
  } else {
    # 如果没有 NA 的 index 行，则创建空的 data.table
    ChIPtest_result_na_idx <- ChIPtest_result_na_idx[, `:=`(pvalue_TS_high = numeric(), pvalue_Deql = numeric(), pvalue_Dnun = numeric())]
  }
  # 合并两部分结果
  ChIPtest_result <- rbindlist(list(ChIPtest_result_non_na_idx, ChIPtest_result_na_idx), fill = TRUE)
  ChIPtest_result$chromosome <- 'chr1'
  ChIPtest_result$TSkn_BH <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "BH")
  ChIPtest_result$Dnun_BH <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "BH")
  ChIPtest_result$Deql_BH <- p.adjust(ChIPtest_result$pvalue_Deql, method = "BH")
  ChIPtest_result$TSkn_bon <- p.adjust(ChIPtest_result$pvalue_TS_high, method = "bonferroni")
  ChIPtest_result$Dnun_bon <- p.adjust(ChIPtest_result$pvalue_Dnun, method = "bonferroni")
  ChIPtest_result$Deql_bon <- p.adjust(ChIPtest_result$pvalue_Deql, method = "bonferroni")
  
  # --- Select most significant p-value for overlapping true sites ---
  dt_ChIPtest_result <- data.table::setDT(as.data.frame(ChIPtest_result))
  
  ChIPtest_result_0 <- dt_ChIPtest_result[is.na(sim_start)]
  ChIPtest_result_1 <- dt_ChIPtest_result[!is.na(sim_start), .(
    index = index[1], start = start[1], end = end[1], sim_end = sim_end[1],
    overlap = overlap[1], chromosome = chromosome[1],
    pvalue_TS_high = min(pvalue_TS_high, na.rm = TRUE),
    pvalue_Deql = min(pvalue_Deql, na.rm = TRUE),
    pvalue_Dnun = min(pvalue_Dnun, na.rm = TRUE),
    TSkn_BH = min(TSkn_BH, na.rm = TRUE), Dnun_BH = min(Dnun_BH, na.rm = TRUE),
    Deql_BH = min(Deql_BH, na.rm = TRUE), TSkn_bon = min(TSkn_bon, na.rm = TRUE),
    Dnun_bon = min(Dnun_bon, na.rm = TRUE), Deql_bon = min(Deql_bon, na.rm = TRUE)
  ), by = sim_start]
  
  if(nrow(ChIPtest_result_0) > 0 && nrow(ChIPtest_result_1) > 0) {
    data.table::setcolorder(ChIPtest_result_1, names(ChIPtest_result_0))
    ChIPtest_selected <- data.table::rbindlist(list(ChIPtest_result_1, ChIPtest_result_0), use.names = TRUE, fill = TRUE)
  } else if (nrow(ChIPtest_result_1) > 0) {
    ChIPtest_selected <- ChIPtest_result_1
  } else {
    ChIPtest_selected <- ChIPtest_result_0
  }
  ChIPtest_selected <- as.data.frame(ChIPtest_selected)
  
  
  # --- Calculate Performance Metrics ---
  cat("Calculating final performance metrics...\n")
  if (!exists("calculate_fdr_recall", mode = "function")) {
    stop("Function 'calculate_fdr_recall' not found.")
  }
  
  metrics_list <- list()
  models <- c("TSkn", "Deql", "Dnun")
  controls <- c("BH", "bon")
  
  if(nrow(ChIPtest_selected) == 0 || !"overlap" %in% names(ChIPtest_selected)) {
    warning("ChIPtest_selected is empty or missing 'overlap' column. Cannot calculate performance metrics.")
    fdr_recall_all_df <- data.frame(Threshold=NA, FDR=NA, Recall=NA, model=NA, control=NA, auc_value=NA)[-1,]
  } else {
    for (model in models) {
      for (control in controls) {
        p_col_name <- paste0(model, "_", control)
        if (!p_col_name %in% names(ChIPtest_selected)) {
          warning(paste("P-value column", p_col_name, "not found in results. Skipping."))
          next
        }
        
        valid_for_roc <- !is.na(ChIPtest_selected[[p_col_name]]) & !is.na(ChIPtest_selected$overlap)
        if(sum(valid_for_roc) < 2 || length(unique(ChIPtest_selected$overlap[valid_for_roc])) < 2) {
          auc_val <- NA
          warning(paste("Not enough distinct data points for ROC calculation for", p_col_name, ". AUC will be NA."))
        } else {
          pred_obj <- ROCR::prediction(1 - ChIPtest_selected[[p_col_name]][valid_for_roc], ChIPtest_selected$overlap[valid_for_roc])
          auc_val <- ROCR::performance(pred_obj, "auc")@y.values[[1]]
        }
        
        temp_metrics <- calculate_fdr_recall(ChIPtest_selected[[p_col_name]], ChIPtest_selected$overlap, thresholds = c(0.05))
        temp_metrics_df <- as.data.frame(temp_metrics)
        temp_metrics_df$model <- model
        temp_metrics_df$control <- control
        temp_metrics_df$auc_value <- auc_val
        metrics_list[[paste0(model, "_", control)]] <- temp_metrics_df
      }
    }
    fdr_recall_all_df <- do.call(rbind, metrics_list)
  }
  
  # --- 6. 保存结果 ---
  message("Saving ChIPtest results...")
  # 文件将保存到 current_analysis_output_dir (当前工作目录)
  write.csv(ChIPtest_selected, file = "results_macs3_ChIPtest.csv", row.names = FALSE)
  write.csv(fdr_recall_all_df, file = "matrics_macs3_ChIPtest.csv", row.names = FALSE)
  message("ChIPtest results saved.")
  
  # --- 7. 返回 ---
  message("MACS3 + ChIPtest pipeline run complete for run_tag: ", run_tag)
  return(
    list(
      results_ChIPtest_df = ChIPtest_selected,
      metrics_ChIPtest_df = fdr_recall_all_df,
      output_path = current_analysis_output_dir,
      # 可以添加其他需要返回的对象，例如 finetune_results
      finetune_params = list(quantile_TSkn_min=quantile_TSkn_min, band_TSkn_min=band_TSkn_min,
                             quantile_Deql_min=quantile_Deql_min, quantile_Dnun_min=quantile_Dnun_min)
    )
  )
}





























































