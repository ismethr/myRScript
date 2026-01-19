# ==============================================================================
# 0. 库加载与环境配置
# ==============================================================================
# 如果没有安装，请取消下面这行的注释并运行一次
# install.packages(c("tidyverse", "trend", "future.apply", "matrixStats"))

library(tidyverse)    # 数据清洗
library(trend)        # MK 检验
library(future.apply) # 并行计算
library(matrixStats)  # 矩阵统计

# --- 配置区域 (macOS 特有修改) ---

# 【关键修改 1】路径格式
# 请在 Finder 中找到文件夹，右键点击 -> 按住 Option 键 -> 选择“拷贝为路径名称”
# 然后粘贴到下面。注意 macOS 路径以 /Users/ 开头
DATA_FOLDER <- "/Users/wendao/Library/CloudStorage/OneDrive-个人/桌面/论文/02.CFI/原始计算结果/格网尺度"

# 文件名模板 (保持不变)
FILENAME_TEMPLATE <- "grid7_30_landscape_metrics_summary_parallel_%s_robustnorm.csv"

# 年份范围
START_YEAR <- 1990
END_YEAR <- 2025
YEARS <- as.character(START_YEAR:END_YEAR)

# 并行核心数
# M4 芯片通常有 8-10 核以上，建议留 2 个核心给系统
# 您可以在终端输入 sysctl -n hw.ncpu 查看总核数
NB_WORKERS <- 6 

# ==============================================================================
# 第一步：数据清洗与整合 (Data Loading & Reshaping)
# ==============================================================================
cat("\n[Step 1] 开始读取并整合每年的 CSV 文件...\n")

read_year_data <- function(year) {
  file_path <- file.path(DATA_FOLDER, sprintf(FILENAME_TEMPLATE, year))
  
  if (file.exists(file_path)) {
    # 只读取 GRID_ID 和 CFI
    # show_col_types = FALSE 用于静默读取，不刷屏
    df <- read_csv(file_path, col_select = c("GRID_ID", "CFI"), show_col_types = FALSE) %>%
      rename(!!year := CFI) # 动态重命名列
    return(df)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)
  }
}

# 批量读取
list_of_dfs <- lapply(YEARS, read_year_data)
list_of_dfs <- list_of_dfs[!sapply(list_of_dfs, is.null)]

# 合并所有年份
combined_df <- list_of_dfs %>%
  reduce(full_join, by = "GRID_ID")

# ------------------------------------------------------------------
# 【新增功能】保存整合后的中间数据
# ------------------------------------------------------------------
raw_output_filename <- "H3_L7_Integrated_Raw_1990_2025.csv"

cat(sprintf("\n[Step 1.5] 正在导出整合后的原始宽表数据至: %s ...\n", raw_output_filename))
cat("这可能需要一点时间，请稍候...\n")

write_csv(combined_df, raw_output_filename)

cat(">>> 导出成功！您现在可以打开该文件检查：\n")
cat("1. 是否包含从 1990 到 2025 的所有年份列\n")
cat("2. GRID_ID 是否正确\n")
cat("3. 数据是否按年份对齐 (没有发生错位)\n")
# ------------------------------------------------------------------

# --- CHECKPOINT 1 ---
cat("\n>>> CHECKPOINT 1: 数据整合完成 <<<\n")
cat("总行数 (格网数):", nrow(combined_df), "\n")
print(head(combined_df))

# ==============================================================================
# 第二步：定义分析函数 (Sen + MK)
# ==============================================================================
calc_trend_stats <- function(x) {
  # 移除 NA
  valid_data <- x[!is.na(x)]
  
  if (length(valid_data) < 10) {
    return(c(Sen_Slope = NA, P_Value = NA, Trend_Type = "Insufficient Data"))
  }
  
  tryCatch({
    # MK 检验
    mk_res <- trend::mk.test(valid_data)
    p_value <- mk_res$p.value
    
    # Sen's Slope
    sen_res <- trend::sens.slope(valid_data)
    slope <- sen_res$estimates
    
    # 趋势分类
    trend_type <- "Stable"
    if (p_value < 0.05) { # 显著
      if (slope > 0) {
        trend_type <- "Significant Deterioration" 
      } else if (slope < 0) {
        trend_type <- "Significant Improvement"   
      }
    } else { # 不显著
      if (slope > 0) {
        trend_type <- "Slight Deterioration"
      } else if (slope < 0) {
        trend_type <- "Slight Improvement"
      }
    }
    
    return(c(Sen_Slope = as.numeric(slope), P_Value = as.numeric(p_value), Trend_Type = trend_type))
    
  }, error = function(e) {
    return(c(Sen_Slope = NA, P_Value = NA, Trend_Type = "Error"))
  })
}

# ==============================================================================
# 第三步：并行计算趋势 (macOS 优化)
# ==============================================================================
cat("\n[Step 2] 正在进行并行趋势分析...\n")

# 准备矩阵 (显式按年份排序)
data_matrix <- as.matrix(combined_df[, YEARS])

# --- CHECKPOINT 2 ---
cat("\n>>> CHECKPOINT 2: 矩阵准备就绪 <<<\n")
cat("验证年份顺序 (前5列): ", paste(colnames(data_matrix)[1:5], collapse=", "), "\n")

# 【关键修改 2】设置并行计划
# 选项 A: 如果您使用 RStudio，请保持 'multisession' (最稳妥)
# 选项 B: 如果您使用 VS Code 或 Terminal 运行 R，可以改为 'multicore' (速度更快，内存更省)
plan(multisession, workers = NB_WORKERS) 

# 开始计算
trend_results_matrix <- future_apply(data_matrix, 1, calc_trend_stats)
trend_results_df <- as.data.frame(t(trend_results_matrix))

# 类型转换
trend_results_df$Sen_Slope <- as.numeric(trend_results_df$Sen_Slope)
trend_results_df$P_Value <- as.numeric(trend_results_df$P_Value)

# --- CHECKPOINT 3 ---
cat("\n>>> CHECKPOINT 3: 趋势分析完成 <<<\n")
print(head(trend_results_df))

# ==============================================================================
# 第四步：计算稳定性 (CV) & 合并
# ==============================================================================
cat("\n[Step 3] 计算变异系数 (CV)...\n")

row_means <- rowMeans(data_matrix, na.rm = TRUE)
row_sds <- matrixStats::rowSds(data_matrix, na.rm = TRUE)
cv_values <- row_sds / row_means
cv_values[row_means == 0] <- NA

final_df <- bind_cols(
  GRID_ID = combined_df$GRID_ID,
  trend_results_df,
  CV = cv_values,
  combined_df[, YEARS]
)

# ==============================================================================
# 第五步：保存结果
# ==============================================================================
output_filename <- "H3_L7_Fragmentation_Trend_Analysis_Mac.csv"
cat(sprintf("\n[Step 4] 保存结果至: %s/%s \n", getwd(), output_filename))

write_csv(final_df, output_filename)

cat("全部完成！请检查输出文件。\n")