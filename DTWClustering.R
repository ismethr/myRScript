# 安装必要的包 (如果还没安装)
# install.packages(c("dtwclust", "tidyverse", "imputeTS", "doParallel"))

library(tidyverse)
library(dtwclust)   # DTW 聚类核心包
library(imputeTS)   # 用于填补时间序列缺失值
library(doParallel) # 并行计算


# 配置路径
INPUT_FILE <- "/Users/wendao/Library/CloudStorage/SynologyDrive-Mac/论文/02.CFI/格网尺度分析/H3_L7_Fragmentation_Trend_Analysis_Mac.csv"

# 读取数据
raw_df <- read_csv(INPUT_FILE, show_col_types = FALSE)

# 1. 提取数值矩阵 (只取年份列)
# 假设年份列是第 2 列到最后一列 (根据您的文件结构调整)
# 或者使用名字匹配
year_cols <- grep("^19|^20", names(raw_df), value = TRUE)
data_matrix <- as.matrix(raw_df[, year_cols])

# 2. 缺失值处理 (DTW 极其讨厌 NA，必须处理)
# 使用线性插值填补中间的 NA，两头缺失用最近值填充
#这一步可能需要几十秒
cat("正在处理缺失值...\n")
data_matrix <- na_interpolation(data_matrix, option = "linear")

# 3. 数据归一化 (Z-Score) —— 【重要决策】
# 选项 A: 使用原始值 (Raw)。可以区分 "一直很高" 和 "一直很低"。
# 选项 B: 使用 Z-Score。只看 "形状" (上升/下降)，忽略绝对高低。
# ---
# 建议：由于您想识别 "Class 3: 稳定低值型"，绝对值很重要，建议【不】进行 Z-Score 标准化，或者使用 Min-Max。
# 但为了计算速度和收敛，通常建议归一化。
# 这里我们采用 【Z-Score 标准化】，然后在最后画图时画回原始值，或者在解释时结合原始均值。
# 如果您非常看重绝对值 (0.8 vs 0.1)，请注释掉下面这行：
data_matrix_norm <- zscore(data_matrix) 

# 为了测试代码，建议先用 5000 个数据跑通，再跑全量
# sample_idx <- sample(nrow(data_matrix_norm), 5000)
# test_data <- data_matrix_norm[sample_idx, ]


# 设置聚类簇数 (您预想是 4-6 类，我们设为 4)
k_clusters <- 4

cat("开始 DTW 聚类 (这可能需要几分钟到几十分钟)...\n")

# 注册并行核心 (Mac M4芯片可以用多核)
registerDoParallel(cores = 6)

# 执行聚类
# type = "partitional": 类似 K-Means 的划分聚类 (适合大数据量)
# dist = "dtw_basic": 基础 DTW 距离 (比完整 DTW 快很多)
# centroid = "dba": DTW Barycenter Averaging (DTW 专用质心算法，能生成平滑曲线)
hc_dtw <- tsclust(data_matrix_norm, 
                  type = "partitional", 
                  k = k_clusters, 
                  distance = "dtw_basic", 
                  centroid = "dba", 
                  seed = 123,           # 固定随机种子，保证结果可复现
                  trace = TRUE,         # 显示进度
                  args = tsclust_args(dist = list(window.size = 5))) # 窗口限制，加速计算

cat("聚类完成！\n")


# 绘制 4 类质心曲线 (Centroids)
# 这张图就是您要放在论文里的 "Figure: Temporal Patterns"
plot(hc_dtw, type = "centroids") + 
  labs(title = "4 Typical Evolutionary Patterns of Fragmentation",
       x = "Time (1990-2025)", 
       y = "Normalized CFI Value") +
  theme_minimal()

# 如果想看每一类里包含多少个格网
print(table(hc_dtw@cluster))


# 1. 提取聚类标签 (1, 2, 3, 4)
cluster_labels <- hc_dtw@cluster

# 2. 合并结果
result_df <- raw_df %>%
  select(GRID_ID) %>%
  mutate(Cluster_Class = cluster_labels)

# 3. (可选) 计算每一类的一些统计特征，方便你在 ArcGIS 里以此命名
# 比如计算每类的 1990年均值和 2025年均值
summary_stats <- bind_cols(raw_df, Cluster_Class = cluster_labels) %>%
  group_by(Cluster_Class) %>%
  summarise(
    Mean_1990 = mean(`1990`, na.rm = TRUE),
    Mean_2025 = mean(`2025`, na.rm = TRUE),
    Count = n()
  )
print(summary_stats)

# 4. 保存 CSV
output_file <- "H3_L7_DTW_Clustering_Result.csv"
write_csv(result_df, output_file)

cat(sprintf("结果已保存至 %s，请在 ArcGIS Pro 中进行制图。\n", output_file))
