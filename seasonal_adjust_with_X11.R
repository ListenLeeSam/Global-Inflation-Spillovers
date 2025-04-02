# 加载必要的库
library(tseries)   # ADF检验
library(moments)   # 计算偏度、峰度
library(zoo)       # 处理时间序列
library(seasonal)  # 时间序列季节性调整
library(knitr)     # 表格展示
library(xtable)    # LaTeX 输出
library(urca)      # 协整检验
library(ggplot2)   # 可视化
library(gridExtra) # 处理多个图
library(reshape2)  # 转换数据格式
library(changepoint)  # 变点检测
library(classInt)  # 用于 Jenks Natural Breaks
library(dplyr)     # 数据处理

# 读取数据
file_path <- "new_data\\transposed_cpi_data_filtered.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)

# 处理成年月格式
df$Date <- as.yearmon(df$Date, format = "%Y/%m/%d")  
df <- df[order(df$Date), ]  

# 去除NA值
df_clean <- na.omit(df)

# 确保df_clean至少有一行数据
if (nrow(df_clean) == 0) {
  stop("Error: 清理后的数据为空，请检查数据预处理步骤！")
}

# 进行季节性调整
df_adjusted <- df_clean
for (i in 2:ncol(df_clean)) {  # 遍历所有数值列
  if (sum(!is.na(df_clean[, i])) > 0) {  # 确保列中至少有一个非NA值
    ts_data <- ts(df_clean[, i], start = c(as.numeric(format(df_clean$Date[1], "%Y")), 
                                           as.numeric(format(df_clean$Date[1], "%m"))), 
                  frequency = 12)  # 设定时间序列的频率为12（按月数据）
    
    seas_adj <- seas(ts_data, x11 = "")  # 使用 X-11ARIMA-SEATS 方法进行季节性调整
    df_adjusted[, i] <- final(seas_adj)  # 提取季节调整后的数据
  }
}

# 计算环比增长率（月度变化率）
df_mom <- df_adjusted
for (i in 2:ncol(df_adjusted)) {
  if (sum(!is.na(df_adjusted[, i])) > 1) {  # 确保至少有2个非NA值，才能计算diff
    df_mom[, i] <- c(NA, diff(df_adjusted[, i]) / df_adjusted[-nrow(df_adjusted), i] * 100)
  }
}

save_path <- "new_data\\final_cpi_data_filtered.csv"
write.csv(df_mom, save_path, row.names = FALSE)

# =========================================================== 描述性统计 ===========================================
# 进行ADF单位根检验
non_stationary_results <- data.frame(Country = character(), P_Value = numeric(), stringsAsFactors = FALSE)

for (col in colnames(df_mom)[-1]) {  
  if (sum(!is.na(df_mom[[col]])) > 0) {  
    adf_test <- adf.test(na.omit(df_mom[[col]])) 
    p_value <- adf_test$p.value 
    
    if (p_value > 0.05) { 
      non_stationary_results <- rbind(non_stationary_results, data.frame(Country = col, P_Value = p_value))
    }
  }
}

# 打印非平稳国家及其ADF检验p值
print("非平稳国家的 ADF 检验 p 值：")
print(non_stationary_results)

# 以表格形式展示
kable(non_stationary_results)

desc_stats <- data.frame(Country = character(),
                         Mean = numeric(),
                         Median = numeric(),
                         Std_Dev = numeric(),
                         Skewness = numeric(),
                         Kurtosis = numeric(),
                         Max = numeric(),
                         Min = numeric(),
                         JB_PValue = numeric(),
                         stringsAsFactors = FALSE)

for (col in colnames(df_mom)[-1]) {  
  data <- as.vector(na.omit(df_mom[[col]])) 
  if (length(data) > 2) { 
    jb_test <- jarque.test(data)  
    desc_stats <- rbind(desc_stats, data.frame(
      Country = col,
      Mean = round(mean(data), 2),
      Median = round(median(data), 2),
      Std_Dev = round(sd(data), 2),
      Skewness = round(skewness(data), 2),
      Kurtosis = round(kurtosis(data), 2),
      Max = round(max(data), 2),
      Min = round(min(data), 2),
      JB_PValue = round(jb_test$p.value, 2)  
    ))
  }
}

# 打印描述性统计表格
print("描述性统计：")
print(desc_stats)
kable(desc_stats)

# ==================================================================== 直方图 ========================================================
tex_output <- xtable(desc_stats)
print(tex_output, file = "descriptive_statistics.tex")

# 生成所有国家的直方图
df_long <- melt(df_mom, id.vars = "Date", variable.name = "Country", value.name = "Inflation")
df_long <- na.omit(df_long)  # 删除缺失值

plot_list <- list()
unique_countries <- unique(df_long$Country)  # 获取所有国家名称

for (i in seq_along(unique_countries)) {
  country_name <- unique_countries[i]  # 当前国家名称
  country_data <- subset(df_long, Country == country_name)  # 选择该国家的数据
  
  p <- ggplot(country_data, aes(x = Inflation)) +
    geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(paste0(i, ". ", country_name)) +  # 添加标题（序号+国家名）
    theme_minimal() +
    xlab("Inflation Rate") + ylab("Frequency") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))  # 标题居中
  
  plot_list[[i]] <- p  # 存入列表
}

final_plot <- do.call(grid.arrange, c(plot_list, ncol = 5))

ggsave("new_data\\Inflation_Histograms.png", plot = final_plot, width = 15, height = 10, dpi = 300)

# =================================================================== 协整检验 ============================================
# 选择需要进行协整检验的列
selected_columns <- df_mom[, -1]  # 删除日期列，选择其他国家的通胀数据
selected_columns_clean <- na.omit(selected_columns)  # 删除NA值

# 转换为矩阵格式
data_matrix <- as.matrix(selected_columns_clean)

# 进行Johansen协整检验
johansen_test <- ca.jo(data_matrix, type = "trace", K = 2, spec = "longrun")

# 输出协整检验结果
summary(johansen_test)

# ==================================================================== 变点1 ========================================================
df_summary <- aggregate(Inflation ~ Date, data = df_long, 
                        FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))

df_summary <- do.call(data.frame, df_summary)  # 拆分嵌套列
colnames(df_summary) <- c("Date", "Mean_Inflation", "Min_Inflation", "Max_Inflation")

cpt <- cpt.meanvar(df_summary$Mean_Inflation, method = "PELT")  # PELT 方法检测均值变点
change_points <- cpts(cpt)  # 获取变点索引
change_dates <- df_summary$Date[change_points]  # 获取变点对应的时间点

# 打印变点的时间点
print("均值通胀率的变点时间：")
print(change_dates)

# 绘图
p <- ggplot(df_summary, aes(x = Date)) +
  geom_ribbon(aes(ymin = Min_Inflation, ymax = Max_Inflation), fill = "gray80", alpha = 0.5) +  # 灰色区间
  geom_line(aes(y = Mean_Inflation), color = "black", size = 1) +  
  geom_vline(xintercept = as.numeric(change_dates), linetype = "dashed", color = "red", size = 1) +  # 变点
  labs(y = "通胀率 (%)", x = "时间") +
  theme_classic() +  
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# 保存图片
ggsave("new_data\\Inflation_TimeSeries.png", 
       plot = p, width = 12, height = 6, dpi = 300)




