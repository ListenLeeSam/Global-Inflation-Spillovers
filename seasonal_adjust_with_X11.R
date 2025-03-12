library(tseries)  
library(moments)  
library(zoo)      
library(seasonal) 

# 读取数据
file_path <- "C:/Users/rog/Desktop/Study/统计学习/CPI_en.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)

# 处理成年月格式
df$Date <- as.yearmon(df$Date, format = "%m/%Y")  
df <- df[order(df$Date), ]  # 按日期排序

# 去除NA值
df_clean <- na.omit(df)

# 进行季节性调整
df_adjusted <- df_clean
for (i in 2:ncol(df_clean)) {  # 遍历所有数值列
  ts_data <- ts(df_clean[, i], start = c(as.numeric(format(df_clean$Date[1], "%Y")), 
                                         as.numeric(format(df_clean$Date[1], "%m"))), 
                frequency = 12)  # 设定时间序列的频率为12（按月数据）
  
  seas_adj <- seas(ts_data, x11 = "")  # 使用 X-11ARIMA-SEATS 方法进行季节性调整
  df_adjusted[, i] <- final(seas_adj)  # 提取季节调整后的数据
}

# 计算环比增长率（月度变化率）
df_mom <- df_adjusted
for (i in 2:ncol(df_adjusted)) {
  df_mom[, i] <- c(NA, diff(df_adjusted[, i]) / df_adjusted[-nrow(df_adjusted), i] * 100)
}

# 进行ADF单位根检验
non_stationary_countries <- c()  # 用于存储非平稳国家
for (col in colnames(df_mom)[-1]) {  # 跳过 Date 列
  if (sum(!is.na(df_mom[[col]])) > 0) {  # 确保有数据
    adf_test <- adf.test(na.omit(df_mom[[col]]))
    if (adf_test$p.value > 0.05) {  # p值 > 0.05 说明数据非平稳
      non_stationary_countries <- c(non_stationary_countries, col)
    }
  }
}

# 打印非平稳国家列表
print("非平稳的国家:")
print(non_stationary_countries)