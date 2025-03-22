# 加载必要的包
library(dplyr)

# 读取两个tsv文件
file1 <- read.delim("GCST90016564_buildGRCh37.tsv", header = TRUE, sep = "\t")
file2 <- read.delim("GCST90029028_buildGRCh37.tsv", header = TRUE, sep = "\t")

# 合并两个数据框
combined_data <- bind_rows(file1, file2)

# 检查是否有重复的行
duplicated_rows <- duplicated(combined_data)

# 输出重复数据的情况
if (any(duplicated_rows)) {
  cat("存在重复数据。\n")
  cat("重复数据的行数:", sum(duplicated_rows), "\n")
  cat("重复数据的详细情况:\n")
  print(combined_data[duplicated_rows, ])
} else {
  cat("没有发现重复数据。\n")
}

# 删除重复数据
unique_data <- distinct(combined_data)

# 保存去重后的数据到新的tsv文件
write.table(unique_data, "path/to/output_unique.tsv", sep = "\t", row.names = FALSE)