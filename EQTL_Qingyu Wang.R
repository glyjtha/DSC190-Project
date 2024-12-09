
install.packages("data.table")
#install.packages("R.utils")
library(data.table)
#library(R.utils)
file_path <- "/Users/wqy/Downloads/GCST90454344.tsv"
gwas_data <- fread(file_path, sep = "\t", header = TRUE)
# 查看数据
str(gwas_data)
head(gwas_data)
colnames(gwas_data)
#尝试1kg
library(dplyr) 
install.packages("readr")
library(readr)   
bim_dir <- "/Users/wqy/Downloads/1KG"
colnames(gwas_data)[colnames(gwas_data) == "chromosome"] <- "CHROM"
colnames(gwas_data)[colnames(gwas_data) == "base_pair_location"] <- "POS"
#尝试1kg
bim_data <- data.frame()
for (chr in 1:22) {
  bim_file <- file.path(bim_dir, paste0("1000G_eur_chr", chr, ".bim"))
  if (file.exists(bim_file)) {
    # 读取 .bim 文件
    bim_chr <- read_table2(bim_file, col_names = c("CHROM", "rsID", "cm", "POS", "REF", "ALT"))
    bim_data <- bind_rows(bim_data, bim_chr)
  }
}
head(bim_data)
# match
annotated_gwas <- gwas_data %>%
  inner_join(bim_data, by = c("CHROM", "POS"))
head(annotated_gwas)
write_tsv(annotated_gwas, "/Users/wqy/Downloads/annotated_gwas_data.tsv")
#尝试1kg

significance_threshold <- 7.3  # 对应 p-value < 5e-8
significant_snps <- annotated_gwas %>%
  filter(neg_log_10_p_value >= significance_threshold)
print(dim(significant_snps))  # 显示显著 SNP 的数量
print(head(significant_snps)) 
# rsID
significant_snp_ids <- significant_snps$SNP
write.table(significant_snp_ids, "significant_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#图
library(ggplot2)
# 显著 SNP 分布图
significant_snps <- annotated_gwas[annotated_gwas$`neg_log_10_p_value` > 7.301, ]
ggplot(significant_snps, aes(x = CHROM, y = -log10(10^(-neg_log_10_p_value)))) +
  geom_point(alpha = 0.6) +
  labs(title = "Significant SNP Distribution",
       x = "Chromosome",
       y = "-log10(P-value)") +
  theme_minimal()
#图

#rna-seq
# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)

# Step 1: Load gene-cell matrix
gene_cell_matrix_path <- "/Users/wqy/Downloads/GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv"
gene_cell_matrix <- read.csv(gene_cell_matrix_path, row.names = 1)

# Step 2: Load metadata
metadata_path <- "/Users/wqy/Downloads/GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv"
metadata <- read.csv(metadata_path)

# Preview the data
cat("Gene-cell matrix dimensions:", dim(gene_cell_matrix), "\n")
head(gene_cell_matrix)

cat("Metadata dimensions:", dim(metadata), "\n")
head(metadata)

# Check for missing values
cat("Missing values in gene-cell matrix:", sum(is.na(gene_cell_matrix)), "\n")
cat("Missing values in metadata:", sum(is.na(metadata)), "\n")
#rna-seq

#eQTL
library(data.table)# 用于快速读取大文件
untar("/Users/wqy/Downloads/GTEx_Analysis_v10_eQTL.tar", exdir = "/Users/wqy/Downloads/GTEx_Analysis_v10_eQTL")
library
library(arrow)
# 加载 eQTL 数据
eqtl_data <- read_parquet("/Users/wqy/Downloads/GTEx_Analysis_v10_eQTL_updated/Breast_Mammary_Tissue.v10.eQTLs.signif_pairs.parquet")
# 查看数据结构
head(eqtl_data)
# 筛选目标基因
read_expression_data <- function(eqtl_data, gene_id = "ENSG00000139618") {
  # 筛选目标基因表达信息
  gene_eqtl <- eqtl_data[grepl(gene_id, eqtl_data$gene_id), ]
  
  # 检查结果是否为空
  if (nrow(gene_eqtl) == 0) {
    stop("Target gene not found in the dataset!")
  }
  
  return(gene_eqtl)
}
# 调用函数
gene_eqtl <- read_expression_data(eqtl_data)
# Check the structure and summary of the target gene eQTL
str(gene_eqtl)
summary(gene_eqtl)

# Check the number of significant SNPs
pval_threshold <- 1e-4
significant_eqtl <- gene_eqtl[gene_eqtl$pval_nominal < pval_threshold, ]
print(paste("Number of significant eQTLs:", nrow(significant_eqtl)))
# Visualize p-values
hist(gene_eqtl$pval_nominal, breaks = 50, main = "Distribution of P-values",
     xlab = "P-value", col = "blue")
# 可视化显著 eQTL 的 P 值分布
library(ggplot2)
ggplot(gene_eqtl, aes(x = tss_distance, y = -log10(pval_nominal))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(1e-4), linetype = "dashed", color = "red") +
  labs(title = "eQTL Analysis for Target Gene",
       x = "Distance to TSS (bp)",
       y = "-log10(P-value)") +
  theme_minimal()
# 构建线性模型：SNP 对基因表达的调控
# 检查 af 和 slope 的分布
summary(significant_eqtl$af)
summary(significant_eqtl$slope)
# 移除包含 NA 的行（如果有）
significant_eqtl <- significant_eqtl[!is.na(significant_eqtl$af) & !is.na(significant_eqtl$slope), ]
# 构建线性模型
lm_model <- lm(slope ~ af, data = significant_eqtl)
# 输出模型摘要
summary(lm_model)
library(ggplot2)
# 绘制线性模型拟合图
ggplot(significant_eqtl, aes(x = af, y = slope)) +
  geom_point(color = "blue", alpha = 0.6) +  # 数据点
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # 拟合线
  labs(title = "Linear Model: SNP Regulation on BRCA2 Expression",
       x = "Allele Frequency (af)",
       y = "eQTL Effect Size (slope)") +
  theme_minimal()