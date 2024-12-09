
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

install.packages("glmnet")  # LASSO/Elastic Net
install.packages("data.table")  # 数据操作
install.packages("dplyr")  # 数据处理
library(data.table)
library(glmnet)
library(dplyr)
gwas_data <- gwas_data %>% filter(!is.na(beta))  # 移除缺失值
# 载入基因表达数据 (假设是GTEx数据)
expression_data <- fread("gene_expression_matrix.csv")  # 行是样本，列是基因
brca2_expr <- expression_data[, "ENSG00000139618", with = FALSE]
# 提取 eQTL 信息
eqtl_snps <- target_eqtl$variant_id
eqtl_beta <- target_eqtl$slope

# 从基因型矩阵中提取 eQTL SNP
geno_matrix <- fread("genotype_matrix.csv")  # SNP矩阵：行是样本，列是SNP
geno_eqtl <- geno_matrix[, eqtl_snps, with = FALSE]

# 使用 LASSO 构建预测模型
lasso_model <- cv.glmnet(as.matrix(geno_eqtl), as.numeric(brca2_expr),
                         alpha = 1, family = "gaussian")

# 获得预测系数
coef_lasso <- as.matrix(coef(lasso_model, s = "lambda.min"))
predicted_expr <- as.matrix(geno_eqtl) %*% coef_lasso[-1]
# 合并 GWAS 和预测表达数据
twas_data <- data.frame(
  SNP = rownames(coef_lasso)[-1],
  weight = coef_lasso[-1],
  gwas_beta = gwas_data$beta[match(rownames(coef_lasso)[-1], gwas_data$variant_id)],
  gwas_pval = gwas_data$pval[match(rownames(coef_lasso)[-1], gwas_data$variant_id)]
)

# 计算 TWAS Z-score
twas_data <- twas_data %>%
  mutate(z_score = weight * gwas_beta / sqrt(weight^2 + gwas_beta^2),
         twas_pval = 2 * pnorm(-abs(z_score)))

# 筛选显著 TWAS 结果
significant_twas <- twas_data %>% filter(twas_pval < 0.05)
library(ggplot2)

# 曼哈顿图展示 BRCA2 TWAS 结果
ggplot(twas_data, aes(x = SNP, y = -log10(twas_pval))) +
  geom_point(color = "blue") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "TWAS Results for BRCA2",
       x = "SNP",
       y = "-log10(TWAS p-value)")
#####