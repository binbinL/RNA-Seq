# install.packages("tinyarray")
# install.packages('VennDiagram')
# install.packages("UpSetR")

library(UpSetR)
library(tinyarray)
library(dplyr)
library(ggplot2)
# 定义函数挑选差异基因
deg_filter <- function(df){
  rownames(df)[df$change != "stable"]
}


load("./data/DEG_edgeR.Rdata")
DEG_edgeR$symbol = rownames(DEG_edgeR)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2
P.Value = 0.05
k1 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)


load("./data/DEG_deseq2.Rdata")
DEG_deseq2$symbol = rownames(DEG_deseq2)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2
P.Value = 0.05
k1 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_deseq2$change)

load("./data/DEG_limma_voom.Rdata")
DEG_limma_voom$symbol = rownames(DEG_limma_voom)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2
P.Value = 0.05
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)

# 取交集
all_degs <- intersect(intersect(deg_filter(DEG_deseq2), deg_filter(DEG_edgeR)), deg_filter(DEG_limma_voom))

# 依据三个包得到的差异基因绘制韦恩图
all_degs_venn <- list(DESeq2 = deg_filter(DEG_deseq2), edgeR = deg_filter(DEG_edgeR), limma = deg_filter(DEG_limma_voom))

all_degs_venn <- draw_venn(all_degs_venn, "ALL DEGs")
all_degs_venn

ggsave(filename = "./figure/all_degs_venn.pdf", plot = all_degs_venn, device = "pdf", width = 5, height = 5)


# 读数据
data<-read.delim("data/gene_count.xls",row.names = 1, sep = '\t',check.names = FALSE)

data <- distinct(data, gene_name, .keep_all = T)
dim(data)

# 去重
# library(limma)
# exp_brca <- avereps(exp_brca, exp_brca$gene)

rownames(data)<-data$gene_name
dim(data)
expr<-data[,1:6]

group <- c(rep('control', 3), rep('treat', 3))
group <- factor(group, levels = c("control", "treat"))
table(group)

# 绘制共同差异基因的热图
expr_heatmap <- expr %>% filter(rownames(expr) %in% all_degs)
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(expr_heatmap) 

p1 <- pheatmap(expr_heatmap, show_colnames = T, show_rownames = T,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename="./figure/heatmap_plot_all_degs.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

#upSet图
all_degs_venn <- list(DESeq2 = deg_filter(DEG_deseq2), edgeR = deg_filter(DEG_edgeR), limma = deg_filter(DEG_limma_voom))
listInput<-all_degs_venn
p2<-upset(fromList(listInput), order.by = c("degree", "freq"), decreasing = c(TRUE, FALSE))
      # empty.intersections = "on")
p2
ggsave(filename="./figure/upset_plot_all_degs.pdf", plot = p2, device = "pdf", width = 5, height = 6)
dev.off()
