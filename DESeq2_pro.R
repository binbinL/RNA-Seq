library(dplyr) 

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

library(DESeq2)
# 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
colData <- data.frame(row.names = colnames(expr),group = group)
colData$group <- factor(colData$group, levels = c("control", "treat"))
head(colData)

# DESeq2 要求输入数据是由整数组成的矩阵，且没有经过标准化

# 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
dds <- DESeqDataSetFromMatrix(countData = expr, # 表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息
head(dds)

# 进行差异表达分析
dds <- DESeq(dds)

# 查看结果的名称
resultsNames(dds)

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
DEG_deseq2 <- na.omit(DEG)

# 输出差异表达基因结果的前几行
head(DEG_deseq2)

save(DEG_deseq2, file = './data/DEG_deseq2.Rdata')

# 火山图和热图

load("./data/DEG_deseq2.Rdata")
DEG_deseq2$symbol = rownames(DEG_deseq2)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2
P.Value = 0.05
k1 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_deseq2$change)
#save(DEG_deseq2, file = './data/DEG_deseq2.Rdata')

#基因信息
#library(ggrepel)
#DEG_deseq2$label=ifelse(DEG_deseq2$pvalue < P.Value & abs(DEG_deseq2$log2FoldChange) >= logFC,rownames(DEG_deseq2),"")

# 火山图
p <- ggplot(data = DEG_deseq2, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()+
  geom_text_repel(data = filter(DEG_deseq2, change != "stable") %>% group_by(change) %>% top_n(6, abs(logFC)),
                  mapping = aes(label = symbol),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  max.overlaps =30,
                  show.legend = FALSE)


p

ggsave(filename = "./figure/volcano_plot_deseq2.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
library(pheatmap)
deg_opt <- DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
exp_heatmap <- expr %>% filter(rownames(expr) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_heatmap) 

p1 <- pheatmap(exp_heatmap, show_colnames = T, show_rownames = T,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_deseq2.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()
