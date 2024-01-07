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

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(expr)

# 创建 DGEList 对象
dge <- DGEList(counts = expr, group = group)

# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 使用线性模型进行拟合
fit <- lmFit(v, design)

# 需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con

# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_limma_voom, file = './data/DEG_limma_voom.Rdata')

# 画图图 —— 火山图和热图

load("./data/DEG_limma_voom.Rdata")
DEG_limma_voom$symbol = rownames(DEG_limma_voom)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.05
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)


# 火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()+
  geom_text_repel(data = filter(DEG_limma_voom, change != "stable") %>% group_by(change) %>% top_n(20, abs(logFC)),
                  mapping = aes(label = symbol),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  max.overlaps =30,
                  show.legend = FALSE)
p

ggsave(filename = "./figure/volcano_plot_limma_voom.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
deg_opt <- DEG_limma_voom %>% filter(DEG_limma_voom$change != "stable")
expr_heatmap <- expr %>% filter(rownames(expr) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(expr_heatmap) 

p1 <- pheatmap(expr_heatmap, show_colnames = T, show_rownames = T,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_limma_voom.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()
