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

library(edgeR)

# 创建 DGEList 对象，用于存储基因表达数据和组信息，还是使用原始计数矩阵
d <- DGEList(counts = expr, group = group)

# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2
# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)
table(keep)
# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]
# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)

# 归一化，TMM 方法
# 注意：归一化并不会直接在counts数值上修改，而是归一化系数会被自动存在 d$samples$norm.factors
d <- calcNormFactors(d)
head(d$samples)

# 将归一化后的数据赋值给 dge 变量
dge = d

# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
# edgeR 涉及到差异表达分析的函数有很多：exactTest、glmFit、glmLRT、glmQLFit、glmQLFTest。 
# qCML估计离散度需要搭配 exact test 进行差异表达分析，对应 exactTest 函数。 
# 而其他四个glm都是与GLM模型搭配使用的函数。其中，glmFit 和 glmLRT 函数是配对使用的，
# 用于 likelihood ratio test (似然比检验)，而 glmQLFit和 glmQLFTest则配对用于 quasi-likelihood F test (拟极大似然F检验)。
fit <- glmFit(dge, design)

# 注意这里的 contrast 和 DESeq2 不一样，这里我们只需要输入 c(-1, 1) 即可
# -1 对应 normal，1 对应 tumor
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG <- topTags(lrt, n = nrow(dge))

# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)
head(DEG_edgeR)

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_edgeR, file = './data/DEG_edgeR.Rdata')

#  火山图和热图
load("./data/DEG_edgeR.Rdata")
DEG_edgeR$symbol = rownames(DEG_edgeR)
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 3
P.Value = 0.01
k1 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)

p <- ggplot(data = DEG_edgeR, 
            aes(x = logFC, 
                y = -log10(PValue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()+
  geom_text_repel(data = filter(DEG_edgeR, change != "stable") %>% group_by(change) %>% top_n(16, abs(logFC)),
                  mapping = aes(label = symbol),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  max.overlaps =30,
                  show.legend = FALSE)
p

ggsave(filename = "./figure/volcano_plot_edgeR.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()

# 差异基因热图
deg_opt <- DEG_edgeR %>% filter(DEG_edgeR$change != "stable")
expr_heatmap <- expr %>% filter(rownames(expr) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(expr_heatmap) 

p1 <- pheatmap(expr_heatmap, show_colnames = T, show_rownames = T,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_edgeR.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()
