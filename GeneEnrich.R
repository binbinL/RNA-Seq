
library(tinyarray)
library(dplyr)
#install.packages('BiocManager') 
#BiocManager::install('clusterProfiler')
library(clusterProfiler)

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

# 定义函数挑选差异基因
deg_filter <- function(df){
  rownames(df)[df$change != "stable"]
}

#all_degs <- intersect(intersect(deg_filter(DEG_deseq2), deg_filter(DEG_edgeR)), deg_filter(DEG_limma_voom))
all_degs <- intersect(deg_filter(DEG_limma_voom), deg_filter(DEG_edgeR))


#已有OrgDb的常见物种#
#BiocManager::install("org.Hs.eg.db") 
library(org.Hs.eg.db)
#查看所有可转化类型
keytypes(org.Hs.eg.db) 

entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                      keys = all_degs, #将输入的gene_name列进行数据转换
                      keytype = "SYMBOL", #输入数据的类型
                      column = "ENTREZID",#输出数据的类型
                      multiVals = 'filter')#删除这些重复匹配

entrezid_all  = na.omit(entrezid_all)  #na省略entrezid_all中不是一一对应的数据情况
entrezid_all = data.frame(entrezid_all) #将entrezid_all变成数据框格式
head(entrezid_all)

###GO富集分析###
GO_enrich = enrichGO(gene = entrezid_all[,1], #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", #输入数据的类型
                     ont = "ALL", #可以输入CC/MF/BP/ALL
                     #universe = 背景数据集 # 表示背景基因，无参的物种选择组装出来的全部unigenes作为背景基因；有参背景基因则不需要。
                     pvalueCutoff = 0.05,qvalueCutoff = 0.05, #表示筛选的阈值，阈值设置太严格可导致筛选不到基因。可指定 1 以输出全部
                     readable = T) #是否将基因ID映射到基因名称。
GO_enrich  = data.frame(GO_enrich) #将GO_enrich导成数据框格式
#数据导出#
write.csv(GO_enrich,'./data/GO_enrich.csv') 

###KEGG富集分析###
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], #即待富集的基因列表
                         keyType = "kegg",
                         pAdjustMethod = 'fdr',  #指定p值校正方法
                         organism= "human",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         qvalueCutoff = 0.05, #指定 p 值阈值（可指定 1 以输出全部）
                         pvalueCutoff=0.05) #指定 q 值阈值（可指定 1 以输出全部）
KEGG_enrich  = data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'./data/KEGG_enrich.csv') #数据导出


library(ggplot2)
go_enrich = read.csv("./data/GO_enrich.csv", sep=',' , header=T, row.names = 1)
go_enrich$term <- paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T)
go_enrich_res =  go_enrich %>% group_by(ONTOLOGY) %>% top_n(n=5,wt=p.adjust)

#纵向柱状图#
ggplot(go_enrich_res, 
       aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()


#横向柱状图#
ggplot(go_enrich, 
       aes(x=term,y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  xlab("GO term") + #x轴标签
  ylab("Gene_Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

#气泡图#
ggplot(go_enrich_res,
       aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()


kegg_enrich = read.csv("./data/KEGG_enrich.csv", sep=',' , header=T, row.names = 1)
kegg_enrich$term <- paste(kegg_enrich$ID, kegg_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
kegg_enrich$term <- factor(kegg_enrich$term, levels = kegg_enrich$term,ordered = T)
kegg_enrich_res =  kegg_enrich %>%  top_n(n=12,wt=p.adjust)

#横向柱状图#
# ggplot(kegg_enrich, 
#        aes(x=term,y=Count, fill=ID)) +  #x、y轴定义；根据ONTOLOGY填充颜色
#   geom_bar(stat="identity", width=0.8) +  #柱状图宽度
#   scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
#   facet_grid(.~ID, scale = 'free_x', space = 'free_x')+
#   xlab("Pathways") + #x轴标签
#   ylab("Gene Ratio") +  #y轴标签
#   labs(title = "KEGG Terms Enrich")+ #设置标题
#   theme_bw() + 
#   theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

#气泡图-kegg#
ggplot(kegg_enrich_res,aes(y=term,x=GeneRatio))+
  
  geom_point(aes(size=Count,color=pvalue))+
  
  scale_color_gradient(low = "red",high ="blue")+
  
  labs(color=expression(pvalue,size="Count"),
       
       x="Gene Ratio",y="Pathways",title="KEGG Pathway Enrichment")+
  
  theme_bw()
