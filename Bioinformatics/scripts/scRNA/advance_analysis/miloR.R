# /home/shazhi11/miniforge3/envs/miloR
# 1.加载包和数据
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
seurat_obj = readRDS('/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/sc_GSE182109/NPC_like_CENPH_MCM3_4groupdayu0.rds')
seurat_obj$sample <- seurat_obj$SampleID
seurat_obj <- as.SingleCellExperiment(seurat_obj )

# 2.建立一个图形和邻域
# d：用于KNN细化的约简维数。建议使用同样的方案d用于 KNN 图构建，或通过检查碎石图来选择 PC。
# k：这会影响DA测试的效力，因为我们需要从每个样本中足够多的细胞在邻域内进行代表性采样，以估计重复样本之间的方差。
# 另一方面，k值过大可能导致过度平滑。我们建议首先使用与KNN图构建聚类和UMAP可视化时相同的k值。
# 之后我们将采用一些启发式方法来评估是否应增加k值。

reduc_list = reducedDimNames(seurat_obj )
if("HARMONY" %in% reduc_list){
  reduced.dim = "HARMONY"
}else if("PCA" %in% reduc_list){
  reduced.dim = "PCA"
}else {
  reduced.dim = reduc_list[1]
}

milo.obj <- buildGraph(seurat_obj, k=30, d=30, reduced.dim = reduced.dim)

#prop：用于随机采样初始细胞的比例。我们建议在细胞数量少于3万时使用该参数。对于更大的数据集，使用该参数已足够（且可加快计算速度）。  
# prop=0.1  ，prop=0.05  
#refined：表示是否要使用采样细化算法，还是仅随机选择细胞。默认且推荐的做法是使用细化算法。
milo.obj<- makeNhoods(milo.obj, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims =  reduced.dim )

# plotNhoodSizeHist(milo.obj)
# 查看 nhood 分布，平均面积建议超过5倍N_samples，且在50到100达到峰值

# 3.计算距离，根据实验设计对细胞计数，并执行DA测试。
milo.obj <- countCells(milo.obj, meta.data = as.data.frame(colData(milo.obj)), sample="sample")
milo_design <- data.frame(colData(milo.obj))[,c("sample", "CENPH_MCM3_group")]

rownames(milo_design) <- NULL
milo_design <- distinct(milo_design )
print(milo_design)
print("222")
rownames(milo_design ) <- milo_design$sample

# 计算邻域连通性
# Milo采用了cydar引入的空间FDR修正的改良，我们修正了邻域间重叠的p值。
#具体来说，每个假设检验的P值都被k个最近邻距离的倒数加权。
#要使用该统计量，首先需要存储Milo对象中最近邻之间的距离。
#这由函数完成 （注：这一步是分析流程中最耗时的步骤，对于大型数据集可能只需几分钟。）

milo.obj <- calcNhoodDistance(milo.obj, d=30, reduced.dim = reduced.dim)

# 4.计算差异
da_results <- testNhoods(milo.obj , design = ~ CENPH_MCM3_group, design.df = milo_design , reduced.dim= reduced.dim)
head(da_results)

# 每个邻域计算折叠变化和修正P值，指示发育阶段间是否存在显著丰度差异。
# logFC表示分组间细胞数的log2FC变化
# PValue：在FDR修正前报告P值
# SpatialFDR：报告P值经过多次测试校正，考虑邻域间的重叠

da_results  <- da_results %>%
  arrange(SpatialFDR) 


# 5.结果展示  
# 柱状图和散点图

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
milo.obj <- buildNhoodGraph(milo.obj)

## Plot single-cell UMAP
umap_pl <- plotReducedDim( milo.obj, dimred = "UMAP",
                          colour_by="group", text_by = "celltype",
                          text_size = 3, point_size=0.5) + guides(fill="none")

## Plot neighbourhood graph
#邻域网络图：每个圈代表一个邻域（Nhood），邻域的大小与其中细胞的数量成比例。
#颜色梯度（从蓝色到红色）表示了 logFC（对数折叠变化）：
#  正的 logFC（蓝色）意味着该邻域中的细胞在对照组中富集，
#  负的 logFC（红色）则表明该邻域中的细胞在处理组中富集。
#邻域之间的连线表示邻域之间的相似性或重叠，线条的厚度表示邻域之间的重叠度（overlap size），
#越厚的线表示两个邻域之间的细胞更为相似。 

nh_graph_pl <- plotNhoodGraphDA(milo.obj, da_results,
                                layout="UMAP",alpha=0.7)  #alpha默认0.1

p1 = umap_pl + nh_graph_pl + plot_layout(guides="collect")
ggsave(plot = p1,"umap_nlgraph.pdf",w = 10, h = 6)

da_results <- annotateNhoods(milo.obj, da_results, coldata_col = "CENPH_MCM3_group")
head(da_results)

ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)
da_results$CENPH_MCM3_group <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$CENPH_MCM3_group)
plotDAbeeswarm(da_results, group.by = "CENPH_MCM3_group")

# 6.寻找DA群体的marker
## Add log normalized count to Milo object
milo.obj <- logNormCounts(milo.obj)
print("333")
da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(milo.obj, da_results, 
                                          subset.row = rownames(seurat_obj ), 
                                          aggregate.samples = TRUE, # 建议设置为TRUE
                                          sample_col = "sample")

#在许多情况下，DA邻域位于KNN图的不同区域，将所有重要的 DA 群体分组在一起可能并不理想,
#因为它们可能包括非常不同细胞类型的细胞。对于这种场景,开发者用了一个邻里函数,它使用社区检测根据
#(1)两个邻域之间共享的小区数量将邻域划分为组;
#(2)DA 社区的倍数变化方向;
#(3)折叠变化的差异。

## Run buildNhoodGraph to store nhood adjacency matrix
milo.obj <- buildNhoodGraph(milo.obj)

## Find groups
da_results <- groupNhoods(milo.obj, da_results, 
                          max.lfc.delta = 2,da.fdr = 0.1)#da.fdr默认0，1
head(da_results)

p2 = plotNhoodGroups(milo.obj, da_results, layout="UMAP") 
ggsave(plot = p2,"NhoodGroups.pdf",w = 8, h = 8)
plotDAbeeswarm(da_results, "NhoodGroup",alpha = 0.1) #alpha默认是0.1

## Exclude zero counts genes
keep.rows <- rowSums(logcounts(milo.obj)) != 0
scRNA <- milo.obj[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(scRNA)
hvgs <- getTopHVGs(dec, n=3000)

# 找到每个邻里群的one-vs-all差异基因表达
set.seed(42)
da_results$NhoodGroup <- as.numeric(da_results$NhoodGroup)
scRNA@colData$CENPH_MCM3_group <- as.factor(scRNA@colData$CENPH_MCM3_group)
scRNA@colData$sample <- as.factor(scRNA@colData$sample)
# 只保留包含 ≥ 5 个 nhood 的组
nhood_per_group <- table(da_results$NhoodGroup)
keep_groups <- names(nhood_per_group[nhood_per_group >= 5])
da_results_filtered <- da_results[da_results$NhoodGroup %in% keep_groups, ]

nhood_markers <- findNhoodGroupMarkers(scRNA, da_results_filtered, subset.row = hvgs, 
                                       aggregate.samples = T, sample_col = "sample")

head(nhood_markers)

# 如果明确了相比较的邻域
# nhood_markers <- findNhoodGroupMarkers(scRNA, da_results, subset.row = hvgs,
#                                       subset.nhoods = da_results$NhoodGroup %in% c('11','2'),
#                                       aggregate.samples = TRUE, sample_col = "orig.ident")

markers <- nhood_markers$GeneID[nhood_markers$adj.P.Val_5 < 0.01 & nhood_markers$logFC_5 > 0]
#我们可以用 来可视化邻域中的表达式。plotNhoodExpressionGroups
print("444")
set.seed(42)
plotNhoodExpressionGroups(embryo_milo, da_results, features=intersect(rownames(embryo_milo), markers[1:10]),
                          subset.nhoods = da_results$NhoodGroup %in% c('6','5'), 
                          scale=TRUE,
                          grid.space = "fixed")


#在某些情况下，你可能想检测同一邻域组内不同条件下细胞间的差异表达。你可以用以下方式实现：testDiffExp
#dge_6 <- testDiffExp(embryo_milo, da_results, design = ~ stage, meta.data = data.frame(colData(embryo_milo)),
#                     subset.row = rownames(embryo_milo)[1:5], subset.nhoods=da_results$NhoodGroup=="6")

#dge_6