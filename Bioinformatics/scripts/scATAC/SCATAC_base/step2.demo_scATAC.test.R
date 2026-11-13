#/usr/bin/env Rscript
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(ensembldb)
library(harmony)
library(scDblFinder)
library(BSgenome.Hsapiens.UCSC.hg38)

library(future)
library(future.apply)
plan("multicore", workers = 10)

#0.合并 多样本peaks
# sample_dirs <- c("sample1/outs/", "sample2/outs/", "sample3/outs/")
# sample_names <- c("S1", "S2", "S3")
sample_info = read.delim("sample_info.txt",header=T)
sample_list = unique(sample_info$sample)

## cellranger atac 完整文件读取 peaks
gr_list <- list()
for (i in seq_along(sample_list)) {
    bed <- list.files(paste0(sample_list[i],"/outs"), pattern = "peaks.bed$", full.names = TRUE)
    df <- read.table(bed, col.names = c("chr", "start", "end"),comment.char = "#")
    gr_list[[i]] <- makeGRangesFromDataFrame(df)
}

all_peaks <- Reduce(c, gr_list)
combined_peaks <- reduce(all_peaks)
keep <- as.character(seqnames(combined_peaks)) %in% paste0("chr", c(1:22, "X", "Y", "M"))
combined_peaks <- combined_peaks[keep]

### 过滤过宽或过窄的 peak
combined_peaks <- combined_peaks[width(combined_peaks) < 10000 & width(combined_peaks) > 20]

dir.create("01_QC")
#1.读取样本数据，创建Seurat对象
## 这里需要看情况，有的GEO 给的是 fragments.tsv.gz ,有的给的是矩阵，根据情况对应预处理。
seurat_list <- list()
if(length(sample_list) > 10){
    print("Too Many Cells , Start Low test")
    for (i in sample_list) {
        sample_dirs = paste0(i,"/outs")
        # metadata <- read.csv(
            # file.path(sample_dirs[i], "singlecell.csv"),
            # row.names = 1
        # )
        # if ("NO_BARCODE" %in% rownames(metadata)) {
            # metadata = metadata[!grepl("^NO_BARCODE", rownames(metadata)),]
        # }
        #filter_cells <- readLines(file.path(sample_dirs,"barcodes.tsv.gz"))
        counts <- FeatureMatrix(
            fragments = CreateFragmentObject(file.path(sample_dirs, "fragments.tsv.gz")),
            features = combined_peaks, cells = NULL
            # cells = filter_cells
            #cells = rownames(metadata[metadata$passed_filters > 200, ]) # 官方默认阈值500
        )
        # 单样本可以  counts_h5 <- Read10X_h5(filename = file.path(sample_dirs, "filtered_peak_bc_matrix.h5"))
        chrom_assay <- CreateChromatinAssay(
            counts = counts, min.cells = 10, min.features = 200,
            fragments = file.path(sample_dirs, "fragments.tsv.gz"),
            genome = "hg38"
        )
        obj <- CreateSeuratObject(
            counts = chrom_assay,
            assay = "peaks",
            project = i
        )
        obj$sample <- i
		obj$group <- sample_info[which(sample_info$sample == i),][,'group']
        #obj <- AddMetaData(obj, metadata)
        seurat_list[[i]] <- obj
    }
} else if(length(sample_list) < 10 ){
    print("Not Too Many Cells , Start Fats test")
    # 以下多样本并行慎用，单样本估计在15G内存，太多样本可能会爆内存。
    # 只开一次并行，用于“多样本循环”
    plan("multisession", workers = 4)
    options(future.globals.maxSize = 200 * 1024^3)
    
    # 把每个样本的处理封装成函数
    process_one_sample <- function(sample_id) {
        sample_dirs <- paste0(sample_id, "/outs")
        #filter_cells <- readLines(file.path(sample_dirs, "barcodes.tsv.gz"))
	    
        counts <- FeatureMatrix(
          fragments = CreateFragmentObject( file.path(sample_dirs, "fragments.tsv.gz") ),
          features = combined_peaks, cells = NULL )
            # cells = filter_cells
	    
        chrom_assay <- CreateChromatinAssay( 
          counts    = counts,  min.cells = 10, min.features = 200,
          fragments = file.path(sample_dirs, "fragments.tsv.gz"),
          genome    = "hg38")
	    
        obj <- CreateSeuratObject(
          counts  = chrom_assay,
          assay   = "peaks",
          project = sample_id
        )
        obj$sample <- sample_id
		obj$group <- sample_info[which(sample_info$sample == sample_id),][,'group']
        obj
    }
    # 并行处理所有样本
    seurat_list <- future_lapply(sample_list, process_one_sample)
    names(seurat_list) <- sample_list
    # 跑完记得关掉并行
    plan("sequential")
}
print("Read Raw Sample Files Done")
combined_ob <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_list )

print(table(combined_ob$sample))

## 只保留标准染色体  不然 计算 GC 会报错
gr <- StringToGRanges(rownames(combined_ob))
keep <- as.character(seqnames(gr)) %in% paste0("chr", c(1:22, "X", "Y", "M"))
gr <- gr[keep]
combined_ob <- subset(combined_ob, features = rownames(combined_ob)[keep])
##  过滤小于1%表达的peaks
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
counts <- GetAssayData(combined_ob, assay = "peaks", layer = "counts")
keep <- rowSums(counts > 0) >= ncol(counts) * 0.01
combined_ob <- subset(combined_ob, features = rownames(counts)[keep])

edb <- EnsDb("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/ens.hsa.v113.sql")
annotations <- GetGRangesFromEnsDb(ensdb = edb)

if (!length(grep("chr|Chr", seqlevels(annotations))) > 0) {
    seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
}
tag_info = unique(genome(combined_ob[["peaks"]])[1:5])  # 提取 版本信息，对应上才可以正常运行。
#print()
#if(unique(genome(annotations)) != tag_info){
#    genome(annotations) <- tag_info
#}
genome(annotations) <- tag_info
print(seqlevels(annotations))
Annotation(combined_ob) <- annotations

saveRDS(combined_ob,"01_QC/raw.combined.rds")

#2 计算QC指标
combined_ob <- NucleosomeSignal(object = combined_ob) ## 计算核小体条带分布模式
combined_ob$nucleosome_group <- ifelse(combined_ob$nucleosome_signal > 4, "NS > 4", "NS < 4")
combined_ob$nucleosome_group = as.factor(combined_ob$nucleosome_group)
# pdf("01_QC/fragment_histogram.pdf", width = 10, height = 5)
p = FragmentHistogram(object = combined_ob, group.by = "nucleosome_group")
ggsave(plot = p, "01_QC/fragment_histogram.pdf", w = 12, h = 6)
ggsave(plot = p, "01_QC/fragment_histogram.png", w = 12, h = 6)

# checkpoint1

# dev.off()
# 这张图来自 Signac 官方 PBMC 教程，是用 FragmentDistribution() 画的片段长度分布直方图，按 核小体信号（Nucleosome Signal, NS） 的高低把细胞分成两组对比展示，本质是一张质控（QC）判读图。
# 图里画了什么
# 两张并列直方图：左红（NS < 4，核小体信号低）和右青（NS > 4，核小体信号高）。
# 横轴 Fragment length (bp)：每条被测到的 DNA 片段长度。
# 纵轴 Count：该片段的计数。注意左右纵轴量级差很多（左到 2 万，右仅约 100），说明高质量细胞贡献了绝大部分有效片段，低信号那组本身片段就很少。
# 怎么解读（核心）
# ATAC-seq 里，Tn5 酶只能切开没有核小体遮挡的开放区域，因此一个细胞测到的片段长度会呈现特征性的"阶梯"：
# ~50 bp：subnucleosomal（核小体之间的开放片段）
# ~180 bp：mononucleosomal（单个核小体保护的长度）——这是主峰
# ~360 bp / ~540 bp：di-/tri-nucleosomal（二聚体、三聚体）
# 左图（NS 低）：能在 ~180、~360、~540 bp 处看到清楚的多个峰，即典型的 nucleosome laddering（核小体梯状条带），说明文库构建干净、细胞完整，是高质量细胞。
# 右图（NS 高）：梯状峰消失或很弱，长片段（>1200 bp，对应核小体信号）占比偏高，通常是死亡/破损细胞或 ambient 游离 DNA 污染的特征，应当被过滤。
# 它和你的 QC 流程怎么对应
# Signac 的核小体信号定义为长片段（>1200 bp）占比，值越低越好。
# 这张图就是教你"看图识质量"：分布有清晰梯状峰→留；分布塌成一片、长片段多→弃。
# 实操上一般对 nucleosome_signal 设阈值（教程常用 < 2；图中用 <4/>4 只是示意两组的形态差异），结合 nFrags、TSS enrichment、blacklist ratio 一起过滤。
# 一句话总结：左图是健康细胞的"教科书式"片段分布，右图是该被丢掉的低质量细胞的分布，两图对比就是为了让你判断 nucleosome_signal 的合理过滤线画在哪里。

combined_ob <- TSSEnrichment(object = combined_ob, fast = FALSE , process_n = 4000 ) ## 计算每个细胞的TSS富集分数,FAST 模式会少layers
combined_ob$high.tss <- ifelse(combined_ob$TSS.enrichment > 2, "High", "Low")
#pdf("01_QC/tss_plot.pdf", width = 10, height = 5)
p2 = TSSPlot(combined_ob, group.by = "high.tss") + NoLegend()
ggsave(plot = p2, "01_QC/tss_plot.pdf", w = 12, h = 6)
ggsave(plot = p2, "01_QC/tss_plot.png", w = 12, h = 6)
#dev.off()
# 当fast=T时，仅计算TSS分数，而不存储每个细胞Tn5插入频率的位置矩阵。
# ATAC 数据必须在 TSS 上有峰。因为启动子区（TSS 周围）通常是开放染色质，片段应该在这里富集。峰越尖越高，数据质量越好。
# 图中两组细胞的差别本质上就是按 TSS 富集度（TSS.enrichment 这个 QC 指标）分出来的高低两组
# 这通常是 TSSEnrichment() 算完 QC 后按中位数或阈值划分的结果，用来对比验证指标是否真能区分好坏细胞。
# "Low" 组（比如 TSS.enrichment < 3~5 甚至更低）通常代表低质量细胞或空液滴/死细胞，在下游分析中应该过滤掉。


# 利用密度散点图可视化nCount_peaks与TSS.enrichment之间的关系
pdf("01_QC/density_scatter.pdf", width = 10, height = 8)
DensityScatter(combined_ob,
    x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE,
    quantiles = TRUE
) # 显示5%、10%、90%、95%分位数的提示线
dev.off()
png("01_QC/density_scatter.png", width = 1000, height = 800)
DensityScatter(combined_ob,
    x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE,
    quantiles = TRUE
) # 显示5%、10%、90%、95%分位数的提示线
dev.off()


# 一般认为TSS enrichment score小于4数据质量比较差(这个也可以调整，2、6、8的阈值我也见过)，考虑建库或上游的分析存在问题。
# nCount_peaks数值高则表明该细胞具有较多的开放染色质区域，反映出数据覆盖范围良好，
# 那么nCount_peaks与TSS.enrichment应该呈现一个正相关的关系且在数据的集中区域有一个高亮的”斑点”。
# 数据分布比较随机、没有高亮斑点，则考虑建库或上游的分析存在问题。

#combined_ob$pct_reads_in_peaks <- combined_ob$peak_region_fragments / combined_ob$passed_filters * 100

load("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/ens.hsa.blacklist.rda")
combined_ob$blacklist_ratio <- FractionCountsInRegion(
    object = combined_ob,
    assay = "peaks",
    regions = blacklist_regions)

## 过滤QC ，阈值参考 https://www.sc-best-practices.org/chromatin-accessibility/quality-control/ 和signac
all_filter = c("nCount_peaks", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal", "pct_reads_in_peaks")
real_filter = intersect(all_filter,colnames(combined_ob@meta.data))

p1 = VlnPlot(
    object = combined_ob,
    features = real_filter,
    pt.size = 0.1,
    ncol = length(real_filter))

filter.combined_ob <- subset( x = combined_ob,
    subset = nCount_peaks > 3000 &
        nCount_peaks < 100000 &
        #pct_reads_in_peaks > 40 &
        blacklist_ratio < 0.2 &
        nucleosome_signal < 4 &
        TSS.enrichment > 4)
if("pct_reads_in_peaks" %in% real_filter){
    filter.combined_ob <- subset(  x = combined_ob,
    subset =  pct_reads_in_peaks > 40 )
}

p2 <- VlnPlot(
    object = filter.combined_ob,
    features = real_filter,
    pt.size = 0.1,
    ncol = length(real_filter))

combined_plot <- p1 / p2
ggsave(plot = combined_plot, "01_QC/before_after.QC.pdf", w = 4 + 3 * length(real_filter) , h = 10)
ggsave(plot = combined_plot, "01_QC/before_after.QC.png", w = 4 + 3 * length(real_filter) , h = 10)

#ccc = AddMetaData(filter.combined_ob,meta_df)

## 过滤双胞
sce <- as.SingleCellExperiment(filter.combined_ob)
sce <- scDblFinder(  sce,
    aggregateFeatures = TRUE,   # scATAC 必加
    nfeatures = 25,             # 聚合后的 meta-feature 数，默认即可
    samples = "sample",      # 如果有多样本，必须指定样本列
    #BPPARAM = BiocParallel::MulticoreParam(10)  # 可选：多线程
)

table(sce$scDblFinder.class)

filter.combined_ob$scDblFinder.class <- sce$scDblFinder.class
filter.combined_ob$scDblFinder.score <- sce$scDblFinder.score
write.csv(as.data.frame(filter.combined_ob@meta.data),"01_QC/QC.meta.csv",sep=",",quote=F)

filter.combined_ob <- subset(filter.combined_ob, subset = scDblFinder.class == "singlet")
saveRDS(filter.combined_ob,"01_QC/filter.combined.rmdoublet.rds")
print(table(filter.combined_ob$sample))
# colors = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
			# "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
			# "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
			# "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
			# "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
colors <- c(
  "#2E2300", "#7F6838", "#FFD700", "#FFB6C1", "#DB7093", 
  "#8B0000", "#FF4500", "#FF8C00", "#FFDAB9", "#8B4513", 
  "#D2B48C", "#228B22", "#006400", "#98FB98", "#00FA9A", 
  "#008B8B", "#00FFFF", "#4682B4", "#000080", "#1E90FF", 
  "#9370DB", "#8A2BE2", "#4B0082", "#EE82EE", "#FF1493", 
  "#FF69B4", "#CD5C5C", "#F08080", "#FA8072", "#E9967A", 
  "#9ACD32", "#556B2F", "#6B8E23", "#3CB371", "#20B2AA", 
  "#66CDAA", "#7B68EE", "#6A5ACD", "#483D8B", "#C71585", 
  "#D8BFD8", "#DDA0DD", "#BA55D3", "#9932CC", "#9400D3", 
  "#4169E1", "#6495ED", "#B0C4DE", "#708090", "#2F4F4F"
)

# 3.降维聚类分析
dir.create("02_recluster")
filter.combined_ob <- RunTFIDF(filter.combined_ob)
filter.combined_ob <- FindTopFeatures(filter.combined_ob, min.cutoff = "q0")
filter.combined_ob <- RunSVD(filter.combined_ob)
## 通常在进行完线性降维后，得到的第一个LSI所代表的通常为技术差异，而非生物学差异，我们一般将第一个LSI去除。
## 可以使用DepthCor函数评估每个LSI成分与测序深度之间的相关性
pdf("02_recluster/depthCor.pdf",w = 10,h = 7)
DepthCor(filter.combined_ob)
dev.off()
png("02_recluster/depthCor.png",w = 1000,h = 700)
DepthCor(filter.combined_ob)
dev.off()

pdf("02_recluster/density_scatter.pdf", width = 10, height = 8)
DensityScatter(filter.combined_ob,
    x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE,
    quantiles = TRUE
) # 显示5%、10%、90%、95%分位数的提示线
dev.off()
png("02_recluster/density_scatter.png", width = 1000, height = 800)
DensityScatter(filter.combined_ob,
    x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE,
    quantiles = TRUE
) # 显示5%、10%、90%、95%分位数的提示线
dev.off()

# 从这个图上可以发现，第一个LSI与细胞的总计数counts 之间存在非常强的相关性，因此，我们将去除这个组分后再执行下游步骤
batch_var <- "sample" # 去批次选择的列信息
if (length(unique(filter.combined_ob@meta.data[[batch_var]])) > 1) {
    filter.combined_ob <- RunHarmony(  filter.combined_ob, group.by.vars = batch_var, reduction.use= "lsi",project.dim = FALSE)
    reduct <- "harmony"; 
    rdims <- seq_len(ncol(Embeddings(filter.combined_ob, reduct)))
} else {
    reduct <- "lsi"; 
    rdims <-  2:30 
}

filter.combined_ob <- RunUMAP(object = filter.combined_ob, reduction = reduct, dims = rdims)
filter.combined_ob <- RunTSNE(object = filter.combined_ob, reduction = reduct, dims = rdims)
filter.combined_ob <- FindNeighbors(object = filter.combined_ob, reduction = "lsi", dims = rdims)
filter.combined_ob <- FindClusters(object = filter.combined_ob, verbose = FALSE, algorithm = 3)
p1 = DimPlot(object = filter.combined_ob, label = F, cols = colors ) 
ggsave(plot = p1, "02_recluster/umap.pdf", w = 7*1.2, h = 7)
ggsave(plot = p1, "02_recluster/umap.png", w = 7*1.2, h = 7)


p1 <- DimPlot(object = filter.combined_ob, label = F,reduction = "tsne" , cols = colors )
ggsave(plot = p1, "02_recluster/tsne.pdf", w = 7*1.2, h = 7)
ggsave(plot = p1, "02_recluster/tsne.png", w = 7*1.2, h = 7)

#4.创建GeneActivity，函数自动执行创建基因活性表达矩阵这一步骤
gene.activities <- GeneActivity(filter.combined_ob) # 得到活性矩阵
# 将活性矩阵添加到pbmc对象中，作为一个新的assays，并进行标准化。
# add the gene activity matrix to the Seurat object as a new assay and normalize it
filter.combined_ob[["RNA"]] <- CreateAssayObject(counts = gene.activities)

seurat_ob <- NormalizeData(
    object = filter.combined_ob,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(filter.combined_ob$nCount_RNA))

DefaultAssay(seurat_ob) <- "RNA"
saveRDS(seurat_ob,"02_recluster/recluster.activate.rds")

print(seurat_ob)
