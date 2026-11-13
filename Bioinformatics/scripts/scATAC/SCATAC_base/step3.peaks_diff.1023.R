# /home/shazhi11/miniforge3/envs/Signac_atac/lib/R/library
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(ensembldb)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)

library(future)
library(future.apply)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 100 * 1024^3)
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


#1.文件准备和整合注释
sc_ob = readRDS("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/scRNA_GSE278450/Glioma_3B_Fully_Annotated.rds")
seurat_ob = readRDS("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/tsinghua_atac/02_recluster/recluster.activate.rds")

dir.create("03_peaks_diff")
setwd("03_peaks_diff")

sc_ob$celltype <- gsub(" ","-",sc_ob$CellType)
Idents(sc_ob) <- "celltype"
if(!"data" %in% Layers(sc_ob@assays$RNA)){
	sc_ob = NormalizeData(sc_ob, assay = "RNA")
}

DefaultAssay(seurat_ob) <- "RNA"

#2. 识别整合的anchors，使用scRNA-seq的数据作为参照 
transfer.anchors <- FindTransferAnchors(reference = sc_ob, query = seurat_ob, 
								features = VariableFeatures(object = sc_ob), 
								reference.assay = "RNA", 
								query.assay = "RNA", 
								reduction = "cca") 
save(transfer.anchors,file = "./transfer.anchors.rda")

## 使用TransferData函数进行数据转移映射 
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = sc_ob$celltype,
								weight.reduction = seurat_ob[["lsi"]] ,  dims = 2:30) 
## 添加预测信息 
seurat_ob <- AddMetaData(seurat_ob, metadata = celltype.predictions) 
write.csv(as.data.frame(seurat_ob@meta.data),"seurat_ob.meta.csv",quote=F,row.names=T)

p1 <- DimPlot(seurat_ob, group.by = "predicted.id", label = TRUE, repel = TRUE ,cols = colors) + 
		ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop = FALSE) 
p2 <- DimPlot(sc_ob, group.by = "celltype", label = TRUE, repel = TRUE  ,cols = colors) + ggtitle("scRNA-seq cells") + NoLegend() 
p_merge = p1 + p2 
ggsave(plot = p_merge,"scATAC_RNA.celltype.pdf",w = 12,h = 6)
ggsave(plot = p_merge,"scATAC_RNA.celltype.png",w = 12,h = 6)

p3 = VlnPlot(seurat_ob, features = "prediction.score.max", group.by = "predicted.id",cols = colors) 
ggsave(plot = p3,"scATAC_score.celltype.pdf",w = 14,h = 6)
ggsave(plot = p3,"scATAC_score.celltype.png",w = 14,h = 6)
DefaultAssay(seurat_ob) <- "peaks"

#3.celltype 的 peaks差异分析
seurat_ob$celltype <- seurat_ob$predicted.id
Idents(seurat_ob) <- "celltype"
print(table(seurat_ob$celltype))
tumor_list= c("AC-like","MES-like","NPC-like","OPC-like")
normal_list = "Normal-Glial"

peak_list = list()
for(diff_type in tumor_list){
	da_peaks <- FindMarkers(object = seurat_ob, ident.1 = diff_type, ident.2 = normal_list,
	only.pos = F,test.use = "wilcox", min.pct = 0.05 )
	da_peaks$cluster <- diff_type
	da_peaks$peaks_loc <- rownames(da_peaks)
	peak_list[[diff_type]] <- da_peaks
}
merged_list <- do.call(rbind, peak_list)
write.csv(merged_list,"allpeaks.celltype.csv",quote=F,row.names=F)

filter_merged_list = merged_list  %>%
	dplyr::group_by(cluster) %>%
	dplyr::filter(avg_log2FC>= 1 & p_val < 0.05 & pct.1 > 0.1 ) %>%
	dplyr::arrange(p_val,desc(avg_log2FC),p_val)
write.csv(filter_merged_list,"filter_0.05_FC1.peaks.celltype.csv",quote=F,row.names=F)

dir.create("cluster_diffpeaks")
for(diff_type in tumor_list){
	sub_df = filter_merged_list[which(filter_merged_list$cluster == diff_type),]
	peaks_loc = head(sub_df$peaks_loc,n = 3)
	plot1 <- VlnPlot(  object = seurat_ob,  features = peaks_loc,  pt.size = 0.1,idents = c(diff_type,normal_list),ncol = 3 ,cols = colors)
	plot2 <- FeaturePlot(object = seurat_ob,features = peaks_loc, pt.size = 0.1 , ,ncol = 3)
	p_merge = plot1 / plot2
	ggsave(plot = p_merge,paste0("cluster_diffpeaks/",diff_type,".top3.peaks.pdf"),w = 15,h = 8)
	ggsave(plot = p_merge,paste0("cluster_diffpeaks/",diff_type,".top3.peaks.png"),w = 15,h = 8)
}

#4. 关注基因的peaks 差异
# find DA peaks overlapping gene of interest
target_genes = c("RFC3","RBBP4","WDHD1","AURKB","MORC4","CDC45","CENPH","MCM3","CNOT10","PCNA","DSCC1","PSMC4","BST1","SRSF1","ZWILCH","MCM2","TAF13","TMED5","YKT6")
for(diff_type in tumor_list){
	sub_peaks = filter_merged_list[which(filter_merged_list$cluster == diff_type),]$peaks_loc
	new_dir = paste0("cluster_diffpeaks/",diff_type)
	dir.create(new_dir,recursive= T)
	for(gene_hot in target_genes){
		print(gene_hot)
		regions_gene <- LookupGeneCoords(seurat_ob, gene_hot)  #不过滤直接绘图
		## 细胞片段丰度图
		tile_plot <- TilePlot(  object = seurat_ob,  region = regions_gene,  idents = c(diff_type, normal_list),tile.cells = 500)
		ggsave(plot = tile_plot,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot_tile_plot.pdf"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
		ggsave(plot = tile_plot,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot_tile_plot.png"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
		
		regions_highlight <- subsetByOverlaps(StringToGRanges(sub_peaks), LookupGeneCoords(seurat_ob, gene_hot)) # 取交集绘图，可信度更高一些
		if(length(regions_highlight) >= 1){
			p1 = CoveragePlot( object = seurat_ob, region = regions_gene, region.highlight = regions_highlight,idents = c(diff_type,normal_list),
						extend.upstream = 1000,extend.downstream = 1000 )
			ggsave(plot = p1,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot_hilight.pdf"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
			ggsave(plot = p1,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot_hilight.png"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
		} else {
			print("Not Found regions_highlight, Plotting genes ")
			p1 = CoveragePlot( object = seurat_ob, region = gene_hot,idents = c(diff_type,normal_list),
						extend.upstream = 1000,extend.downstream = 1000 )
			ggsave(plot = p1,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot.pdf"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
			ggsave(plot = p1,paste0(new_dir,"/",diff_type,".",gene_hot,".CoveragePlot.png"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
		}
	}
}

#5. motif 富集分析
library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38) ## 染色体号对应
library(chromVAR)
library(ggseqlogo)
set.seed(1234)

# pfm <- readJASPARMatrix("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar_filter.txt", matrixClass = "PFM")  下游分析报错

library(JASPAR2024)
library(TFBSTools)
library(RSQLite)
## 1. 建立与 JASPAR2024 数据库的连接
jaspar_db <- JASPAR2024()
sq_conn <- dbConnect(SQLite(), db(jaspar_db))
## 2. 提取脊椎动物（vertebrates）的 CORE 核心集 PFM
pfm <- getMatrixSet(
  x = sq_conn,
  opts = list(species = 9606,
	collection = "CORE",	   # 核心高质量库
	tax_group = "vertebrates", # 脊椎动物分类
	all_versions = TRUE	   # 仅获取最新版本
  )
)
dbDisconnect(sq_conn)
# 查看获取到的基序对象结构
print(pfm)

# ## 不用数据集，联网
# pfm <- getMatrixSet(
  # JASPAR2020,
  # opts = list(species = 9606, collection = "CORE")  # 人 + CORE 集
# )

# 再不济， ”/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/pfm.hsa.rda“


# ## 只保留标准染色体  不然 计算 GC 会报错
# gr <- StringToGRanges(rownames(seurat_ob))
# keep <- as.character(seqnames(gr)) %in% paste0("chr", c(1:22, "X", "Y", "M"))
# gr <- gr[keep]
# seurat_ob <- subset(seurat_ob, features = rownames(seurat_ob)[keep])

## 至少 5% 细胞里可及的 peak 才保留
# seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
# counts <- GetAssayData(seurat_ob, assay = "peaks", layer = "counts")
# keep <- rowSums(counts > 0) >= ncol(counts) * 0.05
# seurat_ob <- subset(seurat_ob, features = rownames(counts)[keep])

seurat_ob <- AddMotifs( object = seurat_ob, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm )
saveRDS(seurat_ob,"celltype.addmotif.rds")

#6. motif 富集分析
dir.create("motif_enrichment")
for(diff_type in tumor_list){
	new_dir2 = paste0("motif_enrichment/",diff_type,"-vs-",normal_list)
	dir.create(new_dir2,recursive=T)
	
	da_peaks <- FindMarkers(object = seurat_ob, ident.1 = diff_type, ident.2 = normal_list,
	only.pos = T,test.use = "wilcox", min.pct = 0.2,logfc.threshold = 0.25 ) #严卡过滤
	da_peaks$cluster <- diff_type
	da_peaks$peaks_loc <- rownames(da_peaks)
	library(dplyr)
	da_peaks  = da_peaks %>% arrange(desc(avg_log2FC),p_val) 
	write.csv(da_peaks,paste0(new_dir2,"/allpeaks.csv"),quote=F,row.names=F)

	filter_peaks  = da_peaks %>% group_by(cluster) %>% 
			filter(p_val < 0.05 & pct.1 > 0.2 ) %>%
			#filter( abs(avg_log2FC) > 0.25 ) %>%
			arrange(desc(avg_log2FC),p_val) 
	write.csv(filter_peaks,paste0(new_dir2,"/filter-p_0.005-fc_0.25.peaks.csv"),quote=F,row.names=F)
	top.da.peak <- filter_peaks[filter_peaks$p_val < 0.005 & filter_peaks$pct.1 > 0.2, ]$peaks_loc

	open.peaks <- AccessiblePeaks(seurat_ob, idents = c(diff_type, normal_list))
	meta.feature <- GetAssayData(seurat_ob, assay = "peaks", layer = "meta.features")
	peaks.matched <- MatchRegionStats(
		meta.feature = meta.feature[open.peaks, ],
		query.feature = meta.feature[top.da.peak, ],
		n = 50000)
	enriched.motifs <- FindMotifs(  object = seurat_ob, features = top.da.peak,background=peaks.matched)
	
	p = MotifPlot( object = seurat_ob, motifs = head(rownames(enriched.motifs),n = 8))
	ggsave(plot = p,paste0(new_dir2,"/motifs_enriched.top8.pdf"),w = 14,h = 8)
	ggsave(plot = p,paste0(new_dir2,"/motifs_enriched.top8.png"),w = 14,h = 8)
	enriched.motifs_anno = read.delim("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/motif2gene.txt",sep="\t")
	merge_df = merge(enriched.motifs,enriched.motifs_anno,all.x=T,by.x="motif",by.y = "motif")
	write.csv(merge_df,paste0(new_dir2,"/motifs_enriched.csv"),quote=F,row.names=F)
}

####============================================================================================================####
# #5. 计算基序活性
# ## seurat_ob <- RunChromVAR( object = seurat_ob, genome = BSgenome.Hsapiens.UCSC.hg38 )   1.17版本移除该函数
# library(chromVAR)
# library(SummarizedExperiment)
# library(Matrix)

# counts <- GetAssayData(object = seurat_ob, assay = "peaks", layer = "counts")
# motif_matrix <- GetMotifData(object = seurat_ob, assay = "peaks")

# library(BSgenome.Hsapiens.UCSC.hg38)
# genome <- BSgenome.Hsapiens.UCSC.hg38

# se <- SummarizedExperiment( assays = list(counts = counts), rowRanges = granges(seurat_ob))
# se <- keepStandardChromosomes(se, pruning.mode = "coarse")
# se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
# # motif.matrix <- matchMotifs(pfm, subject = rowRanges(se),
							# # genome = BSgenome.Hsapiens.UCSC.hg38)

# bg <- getBackgroundPeaks(object = se)
# dev <- computeDeviations(object = se, annotations = motif_matrix , background_peaks = bg)
# variability <- computeVariability(dev)
# motif_activity <- deviationScores(dev)   # 或 t(dev@assays$data$deviations)
# # deviation scores:行是 motif,列是细胞
# da <- deviationScores(dev)
# seurat_ob[["chromvar"]] <- CreateAssayObject( data = as.matrix(da), min.cells = 0, min.features = 0)

# # 注意:motif 名作为 feature 名
# rownames(seurat_ob[["chromvar"]]) <- rownames(da)
# DefaultAssay(seurat_ob) <- 'chromvar'

# saveRDS(seurat_ob,"celltype.addmotif.chromvar.rds")

library(chromVAR)
library(SummarizedExperiment)
library(Matrix)

counts <- GetAssayData(object = seurat_ob, assay = "peaks", layer = "counts")
keep_peaks <- Matrix::rowSums(counts) > 0
se <- subset(seurat_ob, features = names(keep_peaks)[keep_peaks])

# counts <- GetAssayData(object = se, assay = "peaks", layer = "counts")
motif_matrix <- GetMotifData(object = se, assay = "peaks")

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

se <- SummarizedExperiment( assays = list(counts = GetAssayData(object = se, assay = "peaks", layer = "counts")), 
						rowRanges = granges(se))
se <- keepStandardChromosomes(se, pruning.mode = "coarse")
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)

bg <- getBackgroundPeaks(object = se)
dev <- computeDeviations(object = se, annotations = motif_matrix , background_peaks = bg)
variability <- computeVariability(dev)
motif_activity <- deviationScores(dev)   # 或 t(dev@assays$data$deviations)
seurat_ob[["chromvar"]] <- CreateAssayObject( data = as.matrix(motif_activity), min.cells = 0, min.features = 0)

# 注意:motif 名作为 feature 名
rownames(seurat_ob[["chromvar"]]) <- rownames(motif_activity)
DefaultAssay(seurat_ob) <- 'chromvar'

saveRDS(seurat_ob,"celltype.addmotif.chromvar.rds")

####================================================================================================================####

## motif 基序活性分析
dir.create("motif_activate")
for(diff_type in tumor_list){
	new_dir3 = paste0("motif_activate/",diff_type,"-vs-",normal_list)
	dir.create(new_dir3,recursive=T)
	
	da_peaks <- FindMarkers(object = seurat_ob, ident.1 = diff_type, ident.2 = normal_list,
	only.pos = T,test.use = "wilcox", min.pct = 0.2,logfc.threshold = 0.25 ) #严卡过滤
	da_peaks$cluster <- diff_type
	da_peaks$motif <- rownames(da_peaks)
	library(dplyr)
	da_peaks  = da_peaks %>% arrange(desc(avg_log2FC),p_val)  
	write.csv(da_peaks,paste0(new_dir3,"/allmotif.csv"),quote=F,row.names=F)

	filter_motif  = da_peaks %>% group_by(cluster) %>% 
			  filter(p_val < 0.05 & pct.1 > 0.2 ) %>%
			  #filter( abs(avg_log2FC) > 0.25 ) %>%
			  arrange(desc(avg_log2FC),p_val) 
	write.csv(filter_motif,paste0(new_dir3,"/filter-p_0.005-fc_0.25.motif.csv"),quote=F,row.names=F)
	top.da.motif <- filter_motif[filter_motif$p_val < 0.005 & filter_motif$pct.1 > 0.2, ]$motif

	p = MotifPlot( object = seurat_ob, motifs = head(top.da.motif,n = 8) , assay= "peaks")
	ggsave(plot = p,paste0(new_dir3,"/motifs_loc.top8.pdf"),w = 14,h = 8)
	ggsave(plot = p,paste0(new_dir3,"/motifs_loc.top8.png"),w = 14,h = 8)
	
	for(ii in head(top.da.motif,n = 8)){
		p2 <- FeaturePlot( object = seurat_ob, features = ii, min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1,order = T)
		ggsave(plot = p2,paste0(new_dir3,"/motifs_activate.",ii,".feature.pdf"),w = 7,h = 7)
		ggsave(plot = p2,paste0(new_dir3,"/motifs_activate.",ii,".feature.png"),w = 7,h = 7)
	}

	enriched.motifs_anno = read.delim("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/database/motif2gene.txt",sep="\t")
	merge_df = merge(filter_motif,enriched.motifs_anno,all.x=T,by.x="motif",by.y = "motif")
	write.csv(merge_df,paste0(new_dir3,"/motifs_activate.csv"),quote=F,row.names=F)
}

