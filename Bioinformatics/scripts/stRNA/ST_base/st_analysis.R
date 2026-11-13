# 加载必要的包
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(future)
library(future.apply)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 100 * 1024^3)

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


# 设置 Space Ranger 输出结果的根目录
#0.合并 多样本peaks
# sample_dirs <- c("sample1/outs/", "sample2/outs/", "sample3/outs/")
# sample_names <- c("S1", "S2", "S3")

sample_info = read.delim("sample_info.txt",header=T)
sample_list = unique(sample_info$sample)

seurat_list <- list()
if(length(sample_list) > 10){
    print("Too Many Cells , Start Low test")
    for (i in sample_list) {
        sample_dirs = i
        obj = Load10X_Spatial(data.dir = i, assay = "Spatial")
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
        sample_dirs <- paste0(sample_id)
        obj = Load10X_Spatial(data.dir = sample_dirs, assay = "Spatial")
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

names(combined_ob@images) <- sample_list

#########1.qc##########
# 1. 计算线粒体基因比例
combined_ob[["percent.mt"]] <- PercentageFeatureSet(combined_ob, pattern = "^MT-")
p1 = VlnPlot(combined_ob, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), pt.size = 0) +  plot_layout(nrow = 1)
gbm_obj <- subset(combined_ob, subset = nFeature_Spatial > 200 & nCount_Spatial > 500)
p2 = VlnPlot(gbm_obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), pt.size = 0) +  plot_layout(nrow = 1)

combined_plot <- p1 / p2
ggsave(plot = combined_plot, "before_after.QC.pdf", w = 4 + 3 * 1.5 , h = 10)
ggsave(plot = combined_plot, "before_after.QC.png", w = 4 + 3 * 1.5 , h = 10)


# 2. 标准化 (SCTransform 或 标准 LogNormalize)
# 推荐使用 SCTransform，能够更好地去除技术偏差并稳定后续高变基因的选择
gbm_obj <- SCTransform(gbm_obj, assay = "Spatial", verbose = FALSE)

# 降维聚类
gbm_obj <- RunPCA(gbm_obj, assay = "SCT", verbose = FALSE)
gbm_obj <- FindNeighbors(gbm_obj, reduction = "pca", dims = 1:30)
gbm_obj <- FindClusters(gbm_obj, verbose = FALSE, resolution = 0.6)
gbm_obj <- RunUMAP(gbm_obj, reduction = "pca", dims = 1:30)

color_list = length(unique(gbm_obj$seurat_clusters))
new_colpal = my36colors[1:color_list]
names(new_colpal) = levels(gbm_obj$seurat_clusters)

# 1. 在 UMAP 空间中展示聚类
p1 = DimPlot(gbm_obj, reduction = "umap", label = TRUE) + scale_color_manual(values = new_colpal )
ggsave(plot = p1, "umap.pdf", w = 7, h = 7)
ggsave(plot = p1, "umap.png", w = 7, h = 7)


p2 = SpatialDimPlot(gbm_obj, label = F, label.size =1,ncol = 3) + 
		NoLegend() & 
		scale_fill_manual(values = new_colpal )
ggsave(plot = p2, "Spatial.umap.pdf", w = 12, h = 8)
ggsave(plot = p2, "Spatial.umap.png", w = 12, h = 8)

saveRDS(gbm_obj,"st.merge.rds")

# 寻找空间高变基因（Spatial Variable Genes）
# gbm_obj <- FindSpatiallyVariableFeatures( gbm_obj,   assay = "SCT", 
	# selection.method = "moransi", 
	# escape.matrix = FALSE, 
	# # features = VariableFeatures(gbm_obj)[1:1000]
# )

# target_genes = c("RFC3","RBBP4","WDHD1","AURKB","MORC4","CDC45","CENPH","MCM3","CNOT10","PCNA","DSCC1","PSMC4","BST1","SRSF1","ZWILCH","MCM2","TAF13","TMED5","YKT6")
# target_TF = c("E2F1","SPIB","SPIC","XRN2","MYBL2","TP53","MYBL2","SUZ12")

# filter_gene = intersect(c(target_genes,target_TF),rownames(gbm_obj))
################## 2. 基因可视化 ##############################
target_genes = c("CENPH","MCM3")
filter_gene = intersect(target_genes,rownames(gbm_obj))

gbm_obj = NormalizeData(gbm_obj,assay="Spatial")
DefaultAssay(gbm_obj) <- 'SCT'
for(my_gene in filter_gene){
	print(my_gene)
	for( sample_id in unique(gbm_obj$sample)){
		dir.create(paste0(sample_id,"_Plot"))
		print(sample_id)
		sub_ob = gbm_obj[,gbm_obj@meta.data$sample %in% sample_id]
		p = SpatialFeaturePlot(sub_ob, my_gene, min.cutoff = 0.5,image.alpha= 0.7,pt.size.factor= 1.7)
		ggsave(plot = p,paste0(sample_id,"_Plot/",my_gene,".SpatialFeature.pdf"),w = 10,h = 10)
		ggsave(plot = p,paste0(sample_id,"_Plot/",my_gene,".SpatialFeature.png"),w = 10,h = 10)
	}
}

for(my_gene in filter_gene){
	print(my_gene)
	for( sample_id in unique(gbm_obj$sample)){
		new_dir = paste0("Plot/",sample_id)
		dir.create(new_dir, recursive= T)
		print(sample_id)
		sub_ob = gbm_obj[,gbm_obj@meta.data$sample %in% sample_id]
		p = SpatialFeaturePlot(sub_ob, my_gene, min.cutoff = 0, max.cutoff = 1.5 ,image.alpha= 0.7,pt.size.factor= 1.9)
		ggsave(plot = p,paste0(new_dir,"/",my_gene,".SpatialFeature.pdf"),w = 10,h = 10)
		ggsave(plot = p,paste0(new_dir,"/",my_gene,".SpatialFeature.png"),w = 10,h = 10)
	}
}


DefaultAssay(gbm_obj) <- 'Spatial'
gbm_obj = UpdateSeuratObject(gbm_obj)
gbm_obj = JoinLayers(gbm_obj,assay="Spatial")

sub_exp = gbm_obj@assays$Spatial$counts[target_genes,]
exp_mat = as.data.frame(t(sub_exp))
names(exp_mat) <- paste0(target_genes,"_counts")
merge_ob = AddMetaData(gbm_obj,exp_mat)

merge_ob$State <- "Neither detected"
merge_ob@meta.data[which(merge_ob@meta.data$CENPH_count > 0 & merge_ob@meta.data$MCM3_count == 0),"State"] <- "CENPH only"
merge_ob@meta.data[which(merge_ob@meta.data$CENPH_count == 0 & merge_ob@meta.data$MCM3_count > 0),"State"] <- "MCM3 only"
merge_ob@meta.data[which(merge_ob@meta.data$CENPH_count > 0 & merge_ob@meta.data$MCM3_count > 0),"State"] <- "CENPH & MCM3"


color_list = c("#E64B35","#30B35A","#E5BD22","#CACACA")
names(color_list) <- c("CENPH only","MCM3 only","CENPH & MCM3","Neither detected")

merge_ob@meta.data$State <- factor(merge_ob@meta.data$State,levels = c("CENPH only","MCM3 only","CENPH & MCM3","Neither detected"))

for( sample_id in unique(merge_ob$sample)){
	new_dir = paste0("Plot_State/",sample_id)
	dir.create(new_dir, recursive= T)
	print(sample_id)
	sub_ob = merge_ob[,merge_ob@meta.data$sample %in% sample_id]
	####### 统计数量
	summ_cells <- sub_ob@meta.data[c("sample", "State")] %>% 
            dplyr::group_by(sample,State) %>%
            dplyr::summarize(cell_number = dplyr::n()) %>%
            dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100)
	write.table(as.data.frame(summ_cells), paste0(new_dir,"/4State_Counts_and_Pct.xls"), sep="\t", col.names=T, row.names =F)
	p <- ggplot(summ_cells, aes(x = State, y = cell_number, fill = State)) +
		geom_bar(stat = "identity", width = 0.6, color = NA) +
		geom_text(aes(label = paste0(cell_number, " (", round(freq,digits= 1), "%)")),
					vjust = -0.6,  size = 4,  fontface = "plain"  ) +
		scale_fill_manual( values = color_list) +
		theme( legend.position = "none",
			legend.title = element_blank(),
			axis.text = element_text(size = 11, color = "black"),
			axis.title = element_text(size = 12, face = "bold"),
			plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
			plot.subtitle = element_text(size = 11, hjust = 0.5),
			panel.background = element_blank(),
			panel.grid = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black", linewidth = 0.8)) +
		labs( title = "Spatial | CENPH / MCM3 co-detection",
		subtitle = paste0("Total Spots: ", sum(summ_cells$cell_number)), y = "Number of Spot", x = NULL ) +
		theme(axis.ticks.x = element_blank())
	ggsave(plot = p,paste0(new_dir,"/coexp.Number.pdf"),w = 7,h = 5)
	ggsave(plot = p,paste0(new_dir,"/coexp.Number.png"),w = 7,h = 5)
	
	p = SpatialDimPlot(sub_ob,group.by="State" , cols= , pt.size.factor= 1.7 , image.alpha= 0.7, crop= F) & 
					scale_fill_manual(values = color_list)
	ggsave(plot = p,paste0(new_dir,"/coexp.State.pdf"),w = 7,h = 7)
	ggsave(plot = p,paste0(new_dir,"/coexp.State.png"),w = 7,h = 7)
}

summ_cells <- merge_ob@meta.data[c("State")] %>% 
            dplyr::group_by(State) %>%
            dplyr::summarize(cell_number = dplyr::n()) %>%
            dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100)
write.table(as.data.frame(summ_cells), "all.4State_Counts_and_Pct.xls", sep="\t", col.names=T, row.names =F)












# ### RCTD 注释
# # /home/shazhi11/Lin/Glioma/Data/20260930/GSE182109(scRNA)(Object0_r+nd)/V0_Cell_divert_20260930/5.细胞分群/Glioma_3B_Fully_Annotated.rds

# library(doParallel)
# library(ggplot2)
# library(spacexr)
# # library(SpatialExperiment)
# # library(SummarizedExperiment)


# st_ob = readRDS("")
# # 构建一个对象
# puck <- SpatialRNA(coords, counts, nUMI)
# puck

# barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
# plot_puck_continuous(puck, barcodes, puck@nUMI, 
                   # ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                   # size = 2,
                   # title ='plot of nUMI') 




# sc_ob= readRDS("/home/shazhi11/Lin/Glioma/Data/20260930/GSE182109(scRNA)(Object0_r+nd)/V0_Cell_divert_20260930/5.细胞分群/Glioma_3B_Fully_Annotated.rds")

# reference <- Reference(counts, cell_types, nUMI)

# # 运行
# myRCTD <- create.RCTD(puck, reference, max_cores = 10)
# myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
# saveRDS(myRCTD,'myRCTD_visium_full.rds')