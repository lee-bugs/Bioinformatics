# /home/shazhi11/miniforge3/envs/Signac_atac/lib/R/library
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(ensembldb)
library(harmony)
library(future)
library(future.apply)

plan("multicore", workers = 10)


seurat_ob = readRDS("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/all_cells/03_diffpeaks/celltype.addmotif.chromvar.rds")  ##RDS 包含 assays 有 （RNA,peaks,chromvar）
DefaultAssay(seurat_ob) <- "peaks"

target_genes = c("RFC3","RBBP4","WDHD1","AURKB","MORC4","CDC45","CENPH","MCM3","CNOT10","PCNA","DSCC1","PSMC4","BST1","SRSF1","ZWILCH","MCM2")
target_TF = c("E2F1","SPIB","SPIC","XRN2","MYBL2","TP53","MYBL2","SUZ12")

# 1. TF 共表达分析

for(gene_hot in target_TF){
		print(gene_hot)
		motif.name <- ConvertMotifID(seurat_ob, name = gene_hot)
		if(!is.na(motif.name)){
		    gene_plot <- FeaturePlot(seurat_ob, features = gene_hot, reduction = 'umap',assay = "RNA" ,order = T)
		    motif_plot <- FeaturePlot(seurat_ob, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'umap',order = T)
		    p_merge = gene_plot | motif_plot
		    ggsave(plot = p_merge,paste0(gene_hot,".TF2motif.pdf"),w = 10,h = 5)
		    ggsave(plot = p_merge,paste0(gene_hot,".TF2motif.png"),w = 10,h = 5)
		}
}


# footprinting information for sets of motifs 非常慢
use_genelist= c(target_genes, target_TF)
motif_list = c(as.vector(unlist(seurat_ob@assays$peaks@motifs@motif.names)),names(seurat_ob@assays$peaks@motifs@motif.names))

plot_gene = intersect(use_genelist,motif_list)
library(BSgenome.Hsapiens.UCSC.hg38)
seurat_ob <- Footprint( object = seurat_ob, motif.name = plot_gene, genome = BSgenome.Hsapiens.UCSC.hg38)
for(gene_hot in plot_gene){
    p2 <- PlotFootprint(seurat_ob, features = gene_hot)
    ggsave(plot = p2,paste0(gene_hot,".Footprint.pdf"),w = 8,h = 5)
    ggsave(plot = p2,paste0(gene_hot,".Footprint.png"),w = 8,h = 5)
}

###############================================================== OK ================================================###

# 2. 基因共线性 https://zhuanlan.zhihu.com/p/565182343
library(SeuratWrappers)
library(cicero)   # remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
library(Matrix)
library(monocle3)
DefaultAssay(seurat_ob) = "peaks"
bone.cds <- as.cell_data_set(x = seurat_ob)

expr_matrix <- exprs(bone.cds) 
expr_matrix <- as(expr_matrix, "sparseMatrix")
umap_coords <- reducedDims(bone.cds)$UMAP
cell_meta <- as.data.frame(colData(bone.cds))
gene_meta <- as.data.frame(rowData(bone.cds))

bone.cds.clean <- new_cell_data_set(
  expression_data = expr_matrix,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
reducedDims(bone.cds.clean)$UMAP <- umap_coords
bone.cicero <- make_cicero_cds(bone.cds.clean, reduced_coordinates = reducedDims(bone.cds.clean)$UMAP)

genome_list <- seqlengths(seurat_ob)

# annotations <- Annotation(object = seurat_ob[['peaks']])
# isgene <- annotations$gene_name == target_genes
# isgene <- !is.na(x = isgene) & isgene
# annot.sub <- annotations[isgene]
# annot.sub[which(annot.sub$gene_name == 'E2F1'),]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome_list), "length" = genome_list)

# run cicero
conns <- run_cicero(bone.cicero, genomic_coords = genome.df, sample_num = 500)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat_ob) <- links

# annotations <- Annotation(object = seurat_ob[['peaks']])
# isgene <- annotations$gene_name == target_TF
# isgene <- !is.na(x = isgene) & isgene
# annot.sub <- annotations[isgene]

tumor_list= c("AC-like","MES-like","NPC-like","OPC-like")
normal_list = "Normal-Glial"
for(diff_type in tumor_list){
	for(gene_hot in c(target_genes,target_TF)){
		regions_gene <- LookupGeneCoords(seurat_ob, gene_hot)  #不过滤直接绘图
		p1 = CoveragePlot( object = seurat_ob, region = regions_gene, idents = c(diff_type,normal_list),
					extend.upstream = 1000,extend.downstream = 1000 )
		ggsave(plot = p1,paste0(diff_type,".",gene_hot,".CoveragePlot_Links.pdf"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
		ggsave(plot = p1,paste0(diff_type,".",gene_hot,".CoveragePlot_Links.png"),w = 10,h = 3 + length(unique(seurat_ob$celltype)) * 0.4)
	}
}

plan("sequential")
