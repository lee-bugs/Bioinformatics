# /home/shazhi11/miniforge3/envs/R_Seurat
#Author :lip by 25/12/2024

#===========================0. command ===================
library("optparse")
option_list = list(
    make_option( c("--RDS", "-i"), type = "character", default = NULL,
                 help = "the seurat object saved as R object in RDS format."),
    make_option( c("--output","-o"),type="character", default = "scNMF_result",
                help="the output directory of results.", metavar="character"),
    make_option( c("--clusters","-c"),type="integer", default = 10,
                help="the number of top metaprograms for heatmap."),
    make_option( c("--species","-s"),type="character", default = "hsa",
        help="the species info,default is hsa ,ect hsa ,mmu , rno " ),
    make_option(c("--groupby","-q"),type="character", default = "group",
               help="Which groupby for split rds . e.g. clusters or sample or celltype"),
    make_option( c("--genes_num"),type="integer", default = 300,
                help="the number of each metaprograms genes ."),
    make_option( c("--color"),type="character", default = "viridis",
                help="the palette of heatmap , e.g.  viridis or heatmap."),
    make_option( c("--assays", "-a"), type = "character", default = "RNA",
                 help="assays of your data ")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library("tibble")
library("RColorBrewer")

library(ggpubr)
library(dplyr)
library(ggsci)

#===============================0. function load ===================================
colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
# source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/PSN_pipeline/enrichment.r")

myplotMetaPrograms <- function(mp.res,
                            similarity.cutoff=c(0,1),
                            scale = "none",
                            downsample = 20000,
                            showtree = TRUE,
                            palette = color_pal,
                            annotation_colors = NULL,
                            main = "Clustered Heatmap",
                            show_rownames = FALSE,
                            show_colnames = FALSE,
                            ...) {

  J <- mp.res[["programs.similarity"]]
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]
  labs.order <- labels(as.dendrogram(tree))
  
  #downsample, to avoid overloading the graphics
  if (length(cl_members) > downsample) {
    sid.keep <- downsampleMin(cl_members, size = downsample)
    labs.order <- labs.order[labs.order %in% sid.keep]
    
    cl_members <- cl_members[labs.order]
    J <- J[labs.order, labs.order]
    
    #disable 'tree' for downsampled heatmap
    showtree <- FALSE
  }
  
  cl_names <- names(cl_members)
  cl_members[!is.na(cl_members)] <- paste0("MP",cl_members[!is.na(cl_members)])
  names(cl_members) <- cl_names

  #Recover order of MP clusters
  cluster.order <- unique(cl_members[labs.order])
  nMP <- length(cluster.order)
  
  #Gaps for heatmap
  diffs <- diff(as.integer(as.factor(cl_members))) 
  gaps <- which(diffs != 0)
  
  #Annotation column
  annotation_col <- as.data.frame(cl_members)
  colnames(annotation_col) <- "Metaprogram"
  annotation_col[["Metaprogram"]] <- factor( cl_members, levels = cluster.order )
  print(head(annotation_col))

  #Apply trimming to similarity for plotting  
  J[J<similarity.cutoff[1]] <- similarity.cutoff[1]
  J[J>similarity.cutoff[2]] <- similarity.cutoff[2]
  write.csv(J,paste0(output_dir,"/geneNMF.cluster_heatmap.csv"),quote=F)
  if (!showtree) {
    tree <- FALSE
  }
  
  if(!is.null(annotation_colors)){
    color_df = list(Metaprogram = annotation_colors[1:length(unique(cl_members[labs.order]))])
	names(color_df$Metaprogram) <- gtools::mixedsort(unique(cl_members[labs.order]))
	# annotation_col_tmp <- annotation_col
	# annotation_col_tmp$facet <- rownames(annotation_col)
	# merge_df = merge(annotation_col_tmp , color_df,all.x = T , by = "Metaprogram")
	# ann_colors = list(Metaprogram = merge_df$colorset)
	# names(ann_colors$Metaprogram) <- merge_df$facet
  }
  print(color_df)
  
  ph <- pheatmap::pheatmap(J,
                 scale = scale,
                 color = palette,
                 main = main,
                 cluster_rows = tree,
                 cluster_cols = tree,
                 cutree_rows = nMP,
                 cutree_cols = nMP,
                 gaps_row = gaps,
                 gaps_col = gaps,
				 treeheight_row = 7 + nMP*0.1, treeheight_col = 7 + nMP*0.1,
                 annotation_col = annotation_col,
                 annotation_row = annotation_col,
                 annotation_colors = color_df,
                 annotation_names_col = FALSE,
                 annotation_names_row = FALSE,
                 show_rownames = show_rownames,
                 show_colnames = show_colnames,
                 ...
  )
  return(ph)
}  

# species <- opt$species
# #支持的物种列表，看看放在那里比较合适
# species_mapping <- list(
  # "hsa" = "org.Hs.eg.db", # human
  # "mmu" = "org.Mm.eg.db", # mouse
  # "rno" = "org.Rn.eg.db", # rat
  # "dme" = "org.Dm.eg.db", # fruit_fly
  # "dre" = "org.Dr.eg.db", # zebrafish
  # "ath" = "org.At.tair.db", # Arabidopsis
  # "sce" = "org.Sc.sgd.db", # yeast
  # "cel" = "org.Ce.eg.db", # C.elegans
  # "bta" = "org.Bt.eg.db", # Bovine
  # "mcc" = "org.Mmu.eg.db", # monkey
  # "cfa" = "org.Cf.eg.db", # dog
  # "ssc" = "org.Ss.eg.db", # pig
  # "gga" = "org.Gg.eg.db", # chicken
  # "xla" = "org.Xl.eg.db", # frog
  # "ptr" = "org.Pt.eg.db", # chimpanzee
  # "aga" = "org.Ag.eg.db"  # mosquito
# )

# print(species)
# if (species %in% names(species_mapping)) {
	# species_package <- species_mapping[[species]]
	# suppressMessages(library(species_package, character.only = TRUE , lib.loc = "/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/scanpy/lib/R/library"))
	# species.db <- get(species_package)
# }else{
	# species.db =NULL
# }

#=========================1. input file =================
#1.check file and param
if ( is.null(opt$RDS) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    my_test = readRDS(opt$RDS)
    my_test = UpdateSeuratObject(my_test)
}
print(dim(my_test))

DefaultAssay(my_test) <- opt$assays

if ( file.exists(opt$output)){
    output_dir = opt$output
}else{
    output_dir = opt$output
    dir.create(output_dir,recursive = TRUE)
}
my_group = opt$groupby
my_test = Seurat::SetIdent(my_test,value= my_group )
nclusters = opt$clusters
my_test = NormalizeData(my_test,assay = "RNA")

if(opt$color == "viridis"){
	color_pal = viridis(100, option="A", direction=-1)
}else {
	color_pal = colorRampPalette(c("grey85", "#FEA730", "#B62018"))(n=100)
}

##癌细胞确实在不同患者中具有独特的转录表型，使得技术和生物学变异性难以区分。
##通过 NMF 在个别患者中提取基因程序，因为它绕过了批次效应校正的需要,我们宁愿在基因程序活动水平上整合多个样本。

#=========================2. NMF analysis ===============

my_test.list <- SplitObject(my_test, split.by = my_group)
geneNMF.programs <- multiNMF(my_test.list, assay = "RNA", k = 4:10, min.exp = 0.05, nfeatures = 1000 )

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",       #成对相似性算法，cosine 或 Jaccard 指数
                                        weight.explained = 0.8,  #基因集解释比例，越高基因越多
                                        nMP = nclusters,
                                        min.confidence = 0.5,    #可信度，越高过滤掉越多
										max.genes= opt$genes_num,#模块基因数，默认200,基因太多可能会导致模块数少于预期聚类群数目
										hclust.method = "ward.D2"
										)
##聚类的算法，"ward.D"和"ward.D2"均表示采用ward离差平方和法, "single"表示最短距离法, "complete"表示最长距离法, "average" (= UPGMA)表示类平均法, "median" (= WPGMC) 表示中间距离，"centroid" (= UPGMC)表示重心法，默认值为ward.D2。

save(geneNMF.metaprograms,file = paste0(output_dir,"/geneNMF.metaprograms.rda"))

# plot_length = length(unique(names(geneNMF.metaprograms$metaprograms.genes.weights)))

ph <- myplotMetaPrograms(geneNMF.metaprograms, palette = color_pal,
						annotation_colors = colors,
						similarity.cutoff = c(0.1,1))
ggsave(plot = ph,paste0(output_dir,"/geneNMF.cluster_heatmap.pdf"),w = 10,h = 9)
ggsave(plot = ph,paste0(output_dir,"/geneNMF.cluster_heatmap.png"),w = 10,h = 9)

#=========================3. driver_gene list =============
gene2weight_list = list()
for(myterm in names(geneNMF.metaprograms$metaprograms.genes.weights)){
	MP_weight = as.data.frame(geneNMF.metaprograms$metaprograms.genes.weights[myterm])
	names(MP_weight) <- "weight"
	MP_weight$clusters <- myterm
	MP_weight$gene <- rownames(MP_weight)
	gene2weight_list[[myterm]] <- MP_weight
}

gene2weight <- do.call(rbind,gene2weight_list)
gene2weight = gene2weight[c("gene","clusters","weight")]
# if (!is.null(species.db)) {
	# gene2weight$description = unname(mapIds(x = species.db,
			# keys = gene2weight$gene,
			# column = "GENENAME",
			# keytype = "SYMBOL",
			# multiVals = "first"))
# }

write.table(gene2weight,paste0(output_dir,"/gene2weight_anno.xls"),sep="\t", quote = F , row.names = F , col.names = T)
print(table(gene2weight$clusters))

# #=========================4. compare with msigDB ===========
# library(msigdbr)
# library(fgsea)
# #msigdbr::msigdbr_species() 
# #
# top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
    # runGSEA(program, universe=rownames(my_test), 
            # category = "H", 
			# subcategory = NULL,
			# species="Danio rerio",
			# pval.thr= 0.05
			# )
# })

# top_p_list = list()
# for(myterm in names(top_p)){
	# MP_enrich = as.data.frame(top_p[[myterm]])
	# if(dim(MP_enrich)[1] > 0 ){
		# MP_enrich$clusters <- myterm
		# top_p_list[[myterm]] <- MP_enrich
	# }
# }

# top_p_df <- do.call(rbind,top_p_list)
# for( myterm in  unique(gene2weight$clusters)){
	# cluster_dir = file.path(output_dir,'Enrichment',myterm)
	# if(!file.exists(cluster_dir)){
	# dir.create(cluster_dir,recursive = TRUE)
	# }
	# cluster_markers = subset(gene2weight,clusters==myterm)
	# rownames(cluster_markers)<-cluster_markers$gene
	# if(nrow(cluster_markers)>1){
		# genelist = cluster_markers$gene
		# try(enrichment(species = opt$species, outDir = cluster_dir,geneList = genelist))
	# }
	# # write.table(cluster_markers,paste(cluster_dir,paste0(myterm,"_genes.xls"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)
# }


#=================================4. 基因程序的特征打分 ================================

mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(my_test, features = mp.genes, ncores=4, name = "")

add_list = intersect(unique(c("group","SampleID","seurat_clusters",my_group)),names(seu@meta.data))
for(i in add_list){
	p1= VlnPlot(seu, features=names(mp.genes), group.by = i, pt.size = 0, ncol=3,cols = colors )
	ggsave(plot = p1,paste0(output_dir,"/Ucell.splitby_",i,".pdf"),w = 12 + length(unique(seu[[i]])) * 0.3 ,h = 3 + 1.3 * length(names(mp.genes)))
	ggsave(plot = p1,paste0(output_dir,"/Ucell.splitby_",i,".png"),w = 12 + length(unique(seu[[i]])) * 0.3 ,h = 3 + 1.3 * length(names(mp.genes)),bg = "white")
}


matrix <- seu@meta.data[,names(mp.genes)]

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                         cell.embeddings = dimred,
                                         assay.used = "RNA",
                                         key = "MP_",
                                         global = FALSE)

set.seed(123)
seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")

p2 = DimPlot(seu, reduction = "umap_MP", group.by = my_group,cols = colors) + theme(aspect.ratio = 1)
ggsave(plot = p2,paste0(output_dir,"/DimPlot.by",my_group,".pdf" ))
ggsave(plot = p2,paste0(output_dir,"/DimPlot.by",my_group,".png" ))

saveRDS(seu,paste0(output_dir,"/MPs_umap.Score.rds" ))

library(viridis)
p3 = FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", ncol=4,pt.size = 0.5) &
				scale_color_viridis(option="B") &
				theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

ggsave(plot = p3,paste0(output_dir,"/FeaturePlot.MPs.pdf"), w = 4 + 0.6 * length(names(mp.genes)),h = 4 + 0.6 * length(names(mp.genes))/4)
ggsave(plot = p3,paste0(output_dir,"/FeaturePlot.MPs.png"), w = 4 + 0.6 * length(names(mp.genes)),h = 4 + 0.6 * length(names(mp.genes))/4,bg = "white")