#.libPaths("/home/shazhi11/miniforge3/envs/R_Seurat/lib/R/library")
# Rscript scTenifoldKnk_test.R -i ./GSE182109.Object2_rGBM.rds  -x RFC3 -o ./GSE182109.Object2_rGBM/RFC3_knockdown_QC -s 人
#
library("optparse")
#=command line parameters setting=============================
option_list = list(
    make_option( c( "--input","-i"), type = "character",
        help = "the input rds file."),
    make_option( c("--net" ), type = "integer",  default = 10,
        help = "[Otional]the min net ."),
    make_option( c("--ncells"), type = "double",  default = 500,
        help = "[Otional]the min cells."),
    make_option( c("--ngenes"), type = "double",  default = 200,
        help = "[Otional]the min all genes."),
    make_option( c("--var_n"), type = "double",  default = 5000,
        help = " VariableFeatures gene choose ,default is 8000."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--assay","-a"),type="character", default = "RNA",
        help="the assays info,default is RNA." ),
    make_option( c("--groupby"), type = "character", default = "CellType",
        help = "konckdown gene "),
    make_option( c("--species","-s"),type="character", default = "人",
        help="the species info,  人 小鼠，大鼠 eg , default is NULL." ),
    make_option( c("--gmt_df"),type="character", default = "/home/shazhi11/Lin/Kidney_cancer/linshi/database/GeneSet_gmt/KEGG_hsa_pathways_symbol.gmt",
        help="the species gmt, default is NULL." ),
    make_option( c("--topn","-n"),type="integer", default = 20,
                help="gmt top TERM and end TERM plot.", metavar="character"),
    make_option( c("--ko_analysis"),type="character", default = NULL,
        help="the ko_analysis for second plot" ),
    make_option( c("--nCores"),type="double", default = 8,
        help="reverse the default trajectory if necessary,--reverse True or --reverse False" ),
    make_option( c("--koGene", "-x"), type = "character", default = NULL,
        help = "konckdown gene "),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "The subset of cluster ids used for subtyping.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#####################function
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/color/color.R")
#colors = colorls$"NPG"
colors = c("#3A4892","#ED0000","#008B45","#631879","#00817F","#BA0020","#5F549A",
"#A10055","#7F807F","#A93B52","#9F6D25","#545A62","#50557D","#8A5C4E",
"#9A3C5D","#883A77","#96546A","#626489","#7E4372","#C85113","#407354",
"#5C3B7B","#657167","#AC293F","#774989","#9D3860","#727285","#CC2C31",
"#6E7E35","#5E3F6D","#3C6C7E","#A54137","#83497B","#962866","#8D6B75",
"#51568D","#666D92","#ED7777","#468B68","#6E4979","#418180","#BA5D6D",
"#7D779A","#A1517B","#808080","#A9727D","#9F8662","#5B5E62","#67697D",
"#8A736C")

library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
library(data.table)
set.seed(123)

#1.get information

if ( file.exists(opt$output)){
    output_dir = opt$output
}else{
    output_dir = opt$output
    dir.create(output_dir, recursive = T)
}

data_ob = readRDS( opt$input )
DefaultAssay(data_ob) = opt$assay
data_ob = UpdateSeuratObject(data_ob)


if ( !is.null(opt$which_cells)){
    cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
	print(cluster_list)
    seurat_obj = data_ob[,data_ob@meta.data[,opt$groupby] %in% cluster_list ]
} else {
    seurat_obj = data_ob
}
print(seurat_obj)


var_n = opt$var_n
species = opt$species
topn = opt$topn
# ********目标基因变量 ********
target_gene_name <- opt$koGene     #目标基因的名称(需要修改)# 
#3、加载并预处理 Seurat 数据
# 加载Seurat的对象
# ******** Seurat对象变量 ********
# seurat_obj <- data_ob

# DefaultAssay(seurat_obj) <- "RNA"
genelist = rownames(seurat_obj)
pattern <- "^(AC[0-9]|AP[0-9]|AL[0-9]|C[1-9]+orf|LINC[0-9])"
is_drop <- grep(pattern, genelist,value=T)
filter_gene = setdiff(genelist,is_drop)
# print(filter_gene)
#提取表达矩阵
# ********表达矩阵变量 ********
count_matrix <- GetAssayData(seurat_obj, layer = "counts")
#提取高可变基因
# ********高可变基因相关变量 ********

seurat_obj <- FindVariableFeatures(object=seurat_obj[filter_gene,], selection.method="vst", nfeatures = as.numeric(var_n))
high_variable_genes <- VariableFeatures(seurat_obj)

# target_gene = c("RFC3","RBBP4","WDHD1","AURKB","MORC4","CDC45","CENPH","MCM3","CNOT10","PCNA","DSCC1","PSMC4","ZWILCH","MCM2","YKT6","TMED5","UBE2J2","TAF13","TBCE","TCOF1","UTP3","TMEM53","UNC45A","UXS1","TTC4")
target_gene = c("E2F1","SPIB","SPIC","XRN2","MYBL2","TP53","MYBL2","SUZ12")


# ********分析输入数据变量 ********
term_gene= unique(c(target_gene_name,target_gene))
filter_gene = intersect(term_gene,rownames(count_matrix))
print("Maybe some genes not in dataset ,such as below ...")
print(setdiff(term_gene,rownames(count_matrix)))

input_data <- as.data.frame(count_matrix[unique(c(filter_gene, high_variable_genes)),])
print(input_data[1:5,1:5])
if(!is.null(opt$ko_analysis)){
  ko_analysis = normalizePath(opt$ko_analysis)
  print(paste0("start load  scTenifoldKnk ",ko_analysis))
  load(ko_analysis)
} else {
  print("start analysis scTenifoldKnk ... ")
  ko_analysis_result <- scTenifoldKnk(countMatrix = as.matrix(input_data),nCores = opt$nCores ,
          gKO = target_gene_name,      #需要敲除的基因
          qc = FALSE,                #是否进行QC
          qc_maxMTratio = 0.25,        #mt的阈值
          qc_minLibSize = opt$ngenes,         #文库阈值(细胞测到的基因总数)
          nc_nNet = opt$net,              #子网络数量
          nc_nCells = opt$ncells            #每个网络中随机抽取的细胞数
  )
  save(ko_analysis_result,file = paste0(output_dir,"/ko_analysis_result.rda"))
}


#5、提取显著差异结果并保存

#输出差异分析的结果
# ********差异调控结果相关变量 ********
diff_regulation_df <- ko_analysis_result$diffRegulation
diff_regulation_df <- diff_regulation_df[which(diff_regulation_df$gene != target_gene_name), ]
write.table(diff_regulation_df,paste0(output_dir,"/diff_all.scTenifoldKnk.xls"),sep="\t",quote=F,row.names=F,col.names=T)
# Gene：基因名；
# Distance：这个基因在虚拟敲除（KO）网络与原始网络中的扰动程度（网络距离）。数值越大，说明该基因受KO影响越大;
# boxcox-transformed distance：为了满足统计假设，对原始 Distance 做 Box-Cox 变换后的值，通常用于后续Z统计量计算，使数据更接近正态分布;
# Z：标准化后的扰动量：衡量该基因的距离相对于其他基因的离群程度。Z 越大（正或负）表示扰动越明显;
# FC：反映虚拟敲除后基因表达或扰动的变化倍数;
# p-value: 原始P值;
# adjusted p-value: 多重假设校正后的P值。

significant_diff_table <- diff_regulation_df[diff_regulation_df$p.adj <= 0.05,]
if ("direction" %in% colnames(diff_regulation_df)) {
  significant_diff_table$Regulated <- significant_diff_table$direction
} else {
  significant_diff_table$Regulated <- "NoSig"
  significant_diff_table$Regulated <- ifelse(significant_diff_table$Z > 0, "Distance_up", "Distance_down")
}

write.table(significant_diff_table, paste0(output_dir,"/sigDiff_padj.scTenifoldKnk.xls"), sep="\t", quote=F, row.names=F)

#6、可视化结果（3 张图搞定差异分析）
#Top20 差异基因柱状图（看哪些基因受影响最大）
###########################柱状图
#绘制柱状图
# ******** top10差异基因变量 ********
top10_genes <- head(significant_diff_table[order(-significant_diff_table$FC), ], 10)
tail10_genes <- head(significant_diff_table[order(significant_diff_table$FC,decreasing = F), ], 10)
top20_diff_genes = rbind(top10_genes,tail10_genes)
top20_diff_genes = unique(top20_diff_genes)
p1=ggplot(top20_diff_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
#输出图形
pdf(file=paste0(output_dir,"/barplot.pdf"), width=6, height=5)
print(p1)
dev.off()

###########################火山图
#准备火山图的数据 
diff_regulation_df$log_p.adj <- -log10(diff_regulation_df$p.adj)
diff_regulation_df$significant <- ifelse(diff_regulation_df$p.adj < 0.05, "Significant", "Not significant")
labeled_significant_genes <- subset(diff_regulation_df, p.adj < 0.05)

# 绘制火山图
y_axis_upper_limit <- quantile(diff_regulation_df$log_p.adj, 0.999, na.rm = TRUE)
p2=ggplot(diff_regulation_df, aes(x=Z, y=log_p.adj, color=significant)) +
  geom_point(alpha=0.7, size=1) +  # 设置点的透明度和大小
  # 定义图形的颜色
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  #在图形中标注基因的名称
  geom_text_repel(data=labeled_significant_genes, aes(label=gene), size=3, max.overlaps=50) +
  labs(title="", x="Z-score", y="-log10(p.adj)") +
  theme_classic() + 
  # 设置纵坐标范围
  coord_cartesian(ylim = c(0, y_axis_upper_limit)) + theme(legend.position = "none")
#输出图形
pdf(file=paste0(output_dir,"/vol.pdf"), width=6, height=5)
print(p2)
dev.off()

###########################饼图
#准备饼图数据
# ********饼图相关数据变量 ********
significant_gene_count <- table(diff_regulation_df$significant)
significant_gene_df <- as.data.frame(significant_gene_count)
colnames(significant_gene_df) <- c("category", "count")
significant_gene_df$percentage <- paste0(round(significant_gene_df$count / sum(significant_gene_df$count) * 100, 1), "%")
# 绘制饼图
p3 <- ggplot(significant_gene_df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label=percentage), position=position_stack(vjust=0.5), size=4) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values=c("Significant"="red", "Not significant"="lightgray")) + 
  labs(title="Proportion of Significant Genes", fill="") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))
# 输出饼图
pdf(file=paste0(output_dir,"/pie.pdf"), width=6, height=5)
print(p3)
dev.off()

pdf(paste0(output_dir,"/koGene.net.pdf"),w = 10,h = 10 )
plotKO( ko_analysis_result, target_gene_name, q = 0.99, annotate = F, nCategories = 20, fdrThreshold = 0.05 )
dev.off()

##### 富集分析
if(!is.null(species)){
  print("start enrichment analysis")

  # all_json = fromJSON(file = "~/../s00/software/5.AdvancedScript/database/ref.json")
  # enrich_db = all_json[[species]][[1]]$desc[1]
  
  sig_df = significant_diff_table
  # sig_df$avg_log2FC <- log2(sig_df$FC)
  sig_df$cluster <- "KO"
   rownames(sig_df) <- NULL
  sig_df = tibble::rownames_to_column(sig_df,var = "ID")
  write.table(sig_df, paste0(output_dir,"/sigDiff_padj_avg_log2FC.scTenifoldKnk.xls"), sep="\t", quote=F, row.names=F)
  diff_file = paste0(output_dir,"/sigDiff_padj_avg_log2FC.scTenifoldKnk.xls")
  
  # enrich_output = paste0(output_dir,"/enrichment")
  # dir.create(enrich_output , recursive=T )
  
  # print(glue::glue(" source /PERSONALBIO/work/singlecell/s00/software/miniconda3/bin/activate /PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/scanpy && SeuratEnrichment.R --species {species} --allmarkers {diff_file} --trend up --output {enrich_output} --db  {enrich_db}"))
  # system(glue::glue(" source /PERSONALBIO/work/singlecell/s00/software/miniconda3/bin/activate /PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/scanpy && SeuratEnrichment.R --species {species} --allmarkers {diff_file} --trend up --output {enrich_output} --db  {enrich_db}"))
  # print("ehrichment done")

  #diff_all = paste0(output_dir,"/diff_all.scTenifoldKnk.xls")
  if(species %in% c("人","大鼠","小鼠")){
    print("start GSEA analysis")
    if(species =="人"){
      gsea_gmt = "/home/shazhi11/Lin/Kidney_cancer/linshi/database/GeneSet_gmt/CancerHallmarks_integrated.gmt"
    }else {
      gsea_gmt = "/生信共享文件/富集分析参考文件/mh.all.v2025.1.Mm.symbols.gmt"
    }
  }
  if(!is.null(opt$gmt_df)){
    print("load self gmt")
    gsea_gmt = normalizePath(opt$gmt_df)
  }
  gsea_dir =  paste0(output_dir,"/gsea_result")
  dir.create(gsea_dir , recursive=T )
  
  library(fgsea)
  library(dplyr)
  library(clusterProfiler)
  library(ggsci)
  library(pheatmap)
  library(ggpubr)
  
  diff_all = diff_regulation_df
  geneList = diff_all$FC
  names(geneList) = diff_all$gene
  ## 3.排序很重要
  geneList = sort(geneList, decreasing = TRUE)
  db <- read.gmt(gsea_gmt) 
  
  gsea_result <- GSEA(geneList, TERM2GENE = db, verbose=FALSE,minGSSize = 5 ,maxGSSize = 500, pvalueCutoff = 1,eps=1e-10)
    write.table(gsea_result@result,file = paste0(gsea_dir,'/gsea.result.xls'),quote=F,sep="\t",row.names=F)
  save(gsea_result,file = paste0(gsea_dir,'/gsea.result.rda'))
  if(dim(gsea_result@result)[1] > 40 ){
    up_term <- gsea_result@result %>% arrange(desc(NES)) %>% dplyr::top_n(topn,NES)
    down_term <- gsea_result@result %>% arrange(NES) %>% dplyr::top_n(-topn,NES)
    plot_term = unique(c(up_term$ID,down_term$ID))
  } else {
    plot_term = unique(gsea_result@result$ID)
  }
  plot_df = gsea_result@result[(gsea_result@result$ID %in% plot_term),]
  
  PlotGsea = function(gsea_Wnt_NRVC = gsea_Wnt_NRVC,i = i,show_pval = T){
            p6 <- enrichplot::gseaplot2(gsea_Wnt_NRVC,i,
                  #title = gsea_Wnt_NRVC$Description[i],
                  color = "#FA5860", #设置线条颜色
                  base_size = 12, #设置字体大小
                  rel_heights = c(1.5, 0.3, 0.5), # 设置三组图的比例
                  subplots = 1:3, #设置需要显示的图，如c(1,2)
                  pvalue_table = show_pval, #设置显示p值
                  ES_geom = "line" #还有一种'dot'的点图
                  )
            g = sample(gsea_Wnt_NRVC[[gsea_Wnt_NRVC$ID[i]]],15,replace = T)
            p6[[1]] = p6[[1]]+ enrichplot::geom_gsea_gene(g, geom=ggrepel::geom_text_repel)
            return(p6)
  }

    for(pathway in plot_term){
        p6 =PlotGsea(gsea_result,which(gsea_result$ID==pathway)[1],show_pval = F)
        ggsave(plot = p6,paste0(gsea_dir,"/",pathway,".pdf"),width = 8,height = 6)
        ggsave(plot = p6,paste0(gsea_dir,"/",pathway,".png"),width = 8,height = 6)
    }
  print("gsea analysis done")
}