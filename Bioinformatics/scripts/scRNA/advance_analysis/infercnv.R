#加载需要的包和数据
# /home/shazhi11/miniforge3/envs/atac 
library(Seurat)
library(infercnv)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
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
sc_ob = readRDS("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/sc_GSE182109/NPC_like_CENPH_MCM3_4groupdayu0.rds")

sub_ob = sc_ob[,sc_ob@meta.data$CellType %in% c("NPC-like","Normal Glial")]  ## 选择细胞类型
sub_ob$celltype <- gsub(" ","-",sub_ob$CellType)
sub_ob$celltype <- paste(sub_ob$CENPH_MCM3_group)
Idents(sub_ob) <- sub_ob$celltype

ref_group = "CENPH-_MCM3-"

if(dim(sub_ob)[2] > 30000){
	print("cells scRNA more than 3W,downsample ...")
	set.seed(12306)  # 设置随机种子以确保结果可重复
	result <- as.data.frame(sub_ob@meta.data) %>% tibble::rownames_to_column(.,var = "barcode") %>%
	group_by(celltype) %>% slice_sample(n = 4000)
	sub_ob = sub_ob[,result$barcode]
}


#由于CNV分析要用到infercnv包，其CreateInfercnvObject函数需要三个文件创建对象：The raw_counts_matrix；the gene_order_file, contains chromosome, start, and stop position for each gene；The annotations_file, containing the cell name and the cell type classification

#获取counts矩阵
dat = GetAssayData(sub_ob,layers='counts',assay='RNA')
dat[1:4,1:4]

#获取细胞亚群信息
groupinfo=data.frame(v1=colnames(dat),v2= Idents(sub_ob ) )
head(groupinfo)
geneInfor = read.delim("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/scRNA_GSE278450/hsa.geneFile.txt",sep="\t")
head(geneInfor)

#安装染色体和基因位置排序，并将SYMBOL，chr，start，end这几列留下
geneInfor=geneInfor[with(geneInfor,order(chr, start)),]

#去重（这里也可以操作去除性染色体，也可以把染色体排序方式改变）
geneInfor=geneInfor[!duplicated(geneInfor[,1]),] 
length(unique(geneInfor[,1]))
head(geneInfor)

#整理数据
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[intersect( geneInfor[,1], rownames(dat) ),] 
dim(dat)
head(groupinfo)
dat[1:4,1:4]
table(groupinfo$v2)
head(groupinfo)
print("check file done")
# # 为了节约计算机资源，直接抽样 
# kp=sample(1:nrow(groupinfo),2000)
# groupinfo=groupinfo[kp,]
# dat=dat[,kp]
# 如果是真实项目，且计算资源足够,忽略这个抽样的操作

#保存文件
# expFile='expFile.txt'
# write.table(dat,file =expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

#ref_group_names 这里的细胞是正常对照，然后跟其他的细胞比较
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = dat,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = ref_group #这里用 Normal-Glial 做参照
)

infercnv_obj <- infercnv::run(infercnv_obj, output_format = "pdf",
                               cutoff =  0.1, #smart-seq选择1,10X选择0.1
                               out_dir = "infercnv_output", # dir is auto
                               # 是否根据细胞注释文件的分组
                               # 对肿瘤细胞进行分组
                               # 影响read.dendrogram, 如果有多个细胞类型，且设置为TRUE，
                               # 后续的read.dendrogram无法执行
                               cluster_by_groups =  TRUE, #是否根据细胞类型(由细胞注释文件中定义)
                               hclust_method = "ward.D2",# ward.D2 方法进行层次聚类
                               analysis_mode = "subclusters", # 默认是samples，推荐是subclusters
                               denoise = TRUE, # 去噪音
                               HMM = F,  ##特别耗时间,是否要去背景噪音
                               plot_steps = F, #不在每个步骤后生成图形。
                               leiden_resolution = "auto", #可以手动调参数
                               num_threads = 16 # 16线程工作，加快速度，一般电脑不建议哈
                               )

# infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)
saveRDS(infercnv_obj,paste0('infercnv_output/',"inferCNV.rds"))
infercnv::plot_cnv(infercnv_obj, 
                   out_dir = paste0('infercnv_output/'),
                   plot_chr_scale = T, 
                   output_filename = "infercnv_heatmap",output_format = "pdf", 
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))
				   )

infercnv::plot_cnv(infercnv_obj, 
                   out_dir = paste0('infercnv_output/'),
                   plot_chr_scale = T, 
                   output_filename = "infercnv_heatmap",output_format = "png", 
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))
				   )



#  plot 1 获取 inferCNV 校正后的 CNV 信号矩阵 连续变量法（continuous CNV abnormality score）
expr <- infercnv_obj@expr.data
expr2 = expr - 1  # 中心化：对每个值执行 (x - 1)
expr2 = expr2 ^ 2  # 平方：让偏离越大贡献越大，同时避免正负抵消
CNV_score = as.data.frame(colMeans(expr2))  # 对每个细胞取所有基因的均值

colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
CNV_score <- merge(x=CNV_score, y=groupinfo, by.x ="CB", by.y="v1")
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:8] 

# 绘图
p = ggplot(CNV_score, aes(x = v2, y = CNV_score, fill = v2)) +
    geom_violin(trim = FALSE, scale = "width", color = "black") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black") +
    scale_fill_manual(values = color_v) +
    ylim(0, 0.03) +  ## 限制了y的范围，按需调整哦！！
    ylab("CNV abnormality score") +
    xlab("Cell type") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),       # 去掉网格
      panel.border = element_rect(color = "black", fill = NA, size = 1), # 加方框
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"            # 去掉图例
)
ggsave(plot = p,"CNV_score_meadian.vln.pdf",w = 10,h = 6)
ggsave(plot = p,"CNV_score_meadian.vln.png",w = 10,h = 6)
write.csv(CNV_score,"CNV_score.meadian.csv",quote=F,row.names = F)
# 绘图2 =======阈值分级法

if(T){
  tmp1 = expr[,unlist(infercnv_obj@reference_grouped_cell_indices[ref_group])]
  tmp = tmp1
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infercnv_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
   
  cnv_score_table_pts[1:4,1:4]
  str(as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 

  #可视化
  head(cell_scores_CNV) 
  meta = sub_ob@meta.data
  meta$totalCNV = cell_scores_CNV[match(colnames(sub_ob),
                              rownames(cell_scores_CNV)),1] 
  p <- ggplot(meta, aes(x= celltype, 
                   y=totalCNV, 
                   fill=celltype)) +
              geom_boxplot() +
  			scale_fill_manual(values= colors)
  ggsave("CNV_score_Threshold.vln.pdf",plot = p, width = 9,height = 7,dpi = 300)
  ggsave("CNV_score_Threshold.vln.png",plot = p, width = 9,height = 7,dpi = 300)
  write.csv(meta,"CNV_score.Threshold.csv",quote=F,row.names = T)
}