# ==============================================================================
# STAGE 1 | CS-CORE 前置QC + 是否进入阶段2的建议
# ==============================================================================
# 作用：
#   1) 从已注释 Seurat RDS 中提取指定细胞类型；
#   2) 按“每个样本该细胞类型的细胞数”自动过滤样本；
#   3) 输出过滤前/后的样本、目标基因、全基因 Detection QC；
#   4) 判断目标基因是否属于 Detection >=5% 基因集；
#   5) 若目标基因 0<Detection<5%，仍允许进入阶段2，阶段2会强制加入目标基因；
#   6) 保存 Stage1_QC_Result.rds，供阶段2直接读取。
# ==============================================================================

rm(list=ls()); gc(); options(stringsAsFactors=FALSE, warn=1)

# ============================== 只改这里 ======================================
work_dir <- "/home/cy_paperlab/CY/胶质瘤注释/第四批数据集注释尝试（多图形）/目标基因/GSE271959(scRNA)Object1/共表达基因分析"
rds_file <- file.path(work_dir, "Glioma_3B_Fully_Annotated.rds")
target_gene <- "UBE2J2"
cell_type <- "NPC-like"
celltype_col <- "CellType"
sample_col <- "SampleID"
min_cells_per_sample <- 20L   # 该细胞类型中 <20 cells 的样本自动剔除
detection_cutoff_pct <- 5     # 全基因 universe: Detection >=5%
# ============================================================================

safe_tag <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
analysis_dir <- file.path(work_dir, "CSCORE_TwoStage", safe_tag(cell_type), safe_tag(target_gene))
out_dir <- file.path(analysis_dir, "Stage1_QC")
plot_dir <- file.path(out_dir, "Plots")
stage2_dir <- file.path(analysis_dir, "Stage2_Coexpression")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })
save_csv <- function(x, name) write.csv(x, file.path(out_dir,name), row.names=FALSE, fileEncoding="UTF-8")

get_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if (sum(grepl("^counts($|\\.)", Layers(obj[["RNA"]]))) > 1) obj <- JoinLayers(obj, assay="RNA")
  if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA counts layer not found. CS-CORE requires raw UMI counts.")
  x <- LayerData(obj, assay="RNA", layer="counts")
  if (!inherits(x,"sparseMatrix")) x <- Matrix::Matrix(x, sparse=TRUE)
  list(obj=obj, counts=x)
}

sample_target_qc <- function(counts, meta, target, sample_col) {
  smps <- levels(droplevels(factor(meta[[sample_col]])))
  if (!target %in% rownames(counts)) return(data.frame(Sample=smps,Cells=as.integer(table(factor(meta[[sample_col]],levels=smps))),
    Target_Positive_Cells=NA_integer_,Target_Detection_Percent=NA_real_,Target_Mean_UMI=NA_real_,Target_Total_UMI=NA_real_))
  y <- as.numeric(counts[target,,drop=TRUE])
  do.call(rbind,lapply(smps,function(s){i<-meta[[sample_col]]==s; z<-y[i]; data.frame(Sample=s,Cells=sum(i),Target_Positive_Cells=sum(z>0),
    Target_Detection_Percent=mean(z>0)*100,Target_Mean_UMI=mean(z),Target_Total_UMI=sum(z))}))
}

gene_qc <- function(counts, cutoff) {
  depth <- as.numeric(Matrix::colSums(counts)); names(depth) <- colnames(counts)
  if (any(!is.finite(depth)) || any(depth<=0)) stop("Invalid/zero sequencing depth detected.")
  det <- Matrix::rowMeans(counts>0)*100; mean_umi <- Matrix::rowMeans(counts)
  rel <- Matrix::rowMeans(counts %*% Matrix::Diagonal(x=1/pmax(depth,1)))
  ord <- order(rel,decreasing=TRUE,na.last=NA); ranked <- names(rel)[ord]
  rk <- setNames(integer(length(rel)),names(rel)); rk[ranked] <- seq_along(ranked)
  data.frame(Gene=names(rel),Mean_Relative_Expression=as.numeric(rel),Mean_UMI=as.numeric(mean_umi[names(rel)]),
    Detection_Percent=as.numeric(det[names(rel)]),Expression_Rank=as.integer(rk[names(rel)]),
    Pass_Detection=is.finite(det[names(rel)]) & det[names(rel)]>=cutoff,stringsAsFactors=FALSE)
}

# ==============================================================================
# 1. 读取 RDS + 提取目标细胞类型
# ==============================================================================
cat("\n[Stage1-1] Loading RDS and extracting ",cell_type,"...\n",sep="")
if (!file.exists(rds_file)) stop("RDS not found: ",rds_file)
obj <- readRDS(rds_file)
if (!celltype_col %in% colnames(obj@meta.data)) stop("Metadata column not found: ",celltype_col)
if (!sample_col %in% colnames(obj@meta.data)) stop("Metadata column not found: ",sample_col)
if (!"RNA" %in% Assays(obj)) stop("RNA assay not found.")
Idents(obj) <- celltype_col; ct <- subset(obj,idents=cell_type); rm(obj); gc(verbose=FALSE)
if (ncol(ct)==0) stop("No cells found for cell type: ",cell_type)
ret <- get_counts(ct); ct <- ret$obj; counts_before <- ret$counts; rm(ret)
meta_before <- ct@meta.data[colnames(counts_before),,drop=FALSE]
if (anyNA(meta_before[[sample_col]])) stop("NA detected in ",sample_col)
meta_before[[sample_col]] <- droplevels(factor(meta_before[[sample_col]]))

# ==============================================================================
# 2. 样本过滤：只按该细胞类型的细胞数，不看目标基因表达
# ==============================================================================
cat("[Stage1-2] Sample QC: exclude samples with <",min_cells_per_sample," ",cell_type," cells...\n",sep="")
sample_filter <- as.data.frame(table(meta_before[[sample_col]]),stringsAsFactors=FALSE); names(sample_filter)<-c("Sample","Cells")
sample_filter$Cells <- as.integer(sample_filter$Cells)
sample_filter$Decision <- ifelse(sample_filter$Cells<min_cells_per_sample,"Exclude","Keep")
sample_filter$Reason <- ifelse(sample_filter$Decision=="Exclude",paste0("Cells < ",min_cells_per_sample),paste0("Cells >= ",min_cells_per_sample))
save_csv(sample_filter,"01_Sample_Filter_Decision.csv")

keep_samples <- as.character(sample_filter$Sample[sample_filter$Decision=="Keep"])
exclude_samples <- as.character(sample_filter$Sample[sample_filter$Decision=="Exclude"])
keep_cells <- rownames(meta_before)[meta_before[[sample_col]] %in% keep_samples]

qc_target_before <- sample_target_qc(counts_before,meta_before,target_gene,sample_col)
qc_target_before$Stage <- "Before"
qc_target_before$Decision <- sample_filter$Decision[match(qc_target_before$Sample,sample_filter$Sample)]
save_csv(qc_target_before,"02_Target_QC_Before_Filter.csv")

if (length(keep_cells)>0) {
  counts_after <- counts_before[,keep_cells,drop=FALSE]
  meta_after <- meta_before[keep_cells,,drop=FALSE]; meta_after[[sample_col]] <- droplevels(factor(meta_after[[sample_col]]))
  qc_target_after <- sample_target_qc(counts_after,meta_after,target_gene,sample_col); qc_target_after$Stage <- "After"; qc_target_after$Decision <- "Keep"
} else {
  counts_after <- counts_before[,0,drop=FALSE]; meta_after <- meta_before[0,,drop=FALSE]
  qc_target_after <- data.frame(Sample=character(),Cells=integer(),Target_Positive_Cells=integer(),Target_Detection_Percent=numeric(),Target_Mean_UMI=numeric(),Target_Total_UMI=numeric(),Stage=character(),Decision=character())
}
save_csv(qc_target_after,"03_Target_QC_After_Filter.csv")

# ==============================================================================
# 3. 全基因 Detection：过滤前 vs 过滤后
# ==============================================================================
cat("[Stage1-3] Genome-wide Detection QC...\n")
gq_before <- gene_qc(counts_before,detection_cutoff_pct); save_csv(gq_before,"04_Gene_QC_Before_Filter.csv")
gq_after <- if (ncol(counts_after)>0) gene_qc(counts_after,detection_cutoff_pct) else data.frame()
if (nrow(gq_after)>0) save_csv(gq_after,"05_Gene_QC_After_Filter.csv")

target_exists <- target_gene %in% rownames(counts_before)
target_positive_after <- if (target_exists && ncol(counts_after)>0) sum(as.numeric(counts_after[target_gene,,drop=TRUE])>0) else 0L
target_detection_before <- if (target_exists) gq_before$Detection_Percent[match(target_gene,gq_before$Gene)] else NA_real_
target_detection_after <- if (target_exists && nrow(gq_after)>0) gq_after$Detection_Percent[match(target_gene,gq_after$Gene)] else NA_real_
target_rank_after <- if (target_exists && nrow(gq_after)>0) gq_after$Expression_Rank[match(target_gene,gq_after$Gene)] else NA_integer_
genes_detection5 <- if (nrow(gq_after)>0) as.character(gq_after$Gene[gq_after$Pass_Detection]) else character(0)
target_in_detection5 <- target_gene %in% genes_detection5
target_forced <- target_exists && target_positive_after>0 && !target_in_detection5
selected_genes <- unique(c(genes_detection5, if (target_forced) target_gene else character(0)))

# Validation genes = top expressed Detection>=5% genes + target; target may be forced.
validation_n_top <- 100L
if (nrow(gq_after)>0) {
  det5_ranked <- gq_after$Gene[gq_after$Pass_Detection][order(gq_after$Expression_Rank[gq_after$Pass_Detection])]
  validation_genes <- unique(c(head(det5_ranked,validation_n_top), if(target_exists) target_gene else character(0)))
} else validation_genes <- character(0)

# ==============================================================================
# 4. 给出“是否运行阶段2”的意见
# ==============================================================================
retained_n_samples <- length(unique(as.character(meta_after[[sample_col]])))
if (!target_exists) {
  stage2_allowed <- FALSE; advice <- "不建议运行阶段2：目标基因不在 RNA counts 中。"
} else if (retained_n_samples < 2) {
  stage2_allowed <- FALSE; advice <- "不建议运行阶段2：样本过滤后不足2个样本，无法进行当前 SampleID mean-only 校正。"
} else if (target_positive_after==0) {
  stage2_allowed <- FALSE; advice <- "不建议运行阶段2：过滤后目标基因没有任何阳性细胞，无法获得有效共表达。"
} else if (length(selected_genes)<2) {
  stage2_allowed <- FALSE; advice <- "不建议运行阶段2：可用于 CS-CORE 的基因数不足。"
} else if (target_in_detection5) {
  stage2_allowed <- TRUE; advice <- paste0("建议运行阶段2：",target_gene," Detection = ",sprintf("%.2f",target_detection_after),"% >= ",detection_cutoff_pct,"%，正常进入 Detection>=5% 基因集。")
} else {
  stage2_allowed <- TRUE; advice <- paste0("可以运行阶段2，但需谨慎解释：",target_gene," Detection = ",sprintf("%.2f",target_detection_after),"% < ",detection_cutoff_pct,"%。阶段2将只把目标基因强制加入；其他基因仍要求 Detection>=",detection_cutoff_pct,"%。")
}

summary_df <- data.frame(Metric=c("TargetGene","CellType","Cells_Before","Cells_After","Samples_Before","Samples_After","ExcludedSamples",
  "MinCellsPerSample","DetectionCutoffPercent","TargetExists","TargetPositiveCells_After","TargetDetectionPercent_Before","TargetDetectionPercent_After",
  "TargetExpressionRank_After","GenesDetectionGE5_After","TargetInDetectionGE5","TargetWillBeForcedIntoStage2","Stage2Allowed","Advice"),
  Value=c(target_gene,cell_type,ncol(counts_before),ncol(counts_after),nlevels(meta_before[[sample_col]]),retained_n_samples,paste(exclude_samples,collapse=";"),
  min_cells_per_sample,detection_cutoff_pct,target_exists,target_positive_after,target_detection_before,target_detection_after,target_rank_after,length(genes_detection5),
  target_in_detection5,target_forced,stage2_allowed,advice),stringsAsFactors=FALSE)
save_csv(summary_df,"00_Stage1_Summary_and_Advice.csv")

# ==============================================================================
# 5. QC 图：样本过滤 + 目标基因前后 + 全基因 Detection 前后
# ==============================================================================
p1 <- ggplot(sample_filter,aes(x=reorder(Sample,Cells),y=Cells,fill=Decision))+geom_col(color="black",linewidth=.25)+
  geom_hline(yintercept=min_cells_per_sample,linetype="dashed",linewidth=.55)+coord_flip()+
  labs(title=paste0(cell_type," cells per sample"),subtitle=paste0("Samples with <",min_cells_per_sample," cells are excluded"),x=NULL,y="Cell number")+
  scale_fill_manual(values=c(Keep="grey55",Exclude="grey85"))+theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"),plot.subtitle=element_text(hjust=.5))
ggsave(file.path(plot_dir,"01_Sample_Cell_Count_Filter.pdf"),p1,width=7,height=max(4,0.28*nrow(sample_filter)+2.5))

if (target_exists) {
  qcp <- rbind(qc_target_before[,c("Sample","Target_Detection_Percent","Stage")],qc_target_after[,c("Sample","Target_Detection_Percent","Stage")])
  p2 <- ggplot(qcp,aes(Sample,Target_Detection_Percent,fill=Stage))+geom_col(position="dodge",color="black",linewidth=.2)+
    geom_hline(yintercept=detection_cutoff_pct,linetype="dashed",linewidth=.55)+labs(title=paste0(target_gene," detection before vs after sample QC"),x=NULL,y="Detection (%)")+
    theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=45,hjust=1),plot.title=element_text(hjust=.5,face="bold"))
  ggsave(file.path(plot_dir,"02_Target_Detection_Before_After.pdf"),p2,width=8,height=5)
}

if (nrow(gq_after)>0) {
  dd <- rbind(data.frame(Detection_Percent=gq_before$Detection_Percent,Stage="Before"),data.frame(Detection_Percent=gq_after$Detection_Percent,Stage="After"))
  p3 <- ggplot(dd,aes(Detection_Percent,fill=Stage))+geom_histogram(position="identity",alpha=.45,bins=60)+
    geom_vline(xintercept=detection_cutoff_pct,linetype="dashed",linewidth=.55)+labs(title="Genome-wide gene detection before vs after sample QC",x="Detection (%)",y="Gene count")+
    theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"))
  ggsave(file.path(plot_dir,"03_Gene_Detection_Before_After.pdf"),p3,width=7,height=5)

  if (target_exists) {
    rank_df <- gq_after[order(gq_after$Expression_Rank),]
    p4 <- ggplot(rank_df,aes(Expression_Rank,Mean_Relative_Expression))+geom_point(size=.6,alpha=.35)+
      geom_point(data=rank_df[rank_df$Gene==target_gene,,drop=FALSE],size=2.6)+
      labs(title=paste0(target_gene," expression rank after sample QC"),subtitle=paste0("Detection = ",sprintf("%.2f",target_detection_after),"%"),x="Expression rank",y="Mean relative expression")+
      theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5,face="bold"),plot.subtitle=element_text(hjust=.5))
    ggsave(file.path(plot_dir,"04_Target_Expression_Rank_After.pdf"),p4,width=7,height=5)
  }
}

# ==============================================================================
# 6. 保存阶段2唯一需要读取的交接文件
# ==============================================================================
manifest <- list(version="TwoStage_Simple_v1",created=as.character(Sys.time()),stage2_allowed=stage2_allowed,advice=advice,
  rds_file=rds_file,work_dir=work_dir,analysis_dir=analysis_dir,stage1_dir=out_dir,stage2_dir=stage2_dir,
  target_gene=target_gene,cell_type=cell_type,celltype_col=celltype_col,sample_col=sample_col,min_cells_per_sample=min_cells_per_sample,
  detection_cutoff_pct=detection_cutoff_pct,included_cells=colnames(counts_after),retained_samples=keep_samples,excluded_samples=exclude_samples,
  target_exists=target_exists,target_positive_cells=target_positive_after,target_detection_pct=target_detection_after,target_in_detection5=target_in_detection5,
  target_forced=target_forced,genes_detection5=genes_detection5,selected_genes=selected_genes,validation_genes=validation_genes,gene_qc=gq_after,
  sample_adjustment_mode="mean_only",validation_n_top=validation_n_top,validation_est_tol=1e-6,validation_test_tol=1e-4,validation_p_tol=1e-6,
  irls_n_iter=10L,irls_eps=.05,irls_verbose=FALSE,auto_install_cscore=TRUE,fdr_cutoff=.05,strong_positive_estimate=.10,target_omp_threads=12L,blas_threads=12L)
saveRDS(manifest,file.path(out_dir,"Stage1_QC_Result.rds"))

cat("\n======================================================================\n")
cat("STAGE 1 COMPLETED\n")
cat("======================================================================\n")
cat("Target gene                 : ",target_gene,"\n",sep="")
cat("Cell type                   : ",cell_type,"\n",sep="")
cat("Samples before / after      : ",nlevels(meta_before[[sample_col]])," / ",retained_n_samples,"\n",sep="")
cat("Excluded samples            : ",ifelse(length(exclude_samples),paste(exclude_samples,collapse=", "),"None"),"\n",sep="")
cat("Target detection after QC   : ",ifelse(is.finite(target_detection_after),sprintf("%.2f%%",target_detection_after),"NA"),"\n",sep="")
cat("Detection>=5% genes         : ",length(genes_detection5),"\n",sep="")
cat("Target in Detection>=5%     : ",target_in_detection5,"\n",sep="")
cat("Target forced in Stage 2    : ",target_forced,"\n",sep="")
cat("\n>>> 建议：",advice,"\n",sep="")
cat("\nStage 1 result file          : ",file.path(out_dir,"Stage1_QC_Result.rds"),"\n",sep="")
cat("======================================================================\n")
