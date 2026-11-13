library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)


data_ob = readRDS("annotated.rds")
dir.create("1.umap_tsne_summary")
setwd("1.umap_tsne_summary")
data_ob$CellType = factor(data_ob$CellType,levels = sort(levels(data_ob$CellType)))

colors = c("#5BCCB6","#4782B3","#EB7BC0","#FFAA33","#AB80F5","#EACC93","#5C9966","#A0C896","#EB8D8E","#CEAFAF","#FF9B89","#89DACB","#6FA1CA","#EC9DCE","#FFC26D","#FEE376","#C3A4F8","#F0DAB2","#77C785","#BBD8B4","#F1ADAE","#DCC6C6","#FF735A","#FED840","#0272B1","#F6632B","#9B5D00")

###1. umap tsne  分组/不分组 ，样本 summary 占比
p1 = Seurat::DimPlot(object = data_ob,reduction = "umap", label = F, pt.size = 0.9 , group.by = "CellType") + 
      ggplot2::scale_colour_manual( values = colors)
ggsave("umap.pdf", plot = p1,w = 7.5,h = 7)
ggsave("umap.png", plot = p1,w = 7.5,h = 7,dpi=800)

p2 = Seurat::DimPlot(object = data_ob,reduction = "tsne", label = F, pt.size = 0.9, group.by = "CellType") + 
      ggplot2::scale_colour_manual( values = colors)
ggsave("tsne.pdf", plot = p2,w = 7.5,h = 7)
ggsave("tsne.png", plot = p2,w = 7.5,h = 7,dpi=800)


p3 = Seurat::DimPlot(object = data_ob,reduction = "umap", label = F, pt.size = 0.9, group.by = "CellType", split.by= "group") & 
      ggplot2::scale_colour_manual( values = colors)
ggsave("umap.group.pdf", plot = p3,w = 14.5,h = 7)
ggsave("umap.group.png", plot = p3,w = 14.5,h = 7,dpi=800)

p4 = Seurat::DimPlot(object = data_ob,reduction = "tsne", label = F, pt.size = 0.9, group.by = "CellType", split.by= "group") & 
      ggplot2::scale_colour_manual( values = colors)
ggsave("tsne.group.pdf", plot = p4,w = 14.5,h = 7)
ggsave("tsne.group.png", plot = p4,w = 14.5,h = 7,dpi=800)

####--------------------------------- summary --------------------------

DATA <- data_ob@meta.data[,c("SampleID", "CellType")] %>%
                dplyr::group_by( SampleID, CellType) %>%
                dplyr::summarize(cell_number = dplyr::n()) %>%
                dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100)
write.table(as.data.frame(DATA), "sample2CellType.freq.xls",sep="\t", col.names=T, row.names =F)

p1 = ggplot(DATA,aes(SampleID,cell_number,fill=CellType))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) + 
  scale_fill_manual( values = colors)
ggsave("sample2CellType.barplot.pdf",plot = p1, width = 8,height = 10)
ggsave("sample2CellType.barplot.png",plot = p1, width = 8,height = 10,bg = "white")



# DATA <- data_ob@meta.data[,c("group", "CellType")] %>%
                # dplyr::group_by( group, CellType)%>%
                # dplyr::summarize(cell_number = dplyr::n()) %>%
                # dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100)
# write.table(as.data.frame(DATA), "group2CellType.freq.xls",sep="\t", col.names=T, row.names =F)

# p1 = ggplot(DATA,aes(group,cell_number,fill=CellType))+
  # geom_bar(stat="identity",position="fill")+
  # ggtitle("")+
  # theme_bw()+
  # theme(axis.ticks.length=unit(0.5,'cm'))+
  # guides(fill=guide_legend(title=NULL)) + 
  # scale_fill_manual( values = colors)
# ggsave("group2CellType.barplot.pdf",plot = p1, width = 8,height = 10)
# ggsave("group2CellType.barplot.png",plot = p1, width = 8,height = 10,bg = "white")


#------------------------ pie plot
pie_df = as.data.frame(table(data_ob$CellType))
names(pie_df) <- c("CellType","cell_num")
pie_df$labes = paste0(round(pie_df$cell_num / sum(pie_df$cell_num)*100,2) ," %")

p =  ggplot(aes(x = " ", y = cell_num, fill = CellType), data = pie_df) + 
    geom_bar(stat = "identity", width = 1) + 
    coord_polar(theta = "y") + 
    scale_fill_manual( values = colors)
library(scales)
p + theme(axis.text.x=element_blank()) +
  geom_text(aes(y = cell_num/3 + c(0, cumsum(cell_num)[-length(cell_num)]), 
            label = labes), size=5)


p = ggplot(pie_df, aes(x = 0.6, y = cell_num , fill = CellType)) +
    geom_bar(stat = "identity" , position = "stack",color = "white" , width = 1) +
    labs(x = "", y = "", title = "")+
    coord_polar(theta = "y",start = 1, direction = 1, clip = "on")+
    theme_bw()+ggplot2::scale_fill_manual( values = colors) +
    geom_text(stat="identity",color = "black", size=3, check_overlap = TRUE, position=position_stack(vjust = 0.6),
        aes(y = cell_num, x =1.3, label = labes))+
    theme_void()+
    theme(legend.title = element_blank())

ggsave("pie_plot.pdf",plot = p, width = 10,height = 8)
ggsave("pie_plot.png",plot = p, width = 10,height = 8,bg = "white")

##--------------------------------patch plot----------
library(patchwork)

DATA2 <- data_ob@meta.data[,c("SampleID", "CellType")] %>%
                dplyr::group_by( .dots= c("CellType","SampleID")) %>%
                dplyr::summarize(cell_number = dplyr::n()) %>%
                dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100)

new_df = pie_df[order(pie_df$cell_num,decreasing=T),]
new_df$CellType <- factor(new_df$CellType , levels = as.vector(new_df[,'CellType']))

CellType_colors = colors[1:length(levels(data_ob$CellType))]
names(CellType_colors) <- levels(data_ob$CellType)

p11 = ggplot(new_df,aes(x=CellType,y = cell_num,fill = CellType)) + geom_col() + coord_flip() + scale_y_log10()  + theme_bw() +  
    scale_fill_manual( values = CellType_colors) + 
    theme(legend.position = "none") + 
    geom_text(aes(label = cell_num, y = max(new_df$cell_num)) , hjust = -0.4) + 
    theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_blank() ) +
      ylim(0,max(new_df$cell_num) + 1000)

    #+ 
     #geom_text(aes(x = CellType, y = cell_num, label = cell_num))

DATA2$CellType  <- factor(DATA2$CellType , levels = levels(new_df$CellType))

p12 = ggplot(DATA2,aes(CellType,cell_number,fill=group))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+ labs(x = "" , y = "Percent") + 
  scale_fill_manual(values = c("#0272B1","#CA2128")) +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  coord_flip() + 
  theme(axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y =element_blank(),
      axis.ticks.x =element_blank(),
    )

p = wrap_plots(ncol = 2,p11,p12) +  plot_layout(widths = c(3, 2))

ggsave("bar_plot.pdf",plot = p, width = 10,height = 8)
ggsave("bar_plot.png",plot = p, width = 10,height = 8,bg = "white")


#### 拼接 umap 图
colors = c("#5BCCB6","#4782B3","#EB7BC0","#FFAA33","#AB80F5","#EACC93","#5C9966","#A0C896","#EB8D8E","#CEAFAF","#FF9B89","#89DACB","#6FA1CA","#EC9DCE","#FFC26D","#FEE376","#C3A4F8","#F0DAB2","#77C785","#BBD8B4","#F1ADAE","#DCC6C6","#FF735A","#FED840","#0272B1","#F6632B","#9B5D00")

library(gridExtra)
library(grid)
library(cowplot)

p1 =  DimPlot(data_ob,group.by="seurat_clusters",label=T) + scale_color_manual(values= colors) + 
      labs(title = "seurat_clusters") + theme(legend.position = "none")

df2 = read.delim("../anno.xls",sep="\t")[c(2,1)]
df2$new_CellType = paste0(df2$ABV,"\n(",df2$cluster,")") 
# df2$new_CellType[7] <- "Excitatory\n(0,6,7,11,13,\n15,19,34)"

new_df = as.data.frame(data_ob@meta.data)
new_df = tibble::rownames_to_column(new_df,var = "barcodes")
merge_df = merge(new_df,df2[c("ABV","new_CellType")],all.x=T,by.x="CellType",by.y = "ABV")
merge_df = tibble::column_to_rownames(merge_df ,var = "barcodes")

new_ob = AddMetaData(data_ob,merge_df)


p2 =  DimPlot(new_ob,group.by="new_CellType",label=T) + scale_color_manual(values= colors) + 
      labs(title = "CellType") + theme(legend.position = "none")

p3 = plot_grid(p1, p2, align = "h", rel_widths = c(1, 1))
ggsave(plot = p3,"merge_umap.pdf",w = 14,h = 7)
ggsave(plot = p3,"merge_umap.png",w = 14,h = 7)



# p3 = Seurat::DimPlot(object = data_ob,reduction = "umap", label = T, pt.size = 0.9, group.by = "CellType", split.by= "group") + 
      # theme(legend.position = "none") & 
      # ggplot2::scale_colour_manual( values = colors) 
# ggsave("umap.group_label.pdf", plot = p3,w = 14.5,h = 7)
# ggsave("umap.group_label.png", plot = p3,w = 14.5,h = 7,dpi=800)

# p3 = Seurat::DimPlot(object = data_ob,reduction = "tsne", label = T, pt.size = 0.9, group.by = "CellType", split.by= "group") + 
      # theme(legend.position = "none") & 
      # ggplot2::scale_colour_manual( values = colors) 
# ggsave("tsne.group_label.pdf", plot = p3,w = 14.5,h = 7)
# ggsave("tsne.group_label.png", plot = p3,w = 14.5,h = 7,dpi=800)


##marker  和 top 15 基因可视化

allmarker_gene = read.delim("allmarkers.xls",sep="\t")
library(dplyr)
# topn_markers  = allmarker_gene %>%
        # dplyr::group_by(cluster) %>%
        # dplyr::filter(p_val < 0.05)  %>%
        # dplyr::arrange(p_val,desc(avg_log2FC)) %>%
        # dplyr::top_n(15,avg_log2FC)

data_ob = ScaleData(data_ob,features = allmarker_gene$gene)

p = DoHeatmap(data_ob, label = F ,features = allmarker_gene$gene, group.by = "CellType" , assay = "RNA" ,  
    group.colors = colors ) + 
    scale_fill_gradientn(colors = c("Steelblue","white","firebrick"))#设置热图颜色
ggsave(paste0("gene_DoHeatmap.pdf"),plot = p, width = 8 , height = 8 ,limitsize=F)
ggsave(paste0("gene_DoHeatmap.png"),plot = p, width = 8 , height = 8 ,limitsize=F,bg = "white")


# df3 = read.delim("../anno.xls",sep="\t")
# marker_gene = unlist( strsplit(df3$marker,",",perl = T) )


# n_gene = length(unique(marker_gene))
# n_sample = length(unique(data_ob@meta.data$CellType))

# p4 = VlnPlot(data_ob,features= marker_gene,group.by="CellType",stack=T,pt.size = 0.1,flip= T)
# new_df = p4$data
# p5 = ggplot(new_df,aes(x=ident,y = expression,fill = ident)) + 
    # geom_violin(,scale = "width") +
    # geom_jitter(,alpha = 0.5,size = 0.1) + 
    # coord_flip() + 
    # facet_wrap(~feature,nrow = 1,strip.position="top") +
    # scale_fill_manual(values  = colors)  + 
    # theme(strip.background = element_rect(fill = "white"),
      # panel.grid.major =element_blank(), 
      # panel.grid.minor = element_blank(),
      # panel.background = element_blank(),
      # axis.line = element_line(colour = "black") )

# ggsave(paste0("marker_violin.pdf"),plot = p5, width = 4 + 1.2 * n_gene , height = n_sample + 3 ,limitsize=F)
# ggsave(paste0("marker_violin.png"),plot = p5, width = 4 + 1.2 * n_gene , height = n_sample + 3 ,limitsize=F)