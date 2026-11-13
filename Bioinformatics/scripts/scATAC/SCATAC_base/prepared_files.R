#/usr/bin/env Rscript

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
plan("multicore", workers = 4)

## 预处理数据，通过fragments.tsv.gz生成分析流程必须的预置文件和数据
### 1.check tbi index
if ("fragments.tsv.gz.tbi" %in% list.files("./")) {
    print("Found tbi File , Continue...")
} else {
    print("Not Found tbi File , Start Building...")
    system("tabix -p bed fragments.tsv.gz")
    print("Builded tbi File done, Continue...")
}

### 2.check peaks.bed
### 文件不完整，从signac 官网查询到 callpeak 方案 https://stuartlab.org/signac/articles/pbmc_vignette?search-input=CallPeaks
barcodes <- readLines("barcodes.tsv.gz")
fragments <- "fragments.tsv.gz"
total_counts <- CountFragments(fragments)
#cutoff <- 1000 # Change this number depending on your dataset!
#barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
frags <- CreateFragmentObject(path = fragments, cells = barcodes)
peaks <- CallPeaks(frags, macs2.path = "/home/shazhi11/miniforge3/bin/macs3")
mat_peaks <- cbind(
    chr    = as.character(seqnames(peaks)),
    start  = start(peaks),
    end    = end(peaks)
    # width  = width(peaks),
    # strand = as.character(strand(peaks)),
    # mcols(peaks)
)
head_info = read.delim("/home/shazhi11/Lin/Kidney_cancer/linshi/Signac_Test/GSE286415/peaks_header.txt", sep = "\t", header = F)
mat_peaks = as.data.frame(mat_peaks)
names(mat_peaks) = c("V1","V2","V3")
new_peaks = rbind(head_info, (mat_peaks))
write.table(new_peaks,"peaks.bed",sep="\t",row.names = F,col.names = F,quote=F)

# Quantify fragments in each peak
# counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)、
print(list.files("./"))