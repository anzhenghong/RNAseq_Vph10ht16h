##DEG analysis
##BiocManager::install()  ###software installation
library(DESeq2)
library(dplyr)
library(tibble)
library(openxlsx)
####using the openxlsx package in R to convert the count_matrix output by prepDE.py into Vph10T16_count.xlsx
setwd("D:/third/Vph10t16/stringtie/Recent_results/DEG")
a <- read.xlsx("../expression/rawcount/Vph10T16_count.xlsx",rowNames = T)
###Vph10T160d_vs_Vph10T161d_DEG;For the other groups, replace "1d" with the corresponding time point.
b=select(a,'Vph10T160dF11','Vph10T160dF12','Vph10T160dF13','Vph10T161dF11','Vph10T161dF12','Vph10T161dF13')
count <- b[rowMeans(b)>1,]
colData <- data.frame(row.names = colnames(count),
                      condition = factor(c(rep('Vph10T160d',3),rep('Vph10T161d',3)),
                                         levels=c('Vph10T160d','Vph10T161d'))
)
##
dds <- DESeqDataSetFromMatrix(countData = count, colData = colData,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition","Vph10T161d","Vph10T160d"))
out=select(as.data.frame(res),'baseMean','log2FoldChange','pvalue','padj')
out=rownames_to_column(out,'gene_id')
out=filter(out,!is.na(padj))
##
FCcut=2
FDRcut=0.05
diff=filter(out,padj<FDRcut,abs(log2FoldChange)>=abs(log2(FCcut)))
diff$Group <- "NA"
diff$Group[which(diff$log2FoldChange <= -1)] = "down-regulated"
diff$Group[which(diff$log2FoldChange >= 1)] = "up-regulated"
write.xlsx(diff,"Vph10T160d_vs_Vph10T161d_DEG.xlsx")
