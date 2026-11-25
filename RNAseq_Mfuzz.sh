##devtools::install_github("junjunlab/ClusterGVis")
##BiocManager::install()  ###software installation
library(ClusterGVis)
library(openxlsx)
library(eoffice)
setwd("D:/third/Vph10t16/stringtie/Recent_results/4.Cluster")
data <- read.xlsx("Vph10T16_All_DEG_Mean_TPM.xlsx",rowNames = T)
df1 <-as.data.frame(t(apply(data,1,scale)))
colnames(df1) <- colnames(data)
ck <- clusterData(exp = df1,
                  cluster.method = "mfuzz",
                  cluster.num = 5
)
str(ck)
visCluster(
    object = ck,
    plot.type = "line",
    show_row_dend = F
)
results <- as.data.frame(ck$wide.res)
write.xlsx(results,"Cluster_results2.xlsx")
topptx(filename = "Cluster_expression.pptx")
