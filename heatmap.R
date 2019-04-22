#heatmap for differential expression matrix
library(pheatmap)
diff <- 'diffExp.txt'
all_data <- read.table(file=diff,sep="\t",header=T,check.names=F, row.names=1)
hmExp <- log10(all_data+0.001)
hmMat <- as.matrix(hmExp) 
annotation_col = data.frame(Group = factor(rep(c("Normal", "Tumor"), c(6, 120)))) 
rownames(annotation_col) = colnames(all_data) 
pheatmap(hmMat,treeheight_row=100,treeheight_col=20,cluster_cols=FALSE,color=colorRampPalette(c("green","black","red"))(1000),border_color=NA,fontsize=10,fontsize_row=8,fontsize_col=16, annotation_col = annotation_col,show_rownames = F,show_colnames = F,file='heatmap.pdf')
 
