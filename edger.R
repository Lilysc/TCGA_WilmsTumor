#################参数设置
input <- 'mRNA_exp.txt'
group <- '11-01'
foldChange <- 2
padj <- 0.01

library(edgeR)
change_rowname <- function(DATA,name='id'){
  DATA=cbind(id=row.names(DATA), DATA)
  row.names(DATA)=NULL
  return(DATA)
}

data <- read.table(file=input,sep="\t",header=T,check.names=F, row.names=1)
##################extract sample(code:11,01)
normal_symbol <- unlist(strsplit(group,"\\-"))[1]
tumor_symbol <- unlist(strsplit(group,"\\-"))[2]
target_type <- c()
normal_num <- 0
tumor_num <- 0
for(i in colnames(data)){
  sample_type <- unlist(strsplit(i,"\\-"))[4]
  type <- paste0(unlist(strsplit(sample_type,split = "")[1])[1], unlist(strsplit(sample_type,split = "")[1])[2])
  if(type == normal_symbol){
    target_type <- c(target_type,i)
    normal_num  <- normal_num + 1
  }else if(type == tumor_symbol){
    target_type <- c(target_type,i)
    tumor_num <- tumor_num + 1
  }
}
rt <- data[,target_type]

numc <- sapply(rt, is.factor)
rt[numc] <- lapply(rt[numc], function(x) as.numeric(as.character(x)))
rt <- na.omit(rt)
exp <- rt[,1:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
dim(data)
data <- data[rowMeans(data)>0,]
dim(data)

group <- c(rep("normal",normal_num),rep("tumor",tumor_num))
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair=c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et,n=100000)

allDiff <- ordered_tags$table
allDiff <- allDiff[is.na(allDiff$FDR)==FALSE,]
diff <- allDiff
newData <- y$pseudo.counts

diff_save <- change_rowname(diff)
write.table(diff_save,file="edgeOut.xls",sep="\t",quote=F,row.names = FALSE)

diffSig <- diff[(diff$FDR < padj & (diff$logFC > foldChange | diff$logFC < (-foldChange))),]
diffSig_save <-  change_rowname(diffSig)
write.table(diffSig_save,file="diffSig.xls",sep="\t",quote=F,row.names = FALSE)
diffUp <- diff[(diff$FDR < padj & (diff$logFC > foldChange)),]
diffUp_save <-  change_rowname(diffUp)
write.table(diffUp_save,file="up.xls",sep="\t",quote=F,row.names = FALSE)
diffDown <- diff[(diff$FDR < padj & (diff$logFC < (-foldChange))),]
diffDown_save <-  change_rowname(diffDown)
write.table(diffDown_save,file="down.xls",sep="\t",quote=F,row.names = FALSE)

normalizeExp <- rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)
diffExp <-rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

heatmapData <- newData[rownames(diffSig),]
#Volcano
pdf(file="vol.pdf")
xMax <- max(-log10(allDiff$FDR))+1
yMax <- 12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub <- allDiff[(allDiff$FDR < padj & allDiff$logFC > foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub <- allDiff[(allDiff$FDR < padj & allDiff$logFC < (-foldChange)),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

tiff(file="vol.tiff",width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=1200)
xMax <- max(-log10(allDiff$FDR))+1
yMax <- 12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub <- allDiff[(allDiff$FDR < padj & allDiff$logFC > foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub <- allDiff[(allDiff$FDR < padj & allDiff$logFC < (-foldChange)),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
