library(survival)
library('ggplot2')
library('reshape')

picDir="tumor_stage"
dir.create(picDir)

rt=read.table("tumor_surivival.txt",header=T,sep="\t",check.names=F,row.names=1)
df=na.omit(rt)
df[,'fustat'] <- ifelse(df[, 'Vital Status']=='DEAD', 1 ,0)
df <- rename(df,c("Overall Survival Time in Days"="futime"))
df[,"futime"]=df[,"futime"]/365
rt = df

#################
outTab = data.frame()
gene = read.table("import_gene.txt",header=F)
genename = as.vector(gene[,1])[-1]
rt[,'stage'] <- ifelse(rt[, 'Stage'] == ' I' |rt[, 'Stage'] == 'II' ,'I/II' ,'III/IV')
################肿瘤分期，卡方检验
out <- data.frame()
for(i in 1:length(genename)){
  gene1 <- genename[i]
  a <- median(rt[,gene1])
  gene2 <- unlist(strsplit(gene1,"\\|"))[1]
  #将Stage,I和II合并,III和IV合并
  rt[,gene2] <- ifelse(rt[, gene1]>a, 'High' ,'Low')
  chisq_result <- chisq.test(table(rt[,gene2] ,rt[,'stage']))
  out <- rbind(out,data.frame('基因'=gene1,'统计量'=chisq_result$statistic,'p'=chisq_result$p.value))
  #画箱线图
  tiffFile=paste(gene2,"_box_stage.pdf",sep="")
  outTiff=paste(picDir,tiffFile,sep="\\")
  p = ggplot(rt, aes(x=rt[,'Stage'], y=rt[,gene1]),color=rt[,'Stage']) + geom_boxplot(aes(fill=factor(rt[,'Stage']))) + theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) + theme(legend.position="none") + ylab(paste0(gene2,'Expression Level')) + xlab("Tumor Stage")
  ggsave(p, file=outTiff, width=15, height=8)
}
rownames(out) <- NULL
write.csv(out,"stage_chisquare.csv",row.names =F)

################Histology Classification of Primary Tumor
out_1 <- data.frame()
for(i in 1:length(genename)){
  gene1 <- genename[i]
  a <- median(rt[,gene1])
  gene2 <- unlist(strsplit(gene1,"\\|"))[1]
  rt[,gene2] <- ifelse(rt[, gene1]>a, 'High' ,'Low')
  chisq_result <- chisq.test(table(rt[,gene2] ,rt[,'Histology Classification of Primary Tumor']))
  out_1 <- rbind(out_1,data.frame('基因'=gene1,'统计量'=chisq_result$statistic,'p'=chisq_result$p.value))
  #画箱线图
  tiffFile=paste(gene2,"_box_histology.pdf",sep="")
  outTiff=paste(picDir,tiffFile,sep="\\")
  p = ggplot(rt, aes(x=rt[,'Histology Classification of Primary Tumor'], y=rt[,gene1]),color=rt[,'Histology Classification of Primary Tumor']) + geom_boxplot(aes(fill=factor(rt[,'Histology Classification of Primary Tumor']))) + theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) + theme(legend.position="none") + ylab(paste0(gene2,'Expression Level')) + xlab("Histology Classification of Primary Tumor")
  ggsave(p, file=outTiff, width=15, height=8)
}
rownames(out_1) <- NULL
write.csv(out_1,"histology_chisquare.csv",row.names =F)

################Gender Classification of Primary Tumor
out_1 <- data.frame()
for(i in 1:length(genename)){
  gene1 <- genename[i]
  a <- median(rt[,gene1])
  gene2 <- unlist(strsplit(gene1,"\\|"))[1]
  rt[,gene2] <- ifelse(rt[, gene1]>a, 'High' ,'Low')
  chisq_result <- chisq.test(table(rt[,gene2] ,rt[,'Gender']))
  out_1 <- rbind(out_1,data.frame('基因'=gene1,'统计量'=chisq_result$statistic,'p'=chisq_result$p.value))
  #画箱线图
  tiffFile=paste(gene2,"_box_Gender.pdf",sep="")
  outTiff=paste(picDir,tiffFile,sep="\\")
  p = ggplot(rt, aes(x=rt[,'Gender'], y=rt[,gene1]),color=rt[,'Gender']) + geom_boxplot(aes(fill=factor(rt[,'Gender']))) + theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) + theme(legend.position="none") + ylab(paste0(gene2,'Expression Level')) + xlab("Gender")
  ggsave(p, file=outTiff, width=15, height=8)
}
rownames(out_1) <- NULL
write.csv(out_1,"Gender_chisquare.csv",row.names =F)

######################堆叠柱状图
target_gene <- 'ZBTB4|ENSG00000174282|protein_coding'
gene_symbol <- unlist(strsplit(target_gene,"\\|"))[1]
#stage
cross <- table(rt[,gene_symbol] ,rt[,'stage'])
prop_cross <- data.frame(prop.table(cross,2))
prop_cross$Freq <- prop_cross$Freq*100
colnames(prop_cross)[1] <- 'Expression level'
tiffFile=paste(gene_symbol,"_stack_stage.pdf",sep="")
outpdf=paste(picDir,tiffFile,sep="\\")
p1=ggplot(prop_cross,aes(x=prop_cross[,'Var2'], y=prop_cross[,'Freq'],fill=prop_cross[,'Expression level']))+geom_bar(stat='identity',position="stack")+xlab(expression(bold(paste("Clinical Staging(", italic("P"),'<0.05)')))) + ylab("Percentage")+theme(plot.title = element_text(size=20,face = "bold"),axis.text.y=element_text(size=15),axis.text.x=element_text(size=10),axis.title.y=element_text(size=15,face="bold"),axis.title.x=element_text(size=15,face="bold"))+ labs(fill=expression(paste("ZBTB4"," Expression Level")))+theme(legend.position="top")
ggsave(p1, file=outpdf, width=8, height=8)

tiffFile=paste(gene_symbol,"_stack_stage.png",sep="")
outpdf=paste(picDir,tiffFile,sep="\\")
tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=600)
p1=ggplot(prop_cross,aes(x=prop_cross[,'Var2'], y=prop_cross[,'Freq'],fill=prop_cross[,'Expression level']))+geom_bar(stat='identity',position="stack")+xlab(expression(bold(paste("Clinical Staging(", italic("P"),'<0.05)')))) + ylab("Percentage")+theme(plot.title = element_text(size=20,face = "bold"),axis.text.y=element_text(size=15),axis.text.x=element_text(size=10),axis.title.y=element_text(size=15,face="bold"),axis.title.x=element_text(size=15,face="bold"))+ labs(fill=expression(paste("ZBTB4"," Expression Level")))+theme(legend.position="top")
ggsave(p1, file=outpdf, width=8, height=8)

#HF
cross <- table(rt[,gene_symbol] ,rt[,'Histology Classification of Primary Tumor'])
prop_cross <- data.frame(prop.table(cross,2))
prop_cross$Freq <- prop_cross$Freq*100
colnames(prop_cross)[1] <- 'Expression level'
tiffFile=paste(gene_symbol,"_stack_HF.pdf",sep="")
outpdf=paste(picDir,tiffFile,sep="\\")
p1=ggplot(prop_cross,aes(x=prop_cross[,'Var2'], y=prop_cross[,'Freq'],fill=prop_cross[,'Expression level']))+geom_bar(stat='identity',position="stack")+xlab(expression(bold(paste("Histology Classification(", italic("P"),'<0.05)')))) + ylab("Percentage")+theme(plot.title = element_text(size=20,face = "bold"),axis.text.y=element_text(size=15),axis.text.x=element_text(size=10),axis.title.y=element_text(size=15,face="bold"),axis.title.x=element_text(size=15,face="bold"))+ labs(fill=expression(paste("ZBTB4"," Expression Level")))+theme(legend.position="top")
ggsave(p1, file=outpdf, width=8, height=8)

tiffFile=paste(gene_symbol,"_stack_HF.png",sep="")
outpdf=paste(picDir,tiffFile,sep="\\")
p1=ggplot(prop_cross,aes(x=prop_cross[,'Var2'], y=prop_cross[,'Freq'],fill=prop_cross[,'Expression level']))+geom_bar(stat='identity',position="stack")+xlab(expression(bold(paste("Histology Classification(", italic("P"),'<0.05)')))) + ylab("Percentage")+theme(plot.title = element_text(size=20,face = "bold"),axis.text.y=element_text(size=15),axis.text.x=element_text(size=10),axis.title.y=element_text(size=15,face="bold"),axis.title.x=element_text(size=15,face="bold"))+ labs(fill=expression(bold(paste("ZBTB4"," Expression Level"))))+theme(legend.position="top")
ggsave(p1, file=outpdf, width=8, height=8)
