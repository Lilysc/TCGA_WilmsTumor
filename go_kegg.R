library("org.Hs.eg.db")
library(clusterProfiler)

inputfile="input.txt"

gene_symbol=read.table(inputfile,sep="\t",check.names=F,header=F)
gene_name=as.vector(gene_symbol[,1])
geneID <- mget(gene_name, org.Hs.egSYMBOL2EG, ifnotfound=NA)
geneID <- as.character(geneID)
data=cbind(symbol=gene_symbol,entrezID=geneID)

#GO
go <- enrichGO(gene=data$V1, OrgDb = org.Hs.eg.db, ont='ALL',pvalueCutoff = 0.05, keyType = 'SYMBOL')


write.csv(go,"go.csv",row.names =F)
pdf('go_barplot.pdf',width=13)
barplot(go, showCategory=20,ylabel='a')

dev.off()

pdf('go_dot.pdf',width=13)
dotplot(go,showCategory=20)
dev.off()

##KEGG
kegg <- enrichKEGG(gene = data$entrezID,organism ="human",pvalueCutoff = 0.05)

kegg_sym =c()
for(i in 1:length(kegg$geneID)){
  a = unlist(strsplit(kegg$geneID[i],"\\/"))
  kegg_s=''
  for(j in a){
    for(z in 1:nrow(data)){
      if(j == data[z,2]){
        kegg_s = paste0(kegg_s,data[z,1],sep=",")
      }
    }
  }
  kegg_sym[i] = substr(kegg_s, 1,nchar(kegg_s)-1)
}
kegg = cbind(data.frame(kegg),'geneSymbol'=kegg_sym)
kegg = kegg[,-8]
write.csv(kegg,"KEGG.csv",row.names =F)
pdf('kegg_barplot.pdf')
barplot(kegg, showCategory=20)
dev.off()
pdf('kegg_dot.pdf')
dotplot(kegg,showCategory=20)
dev.off()


tiff(file='go_barplot.tiff',width = 25,height = 15,units ="cm",compression="lzw",bg="white",res=1800)
barplot(go, showCategory=20)
dev.off()

tiff(file='go_dot.tiff',width = 25,height = 15,units ="cm",compression="lzw",bg="white",res=1800)
dotplot(go,showCategory=20)
dev.off()

tiff(file='kegg_barplot.tiff',width = 20,height = 15,units ="cm",compression="lzw",bg="white",res=1800)
barplot(kegg, showCategory=20)
dev.off()
tiff(file='kegg_dot.tiff',width = 17,height = 15,units ="cm",compression="lzw",bg="white",res=1800)
dotplot(kegg,showCategory=20)
dev.off()


display_number = c(9, 10, 15)
ego_MF <- enrichGO(gene=data$V1, OrgDb = org.Hs.eg.db, ont='MF',pvalueCutoff = 0.05, keyType = 'SYMBOL')
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ] 
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),] 

ego_CC <- enrichGO(gene=data$V1, OrgDb = org.Hs.eg.db, ont='CC',pvalueCutoff = 0.05, keyType = 'SYMBOL')
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ] 
ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),] 

ego_BP <- enrichGO(gene=data$V1, OrgDb = org.Hs.eg.db, ont='BP',pvalueCutoff = 0.05, keyType = 'SYMBOL')
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ]) 
ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),] 

go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID), Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description), GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count), type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]), rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process"))) 
## numbers as data on x axis 
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms 
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],collapse=" "), "...", sep="")
    return(x)
  }else{
    return(x)
  }
}
labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df)) 
## colors for bar // green, blue, orange 
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5") 
library(ggplot2) 
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) + geom_bar(stat="identity", width=0.8) + coord_flip() + scale_fill_manual(values = CPCOLS) + theme_bw() + scale_x_discrete(labels=labels) + xlab("GO term") + theme(axis.text=element_text(face = "bold", color="gray50")) + labs(title = "The Most Enriched GO Terms") 
pdf("go_enrichment.pdf",width=9,height=8) 
p
dev.off() 
svg("go_enrichment.svg") 
p 
dev.off()

