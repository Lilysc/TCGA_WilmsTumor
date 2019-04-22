library("survival")
library("KMsurv")
library('reshape')

picDir="picture"
dir.create(picDir)

rt=read.table("tumor_surivival.txt",header=T,sep="\t",check.names=F,row.names=1)
df=na.omit(rt)
df[,'fustat'] <- ifelse(df[, 'Vital Status']=='DEAD', 1 ,0)
df <- rename(df,c("Overall Survival Time in Days"="futime"))
df[,"futime"]=df[,"futime"]/365
rt = df
outTab = data.frame()
gene = c()
for(i in colnames(rt[,19:(ncol(rt)-1)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  med=median(rt[,i])
  if(med!=0){
    a=rt[,i]>med
    rt1=rt[a,]
    b=setdiff(rownames(rt),rownames(rt1))
    rt2=rt[b,]
    n1=nrow(rt1)
    n2=nrow(rt2)
    surTab1=summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
    surTab2=summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
    medianTab1=surTab1$table
    medianTab2=surTab2$table
    diff=survdiff(Surv(futime, fustat) ~a,data = rt)
    fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    outTab=rbind(outTab,cbind(gene=i,coxSummary$coefficients,coxSummary$conf.int,KM=pValue,
                              H_med=medianTab1["median"],H_0.95LCL=medianTab1["0.95LCL"],H_0.95UCL=medianTab1["0.95UCL"],
                              L_med=medianTab2["median"],L_0.95LCL=medianTab2["0.95LCL"],L_0.95UCL=medianTab2["0.95UCL"]))
    pval=0
    if(pValue<0.05){
      pval=signif(pValue,4)
      pval=format(pval, scientific = TRUE)
    }else{
      pval=round(pValue,3)
    }
    if(pValue<0.05){
      geneName=unlist(strsplit(i,"\\|"))[1]
      gene = c(gene,i)
      tiffFile=paste(geneName,".survival.tiff",sep="")
      outTiff=paste(picDir,tiffFile,sep="\\")
      tiff(file=outTiff,width = 15,height = 15,units ="cm",compression="lzw",bg="white",res=1800)
      plot(fit, col=c("blue","red"), xlab="Time (years)", ylab="Overall survival",
           main=substitute({geneName1} *"("*{italic(P) == pval1}*")", list(geneName1=geneName,pval1=pval)),mark.time=T,ylim=c(0,1.1),
           lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
      legend("topright", c(paste("Low expression"), 
                           paste("High expression")), 
             col=c("blue","red"), bty="n", lwd = 2, cex=0.8)
      dev.off()
    }
  }
}
write.csv(gene,"important_gene.txt",row.names =F)
