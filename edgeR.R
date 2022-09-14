######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

logFoldChange=2
padj=0.05

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("edgeR")

rt=read.table("mRNA.txt",sep="\t",header=T,check.names=F)                           #�ĳ��Լ����ļ���
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

#group=c("normal","tumor","tumor","normal","tumor")
group=c(rep("normal",4),rep("tumor",178))                         #���հ�֢��������Ʒ��Ŀ�޸�
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>logFoldChange | diff$logFC<(-logFoldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>logFoldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-logFoldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #������л���У����ı���ֵ��normalizeExp.txt��
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)         #����������У����ı���ֵ��diffmiRNAExp.txt��

#volcano
tiff(file="vol.tiff",
       width = 12,            #ͼƬ�Ŀ���
       height =12,            #ͼƬ�ĸ߶�
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
xMax=8
yMax=10
plot(allDiff$logFC, -log10(allDiff$FDR), xlab="logFC",ylab="-log10(FDR)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>logFoldChange,]
points(diffSub$logFC, -log10(diffSub$FDR), pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-logFoldChange),]
points(diffSub$logFC, -log10(diffSub$FDR), pch=20, col="green",cex=0.4)
abline(v=0,lty=2,lwd=3)
dev.off()

######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056