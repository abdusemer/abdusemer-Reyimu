######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("survivalROC")

library(survivalROC)
setwd("C:\\Users\\lexb4\\Desktop\\����\\TCGAcox\\15.ROC")
rt=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)
tiff(file="ROC.tiff",
       width = 14,            #ͼƬ�Ŀ���
       height =14,            #ͼƬ�ĸ߶�
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
      predict.time =3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
  xlab="False positive rate", ylab="True positive rate",
  main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056