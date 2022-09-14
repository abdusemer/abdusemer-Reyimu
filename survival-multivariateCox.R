######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)
setwd("C:\\Users\\lexb4\\Desktop\\资料\\TCGAcox\\13.multiCox")
rt=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt[,"futime"]=rt[,"futime"]/365
cox <- coxph(Surv(futime, fustat) ~ ., data = rt)
cox=step(cox,direction = "both")
riskScore=predict(cox,type="risk",newdata=rt)
summary=summary(cox)
coxGene=rownames(summary$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)
write.table(cbind(id=coxGene,summary$coefficients),
    file="coxResult.xls",
    sep="\t",
    quote=F,
    row.names=F)
cox

tiff(file="forest.tiff",
       width = 18,            #图片的宽度
       height =12,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
ggforest(cox,
         data=rt,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056