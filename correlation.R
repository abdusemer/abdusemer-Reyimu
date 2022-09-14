###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######������ѧ��: http://www.biowolf.cn/
######�������䣺2749657388@qq.com
######����΢��: 18520221056

#install.packages("corrplot")

library(corrplot)
setwd("")      #���ù���Ŀ¼
rt=read.table(".txt",sep="\t",header=T,row.names=1,check.names=F)
rt=t(rt)
res1 <- cor.mtest(rt, conf.level = 0.95)

pdf("correlation.pdf",height=8,width=8)              #����ͼƬ���ļ�����
corrplot(corr=cor(rt),
         method = "circle",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         p.mat = res1$p,
         sig.level = 0.001,
         insig = "pch",
         number.cex = 1,
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50),
         )
dev.off()

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######������ѧ��: http://www.biowolf.cn/
######�������䣺2749657388@qq.com
######����΢��: 18520221056