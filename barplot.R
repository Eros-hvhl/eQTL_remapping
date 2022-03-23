library("RColorBrewer")
display.brewer.all()



snpstats<- data.frame(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4),c(1,2,3,4))
# binmapped, untrimmed, trimmed, (imputed?)
  # numbersnips -> 0 calls -> 0.5 calls -> 1/-1 calls.
dfbm<-c(1059,38554,19769)
dfut<-c(38554,)

# basic 1 vec stats.
ggplot(data.frame(dfbm),aes(seq_along(dfbm),dfbm))+geom_bar(stat="identity")
# 


GT.BAY<- parent.GT[,grep("BAY",colnames(parent.GT))]
GT.SHA<- parent.GT[,grep("SHA",colnames(parent.GT))]

BAY.allele<- apply(GT.BAY,1,getmode)
SHA.allele<-apply(GT.SHA,1,getmode)


sb<-c(sum(binmap==-3),sum(binmap==0),sum(binmap==1),sum(binmap==-1))
tr<-c(sum(mapfile==0),sum(mapfile==0.5),sum(mapfile==1),sum(mapfile==-1))
clm<-c(sum(cleanmap==0),sum(cleanmap==0.5),sum(cleanmap==1),sum(cleanmap==-1))
sdataf<-data.frame(sb,tr,utr)
rownames(sdataf)<-c('no call','heterozyhous','Bay','Sha')
sdataf<-as.matrix(sdataf)

sdataf2<-as.matrix(data.frame(sb,tr,clm))

data_percentage <- apply(sdataf2, 2, function(x){x*100/sum(x,na.rm=T)})
colnames(data_percentage)<-c("binmapped","trimmed","imputed")
rownames(data_percentage)<-c('no call','heterozygous','Bay','Sha')

ggplot(data.frame(dfbm),aes(seq_along(dfbm),dfbm))+geom_bar(stat="identity")

barplot(dfbm,col = brewer.pal(n = 3, name = "Blues"))

par(xpd=NA)
x<-barplot(dfbm,yaxp=c(0, 40000, 8),col = brewer.pal(n = 3, name = "Blues"),
           ylab="Number of SNPs",
           font.axis=2, 
           beside=F, 
           names.arg=c("binmapped","untrimmed","trimmed"),
           font.lab=2,cex.names = 1.5,cex.lab=1.5)
text(x,dfbm+1500,labels=as.character(dfbm),font = 2)


# STACKED BARPLOT (in percentages)

par(xpd=NA)
par(mar=c(5,6,4,1)+.1)
layout(matrix(c(1,2), nrow = 1), widths = c(0.7, 0.2))
barplot(data_percentage, 
        col=brewer.pal(n = 4, name = "Blues"), 
        border="white",
        ylab="percentage of SNPs",
        font.axis=2, 
        beside=F, 
        font.lab=2,cex.names = 1.5,cex.lab=1.5)
legend(3.9, 98, rownames(data_percentage),brewer.pal(n = 4, name = "Blues"))







plotDat1 <- reshape::melt(mapfile)
plotDat1<-plotDat1[order(plotDat1$X1,decreasing = F),]
idx<-rep(1:length(unique(plotDat1$X1)),each=160) 
plotDat<-cbind(plotDat1,idx)

kl<-c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue")
x1<-ggplot(plotDat, aes(idx, X2,fill=as.factor(value))) +
  geom_tile() +
  ggtitle("test")+ 
  scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
  theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))

x2<-ggplot(plotDat, aes(X1, X2,fill=as.factor(value))) +
  geom_tile() +
  ggtitle("test")+ 
  scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
  theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))