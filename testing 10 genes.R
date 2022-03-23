## ---------------------------
## Script name: 10 random genes.
## Date Created: 13-06-2021
## ---------------------------
library(ggplot2)
library(vcfR)
library(cowplot)
library(plotly)
library(ggpubr)
library(reshape)
workdir<-"C:/Users/Eros/Desktop/R_WORKING_DIR/"
setwd(workdir)
library(ggrepel)
library(reshape)
source("FUNCTIONS.R")
##----------------------------
binmap<- read.table("http://www.bioinformatics.nl/~nijve002/Genotypes_Bay_Sha_100Kbin.txt")
binmarkers<- rownames(binmap)
cleanmap <- imputesnps(mapfile,5)
set.seed(158)
geneselection<-sample(1:sum(localgenes$LOD_score>80),10)
o.localgenes<-localgenes[order(localgenes$LOD_score,decreasing = T),]
gnames<-o.localgenes[geneselection,1]
# could have been a apply, but lazy.

g2<-placeholder(gnames[2],cleanmap,snplist,0)               # <- int? 22232 rn "AT3G52680"
g3<-placeholder(gnames[3],cleanmap,snplist,0)
g4<-placeholder(gnames[4],cleanmap,snplist,0)
g5<-placeholder(gnames[5],cleanmap,snplist,0)
g6<-placeholder(gnames[6],cleanmap,snplist,0)
g7<-placeholder(gnames[7],cleanmap,snplist,0)
g8<-placeholder(gnames[8],cleanmap,snplist,0)
g9<-placeholder(gnames[9],cleanmap,snplist,0)  
g10<-placeholder(gnames[10],cleanmap,snplist,0)
 
g2.2<-placeholder(gnames[2],cleanmap,snplist,50000)  #19526968-19529062 




mapfile2<-cleanmap[,-69]
c2<-c_norm_log[,-69]
c2<-c2[,-trm]
mapfile2<-mapfile2[,-trm]
calc_LOD_score.t <- function(marker,gene,mfile) {
  as.numeric(anova(lm(unlist(c2[gene,])[mfile[marker,]*mfile[marker,] == 1] ~ unlist(mfile[marker,])[mfile[marker,]*mfile[marker,] == 1]))[1,5])
}

test1<- mapfile[match(22207,rownames(mapfile)),]

-log10(calc_LOD_score.t(match(22207,rownames(mapfile2)),10868,mapfile2))



 ub<-closest(g2.2[[3]][,2],19529062)
 lb<-closest(g2.2[[3]][,2],19526968)
ub == lb

# g2.2 [[8]]
# g2.2 [[21]]


# ( length(unique(vec)>2 ) binair


#subset local genes

lgenes<- genes[sort(match(localgenes$gene_name,genes))]
snpongene<-sapply(lgenes,chuck)

#find high LOD/snp on gene.

hllg<-localgenes

sg<-snpongene[-grep("FALSE|No snps",snpongene)]

gx2<-placeholder("AT1G76240",cleanmap,snplist,50000)  

# TOTAL of 113 local genes with SNP's on the gene using the cleaned snpmap. allemaal een enkele smh



gx<-placeholder(genes[1],cleanmap,snplist,0)


lui <- function(x) {
  placeholder(localgenenames[x],cleanmap,snplist,0)[[1]][[2]]
}


localgenenames<-o.localgenes[,1]


#-=---------------------------------------------------------
#
#                      DEEL 2
#
#-=---------------------------------------------------------

x1+geom_label_repel(aes(label=colnames(c_norm_log)))


#  "AT4G20480"
int.1<-placeholder(localgenenames[3],cleanmap,snplist,0)
int.1.1<-placeholder(localgenenames[3],mapfile,snplist,0)
bad(int.1,1,6)

# "AT1G59700" (what is going on?)
int.2<-placeholder(localgenenames[34],cleanmap,snplist,0) 
int.2.1<-placeholder(localgenenames[34],mapfile,snplist,0) 

#  "AT4G22285"
int.3<-placeholder(localgenenames[48],cleanmap,snplist,0) 
int.3.1<-placeholder(localgenenames[48],mapfile,snplist,0) 

#  "AT5G56920"
int.4<-placeholder(localgenenames[49],cleanmap,snplist,0) 
int.4.1<-placeholder(localgenenames[49],mapfile,snplist,0) 


# "AT2G36270"  	15204659 - 15207582

# use ph2(placeholder2)
qtl.serin<-readRDS("eQTL_list.rds")
match("AT2G36270",genes)
plot(qtl.serin[[15297]])
mindx<-binmarkers[order(order(qtl.serin[[15297]],decreasing = T))[1]]

sebm<-"RSM_5_18.85"
tg1<-ph2(mindx,cleanmap,snplist,0,"AT2G36270")
tg1.s<-ph2(sebm,cleanmap,snplist,0,"AT2G36270")



#figuurtjes maken "AT3G52680"
lodp<-eQTL_list[[10868]]
plot(lodp)

chromie<-gsub("RSM_|_.*","",binmarkers)

data<-data.frame(lodp,chromie,binmarkers)


ggplot(data,aes(1:1059,lodp))+geom_line()+geom_point(size=0.9)+
  facet_grid(. ~ paste("chromosome",chromie),scales="free")+
  xlab("Binmarker")+
  ylab("LOD-score")+
  ggtitle("AT3G52680")+
  theme(axis.text=element_text(size=10),
  axis.title=element_text(size=18,face="bold"),strip.text.x = element_text(size = 15,face = "bold"))


dc3<-data[data$chromie==3,]

#"RSM_3_19.55" = 614

ggplot(dc3,aes(as.numeric(rownames(dc3)[1]):as.numeric(rownames(dc3)[207]),lodp))+geom_line()+geom_point(size=0.9)+
  xlab("Binmarker")+
  ylab("LOD-score")+
  ggtitle("AT3G52680       chromosome 3")+
  geom_segment(aes(x=614,y=12,xend=614,yend=11),size=1.3,arrow = arrow(length=unit(0.02,"npc")))+
  annotate("text", x=614, y=12.3, label= "A")+
  geom_segment(aes(x=620,y=12.5,xend=620,yend=11.5),size=1.3,arrow = arrow(length=unit(0.02,"npc")))+
  annotate("text", x=620, y=12.8, label= "B")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"))


#ph3 is temp with only + range
placeholder3()

g2c<-placeholder3(gnames[2],cleanmap,snplist,100000)
g2c+ggtitle("huh")
g2d<-placeholder3(gnames[2],mapfile,snplist,100000) 

# mapfiles, imputed and unimputed side by side
pa<-g2d[[1]][[2]]+ rremove("legend")+  theme(legend.position = "none",
                                             panel.grid = element_blank(),
                                             axis.title = element_blank(),
                                             axis.text = element_blank(),
                                             axis.ticks = element_blank(),
                                             panel.background = element_blank()) 
pb<-g2c[[1]][[2]]+ rremove("legend") +  theme(legend.position = "none",
                                              panel.grid = element_blank(),
                                              axis.title = element_blank(),
                                              axis.text = element_blank(),
                                              axis.ticks = element_blank(),
                                              panel.background = element_blank()) 


plot_grid(pa,pb,ncol = 1,
          labels='AUTO', label_size = 12,
          label_x = 0, label_y = 0,
          hjust = -0.5, vjust = -0.5)




#combination plot LOD scores 2 bins. 

combdf<-data.frame(g2d[[4]][,1],g2c[[4]][,1])
colnames(combdf)<-c("raw","imputed")

colors<-c("raw"="red","imputed"="green")
g <- ggplot(combdf, aes(1:64))+
geom_line(aes(y=raw, color="raw"),size=1)+geom_point(aes(y=raw,colour= "raw"),size=1.5)+
geom_line(aes(y=imputed,color="imputed"),size=1) +geom_point(aes(y=imputed,colour= "imputed"),size=1.5)+ylab("LOD score")+
  xlab("marker")+ggtitle("LOD scores of raw and imputed map on chromosome 3 19500000:19600000")+
  labs(color = "Map type") +
  scale_color_manual(values = colors)



colnames(combdf)<-c("raw","imputed")

colors<-c("raw"="red","imputed"="green")
ggplot(combdf, aes(1:64))+ geom_line(aes(y=imputed),size=1)+geom_point(aes(y=imputed),size=1.5)+ylab("LOD score")+
xlab("marker")+ggtitle("LOD scores of imputed map on chromosome 3 19500000:19600000")



exgen<-data.frame(eQTL_list[[10868]],c(1:1059))
colnames(exgen)<-c("a","b")

ggplot(exgen,aes(b,a))+geom_line()


geom_segment(aes(x=620,y=12.5,xend=620,yend=11.5),size=1.3,arrow = arrow(length=unit(0.02,"npc")))


pd2 <- reshape::melt(as.matrix(binmap))


pd2<-pd2[order(pd2$X1,decreasing = F),]
idx<-rep(1:length(unique(pd2$X1)),each=160) 
pd2<-cbind(pd2,idx)
kl<-c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue")

mapplot<-ggplot(pd2, aes(idx, X2,fill=as.factor(value))) +
  geom_tile() +
  ggtitle(paste(gene))+ 
  scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
  theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))







