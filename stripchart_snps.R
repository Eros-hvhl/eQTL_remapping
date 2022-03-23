## ---------------------------
## Script name: Plot underlying data in depth snps.
## Date Created: 9-05-2021
## ---------------------------
library(vcfR)
workdir<-"C:/Users/Eros/Desktop/R_WORKING_DIR/"
setwd(workdir)

snplist<-read.table("snpframe.txt", header = T)
counts<-read.table("http://www.bioinformatics.nl/~nijve002/Serin_etal_counts.txt")
counts$RIL72<-NULL
mapfile<-read.table("temp_mapfile.txt", header = T)
mapfile[,93]<-NULL
genes<-rownames(counts)


normalize<-function(x){
  x/sum(x)*1e6
}

c_norm <- apply(counts,2,normalize)
c_norm_log = apply(c_norm + 0.1,2,log2)



calc_LOD_score <- function(marker,gene,mfile) {
  as.numeric(anova(lm(unlist(c_norm_log[gene,])[mfile[marker,] != 0 | 0.5 ] ~ unlist(mfile[marker,])[mfile[marker, ] != 0 | 0.5]))[1,5])
}

sample.data<-read.table("C:/Users/Eros/Desktop/samples_athal.txt")
ril.sample.data<-sample.data[grep("RIL",sample.data$V2),]
ril.sample.data<-ril.sample.data[-68,]

order.rils<-order(order(ril.sample.data[,1]))
mapfile<-mapfile[,order(colnames(mapfile))]
mapfile<-mapfile[,order.rils]

#remove snp with +50% 0/0.5 calls
n.badcall<-apply(mapfile, 1, function(x) sum(x==0,x==0.5))
clean.mapfile<-mapfile[n.badcall<80,]
clean.snplist<-snplist[n.badcall<80,]
rownames(clean.snplist)<-c(1:nrow(clean.snplist))

#usefull: marker list.

markers<-paste(clean.snplist[,1],clean.snplist[,2])

#=====================================================
closest.marker <- function(chrom,pos) {
  x<-clean.snplist[clean.snplist[,1] == paste(chrom,"",sep = ""),]
  x[which.min(abs(x[,2]-pos)),]
}

rownames(clean.mapfile)<-c(1:nrow(clean.mapfile))
###########################
#                         # 
#        AT5G11200        #  RSM_5_3.55 5 	AGI 	nuc_sequence 	3567175 - 3571015 bp
#                         #  
###########################

# closest marker to center of gene.

markers[as.numeric(rownames(closest.marker(5,mean(3567175,3571015))))]

t.l<-clean.mapfile[as.numeric(rownames(closest.marker(5,mean(3567175,3571015)))),]
c.mapdf<-t(rbind(as.numeric(t.l),counts[32959,]))


colnames(c.mapdf)<-c("allele", "expr")

stripchart(expr ~ allele, 
           data = c.mapdf, 
           vertical = TRUE,
           method = "jitter", 
           xlab = "allele",
           main=markers[as.numeric(rownames(closest.marker(5,mean(3567175,3571015))))])

#highest LOD-score. 12748 

t.l<-clean.mapfile[12748,]
c.mapdf<-t(rbind(as.numeric(t.l),counts[32959,]))
colnames(c.mapdf)<-c("allele", "expr")


stripchart(expr ~ allele, 
           data = c.mapdf, 
           vertical = TRUE,
           method = "jitter", 
           xlab = "allele",
           main=markers[12748])
















