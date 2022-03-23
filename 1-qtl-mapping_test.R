##############################
##### seed eQTL analysis #####
##############################

#### prepare the script working environment ####
remove(list = ls())
gc()
set.seed(1000)  

# Set working directory ####
work.dir <- getwd()
setwd(work.dir)

# library
library(ggplot2)
library(dplyr)

# load function
setwd(paste0(work.dir, '/eQTL_pipeline-master-R/R'))
for(i in 1:length(dir())){
  source(dir()[i])
} # read function from Mark
setwd(work.dir)

map.per.marker <- function(trait, marker) {
  model <- lm(terms(trait ~ marker, keep.order = FALSE))
  summ <- summary(model)
  pval <- summ$coefficients[2, 4]
  lod <- -log10(pval)
  eff <- summ$coefficients[2, 1]
  output <- c(lod, eff)
  return(output)
} # map the QTL at a marker location
map.all.markers <- function(trait, markers) {
  eff.out <- rep(NA, nrow(markers))
  pval.out <- rep(NA, nrow(markers))
  for (i in 1:nrow(markers)) {
    if(i == 1) {
      out.tmp <- map.per.marker(trait, markers[i, ])
    }
    if( i != 1 & sum(abs(as.numeric(markers[i-1,])  - as.numeric(markers[i,])), na.rm = T) != 0 ) {
      out.tmp <- map.per.marker(trait, markers[i, ])
    }
    if( i != 1 & sum(abs(as.numeric(markers[i-1,]) - as.numeric(markers[i,])), na.rm = T) == 0 ) {
      out.tmp <- out.tmp
    }
    pval.out[i] <- out.tmp[1]
    eff.out[i] <- out.tmp[2]
    output.lod <- cbind(pval.out, eff.out)
    colnames(output.lod) <- c('LOD', 'Eff')
  }
  return(output.lod)
} # map the QTL using genome wide markers
write.EleQTL <- function(map1.output,filename){
  
  selector <- cbind(trait = rownames(map1.output$LOD), pval = apply(map1.output$LOD,1,max,na.rm=T)) %>%
    data.frame()
  
  rownames(selector) <- NULL
  
  lod <- map1.output$LOD
  lod <- lod[rownames(lod) %in% selector[,1],]
  rownames(lod) <- selector$trait
  colnames(lod) <- map1.output$Marker[,1]
  
  eff <- map1.output$Effect
  eff <- eff[rownames(eff) %in% selector[,1],]
  rownames(eff) <- selector$trait
  colnames(eff) <- map1.output$Marker[,1]
  
  lod.eff <- lod*sign(eff)
  
  dat <- map1.output$Trait
  dat <- dat[rownames(dat) %in% selector[,1],]
  rownames(dat) <- selector$trait
  colnames(dat) <- colnames(map1.output$Map)
  
  map <- map1.output$Map
  rownames(map) <- map1.output$Marker[,1]
  
  marker <- map1.output$Marker
  
  write.table(lod,file=paste(filename,"_lod.txt",sep=""),sep="\t",quote=F)
  write.table(eff,file=paste(filename,"_eff.txt",sep=""),sep="\t",quote=F)
  write.table(lod.eff,file=paste(filename,"_lodxeff.txt",sep=""),sep="\t",quote=F)
  write.table(marker,file=paste(filename,"_marker.txt",sep=""),sep="\t",quote=F)
  write.table(dat,file=paste(filename,"_data.txt",sep=""),sep="\t",quote=F)
  write.table(map,file=paste(filename,"_map.txt",sep=""),sep="\t",quote=F)                         
} # function to convert mapping result to tables

# load data
map <- read.table("http://www.bioinformatics.nl/~nijve002/Genotypes_Bay_Sha_100Kbin.txt")

#traits enzo

counts <- read.table("http://www.bioinformatics.nl/~nijve002/Serin_etal_counts.txt")
counts$RIL72<-NULL

normalize<-function(x){
  x/sum(x)*1e6
}

c_norm <- apply(counts,2,normalize)
c_norm_log = apply(c_norm + 0.1,2,log2)


marker <- as.matrix(read.csv('marker.csv', row.names = 1))
map <- as.matrix(read.csv('genetic-map.csv'))
traits <- as.matrix(read.csv('phenotype-matrix.csv', row.names = 1))


