## ---------------------------
## Script name: Source Functions
## Date Created: 11-05-2021
## ---------------------------

workdir<-"C:/Users/Eros/Desktop/R_WORKING_DIR/"
setwd(workdir)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

normalize<-function(x){
  x/sum(x)*1e6
}


calc_LOD_score <- function(marker,gene,mfile) {
  as.numeric(anova(lm(unlist(c_norm_log[gene,])[mfile[marker,]*mfile[marker,] == 1] ~ unlist(mfile[marker,])[mfile[marker,]*mfile[marker,] == 1]))[1,5])
}

closest.marker <- function(chrom,pos) {
  x<-clean.snplist[clean.snplist[,1] == paste(chrom,"",sep = ""),]
  x[which.min(abs(x[,2]-pos)),]
}

getgenesnp <- function(gene) {
  g.row<-match(gene,genelocations[,1])
  if (is.na(g.row)) {
    stop("no match")
  }else{
    g.data<-genelocations[g.row,]
    g.snplist<-snplist[snplist$CHROM == g.data[,2] & snplist$POS >= g.data[,3]-1e5 & snplist$POS <= g.data[,4]+1e5,1:2]
    g.map.subset<-mapfile[rownames(g.snplist),]
    list(g.data,g.snplist,g.map.subset)
  }
}

scrubsnps <- function(table,window) {
  for (x in 1:ncol(table)) {
    tmpv<-table[,x]
    for(i in 1:length(tmpv)){
      u<-ifelse(i + window<=length(tmpv),i + window,length(tmpv))
      l<-ifelse(i - window>0,i-window,1)
      tmpv[i]<-getmode(tmpv[l:u])
    }
    table[,x]<-tmpv
  } 
  return(table)
}   


# NEEDS IMPROVING 
# Funny lines . Improve imputing (dont compare parents if 0.5).


imputesnps <- function(table,window) {
  for (x in 1:ncol(table)) {
    tmpv<-table[,x]
    for(i in 1:length(tmpv)){
      if ((tmpv[i])^2!=1){
        if (ifelse(i>1,tmpv[i-1],1) == ifelse(i<length(tmpv),tmpv[i+1],length(tmpv))) {
          tmpv[i]<-tmpv[i+1]
        }else{
          u<-ifelse(i + window<=length(tmpv),i + window,length(tmpv))          ###### look
          l<-ifelse(i - window>0,i-window,1)
          tmpv[i]<-getmode(tmpv[l:u])
        }
      }else{
        tmpv[i]<-ifelse(ifelse(i>1,tmpv[i-1],1) == ifelse(i<length(tmpv),tmpv[i+1],length(tmpv)),tmpv[i+1],tmpv[i])
      }
      tmpv[i]<-tmpv[i]
    }
    table[,x]<-tmpv
  } 
  return(table)
}

#files?\



genelocations<- read.table("GeneLoc.txt")
data<-readRDS("snpdata.rds")
snplist<-data[[2]]
mapfile<-data[[1]]
counts<-read.table("http://www.bioinformatics.nl/~nijve002/Serin_etal_counts.txt")
counts$RIL72<-NULL
genes<-rownames(counts)
fullmap<-readRDS("snpdata.rds")[[1]]
localgenes<-read.csv("Serin_etal_2017_local.csv")
c_norm <- apply(counts,2,normalize)
c_norm_log = apply(c_norm + 0.1,2,log2)
spls<-read.table("samples_athal.txt")
snames<-gsub("At_|_..","",spls$V2[grep("...RIL.*_..",spls[,2])]) #"RIL[0-9]*"
snames<-snames[-grep("RIL72",snames)]

tempmapfile<-imputesnps(mapfile,5)

placeholder <- function(gene,mapfile,snplist,ex.range) {
  if (is.na(match(gene,localgenes$gene_name))) {
    stop("not a local gene")}
  if (nrow(mapfile)!=nrow(snplist)) {
    stop("incompatible map/snplist")}   
  x<-localgenes[match(gene,localgenes$gene_name),2]
  chr<-gsub("RSM_|_.*","",x)
  bp<-as.numeric(gsub(".*_","",x))*1e6
  gss<- mapfile[snplist$CHROM==chr & snplist$POS >= bp-(5e4+ex.range) & snplist$POS <= bp + (5e4+ex.range), ]
  if (nrow(gss)==0) {
    stop("No SNP's in this interval")
  }
  print(paste(nrow(gss),"snps found",(5e4+ex.range),"basepairs around",x))
  #-object 1 3 plots + mapfile.
  list<-list(list())
  colnames(gss)<- snames
  plotDat <- reshape::melt(gss)
  plotDat<-plotDat[order(plotDat$X1,decreasing = F),]
  idx<-rep(1:length(unique(plotDat$X1)),each=160) 
  plotDat<-cbind(plotDat,idx)
  kl<-c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue")
  mapplot<-ggplot(plotDat, aes(idx, X2,fill=as.factor(value))) +
    geom_tile() +
    ggtitle(paste(gene))+ 
    scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
    theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))
  
  list[[1]][[2]]<-mapplot
  
  lodscores<-sapply(1:nrow(gss), function(x) -log10(calc_LOD_score(x,match(gene,genes),gss)))   # Lod Score, only peak                    #plot  
  tdf<-data.frame(lodscores,1:length(lodscores))
  colnames(tdf)<-c("A","B")
  xplot<-ggplot(data = tdf, aes(x=B,y=A))+
    geom_point()
  list[[1]][[3]]<-xplot
  
  xplot <- xplot + rremove("legend") +  theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
  
  combiplot<-plot_grid(xplot,NULL,mapplot, ncol = 1, align = "h", rel_widths = c(1,-0.2, 1), rel_heights = c(1,-0.15, 1),axis = "l")
  list[[1]][[1]]<-combiplot
  #=================== # creating and adding to list in loop is very slow, fix if possible !
  
  templist<-vector(mode="list",length = nrow(gss))
  for (x in 1:nrow(gss)) {
    snpvec<-gss[x,]
    tbl<-data.frame(as.numeric(snpvec),as.numeric(c_norm_log[match(gene,genes),]))
    colnames(tbl)<-c("allele", "expr")
    plot<-ggplot(data = tbl,aes(allele, expr, color = as.factor(allele)))+
      geom_jitter(width = 0.05)+
      scale_color_manual(name= "allele",values = c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue"))+
      ggtitle(paste(snplist[match(rownames(gss)[x],rownames(snplist)),1],"-",snplist[match(rownames(gss)[x],rownames(snplist)),2]))+
      theme(legend.position = "top")+
      scale_fill_discrete(labels = c("bay", "no data", "heterozygous","sha"))
    
    templist[[x]]<-plot
  }
  snpdata<-snplist[match(rownames(gss),rownames(snplist)),]
  list[[2]]<-templist
  list[[3]]<-snpdata
  list[[4]]<-tdf
  return(list)
}

bad <- function(obj,lb,ub) {
  pl<-obj[[2]]
  grid<-plot_grid(plotlist=pl[lb:ub])
  return(grid)  
}


sneed <- function(obj,rng) {
  pl<-obj[[2]]
  x<-order(obj[[4]][,1],decreasing = T)[1]
  ub<-ifelse(x+rng>=length(obj[[2]]),x+rng,length(obj[[2]]))
  lb<-ifelse(x-rng>=0,x-rng,1)
  grid<-plot_grid(plotlist=pl[lb:ub])
  return(grid)  
}



# --- 74 79.8334   -> 90.84
# ---71.98259
# 74/73/93 -- 105.8763


#uses cleanmap by default
chuck <- function(gene) {
  z<-localgenes[match(gene,localgenes[,1]),2]
  if (is.na(z)) {
    warning("not a local gene")
    return("Hoe?")
  }
  x<-genelocations[match(gene,genelocations[,1]),3:4]
  chr<-gsub("RSM_|_.*","",z)
  bp<-as.numeric(gsub(".*_","",z))*1e6
  gss<- cleanmap[snplist$CHROM==chr & snplist$POS >= bp-(5e4) & snplist$POS <= bp + (5e4), ]
  if (all(is.na(gss))) {
    warning("no snps in bin")
    return("No snps")
  }
  ubg<-x[2]
  lbg<-x[1]
  if (is.vector(gss)) {
    return(rownames(cleanmap)[snplist$CHROM==chr & snplist$POS >= bp-(5e4) & snplist$POS <= bp + (5e4)])
  }
  sig<-gss[gss[,2]>=lbg & gss[,2]<=ubg]
  if (all(is.na(sig))) {
    return(FALSE)
  }else{
    return(TRUE)
  }
}



closest <- function(vec,number) {
  print(paste0(vec[which.min(abs(vec-number))]," is closest to ",number))
  print(number - vec[which.min(abs(vec-number))])
  which.min(abs(vec-number))
}


#if(missing(mapfile)) {
#  mapfile<- mapfile}
#if(missing(snplist)) {
#  snplist <- snplist}
#if(missing(snplist)) {
#  snplist <- snplist}
#if(missing(ex.range)) {
#  ex.range <- 0}


ph2 <- function(marker,mapfile,snplist,ex.range,gene) {
  chr<-gsub("RSM_|_.*","",marker)
  bp<-as.numeric(gsub(".*_","",marker))*1e6
  gss<- mapfile[snplist$CHROM==chr & snplist$POS >= bp-(5e4+ex.range) & snplist$POS <= bp + (5e4+ex.range), ]
  if (nrow(gss)==0) {
    stop("No SNP's in this interval")
  }
  print(paste(nrow(gss),"snps found",(5e4+ex.range),"basepairs around",marker))
  #-object 1 3 plots + mapfile.
  list<-list(list())
  colnames(gss)<- snames
  plotDat <- reshape::melt(gss)
  plotDat<-plotDat[order(plotDat$X1,decreasing = F),]
  idx<-rep(1:length(unique(plotDat$X1)),each=160) 
  plotDat<-cbind(plotDat,idx)
  kl<-c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue")
  mapplot<-ggplot(plotDat, aes(idx, X2,fill=as.factor(value))) +
    geom_tile() +
    ggtitle(paste(gene))+ 
    scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
    theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))
  
  list[[1]][[2]]<-mapplot
  
  lodscores<-sapply(1:nrow(gss), function(x) -log10(calc_LOD_score(x,match(gene,genes),gss)))   # Lod Score, only peak                    #plot  
  tdf<-data.frame(lodscores,1:length(lodscores))
  colnames(tdf)<-c("A","B")
  xplot<-ggplot(data = tdf, aes(x=B,y=A))+
    geom_point()
  list[[1]][[3]]<-xplot
  
  xplot <- xplot + rremove("legend") +  theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
  
  combiplot<-plot_grid(xplot,NULL,mapplot, ncol = 1, align = "h", rel_widths = c(1,-0.2, 1), rel_heights = c(1,-0.15, 1),axis = "l")
  list[[1]][[1]]<-combiplot
  #=================== # creating and adding to list in loop is very slow, fix if possible !
  
  templist<-vector(mode="list",length = nrow(gss))
  for (x in 1:nrow(gss)) {
    snpvec<-gss[x,]
    tbl<-data.frame(as.numeric(snpvec),as.numeric(c_norm_log[match(gene,genes),]))
    colnames(tbl)<-c("allele", "expr")
    plot<-ggplot(data = tbl,aes(allele, expr, color = as.factor(allele)))+
      geom_jitter(width = 0.05)+
      scale_color_manual(name= "allele",values = c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue"))+
      ggtitle(paste(snplist[match(rownames(gss)[x],rownames(snplist)),1],"-",snplist[match(rownames(gss)[x],rownames(snplist)),2]))+
      theme(legend.position = "top")+
      scale_fill_discrete(labels = c("bay", "no data", "heterozygous","sha"))
    
    templist[[x]]<-plot
  }
  snpdata<-snplist[match(rownames(gss),rownames(snplist)),]
  list[[2]]<-templist
  list[[3]]<-snpdata
  list[[4]]<-tdf
  return(list)
}




placeholder3 <- function(gene,mapfile,snplist,ex.range) {
  if (is.na(match(gene,localgenes$gene_name))) {
    stop("not a local gene")}
  if (nrow(mapfile)!=nrow(snplist)) {
    stop("incompatible map/snplist")}   
  x<-localgenes[match(gene,localgenes$gene_name),2]
  chr<-gsub("RSM_|_.*","",x)
  bp<-as.numeric(gsub(".*_","",x))*1e6
  gss<- mapfile[snplist$CHROM==chr & snplist$POS >= bp-(5e4) & snplist$POS <= bp + (5e4+ex.range), ]
  if (nrow(gss)==0) {
    stop("No SNP's in this interval")
  }
  print(paste(nrow(gss),"snps found",(5e4),"basepairs around",x))
  #-object 1 3 plots + mapfile.
  list<-list(list())
  colnames(gss)<- snames
  plotDat <- reshape::melt(gss)
  plotDat<-plotDat[order(plotDat$X1,decreasing = F),]
  idx<-rep(1:length(unique(plotDat$X1)),each=160) 
  plotDat<-cbind(plotDat,idx)
  kl<-c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue")
  mapplot<-ggplot(plotDat, aes(idx, X2,fill=as.factor(value))) +
    geom_tile() +
    ggtitle(paste(gene))+ 
    scale_fill_manual(values = kl)+ theme_grey(base_size = 9)+
    theme(legend.position = "bottom")+ guides(fill=guide_legend(title="allele"))
  
  list[[1]][[2]]<-mapplot
  
  lodscores<-sapply(1:nrow(gss), function(x) -log10(calc_LOD_score(x,match(gene,genes),gss)))   # Lod Score, only peak                    #plot  
  tdf<-data.frame(lodscores,1:length(lodscores))
  colnames(tdf)<-c("A","B")
  xplot<-ggplot(data = tdf, aes(x=B,y=A))+
    geom_point()
  list[[1]][[3]]<-xplot
  
  xplot <- xplot + rremove("legend") +  theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())
  
  combiplot<-plot_grid(xplot,NULL,mapplot, ncol = 1, align = "h", rel_widths = c(1,-0.2, 1), rel_heights = c(1,-0.15, 1),axis = "l")
  list[[1]][[1]]<-combiplot
  #=================== # creating and adding to list in loop is very slow, fix if possible !
  
  templist<-vector(mode="list",length = nrow(gss))
  for (x in 1:nrow(gss)) {
    snpvec<-gss[x,]
    tbl<-data.frame(as.numeric(snpvec),as.numeric(c_norm_log[match(gene,genes),]))
    colnames(tbl)<-c("allele", "expr")
    plot<-ggplot(data = tbl,aes(allele, expr, color = as.factor(allele)))+
      geom_jitter(width = 0.05)+
      scale_color_manual(name= "allele",values = c("-1" = "red", "0" = "gray", "0.5" = "black", "1" = "blue"))+
      ggtitle(paste(snplist[match(rownames(gss)[x],rownames(snplist)),1],"-",snplist[match(rownames(gss)[x],rownames(snplist)),2]))+
      theme(legend.position = "top")+
      scale_fill_discrete(labels = c("bay", "no data", "heterozygous","sha"))
    
    templist[[x]]<-plot
  }
  snpdata<-snplist[match(rownames(gss),rownames(snplist)),]
  list[[2]]<-templist
  list[[3]]<-snpdata
  list[[4]]<-tdf
  return(list)
}




























