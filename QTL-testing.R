map = read.table("http://www.bioinformatics.nl/~nijve002/Genotypes_Bay_Sha_100Kbin.txt")
counts = read.table("http://www.bioinformatics.nl/~nijve002/Serin_etal_counts.txt")

# RIL72 mist in counts.
counts$RIL72<-NULL


genes<-rownames(counts)
markers<-rownames(map)


     
#Normaliseren tot CPM
normalize<-function(x){
  x/sum(x)*1e6
}

c_norm <- apply(counts,2,normalize)
c_norm_log = apply(c_norm + 0.1,2,log2)

gene= "AT5G11170"

LOD_scores = c();
  for (m in markers) {
    expression = unlist(c_norm_log[gene,])
    alleles = unlist(map[1,])
    pval = anova(lm(expression ~ alleles))[1,5];
    LOD_scores = c(LOD_scores,-log10(pval));

  }

LOD_profile = data.frame(LOD_scores,markers)

plot(seq(length(LOD_scores)),LOD_scores)
