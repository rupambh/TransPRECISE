rm(list=ls())
gc()
library(fastcluster)
source("NetworkPrediction.R")
set.seed(101)
pw.array = c("Apoptosis","Breast reactive","Cell cycle","Core reactive","DNA damage response","EMT","PI3K/AKT","RAS/MAPK","RTK","TSC/mTOR","Hormone receptor","Hormone signaling (Breast)")

load("Scores.rda")
samplename = FullAdd[,1]
tumortype = as.character(FullAdd[,2])
scoreMat = matrix(as.numeric(unlist(FullAdd[,-(1:2)])),ncol=ncol(FullAdd)-2)
signMat = matrix(as.numeric(unlist(FullMax[,-(1:2)])),ncol=ncol(FullMax)-2)
pwnames = colnames(FullAdd)[-c(1,2)]
pw.num.nodes = c(9,6,7,5,11,7,7,9,5,5,3,3)
clineages = unique(tumortype)[32:47]
datype = rep("TCGA",length=length(tumortype))
datype[tumortype%in%clineages] = "MCLP"

w = match(pw.array,pwnames)
scoreMat = scoreMat[,w]
signMat = signMat[,w]
pwnames=pwnames[w]
scaled.scoreMat = t(apply(scoreMat,1,function(x)x/pw.num.nodes))
colnames(scaled.scoreMat) = pw.array

###### Optimal Clusters
#load("../data/fit_Gap_scaled_AllScores.RData")
#a = fit$Tab[,3]
#b = fit$Tab[,3] - fit$Tab[,4]
#b = c(b[-1],0)
#as.matrix(a>b) ## maybe 29 clusters?

hclust.cut <- function(x,k) {
    list(cluster=cutree(hclust(dist(x)),k))
}
opt.n = 29 # optimal number of clusters
hc = hclust.cut(x=scaled.scoreMat,k=opt.n)
table(hc$cluster,datype)
#######
source("heatmap.3.R")
library(gplots)
library(RColorBrewer)
library(colorRamps)

#### Sorting the samples within cluster ####
mmat = numeric(0)
cltype = numeric(0)
ttype = numeric(0)
dtype = numeric(0)
for (cc in 1:opt.n) {
    mat = scaled.scoreMat[hc$cluster==cc,]
    mhc = hclust(dist(mat))
    mmat = rbind(mmat,mat[mhc$order,])
    cltype = c(cltype,rep(cc,nrow(mat)))
    ttype = c(ttype,(tumortype[hc$cluster==cc])[mhc$order])
    dtype = c(dtype,(datype[hc$cluster==cc])[mhc$order])

}



########################### Draw Heatmap
library(pheatmap)
library(RColorBrewer)
library(colorRamps)

#cols = c("#00FFFF",c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"),primary.colors(3))[sample(1:30)])
## Tumor type colors
unilineage = unique(tumortype[datype=="MCLP"])
unitumor = unique(tumortype[datype=="TCGA"])

# cols for unitumor #
cols = c("#00FFFF",c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"),primary.colors(10)[c(8,9,5)]))
cols[31] = "chartreuse"
cols[14] = "aquamarine"
cols[18] = "darkcyan"
cols[17] = "grey"
cols[23] = "darkred"
cols[20] = "darkkhaki"
cols[25] = "OrangeRed"

# cols for unilineage (matching color with tumor type)#
lcols = rep(NA,length(unilineage))
lcols[unilineage=="bladder"] = cols[unitumor=="BLCA"]
lcols[unilineage=="blood"] = cols[unitumor=="DLBC"]
lcols[unilineage=="bone"] = cols[unitumor=="SARC"]
lcols[unilineage=="brain"] = cols[unitumor=="GBM"]
lcols[unilineage=="breast"] = cols[unitumor=="BRCA"]
lcols[unilineage=="colon"] = cols[unitumor=="CORE"]
lcols[unilineage=="head and neck"] = cols[unitumor=="HNSC"]
lcols[unilineage=="kidney"] = cols[unitumor=="KIRC"]
lcols[unilineage=="liver"] = cols[unitumor=="LIHC"]
lcols[unilineage=="lung"] = cols[unitumor=="LUSC"]
lcols[unilineage=="ovary"] = cols[unitumor=="OV"]
lcols[unilineage=="pancreas"] = cols[unitumor=="PAAD"]
lcols[unilineage=="sarcoma"] = cols[unitumor=="SARC"]
lcols[unilineage=="skin"] = cols[unitumor=="SKCM"]
lcols[unilineage=="stomach-oesophagus"] = cols[unitumor=="STAD"]
lcols[unilineage=="uterus"] = cols[unitumor=="UCEC"]

unitumorlineage = c(unitumor,unilineage)
unicols = c(cols,lcols)

# Tumor type colors
patientcolors1 = rep(NA,length(ttype))
for (i in 1:length(unitumorlineage)) {
    tt = unitumorlineage[i]
    patientcolors1[ttype==tt] = unicols[i]
}

# Data type colors
patientcolors2 = rep("white",length(dtype))
patientcolors2[dtype=="MCLP"] = "black"

# Cluster colors #
patientcolors3 = rep(NA,length(ttype))
for (i in 1:opt.n) {
    patientcolors3[cltype==i] = unicols[i]
}

rowsidecolors = rbind(patientcolors1,patientcolors1,patientcolors1,patientcolors2,patientcolors3,patientcolors3)
rownames(rowsidecolors) = c("Lineage1","Lineage2","Lineage3","Datatype","Cluster1","Cluster2")


pdf("heatmap_score_TCGA_MCLP.pdf",width=8,height=20)
par(mar=c(0,0,0,0))
fit = heatmap.3(mmat[8354:1,],dendrogram-"none",trace="non",col=brewer.pal(9,"YlOrRd")[c(1,3,5,7,8,9)],labCol="",labRow="",key=T,Colv=F,Rowv=F
,RowSideColors=rowsidecolors[,8354:1])
#fit = heatmap.3(mmat[6844:1,],Colv=F,Rowv=F,dendrogram="none",trace="none",RowSideColors=rowsidecolors[,6844:1],labRow = "",col=brewer.pal(9,"YlOrRd")[c(1,3,7,8,9)],labCol="",key=F)
#legend(-0.03,0.7,legend=unitumor,col=cols,pch=15,bty="n",cex=0.8)
#legend(-0.03,0.3,legend=clnames,col=cols.clusters,pch=15,bty="n",cex=0.8)
dev.off()


png("heatmap_score_TCGA_MCLP.png",width=1000,height=3000,res=100)
par(mar=c(0,0,0,0))
fit = heatmap.3(mmat[8354:1,],dendrogram-"none",trace="non",col=brewer.pal(9,"YlOrRd")[c(1,3,5,7,8,9)],labCol="",labRow="",key=T,Colv=F,Rowv=F
,RowSideColors=rowsidecolors[,8354:1])
#fit = heatmap.3(mmat[6844:1,],Colv=F,Rowv=F,dendrogram="none",trace="none",RowSideColors=rowsidecolors[,6844:1],labRow = "",col=brewer.pal(9,"YlOrRd")[c(1,3,7,8,9)],labCol="",key=F)
#legend(-0.03,0.7,legend=unitumor,col=cols,pch=15,bty="n",cex=0.8)
#legend(-0.03,0.3,legend=clnames,col=cols.clusters,pch=15,bty="n",cex=0.8)
dev.off()

pdf("heatmap_score_legend.pdf",width=6,height=10)
par(mar=c(0,0,0,0))
plot(rep(0.1,47),47:1,pch=15,cex=2,col=unicols,xlab="",ylab="",xlim=c(0,8),axes=F,frame.plot=F)
text(rep(1.5,47),47:1,unitumorlineage,cex=1)
dev.off()


# Cross tab

ta = table(ttype,cltype)
w = match(unitumorlineage,rownames(ta))
cbind(unitumorlineage,rownames(ta)[w])
ta = ta[w,]
write.table(ta,file="table_cancer_cluster.csv",sep=",",quote=F,col.names=T,row.names=T)
