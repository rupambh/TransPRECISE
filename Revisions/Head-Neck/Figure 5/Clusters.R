library(fastcluster)
load("Combined.rda")
set.seed(69)

pw.array = c("Apoptosis","Breast reactive","Cell cycle","Core reactive","DNA damage response","EMT","PI3K/AKT","RAS/MAPK","RTK","TSC/mTOR","Hormone receptor","Hormone signaling (Breast)")

samplename = FullAdd[,1]
tumortype = as.character(FullAdd[,2])
scoreMat = matrix(as.numeric(unlist(FullAdd[,-(1:2)])),ncol=ncol(FullAdd)-2)
signMat = matrix(as.numeric(unlist(FullMax[,-(1:2)])),ncol=ncol(FullMax)-2)
pwnames = colnames(FullAdd)[-c(1,2)]

clineages = unique(tumortype)[32:47]
pw.num.nodes = c(9,6,7,5,11,7,7,9,5,5,3,3)
datype = rep("TCGA",length=length(tumortype))
datype[tumortype%in%clineages] = "MCLP"

w = match(pw.array,pwnames)
scoreMat = scoreMat[,w]
signMat = signMat[,w]
pwnames=pwnames[w]

scaled.scoreMat = t(apply(scoreMat,1,function(x)x/pw.num.nodes))
colnames(scaled.scoreMat) = pw.array

hclust.cut <- function(x,k) {
    list(cluster=cutree(hclust(dist(x)),k))
}

opt.n = 29 # optimal number of clusters
hc = hclust.cut(x=scaled.scoreMat,k=opt.n)
table(hc$cluster,datype)

#### Sorting the samples within cluster ####

mmat = numeric(0)
sname = character(0)
dtype = numeric(0)
cltype = numeric(0)
ttype = numeric(0)

for (cc in 1:opt.n) {
    mat = scaled.scoreMat[hc$cluster==cc,]
    mhc = hclust(dist(mat))
    
    mmat = rbind(mmat,mat[mhc$order,])
    
    cltype = c(cltype,rep(cc,nrow(mat)))
    ttype = c(ttype,(tumortype[hc$cluster==cc])[mhc$order])
    sname = c(sname,(as.character(samplename)[hc$cluster==cc])[mhc$order])
    dtype = c(dtype,(datype[hc$cluster==cc])[mhc$order])
}

Clustered = data.frame(Sample=sname,Source=dtype,Cancer=ttype,Cluster=cltype)
HNSubset = Clustered[Clustered$Cluster%in%c(3,4,14,28,29)&Clustered$Cancer%in%c(TumorTypes,"head and neck"),]
ScoreHN = mmat[as.numeric(rownames(HNSubset)),]
rownames(ScoreHN) = HNSubset$Sample

HNSubset$Cluster=as.factor(HNSubset$Cluster)
HNSubset$Cluster=plyr::revalue(HNSubset$Cluster,c("3"="A","4"="B","14"="C","28"="D","29"="E"))
Annotation=data.frame(Cluster=as.character(HNSubset$Cluster),Cancer=as.character(HNSubset$Cancer),Source=HNSubset$Source)
rownames(Annotation)=HNSubset$Sample

Colours.Cluster=colorRampPalette(grDevices::rainbow(5))(5)
names(Colours.Cluster)=c("A","B","C","D","E")

Colours.Cancer=colorRampPalette(grDevices::rainbow(length(unique(HNSubset$Cancer))))(length(unique(HNSubset$Cancer)))
names(Colours.Cancer)=unique(HNSubset$Cancer)

Colours.Source=c("white","black")
names(Colours.Source)=c("TCGA","MCLP")

Colours=list(Cluster=Colours.Cluster,Cancer=Colours.Cancer,Source=Colours.Source)

pheatmap::pheatmap(ScoreHN,cluster_rows=FALSE,cluster_cols=FALSE,color=gplots::colorpanel(60,low="white",high="red"),annotation_row=Annotation,show_rownames=F,annotation_colors=Colours)