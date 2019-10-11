load("Correlations.rda")

High=function(A,B,cutoff)
{
  patients=which(FullAdd$Type==A)
  celllines=which(FullAdd$Type==B)
  
  Spearman.subset=as.vector(Spearman$r[patients,celllines])
  Pearson.subset=as.vector(Pearson$r[patients,celllines])
  
  High.summary=c(100*mean(abs(Spearman.subset)>cutoff),100*mean(abs(Pearson.subset)>cutoff))
  return(High.summary)
}

cutoff=90

Results.Spearman=Results.Pearson=matrix(0,length(Tumors),length(Lineages))

for(i in 1:length(Tumors))
{
  for(j in 1:length(Lineages))
  {
    Store=High(Tumors[i],Lineages[j],0.9)
    
    Results.Spearman[i,j]=round(Store[1])
    Results.Pearson[i,j]=round(Store[2])
    
    Results.Spearman[i,j]=Results.Spearman[i,j]*(Results.Spearman[i,j]>cutoff)
    Results.Pearson[i,j]=Results.Pearson[i,j]*(Results.Pearson[i,j]>cutoff)
  }
}

rownames(Results.Spearman)=rownames(Results.Pearson)=Tumors
colnames(Results.Spearman)=colnames(Results.Pearson)=Lineages
colnames(Results.Spearman)[c(7,15)]=colnames(Results.Pearson)[c(7,15)]=c("head-neck","stomach")
write.table(Results.Spearman,file=paste0("Spearman ",eval(cutoff),".txt"),sep="\t",row.names=TRUE,col.names=TRUE)
write.table(Results.Pearson,file=paste0("Pearson ",eval(cutoff),".txt"),sep="\t",row.names=TRUE,col.names=TRUE)