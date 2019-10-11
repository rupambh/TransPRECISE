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

Results.Spearman=Results.Pearson=matrix(0,length(Tumors),length(Lineages))

for(i in 1:length(Tumors))
{
  for(j in 1:length(Lineages))
  {
    Store=High(Tumors[i],Lineages[j],0.9)
    Results.Spearman[i,j]=round(Store[1],2)
    Results.Pearson[i,j]=round(Store[2],2)
  }
}

rownames(Results.Spearman)=rownames(Results.Pearson)=Tumors
colnames(Results.Spearman)=colnames(Results.Pearson)=Lineages
write.csv(Results.Spearman,file="Spearman.csv")
write.csv(Results.Pearson,file="Pearson.csv")