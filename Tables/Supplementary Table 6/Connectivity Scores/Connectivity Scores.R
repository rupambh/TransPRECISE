load("Graphs.rda")
CScores=matrix(0,nrow(MCLPGraphs)+nrow(TCGAGraphs),length(Pathways))
rownames(CScores)=c(Tumors,Lineages)
colnames(CScores)=Pathways
Cutoff=0.3

for(i in 1:nrow(CScores))
{
  for(j in 1:ncol(CScores))
  {
    if(i<=31)
    {
      CScores[i,j]=sum(TCGAGraphs[i,j][[1]]>Cutoff)/(nrow(TCGAGraphs[i,j][[1]])*(nrow(TCGAGraphs[i,j][[1]])-1))
    }
    else
    {
      CScores[i,j]=sum(MCLPGraphs[i-31,j][[1]]>Cutoff)/(nrow(MCLPGraphs[i-31,j][[1]])*(nrow(MCLPGraphs[i-31,j][[1]])-1))
    }
  }
}

DScores=rbind(apply(CScores,2,sd),apply(CScores[1:31,],2,sd),apply(CScores[32:47,],2,sd))
rownames(DScores)=c("All","TCGA","MCLP")

save.image("Scores.rda")