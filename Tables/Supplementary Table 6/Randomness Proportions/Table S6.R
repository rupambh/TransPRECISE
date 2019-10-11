load("Scores.rda")
load("Cluster Output.rda")

PValue=function(x,y)
{
  if(x==1)
  {
    return(0)
  }
  else
  {
    return(sum(x<=y)/length(y))
  }
}

PValues=matrix(0,47,12)
{
  for(i in 1:47)
  {
    for(j in 1:12)
    {
      PValues[i,j]=PValue(CScores[i,j],CSFull[i,j][[1]])
    }
  }
}

Print=matrix(0,47,12)
rownames(Print)=c(Tumors,Lineages)
colnames(Print)=Pathways
for(i in 1:47)
{
  for(j in 1:12)
  {
    Print[i,j]=paste0(round(CScores[i,j],2)," (",round(PValues[i,j],3),")")
  }
}

write.csv(Print,file="Table 6.csv")