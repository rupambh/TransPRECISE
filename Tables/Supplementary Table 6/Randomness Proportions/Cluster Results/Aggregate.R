CSFull=matrix(list(0),47,12)

for(i in 1:47)
{
  for(j in 1:12)
  {
    for(k in 1:1000)
    {
      load(paste0("output",k,".rda"))
      CSFull[i,j][[1]]=c(CSFull[i,j][[1]],CSDists[i,j][[1]])
    }
    
    CSFull[i,j][[1]]=CSFull[i,j][[1]][-1]
    CSFull[i,j][[1]]=CSFull[i,j][[1]]/2
  }
}

save(CSFull,file="Cluster Output.rda")