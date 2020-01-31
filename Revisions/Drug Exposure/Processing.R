load("Results.rda")

AUCs=Matches=matrix(0,nrow(Lineages),length(Drugs.Both))

GetAUC=function(i,j,grids)
{
  YTrue=True[i,j][[1]]
  YHat=as.numeric(Estimates[i,j][[1]]>0.5)
  
  evals=seq.default(0,1,length.out=grids)
  
  if(sum(YTrue==1)==0)
  {
    TPR=rep(0,length(evals))
  } else {
    TPR=unlist(lapply(evals,function(x){return(sum((YTrue==1)*(YHat>x))/sum(YTrue==1))}))
  }
  
  if(sum(YTrue==0)==0)
  {
    FPR=rep(0,length(evals))
  } else {
    FPR=unlist(lapply(evals,function(x){return(sum((YTrue==0)*(YHat>x))/sum(YTrue==0))}))
  }
  
  return(MESS::auc(FPR,TPR))
}

for(i in 1:nrow(Lineages))
{
  for(j in 1:length(Drugs.Both))
  {
    if(SamplesUsed[i,j]>=10&Errors[i,j]==0)
    {
      Result=try(GetAUC(i,j,1000))
      
      if(isTRUE(class(Result)=="try-error"))
      {
        Errors[i,j]=1
        print(c(i,j))
        next
      }
      else
      {
        AUCs[i,j]=Result
        print(c(i,j))
      }
    }
  }
}

for(i in 1:nrow(Lineages))
{
  for(j in 1:length(Drugs.Both))
  {
    if(SamplesUsed[i,j]>=10&Errors[i,j]==0)
    {
      YTrue=True[i,j][[1]]
      YHat=as.numeric(Estimates[i,j][[1]]>0.99)
      Matches[i,j]=sum(YHat==YTrue)/length(YTrue)
    }
  }
}