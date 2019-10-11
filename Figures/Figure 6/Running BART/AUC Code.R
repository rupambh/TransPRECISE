library(MESS)
load("Results.rda")
AUCs=matrix(0,nrow(Errors),ncol(Errors))
dimnames(AUCs)=dimnames(Errors)
GetAUC=function(i,j,grids)
{
  YTrue=True[i,j][[1]]
  YHat=Estimates[i,j][[1]]
  YTrue=unlist(lapply(YTrue,function(x){ifelse(x=="resistant",1,0)}))
  evals=seq.default(0,1,length.out=grids)
  
  TPR=unlist(lapply(evals,function(x){return(sum((YTrue==1)*(YHat>x))/sum(YTrue==1))}))
  FPR=unlist(lapply(evals,function(x){return(sum((YTrue==0)*(YHat>x))/sum(YTrue==0))}))
  return(auc(FPR,TPR,type="spline"))
}

for(i in 1:nrow(Errors))
{
  for(j in 1:ncol(Errors))
  {
    if(Errors[i,j]==0)
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