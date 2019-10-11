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

MissPercentage=MissPercentage[-1]
options(stringsAsFactors=FALSE)

Top=function(x)
{
  y=rep(0,length(x))
  y[which.max(x)]=1
  return(y)
}

Inputs=data.frame(Pathways)

suppressWarnings(
  for(i in 1:length(Cancers))
  {
    Input=data.frame(Pathways)
    
    for(j in 1:length(MissPercentage))
    {
      if(Errors[i,j]==0 & MissPercentage[j]<0.1 & is.na(AUCs[i,j])==0)
      {
        if(AUCs[i,j]>0.85)
        {
          Input=data.frame(Input,log1p(Inclusions[i,j][[1]]))
          colnames(Input)[ncol(Input)]=colnames(Inclusions)[j]  
        }
      }
    }
    
    Input=Input[,-1]
    Input=data.matrix(Input)
    Input=Input[rowSums(Input)!=0,]
    
    Input.final=apply(Input,2,Top)
    dimnames(Input.final)=dimnames(Input)
    Input.final=Input.final[rowMeans(Input.final)>0,]
    
    if(nrow(Input.final)!=0)
    {
      Limit=range(rowMeans(Input.final))
      
      Limit[1]=0.95*Limit[1]
      Limit[2]=1.05*Limit[2]
      
      Removed=setdiff(Pathways,rownames(Input.final))
      if(length(Removed)>0)
      {
        Input.final=rbind(Input.final,matrix(0,length(Removed),ncol(Input.final)))
        rownames(Input.final)[(13-length(Removed)):12]=Removed
      }
      
      Plotted=rowMeans(Input.final)
      Plotted=Plotted[order(names(Plotted))]
      Inputs=data.frame(Inputs,Plotted)
    }
  })

Inputs=Inputs[,-1]
rownames(Inputs)=Pathways
colnames(Inputs)=Cancers

Inputs=as.data.frame(t(Inputs))
Inputs.log=log1p(Inputs)
library(fmsb)

Colours=c(rainbow(8)[-3],"black")
radarchart(Inputs.log,axistype=1,maxmin=FALSE,pcol=Colours)
legend(x=1.5,y=0.75,legend=rownames(Inputs),pch=20,col=Colours,cex=0.9,pt.cex=3,y.intersp=0.45,bty="n")