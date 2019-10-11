load("Scores.rda")

Drugs=read.delim("Drugs.txt")
Available=intersect(Drugs$UID,MCLP$ID)
MCLP=MCLP[MCLP$ID%in%Available,]

Drugs=Drugs[Drugs$UID%in%Available,]
Drugs=data.frame(MCLP$Type[order(MCLP$ID)],Drugs)
Drugs=Drugs[order(Drugs[,1]),]
colnames(Drugs)[1]="Type"

Missing=is.na(Drugs[,-1])
rownames(Missing)=Drugs$UID
MissPercentage=sort(colMeans(Missing))
SampleSizes=sort(summary(MCLP$Type),decreasing=TRUE)
Cancers=names(SampleSizes[SampleSizes>10])
Cancers=sort(Cancers)

require(bartMachine)

TopPathways=function(i,j)
{
  y=Drugs[Drugs$Type==Cancers[i],j+2]
  x=MCLP[MCLP$Type==Cancers[i],-(1:2)]
  x=x[is.na(y)==FALSE,]
  y=y[is.na(y)==FALSE]
  
  ModelFitted=bartMachine(x,y)
  Proportions=get_var_props_over_chain(ModelFitted)
  Probabilities=predict(ModelFitted,x)
  
  return(list(y,Probabilities,Proportions,length(y)))
}

Errors=SamplesUsed=matrix(0,length(Cancers),ncol(Drugs)-2)
True=Estimates=Inclusions=matrix(list(),length(Cancers),ncol(Drugs)-2)
colnames(True)=colnames(Estimates)=colnames(Inclusions)=colnames(SamplesUsed)=colnames(Errors)=colnames(Drugs)[-(1:2)]
rownames(True)=rownames(Errors)=rownames(Inclusions)=rownames(Estimates)=rownames(SamplesUsed)=Cancers

for(i in 1:length(Cancers))
{
  for(j in 1:(ncol(Drugs)-2))
  {
    Results=try(TopPathways(i,j))
    
    if(isTRUE(class(Results)=="try-error"))
    {
      Errors[i,j]=1
      next
    }
    else
    {
      True[i,j][[1]]=Results[[1]]
      Estimates[i,j][[1]]=Results[[2]]
      Inclusions[i,j][[1]]=Results[[3]]
      SamplesUsed[i,j]=Results[[4]]
    }
    
    cat("\014")
    
    print(c(i,j))
  }
}

save.image(file="Results.rda")
rm(list=ls())