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

Exposure=read.csv("Exposure.csv")
Drugs.CL=tolower(sort(colnames(Drugs)[-(1:2)]))
Drugs.PT=tolower(sort(colnames(Exposure)[-1]))
Drugs.Both=intersect(Drugs.CL,Drugs.PT)

i=1

while(i<=ncol(Drugs)-2)
{
  if(tolower(colnames(Drugs)[i+2])%in%Drugs.Both)
  {
    i=i+1
  }
  else
  {
    Drugs=Drugs[,-(i+2)]
  }
}

i=1

while(i<=ncol(Exposure)-1)
{
  if(tolower(colnames(Exposure)[i+1])%in%Drugs.Both)
  {
    i=i+1
  }
  else
  {
    Exposure=Exposure[,-(i+1)]
  }
}

TCGA$ShortID=substr(as.character(TCGA$ID),1,12)
CommonSamples=intersect(TCGA$ShortID,Exposure$PatientNames)
Exposure=Exposure[Exposure$PatientNames%in%CommonSamples,]
TCGA=TCGA[TCGA$ShortID%in%CommonSamples,]
TCGA=TCGA[!duplicated(TCGA$ShortID),]

TCGA=TCGA[order(TCGA$ShortID),]
Exposure=Exposure[order(Exposure$PatientNames),]
Exposure$Type=TCGA$Type

Lineages=readxl::read_xlsx("Lineages.xlsx")

require(bartMachine)

PredictMatch=function(i,j)
{
  PTType=Lineages$TCGA[i]
  CLType=Lineages$Tissue[i]
  DrugName=Drugs.Both[j]
  
  y=Drugs[Drugs$Type==CLType,tolower(colnames(Drugs))==DrugName]
  x=MCLP[MCLP$Type==CLType,-(1:2)]
  x=x[is.na(y)==FALSE,]
  y=y[is.na(y)==FALSE]
  
  n=length(y)
  
  yTrue=Exposure[Exposure$Type==PTType,tolower(colnames(Exposure))==DrugName]
  xPredict=TCGA[TCGA$Type==PTType,3:14]
  
  ModelFitted=bartMachine(x,y)
  yHat=predict(ModelFitted,xPredict)
  return(list(yTrue,yHat,n))
}

Errors=SamplesUsed=matrix(0,nrow(Lineages),length(Drugs.Both))
True=Estimates=matrix(list(),nrow(Lineages),length(Drugs.Both))

colnames(True)=colnames(Estimates)=colnames(SamplesUsed)=colnames(Errors)=Drugs.Both
rownames(True)=rownames(Errors)=rownames(Estimates)=rownames(SamplesUsed)=Lineages$TCGA

for(i in 1:nrow(Lineages))
{
  for(j in 1:length(Drugs.Both))
  {
    Results=try(PredictMatch(i,j))
    
    if(isTRUE(class(Results)=="try-error"))
    {
      Errors[i,j]=1
      next
    }
    else
    {
      True[i,j][[1]]=Results[[1]]
      Estimates[i,j][[1]]=Results[[2]]
      SamplesUsed[i,j]=Results[[3]]
    }
    
    cat("\014")
    
    print(c(i,j))
  }
}

save.image(file="Results.rda")
rm(list=ls())