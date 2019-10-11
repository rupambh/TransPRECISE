require(mgcv)
load("RPPA.rda")
AllFinal=rbind(TCGAFinal,MCLPFinal)
Combined=list()

for(i in 1:nrow(AllFinal))
{
  FullData=do.call(cbind,AllFinal[i,])
  FullData=t(uniquecombs(t(FullData)))
  
  Details=FullData[,1:2]
  Data=matrix(as.numeric(FullData[,-(1:2)]),nrow(FullData),ncol(FullData)-2)
  NewData=data.frame(Details,Data)
  
  colnames(NewData)=colnames(FullData)
  Combined[[i]]=NewData
}

names(Combined)=c(Tumors,Lineages)

save(Combined,PathwayProteins,Lineages,Pathways,Tumors,file="RPPA.rda")
rm(list=ls())