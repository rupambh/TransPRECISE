load("Scores.rda")

library(readxl)
Conversion=read_xlsx("Tissues.xlsx")
Conversion=Conversion[Conversion$CellLine_tissue%in%Lineages,]
Conversion=Conversion[Conversion$PanCan%in%Tumors,]

Boxes=matrix(0,nrow(Conversion),ncol(CScores))
for(i in 1:nrow(Conversion))
{
  for(j in 1:ncol(CScores))
  {
    Boxes[i,j]=(CScores[Conversion$CellLine_tissue[i],j]>2/3)*(CScores[Conversion$PanCan[i],j]>2/3)
  }
}
Boxes=data.frame(Conversion,Boxes)

colnames(Boxes)=c("MCLP","TCGA",Pathways)
write.csv(Boxes,file="High Pairs.csv")