load("All C & D Scores.rda")
dimnames(PValues)=dimnames(CScores)
Cutoff=0.10

Conversion=expand.grid(Lineages,Tumors)

colnames(Conversion)=c("CellLine_tissue","PanCan")

Conversion=Conversion[Conversion$CellLine_tissue=="head and neck",]

Plotthis=NULL

for(i in 1:nrow(Conversion))
{
  for(j in 1:ncol(PValues))
  {
    Box=(PValues[as.character(Conversion$CellLine_tissue[i]),j]>Cutoff)*(PValues[as.character(Conversion$PanCan[i]),j]>Cutoff)
    
    if(Box==1)
    {
      Plotthis=rbind(Plotthis,data.frame(Conversion[i,1],Conversion[i,2],Pathways[j]))
    }
  }
}

colnames(Plotthis)=c("MCLP","TCGA","Pathway")

Plotthis$MCLP=as.factor(Plotthis$MCLP)
Plotthis$Pathway=as.factor(Plotthis$Pathway)
Plotthis$TCGA=as.factor(Plotthis$TCGA)

Everything=Plotthis

library(flipPlots)
library(htmlwidgets)
library(webshot)

saveWidget(SankeyDiagram(Everything[,c(1,3,2)],font.size=0,label.show.varname=FALSE,link.color="Source"),"Everything.html")
webshot("Everything.html","Everything.png")