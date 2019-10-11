load("Input.rda")

for(i in 1:length(Pathways))
{
  x=PanCanGraphsMCLP[[i]]
  x=x[order(rownames(x)),order(colnames(x))]
  x=pmax(x,t(x))
  
  PanCanGraphsMCLP[[i]]=x
  
  x=PanCanGraphsTCGA[[i]]
  x=x[order(rownames(x)),order(colnames(x))]
  x=pmax(x,t(x))
  
  PanCanGraphsTCGA[[i]]=x
}

for(i in 1:length(Pathways))
{
  Keep=intersect(rownames(PanCanGraphsTCGA[[i]]),rownames(Types[[i]]))
  Keep=intersect(rownames(PanCanGraphsMCLP[[i]]),Keep)
  
  PanCanGraphsMCLP[[i]]=PanCanGraphsMCLP[[i]][rownames(PanCanGraphsMCLP[[i]])%in%Keep,rownames(PanCanGraphsMCLP[[i]])%in%Keep]
  PanCanGraphsTCGA[[i]]=PanCanGraphsTCGA[[i]][rownames(PanCanGraphsTCGA[[i]])%in%Keep,rownames(PanCanGraphsTCGA[[i]])%in%Keep]
  Types[[i]]=Types[[i]][rownames(Types[[i]])%in%Keep,rownames(Types[[i]])%in%Keep]
}

Cutoff=0.5

MCLPOld=MCLPNew=TCGAOld=TCGANew=list()

for(i in 1:length(Pathways))
{
  y=Types[[i]]
  
  MCLPOld[[i]]=MCLPNew[[i]]=""
  
  x=PanCanGraphsMCLP[[i]]
  for(j in 1:(nrow(Types[[i]])-1))
  {
    for(k in (j+1):ncol(Types[[i]]))
    {
      if(x[j,k]>=16*Cutoff)
      {
        if(y[j,k]==1)
        {
          MCLPOld[[i]]=paste(MCLPOld[[i]],paste(paste(rownames(x)[j],colnames(x)[k],sep="-"),paste("(",bquote(.(x[j,k])),")",sep=""),sep=""),sep=", ")
        }
        else
        {
          MCLPNew[[i]]=paste(MCLPNew[[i]],paste(paste(rownames(x)[j],colnames(x)[k],sep="-"),paste("(",bquote(.(x[j,k])),")",sep=""),sep=""),sep=", ")
        }
      }
    }
  }
  names(MCLPOld)[i]=names(MCLPNew)[i]=Pathways[i]
  
  TCGAOld[[i]]=TCGANew[[i]]=""
  
  x=PanCanGraphsTCGA[[i]]
  for(j in 1:(nrow(Types[[i]])-1))
  {
    for(k in (j+1):ncol(Types[[i]]))
    {
      if(x[j,k]>=31*Cutoff)
      {
        if(y[j,k]==1)
        {
          TCGAOld[[i]]=paste(TCGAOld[[i]],paste(paste(rownames(x)[j],colnames(x)[k],sep="-"),paste("(",bquote(.(x[j,k])),")",sep=""),sep=""),sep=", ")
        }
        else
        {
          TCGANew[[i]]=paste(TCGANew[[i]],paste(paste(rownames(x)[j],colnames(x)[k],sep="-"),paste("(",bquote(.(x[j,k])),")",sep=""),sep=""),sep=", ")
        }
      }
    }
  }
  names(TCGAOld)[i]=names(TCGANew)[i]=Pathways[i]
}

library(xlsx)
write.xlsx(MCLPOld,file="Findings.xlsx",sheetName="MCLP Old")
write.xlsx(MCLPNew,file="Findings.xlsx",sheetName="MCLP New",append=TRUE)
write.xlsx(TCGAOld,file="Findings.xlsx",sheetName="TCGA Old",append=TRUE)
write.xlsx(TCGANew,file="Findings.xlsx",sheetName="TCGA New",append=TRUE)
save.image("All Tables.rda")
rm(list=ls())