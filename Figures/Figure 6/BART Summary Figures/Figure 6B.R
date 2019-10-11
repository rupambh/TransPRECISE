library(MESS)
library(dplyr)
library(igraph)

load("Results.rda")
MissPercentage=MissPercentage[-1]
options(stringsAsFactors=FALSE)

Top2=function(x)
{
  A=which.max(x)
  x[which.max(x)]=0
  B=which.max(x)
  
  return(sort(names(x)[c(A,B)]))
}

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
    Input=data.frame(t(apply(Input,2,Top2)))
    Input=Input%>%group_by(X1,X2)%>%summarize(count=n())
    colnames(Input)=c("from","to","weight")
    Input=data.frame(Input)
    
    Removed=setdiff(Pathways,c(Input[,1],Input[,2]))
    Input=rbind(Input,data.frame(from=Removed,to=Removed,weight=0))
    Input=Input[order(Input$to),]
    
    Input=Input[order(Input$from),]
    
    m=data.frame(c(Input[,1],Input[,2]),rep(Input[,3],2))
    m=m%>%group_by(c.Input...1...Input...2..)%>%tally()
    
    net=graph_from_data_frame(Input,directed=FALSE)
    net=delete.edges(net,which(E(net)$weight==0))
    
    V(net)$color="indianred3"
    
    V(net)$size=3*(m$n)
    E(net)$width=E(net)$weight/3
    E(net)$colors="gray"
    
    E(net)$arrow.size=0
    
    ps.options(font="serif")
    bitmap(paste0(gsub("/","-",Cancers[i]),".tiff"),width=3000,height=3000,units='px',res=720,colormodel="cmyk")
    plot(net,layout=layout.circle(net),vertex.label.color="black",vertex.label.cex=1,edge.label=E(net)$weight,edge.color=E(net)$colors,edge.label.color="blue",main=Cancers[i])
    dev.off()
  })