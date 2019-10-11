load("PPI.rda")
load("Graphs.rda")
load("Input.rda")

Lineages[c(7,15)]=c("head-neck","stomach")

require(e1071)
require(igraph)
require(reshape2)
require(matrixcalc)
require(ggplot2)
require(gplots)
require(dils)

#Fill Up Missing Graphs#

for(j in 10:11)
{
  for(i in 1:31)
  {
    TCGAGraphs[i,j][[1]]=TCGAGraphs[i,j][[1]][-(14-j),-(14-j)]
  }
}

c=c(129,131,150,152)

for(i in 1:length(c))
{
  TCGAGraphs[[c[i]]]=matrix(runif(11^2),11,11)
}

c=c(99,100)

for(i in 1:length(c))
{
  MCLPGraphs[[c[i]]]=matrix(runif(3^2),3,3)
}

#Create Fitted Pan-Cancer Graphs#

alpha=0.5

PanCanGraphsTCGA=PanCanGraphsMCLP=list()

for(j in 1:length(Pathways))
{
  PanCanGraphsTCGA[[j]]=matrix(0,nrow(TCGAGraphs[1,j][[1]]),ncol(TCGAGraphs[1,j][[1]]))
  
  for(i in 1:length(Tumors))
  {
    PanCanGraphsTCGA[[j]]=PanCanGraphsTCGA[[j]]+(TCGAGraphs[i,j][[1]]>alpha)
  }
  
  PanCanGraphsMCLP[[j]]=matrix(0,nrow(MCLPGraphs[1,j][[1]]),ncol(MCLPGraphs[1,j][[1]]))
  
  for(i in 1:length(Lineages))
  {
    PanCanGraphsMCLP[[j]]=PanCanGraphsMCLP[[j]]+(MCLPGraphs[i,j][[1]]>alpha)
  }
}

#Process PPI Information#
cutoff=0.5

Types=list()

for(i in 1:length(ppi$ppi))
{
  Pairs=ppi$ppi[[i]]
  Pairs=data.frame(Pairs[,1],Pairs[,2],1*(as.numeric(Pairs[,3])>cutoff*1000))
  PPI=AdjacencyFromEdgelist(Pairs)
  
  rownames(PPI$adjacency)=colnames(PPI$adjacency)=PPI$nodelist
  Types[[i]]=PPI$adjacency
}

#Plot Common Parts of Figure 2#

ps.options(fonts="serif")

bitmap("Legends 1.tiff",width=4500,height=2000,units='px',res=720,colormodel="cmyk")

plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x="top",ncol=2,legend=c("TCPA","MCLP"),fill=c("black","indianred3"),title="Data Source")
legend(x="bottom",ncol=2,legend=c(paste("PPI Score <",eval(cutoff)),paste("PPI Score >",eval(cutoff))),col=rep("black",2),lty=c(2,1),title="Line Type")

dev.off()

bitmap("Legends 2.tiff",width=4500,height=2000,units='px',res=720,colormodel="cmyk")

plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend(x="top",ncol=2,legend=c(paste("Number of Cancer Sites <",eval(round(cutoff*31))),paste("Number of Cancer Sites >",eval(round(cutoff*31)))),col=c("gray","blue"),lty=rep(1,2),title="Line Colour")
legend(x="bottom",ncol=2,legend=c(paste("Number of Cancer Sites <",eval(round(cutoff*16))),paste("Number of Cancer Sites >",eval(round(cutoff*16)))),col=c("gray","blue"),lty=rep(1,2),title="Line Colour")

dev.off()

#Plot Pathway-Specific Parts of Figure 2#

for(i in 1:length(Pathways))
{
  #Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.27/bin/gswin64c.exe")
  
  bitmap(paste(gsub("/","-",Pathways[i])," Heatmap.tiff",sep=""),width=4500,height=2000,units='px',res=720,colormodel="cmyk")
  
  Input=FigureTwoContents[[i]][[1]]
  colnames(Input)=c(Tumors,Lineages)
  heatmap.2(Input,margins=c(8.1,8.1),cexCol=1.62,colCol=c(rep("black",31),rep("indianred3",16)),cexRow=ifelse(i==5||i==9,0.45,0.675),tracecol=NA,density.info="none",dendrogram="none",Rowv=TRUE,Colv=FALSE,col=colorRampPalette(c("white","blue"))(n=299),lhei=c(1,3.6),lwid=c(1.2,4.5))
  
  dev.off()
  
  bitmap(paste(gsub("/","-",Pathways[i])," Network TCGA.tiff",sep=""),width=3000,height=3000,units='px',res=720,colormodel="cmyk")
  
  m=Types[[i]]
  p=graph_from_adjacency_matrix(m,mode="undirected",weighted=TRUE,diag=FALSE)
  
  m=PanCanGraphsTCGA[[i]]
  net=graph_from_adjacency_matrix(m,mode="undirected",weighted=TRUE,diag=FALSE)
  
  V(net)$color="black"
  V(net)$size=rowSums(m)/3.1
  E(net)$width=E(net)$weight/3.1
  E(net)$colors=c("gray","blue")[(E(net)$weight>cutoff*31)+1]
  E(net)$edge.lty=c("dashed","solid")[E(net)%in%E(p)+1]
  E(net)$arrow.size=0
  
  plot(net,layout=layout.star(net),vertex.label.color="white",vertex.label.cex=1.08,edge.label.color="gray30",edge.label=E(net)$weight,edge.color=E(net)$colors,edge.lty=E(net)$edge.lty,cex.main=27)
  
  dev.off()
  
  bitmap(paste(gsub("/","-",Pathways[i])," Network MCLP.tiff",sep=""),width=3000,height=3000,units='px',res=720,colormodel="cmyk")
  
  m=PanCanGraphsMCLP[[i]]
  net=graph_from_adjacency_matrix(m,mode="undirected",weighted=TRUE,diag=FALSE)
  
  V(net)$color="indianred3"
  V(net)$size=rowSums(m)/1.6
  E(net)$width=E(net)$weight/1.6
  E(net)$colors=c("gray","blue")[(E(net)$weight>cutoff*16)+1]
  E(net)$edge.lty=c("dashed","solid")[E(net)%in%E(p)+1]
  E(net)$arrow.size=0
  
  plot(net,layout=layout.star(net),vertex.label.color="black",vertex.label.cex=1.2,edge.label.color="gray30",edge.label=E(net)$weight,edge.color=E(net)$colors,edge.lty=E(net)$edge.lty,cex.main=27)
  
  dev.off()
}