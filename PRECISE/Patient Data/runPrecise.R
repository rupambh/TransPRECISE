runPrecise=function(i,j)
{
  #Load Package#
  library(PRECISE)
  
  #Load Environment for Data#
  load("Input.rda")
  
  #Pathway Headers#
  pw.array=Pathways
  
  #Load Data : Tumor x Pathway#
  RPPAdat=Final[i,j][[1]]
  
  Details=RPPAdat[,1:2]
  RPPAdat=as.matrix(RPPAdat[,-(1:2)])
  p=ncol(RPPAdat)
  
  #Uniform Prior#
  Gmat=matrix(0.5,p,p)
  diag(Gmat)=0
  
  #Final Response Vector + Covariate Matrix#
  dat=getregcovDat(Gmat=Gmat,RPPAdat=RPPAdat)
  
  #Fit Model Now#
  bmsfit=getBMS(dat,Gmat)
  
  #Get Tumor-Specific Pathway Network#
  nodes=names(dat$ylist)
  netfit=getPosteriors(bmsfit$outlist,nodes)
  Network=netfit$G
  
  #Get Betas for Tumor-Specific Pathway Network#
  Coeffs=Network
  
  for(i in 1:nrow(Coeffs))
  {
    Coeffs[i,]=append(bmsfit$outlist[[i]][[1]]$b1mo,0,after=i-1)
  }
  
  #Create Patient-Specific Networks#
  delta=0.5
  psNet=getPRECISE(bmsfit$outlist,bmsfit$pdlist,nodes,delta)
  Scores=cbind(Details,psNet$score.mat)
  
  #Calibrate PRECISE Scores#
  Status=cbind(Details,apply(psNet$score.mat,1,which.max))
  
  #Produce Results#
  Results=list(Scores,Status,Network,Coeffs)
  names(Results)=c("Scores","Status","Network","Coeffs")
  
  return(Results)
}

load("Begin.rda")

Betas=matrix(list(),length(Tumors),length(Pathways))
Scores=matrix(list(),length(Tumors),length(Pathways))
Status=matrix(list(),length(Tumors),length(Pathways))
Graphs=matrix(list(),length(Tumors),length(Pathways))
Errors=matrix(0,length(Tumors),length(Pathways))

for(i in 1:length(Tumors))
{
  for(j in 1:length(Pathways))
  {
    Results=try(runPrecise(i,j))
    
    if(isTRUE(class(Results)=="try-error"))
    {
      Errors[i,j]=1
      next
    }
    else
    {
      Scores[i,j][[1]]=Results[[1]]
      Status[i,j][[1]]=Results[[2]]
      Graphs[i,j][[1]]=Results[[3]]
      Betas[i,j][[1]]=Results[[4]]
    }
    
    dev.off()
    cat("\014")
    print(c(i,j))
  }
}

save.image("Results.rda")