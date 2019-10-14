CSDist=function(i,j,B)
{
    library(PRECISE)
    CScores=rep(0,B)
    load("RPPA.rda")
    
    p=nrow(PathwayProteins[[j]])
    
    for(b in 1:B)
    {
        #Pick Proteins#
        Sample=sample(colnames(Combined[[i]])[-(1:2)],p)
        RPPAdat=Combined[[i]][,intersect(colnames(Combined[[i]]),Sample)]
        Gmat=matrix(0.5,p,p)
        diag(Gmat)=0
        
        #Run PRECISE Model#
        dat=getregcovDat(Gmat=Gmat,RPPAdat=RPPAdat)
        bmsfit=getBMS(dat,Gmat)
        
        #Get C Score#
        nodes=names(dat$ylist)
        netfit=getPosteriors(bmsfit$outlist,nodes)
        CScores[b]=sum(netfit$G>0.3)/(p*(p-1)/2)
        
        cat("\014")
        dev.off()
    }
    
    return(CScores)
}

load("RPPA.rda")

CSDists=matrix(list(),length(Combined),length(Pathways))
Errors=matrix(0,length(Combined),length(Pathways))

for(i in 1:length(Combined))
{
    for(j in 1:length(Pathways))
    {
        Results=try(CSDist(i,j,10))
        
        if(isTRUE(class(Results)=="try-error"))
        {
            Errors[i,j]=1
            next
        }
        else
        {
            CSDists[i,j][[1]]=Results
        }
        
        print(c(i,j))
    }
}

save.image("CS Distributions.rda")
