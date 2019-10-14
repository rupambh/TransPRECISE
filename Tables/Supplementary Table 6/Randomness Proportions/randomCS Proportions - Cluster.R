# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

CSDist=function(i,j)
{
    library(PRECISE)
    load("RPPA.rda")
    CScores=0
    
    p=nrow(PathwayProteins[[j]])
    
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
    CScores=sum(netfit$G>0.3)/(p*(p-1))
    
    return(CScores)
}

load("RPPA.rda")

CSDists=matrix(0,length(Combined),length(Pathways))
Errors=matrix(0,length(Combined),length(Pathways))

for(i in 1:length(Combined))
{
    for(j in 1:length(Pathways))
    {
        Results=try(CSDist(i,j))
        
        if(isTRUE(class(Results)=="try-error"))
        {
            Errors[i,j]=1
            next
        }
        else
        {
            CSDists[i,j]=Results
        }
    }
}

save(CSDists,Errors,file=paste0("output",slurm_arrayid,".rda"))
