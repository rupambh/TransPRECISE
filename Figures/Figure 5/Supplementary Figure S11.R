C4.Pts=c(80,38,22,11)
C4.CLs=c(81,72,48,35,34,33,20)
C4=round((C4.CLs%*%t(C4.Pts))/100)
rownames(C4)=c("ovary","head-neck","skin","lung","kidney","breast","pancreas")
colnames(C4)=c("PAAD","HNSC","BLCA","OV")

write.table(C4,file=paste0("C4.txt"),sep="\t",row.names=TRUE,col.names=TRUE)

C15.CLs=38
C15.Pts=c(60,48,40,30,19,19)
C15=round((C15.CLs%*%t(C15.Pts))/100)
colnames(C15)=c("HNSC","CESC","ESCA","SARC","KIRP","GBM")
rownames(C15)=c("stomach-oesophagus")

write.table(C15,file=paste0("C15.txt"),sep="\t",row.names=TRUE,col.names=TRUE)