#metiendo datos
counts.all<-read.table(file="/home/sandra/Documents/BGI/statistics_tomato/RPKMstomatoRNaseq.txt",header=T, sep="\t",row.names=1)

#Factores
stages<-factor(c("fdev1","fdev1","fdev2","fdev2","fdev3","fdev3","fMG","fMG","fB","fB","fB10","fB10"))
d.all<- DGEList(counts.all,group=stages)#toda la matrix completa


pdf(file="edgeR_tomatoMDSplot.pdf")
plotMDS(d.all, main = "Tomato Maturation RNAseq", xlim = c(-8, 8), ylim=c(-4,4))
dev.off()### checar para ejemplo del manual

design <- model.matrix(~0+ stages)
#d.all <- estimateGLMCommonDisp(d.all, design, verbose=TRUE) 
#d.all<- estimateGLMTagwiseDisp(d.all, design)
colnames(design)<- levels(stages)
design

#haciendo contrastes
contVector <- c(
"fdev2vsfdev1"="fdev2-fdev1",#contro fruto pequeño
"fdev3vsfdev1"="fdev3-fdev1",
"fMGvsfdev1"="fMG-fdev1",
"fBvsfdev1"="fB-fdev1",
"fB10vsfdev1"="fB10-fdev1",
"fBvsfMG"="fB-fMG",
"fB10vsfMG"="fB10-fMG")

contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)


#Haciendo el test
glmfit.all <- glmFit(d.all, design,dispersion=0.1)

for (cont in colnames(contMatrix)){
print(cont)
lrt.all<-glmLRT(glmfit.all, contrast=contMatrix[,cont])
print(topTags(lrt.all))
print(summary(decideTestsDGE(lrt.all, p.value=0.05)))
top<-topTags(lrt.all, n=nrow(d.all$counts))
write.table(top$table[,c(1,5)], file=paste(cont, ".txt", sep=""),sep="\t",col.names=TRUE)
}



