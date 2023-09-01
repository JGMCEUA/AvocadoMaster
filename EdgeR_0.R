#source("/home/xocanol/Documents/Analisis_transcrip/salida_R/variacion_biologica/EdgeR_0.R") ubicacion del script

library(edgeR)

counts.all<-read.table(file="/home/xocanol/Documents/Analisis_transcrip/counts/variacion_bio.txt",header=T, sep="\t",row.names=1)
#head(counts.all) # ver datos

#groups
samples<-factor(c(rep(c("IA.0","IC.0","QIA.0","QIC.0"),2)))

#Analisis completo
y<- DGEList(counts.all,group=samples) #toda la matrix completa
y<- y[rowSums(y$counts) >= 3, ]
y<- calcNormFactors(y) #Normalizacion
#y$samples #ver grupos y norm.factores

#graficaMDS
pdf(file="/home/xocanol/Documents/Analisis_transcrip/salida_R/variacion_biologica/plotMDS_edeR_0.pdf")
plotMDS(y, main = "Distribución", xlim = c(-3, 3))
dev.off()


#design <- model.matrix(~0+ samples) #los nombre sampleIA.0...etc

#colnames(design)<- levels(samples) #cambiar el nombre (IA.0...etc)

#diseño
design <- model.matrix(~0+ samples)
y<- estimateGLMCommonDisp(y, design, verbose=TRUE)
y<-estimateGLMTrendedDisp(y, design)
y<- estimateGLMTagwiseDisp(y, design)

pdf(file="BCV_Desiccation.pdf")
plotBCV(y, main ="BCV")
dev.off()

colnames(design)<- levels(samples)
design


#design: a los tratamientos se les resta los testigos

contVector <- c("QvsnQ"="((QIC.0 + QIA.0) - (IC.0 + IA.0))",
"QIAvsIA"= "(QIA.0 - IA.0)", 
"QICvsIC0"= "(QIC.0 - IC.0)", 
"ICvsIA_0"= "(IC.0 - IA.0)"
)
contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)

#test
glmfit.all <- glmFit(y, design,dispersion=0.05) 



for (cont in colnames(contMatrix)){
print(cont)
lrt.all<-glmLRT(glmfit.all, contrast=contMatrix[,cont])
print(topTags(lrt.all))
print(summary(decideTestsDGE(lrt.all, p.value=0.05)))
top<-topTags(lrt.all, n=nrow(y$counts))
write.table (top$table[,c(1,5)], file=paste(cont, ".txt", sep=""),sep="\t",col.names=TRUE)
}



