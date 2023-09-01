#source("/home/xocanol/Desktop/Diario/Diario_EdgeR/EdgeR_5.R") ubicacion del script

library(edgeR)

counts.all<-read.table(file="/home/xocanol/samples/conteo/muestras_todas.txt",header=T, sep="\t",row.names=1)
#head(counts.all) # ver datos

#groups
samples<-factor(c("IA.0","IA.1","IA.69","IA.24","IC.0","IC.1","IC.69","IC.24","QIA.0","QIA.1","QIA.69","QIA.24","QIC.0","QIC.1","QIC.69","QIC.24"))

#Analisis completo
y<- DGEList(counts.all,group=samples) #toda la matrix completa
y<- y[rowSums(y$counts) >= 3, ]
y<- calcNormFactors(y) #Normalizacion
#y$samples #ver grupos y norm.factores

#graficaMDS
pdf(file="/home/xocanol/samples/arch_R/plotMDS_edeR_5.pdf")
plotMDS(y, main = "Distribución", xlim = c(-2, 2))
dev.off()


design <- model.matrix(~0+ samples) #los nombre sampleIA.0...etc

colnames(design)<- levels(samples) #cambiar el nombre (IA.0...etc)

#design: a los tratamientos se les resta los testigos

contVector <- c("QvsnQ1"="(QIC.24+QIC.0+QIC.1+QIC.69)+(QIA.24+QIA.0+QIA.1+QIA.69)-(IC.24+IC.0+IC.1+IC.69)-(IA.24+IA.0+IA.1+IA.69)",
#"QvsnQ2"="(QIA.0 - IA.0) + (QIA.1 - IA.1) + (QIA.69 - IA.69) + (QIA.24 - IA.24) + (QIC.0 - IC.0) + (QIC.1 - IC.1) + (QIC.69 - IC.69) + (QIC.24 - IC.24) + (IC.1 - IA.0) + (IC.69 - IA.0) + (IC.24 - IA.0)"
"QIAvsIA_0"= "(QIA.0 - IA.0)", #COMPARACIONES POR PARES QUITOSANO NO INOCULADO VS TESTIGO
"QIAvsIA_69"= "(QIA.69 - IA.69)",
"QIAvsIA_24"= "(QIA.24 - IA.24)",
"QIAvsIA_1"= "(QIA.1 - IA.1)",
"QICvsIC_0"= "(QIC.0 - IC.0)", #COMPARACIONES POR PARES QUITOSANO, INOCULADO VS TESTIGO
"QICvsIC_69"= "(QIC.69 - IC.69)",
"QICvsIC_24"= "(QIC.24 - IC.24)",
"QICvsIC_1"= "(QIC.1 - IC.1)",
"ICvsIA_0"= "(IC.0 - IA.0)", # COMPARACIONES POR PARES NO QUITOSANO, INOCULADO VS TESTIGO
"ICvsIA_1"= "(IC.1 - IA.1)",
"ICvsIA_69"= "(IC.69 - IA.69)",
"ICvsIA_24"= "(IC.24 - IA.24)"
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



