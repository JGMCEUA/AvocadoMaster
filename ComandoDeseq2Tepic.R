#Algunas Librerías utilizadas
library(DESeq2)
library(ggplot2)
library(scales)
library(viridis)
library(reshape2)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(vsn)
library(dplyr)
library(tximport)

###############################################################################
#
#                               T0 DGE
#
###############################################################################
#Cargar las direcciones de trabajo
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results","F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3")
#Importar los datos
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
#Determinar las condiciones del estudio
condition <- data.frame(condition = factor(c(rep("FC",3), rep("FF",3))))
#Filtro de contigs con conteo 0 o con longitud 0
##Cargar tablas de la lista txi.rsem
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
LengthData <- as.data.frame(txi.rsem$length, row.names = rownames(txi.rsem$length))
rownames(condition) = colnames(txi.rsem$counts)
CountData <- CountData[!(CountData$FC0_1 == 0 & CountData$FC0_2 == 0 & CountData$FC0_3 == 0 & CountData$FF0_1 == 0 & CountData$FF0_2 == 0 & CountData$FF0_3 == 0), ]
CountData$gene_id <- row.names(CountData)
LengthData <- LengthData[!(LengthData$FC0_1 == 0 | LengthData$FC0_2 == 0 | LengthData$FC0_3 == 0 | LengthData$FF0_1 == 0 | LengthData$FF0_2 == 0 | LengthData$FF0_3 == 0) ,]
LengthData$gene_id <- row.names(LengthData)
## Filtro realizado a la lista txi.rsem
DF <- inner_join(CountData,LengthData,"gene_id")
txi.rsem$abundance <- txi.rsem$abundance[DF$gene_id,]
txi.rsem$counts <- txi.rsem$counts[DF$gene_id,]
txi.rsem$length <- txi.rsem$length[DF$gene_id,]
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromTximport(txi.rsem,condition, ~ condition)
### Analisis Paramétrico de deseq, tomando la distribución local de las librerías
dds <- DESeq(dds , fitType='local')
#AnÃ¡lisis de expresion dierencial
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
res05 <- results(dds, alpha=0.05)
summary(res05)

################################################################################
###Filtrar genes que no obtuvieron valores significativos en el GED con distribución local
Nares <- na.omit(res05, padj)
Nares05 <- Nares[Nares$padj<=0.05,]
txi.rsem$abundance <- txi.rsem$abundance[rownames(Nares05), ]
txi.rsem$counts <- txi.rsem$counts[rownames(Nares05),]
txi.rsem$length <- txi.rsem$length[rownames(Nares05),]
### Realizar DGE con distribución paramétrica
dds <- DESeqDataSetFromTximport(txi.rsem,condition, ~ condition)
dds <- DESeq(dds)
####Resultados
res <- results(dds, c("condition","FF","FC"))
resOrdered <- res[order(res$padj),]
res05 <- results(dds, alpha=0.05)
summary(res05)

