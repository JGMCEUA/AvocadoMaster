library(readr)
library(tidyverse)
library(dslabs)
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

setwd("F:/ETAPA 2/DEG/AnalisisR")
hola <- dist

T0 <- as.data.frame(read.csv("res05UD-T0.csv"))
T02 <- T0[abs( T0$logFC) > 2,]
#T1 <- as.data.frame(read.csv("res05UD-T1.csv"))
T24 <- as.data.frame(read.csv("res05UD-T24.csv"))
T242 <- T0[abs( T0$logFC) > 2,]
T72 <- as.data.frame(read.csv("res05UD-T72.csv"))
T722 <- T0[abs( T0$logFC) > 2,]

ED1 <- full_join(T02,T242, "genes")
ED2 <- full_join(ED1,T722, "genes")
write.csv(ED2,"ConteosUD2-OK.csv") ##Modificar archivo desde excel para que solo contenga los resutlados del LogFC primera columna debe ser Gene
Res <- as.data.frame(read.csv("ConteosUD2-OK.csv"))
condition <- data.frame(condition = factor(c("1","24","72")))
################################################################################
#                                     HeatMap
################################################################################
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(Res, id.vars = c("genes"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, genes ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$genes
ddsVSTMatrix$genes <- NULL
#####
#Calcular distancia en ambas dimensiones de la matriz
#####
distanceGene <- dist(ddsVSTMatrix)
distanceSample <- dist(t(ddsVSTMatrix))
#####
#Agrupar de acuerdo a los calculos de la distancia
#####
clusterGene <- hclust(distanceGene, method = "average")
clusterSample <- hclust(distanceSample, method = "average")
#Construir un dendograma para las muestras
geneModel <- as.dendrogram(clusterGene)
sampleModel <- as.dendrogram(clusterSample)

geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))

sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x=x, y=y, xend= xend, yend=yend)) + theme_dendro()
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()
#####
#Reasignar factores para ggplot2
#####
ddsVST$variable <- factor(ddsVST$variable, levels = clusterSample$labels[clusterSample$order])
ddsVST$Gene <- factor(ddsVST$Gene, levels = clusterGene$labels[clusterGene$order])
#####
#Modificar el objeto ggplot
#####
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster()  + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand = c(.0085, 0.0085)) + scale_y_continuous(expand = c(0,0))
heatmap_1 <- heatmap + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))
#####
#Convertir ambos objetos grid a grobs
#####
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))
heatmapGrob <- ggplotGrob(heatmap_1)
#####
#checar el ancho de cada grob
#####
#sampleDendrogramGrob$widths
#heatmapGrob$widths
#length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
#####
#A?adir en las columnas perdidas
#####
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7],6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8],7)
#####
#Garantizar que el ancho de los dos grob es el mismo
#####
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

###############################################################################
#                     CARGAR DIFERENCIA DE MUESTRA
###############################################################################
#####
#cargar dise?o experimental
#####
sampleData_v2 <- data.frame( colnames(Res[-1]),condition)

colnames(sampleData_v2) <- c("id","condition")
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B","Red","Blue")
sampleClinical <- ggplot(sampleData_v2, aes(x=id, y=1, fill= condition)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tratamiento", values=colours) + theme_void()
#####
# Convertir la gr?fica del dise?o experimental a un grob
#####
sampleClinicalGrob <- ggplotGrob(sampleClinical)
#####
# Asegurarse que el ancho de todos los grobs son los mismos
#####
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)
################################################################################
# Concatenaci?n Final
################################################################################
blankPanel <- grid.rect(gp=gpar(col="white"))
#####
# Concatenar todas las gr?ficas juntas
#####
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))
#####
# Dibujar gr?fica final
#####
grid.draw(finalGrob_v2)
################################################################################

