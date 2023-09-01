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

###############################################################################
#
#                               T0 DGE
#
###############################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results","F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)

condition <- data.frame(condition = factor(c(rep("FC",3), rep("FF",3))))
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
rownames(condition) = colnames(CountData)
coldata <- data.frame( colnames(CountData),condition)
colnames(coldata) <- c("id","condition")
CountData <- CountData[, coldata$id]
CountData2 <- ceiling(CountData)
CountData2 <- CountData2[CountData2$FC0_1 != 0 & CountData2$FC0_2 != 0 & CountData2$FC0_3 != 0 & CountData2$FF0_1 != 0 & CountData2$FF0_2 != 0 & CountData2$FF0_3 != 0 ,]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(coldata$idn %in% colnames(CountData2))
all(coldata$idn == colnames(CountData2))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = CountData2, colData = condition, design = ~ condition)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
#Análisis de expresión
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "GEDT0.csv")
# Gráfica de dispersión
png("qc-dispersionsT0.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T0")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj <= 0.05, na.rm=TRUE)
#####
#Imprimir resultados del segundo filtro
#####
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
#####
# crear un MA plot de los resultados 
#####
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT0.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05T0.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
ddsVST <- vst(dds)
ddsVST
#####
#Conversi?n de los objetos DEseq tranformados a data.frame
#####
ddsVST <- assay(ddsVST)
ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
#####
#Filtro de los genes significativos y con cierto valor de LFC
#####
sigGenes <- rownames(res05DF[res05DF$padj <= 0.05 & abs(res05DF$log2FoldChange) > 2,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
#####
#Comparaci?n de datos en ancho y largo
#####
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
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
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
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
sampleData_v2 <- coldata
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B")
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

###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange >= 0))
head(up_reg)
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))
head(down_reg)

DGE_T0 <- row.names(subset(res, res$padj <= 0.05))
write.table( DGE_T0, file = "DGE_T0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_regT0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_regT0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3, colors= "black")

###############################################################################
#
#                               T1 DGE
#
###############################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","FH1_3_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FC1_3","FF1_1","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
condition <- data.frame(condition = factor(c(rep("FC",3), rep("FF",3))))
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
rownames(condition) = colnames(CountData)
coldata <- data.frame( colnames(CountData),condition)
colnames(coldata) <- c("id","condition")
CountData <- CountData[, coldata$id]
CountData2 <- ceiling(CountData)
CountData2 <- CountData2[CountData2$FC1_1 != 0 & CountData2$FC1_2 != 0 & CountData2$FC1_3 != 0 & CountData2$FF1_1 != 0 & CountData2$FF1_2 != 0 & CountData2$FF1_3 != 0 ,]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(coldata$idn %in% colnames(CountData2))
all(coldata$idn == colnames(CountData2))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = CountData2, colData = condition, design = ~ condition)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
#Análisis de expresión
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "GEDT1.csv")
# Gráfica de dispersión
png("qc-dispersionsT1.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T1")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj <= 0.05, na.rm=TRUE)
#####
#Imprimir resultados del segundo filtro
#####
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
#####
# crear un MA plot de los resultados 
#####
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT1.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05T1.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
ddsVST <- vst(dds)
#####
#Conversi?n de los objetos DEseq tranformados a data.frame
#####
ddsVST <- assay(ddsVST)
ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
#####
#Filtro de los genes significativos y con cierto valor de LFC
#####
sigGenes <- rownames(res05DF[res05DF$padj <= 0.05 & abs(res05DF$log2FoldChange) > 2,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
#####
#Comparaci?n de datos en ancho y largo
#####
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
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
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
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
sampleData_v2 <- coldata
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B")
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

###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange >= 0))
head(up_reg)
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj <= 0.05))
write.table( DGE_T1, file = "DGE_T1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_regT1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_regT1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3, colors= "black")


###############################################################################
#
#                               T1-1 DGE
#
###############################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FF1_1","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
condition <- data.frame(condition = factor(c(rep("FC",2), rep("FF",3))))
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
rownames(condition) = colnames(CountData)
coldata <- data.frame( colnames(CountData),condition)
colnames(coldata) <- c("id","condition")
CountData <- CountData[, coldata$id]
CountData2 <- ceiling(CountData)
CountData2 <- CountData2[CountData2$FC1_1 != 0 & CountData2$FC1_2 != 0 & CountData2$FF1_1 != 0 & CountData2$FF1_2 != 0 & CountData2$FF1_3 != 0 ,]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(coldata$idn %in% colnames(CountData2))
all(coldata$idn == colnames(CountData2))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = CountData2, colData = condition, design = ~ condition)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
#Análisis de expresión
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "GEDT1P.csv")
# Gráfica de dispersión
png("qc-dispersionsT1P.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T1")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj <= 0.05, na.rm=TRUE)
#####
#Imprimir resultados del segundo filtro
#####
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
#####
# crear un MA plot de los resultados 
#####
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT1P.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05T1P.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
ddsVST <- vst(dds)
#####
#Conversi?n de los objetos DEseq tranformados a data.frame
#####
ddsVST <- assay(ddsVST)
ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
#####
#Filtro de los genes significativos y con cierto valor de LFC
#####
sigGenes <- rownames(res05DF[res05DF$padj <= 0.05 & abs(res05DF$log2FoldChange) > 2,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
#####
#Comparaci?n de datos en ancho y largo
#####
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
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
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
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
sampleData_v2 <- coldata
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B")
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

###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange >= 0))
head(up_reg)
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj <= 0.05))
write.table( DGE_T1, file = "DGE_T1P.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_regT1P.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_regT1P.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3, colors= "black")


###############################################################################
#
#                               T1-2 DGE
#
###############################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
condition <- data.frame(condition = factor(c(rep("FC",2), rep("FF",2))))
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
rownames(condition) = colnames(CountData)
coldata <- data.frame( colnames(CountData),condition)
colnames(coldata) <- c("id","condition")
CountData <- CountData[, coldata$id]
CountData2 <- ceiling(CountData)
CountData2 <- CountData2[CountData2$FC1_1 != 0 & CountData2$FC1_2 != 0 & CountData2$FF1_2 != 0 & CountData2$FF1_3 != 0 ,]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(coldata$idn %in% colnames(CountData2))
all(coldata$idn == colnames(CountData2))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = CountData2, colData = condition, design = ~ condition)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
#Análisis de expresión
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "GEDT1PP.csv")
# Gráfica de dispersión
png("qc-dispersionsT1P.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T1PP")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj <= 0.05, na.rm=TRUE)
#####
#Imprimir resultados del segundo filtro
#####
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
#####
# crear un MA plot de los resultados 
#####
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT1PP.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05T1PP.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
ddsVST <- vst(dds)
#####
#Conversi?n de los objetos DEseq tranformados a data.frame
#####
ddsVST <- assay(ddsVST)
ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
#####
#Filtro de los genes significativos y con cierto valor de LFC
#####
sigGenes <- rownames(res05DF[res05DF$padj <= 0.05 & abs(res05DF$log2FoldChange) > 2,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
#####
#Comparaci?n de datos en ancho y largo
#####
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
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
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
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
sampleData_v2 <- coldata
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B")
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

###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange >= 0))
head(up_reg)
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj <= 0.05))
write.table( DGE_T1, file = "DGE_T1PP.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_regT1PP.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_regT1PP.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3, colors= "black")

###############################################################################
#
#                               T24 DGE
#
###############################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results","F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC24_1","FC24_2","FC24_3","FF24_1","FF24_2","FF24_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
class <- c( rep("FC",3), rep("FF",3))
condition <- data.frame(condition = factor(c(rep("FC",3), rep("FF",3))))
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
rownames(condition) = colnames(CountData)
coldata <- data.frame( colnames(CountData),condition)
colnames(coldata) <- c("id","condition")
CountData <- CountData[, coldata$id]
CountData2 <- ceiling(CountData)
CountData2 <- CountData2[CountData2$FC24_1 != 0 & CountData2$FC24_2 != 0 & CountData2$FC24_3 != 0 & CountData2$FF24_1 != 0 & CountData2$FF24_2 != 0 & CountData2$FF24_3 != 0 ,]
##Checar si todas las ID en colData está también en CountData e indica el orden
#all(coldata$idn %in% colnames(CountData2))
#all(coldata$idn == colnames(CountData2))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = CountData2, colData = condition, design = ~ condition)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
#Análisis de expresión
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
resOrdered <- res[order(res$padj),]
write.csv(resOrdered,file = "GEDT24.csv")
# Gráfica de dispersión
png("qc-dispersionsT24.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T24")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj <= 0.05, na.rm=TRUE)
#####
#Imprimir resultados del segundo filtro
#####
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
#####
# crear un MA plot de los resultados 
#####
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT24.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05FCT24.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
ddsVST <- vst(dds)
ddsVST
#####
#Conversi?n de los objetos DEseq tranformados a data.frame
#####
ddsVST <- assay(ddsVST)
ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
#####
#Filtro de los genes significativos y con cierto valor de LFC
#####
sigGenes <- rownames(res05DF[res05DF$padj <= 0.05 & abs(res05DF$log2FoldChange) > 2,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
#####
#Comparaci?n de datos en ancho y largo
#####
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
#####
#Ajustar los datos de acuerdo al formato long
#####
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
################################################################################
##Agrupar resultados del HeatMap
################################################################################
#####
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
#####
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
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
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
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
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
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
sampleData_v2 <- coldata
#####
# Re ordenar el dise?o experimental
#####
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])
#####
# Construir una gr?fica para mostrar los tratamientos
#####
colours <- c("#743B8B", "#8B743B")
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

###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange >= 0))
head(up_reg)
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))
head(down_reg)

DGE_T24 <- row.names(subset(res, res$padj <= 0.05))
write.table( DGE_T24, file = "DGE_T24.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_regT24.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_regT24.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"

###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3, colors= "black")

