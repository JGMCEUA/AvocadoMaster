library(DESeq2)
library(ggplot2)
library(scales)
library(viridis)
library(reshape2)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
###############################################################################
#
#                               TTOTAL DGE
#
###############################################################################
#setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RTtotal")
setwd("E:/R/DEG/RTtotal")
condition <- factor(c( rep("FC0",3),rep("FC1",3),rep("FC24",3),rep("FC72",3),rep("FF0",3),rep("FF1",3),rep("FF24",3),rep("FF72",3)))
time <- factor(c(rep("0",3),rep("1",3),rep("24",3),rep("72",3),rep("0",3),rep("1",3),rep("24",3),rep("72",3)))
condition = relevel(condition, ref = "FC0")
countDataT1 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))

#coldata <- data.frame( row.names = colnames(countDataT1),condition)
coldata <- data.frame( colnames(countDataT1),condition, time)
colnames(coldata) <- c("id","condition","time")
head(coldata)
countDataT1 <- countDataT1[, coldata$id]
##Checar si todas las ID en colData estÃ¡ tambiÃ©n en CountData e indica el orden
all(coldata$idn %in% colnames(countDataT1))
all(coldata$idn == colnames(countDataT1))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = countDataT1, colData = coldata, design = ~ condition)
###AnÃ¡lisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
head(dds)
#AnÃ¡lisis de expresiÃ³n
res <- results(dds)
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "GED.csv")
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
### crear un MA plot de los resultados 
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05DF$significant <- ifelse(res05DF$padj < 0.05, "Significant", NA)
Holakhace <- subset(res05DF, subset= res05DF$significant=="Significant")
write.csv(Holakhace, file = "res05DF.csv")

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=significant))+
  geom_point(size=1) + scale_y_continuous(limits = c(-15,10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) +
  labs(x="medias de conteoss normalizados", y= "Log Fold Change") +
  scale_colour_manual(name="q-value", values = ("Significant"="red"),
                      na.value = "grey50") + theme_bw()

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-15, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

sampleData_v2 <- coldata

###HeatMap
ddsVST <- vst(dds)
ddsVST <- assay(ddsVST)
ddsVST<- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
head(ddsVST)

sigGenes <- rownames(res05DF[res05DF$padj <= .05 &
                               abs(res05DF$log2FoldChange) > 3,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]
ddsVST_wide <- ddsVST
ddsVST_long <- melt(ddsVST, id.vars = c("Gene"))
head(ddsVST_wide)
head(ddsVST_long)
ddsVST <- melt(ddsVST, id.vars = c("Gene"))
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
heatmap
###############################################################################
##Agrupar resultados del HeatMap
###############################################################################
#Convertir los datos de los genes significantes a una matriz para el agrupamiento
ddsVSTMatrix <- dcast(ddsVST, Gene ~ variable)
rownames(ddsVSTMatrix) <- ddsVSTMatrix$Gene
ddsVSTMatrix$Gene <- NULL
#Calcular distancia en ambas dimensiones de la matriz
distanceGene <- dist(ddsVSTMatrix)
distanceSample <- dist(t(ddsVSTMatrix))
#Agrupar de acuerdo a los calculos de la distancia
ClusterGene <- hclust(distanceGene, method = "average")
clusterSample <- hclust(distanceSample, method = "average")
#Construir un dendograma para las muestras
sampleModel <- as.dendrogram(clusterSample)
geneModel <- as.dendrogram(ClusterGene)

sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x=x, y=y, xend= xend, yend=yend)) + theme_dendro()
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()
#Reasignar factores para ggplot2
ddsVST$variable <- factor(ddsVST$variable, levels = clusterSample$labels[clusterSample$order])
ddsVST$Gene <- factor(ddsVST$Gene, levels = ClusterGene$labels[ClusterGene$order])

heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) +  geom_raster() + scale_fill_viridis(trans="sqrt") +  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
#heatmap
#Combina el dendogram y el Heatmap
#grid.arrange(sampleDendrogram, heatmap, ncol=1, heights = c(1,5))
###############################################################################
##Ajustar datos
#Modificar el objeto ggplot
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand = c(.0085, 0.0085)) + scale_y_continuous(expand = c(0,0))
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))
heatmap_1 <- heatmap + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))


#Convertir ambos objetos grid a grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
geneDendrogramGrob <- ggplotGrob(geneDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)
#checar el ancho de cada grob
#sampleDendrogramGrob$widths
#heatmapGrob$widths
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)

#Añadir en las columnas perdidas
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7],6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8],7)
#Garantizar que el ancho de los dos grob es el mismo
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

###############################################################################
#                     CARGAR DIFERENCIA DE MUESTRA
###############################################################################
# Re-order the sample data to match the clustering we did
sampleData_v2$id <- factor(sampleData_v2$id, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B")
sampleClinical <- ggplot(sampleData_v2, aes(x=id, y=1, fill= condition)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tratamiento", values=colours) + theme_void()

# Convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

###############################################################################

sampleDendrogramGrob$widths
heatmapGrob$widths
sampleClinicalGrob$widths


blankPanel <- grid.rect(gp=gpar(col="white"))
# arrange all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))

# draw the final result
grid.draw(finalGrob_v2)


#Promedio de grobs en la gráfica
#finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol = 1, heights = c(2,5))
#Dibujar la gráfica
#grid.draw(finalGrob)




###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange >= 1))
head(up_reg)
down_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange <=- 1))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj < 0.05))
write.table( DGE_T1, file = "DGE_T0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
#guardar tablas de ID
write.table( up_reg, file = "up_reg_T1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table( down_reg, file = "down_reg_T1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
##crear un data frame para realizar un volcano plot que contenga todos los datos de genes diferenciales
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


###PCA
rld<- rlogTransformation(dds, blind=TRUE)
plotPCA(rld, intgroup=c('condition'))
#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=2,aes(label = name), nudge_y = -0.5, colors= "black")

## Volcano plot with "significant" genes labeled
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot_FHvsF_T1", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot_T1.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()


###comando para crear un volcano con paqueteria enhanced
library(EnhancedVolcano)
EnhancedVolcano(res05,
                lab = rownames(res05) ,
                x = 'log2FoldChange',
                y = 'pvalue')



###PCA
rld<- rlogTransformation(dds, blind=TRUE)
plotPCA(rld, intgroup=c('condition'))


#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=3,aes(label = name), nudge_y = -0.5, colour = "black")
