library(DESeq2)
library(ggplot2)
library(scales)
library(viridis)
library(reshape2)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RTtotal")
condition <- factor(c(rep("FC",12), rep("FF",12)))
time <- factor(c(rep(0,3),rep(1,3),rep(24,3),rep(72,3),rep(0,3),rep(1,3),rep(24,3),rep(72,3)))
condition = relevel(condition, ref = "FC")
time = relevel(time, ref = "0")
countDataT1 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
coldata <- data.frame( row.names = colnames(countDataT1), condition,  time)
head(coldata)
coldata
countDataT1 <- countDataT1[, rownames(coldata)]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(rownames(coldata) %in% colnames(countDataT1))
all(rownames(coldata)== colnames(countDataT1))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = countDataT1, colData = coldata, design = ~ condition + time)
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds)
head(dds)
#Análisis de expresión
res <- results(dds)
###Acomodar por su p-value y muestra
(resOrdered <- res[order(res$padj),])
write.csv(resOrdered,file = "PruebaT1.csv")
# Gráfica de dispersión
png("qc-dispersionsT1.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot_T1")
dev.off()
####ajustar los resultados de expresi?n a un pvalue de 0.5
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
### crear un MA plot de los resultados 
plotMA(res05, ylim=c(-5,5))
res05DF <- as.data.frame(res)
head(res05DF)
res05DF$significant <- ifelse(res05DF$padj < .05, "Significant", NA)
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

##Agrupar resultados del HeatMap
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
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) +
  geom_segment(aes(x=x, y=y, xend= xend, yend=yend)) +
  theme_dendro()
#Reasignar factores para ggplot2
ddsVST$variable <- factor(ddsVST$variable,
                          levels = clusterSample$labels[clusterSample$order])
heatmap <- ggplot(ddsVST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
heatmap
#Combina el dendogram y el Heatmap
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights = c(1,5))

##Ajustar datos
#Modificar el objeto ggplot
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand = c(.0085, 0.0085)) +
  scale_y_continuous(expand = c(0,0))
heatmap_1 <- heatmap + scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
#Convertir ambos objetos grid a grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)
#checar el ancho de cada grob
sampleDendrogramGrob$widths
heatmapGrob$widths
#A�adir en las columnas perdidas
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob,
                                        heatmapGrob$widths[7],6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob,
                                        heatmapGrob$widths[8],7)
#Garantizar que el ancho de los dos grob es el mismo
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
#Promedio de grobs en la gr�fica
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol = 1,
                         heights = c(2,5))
#Dibujar la gr�fica
grid.draw(finalGrob)




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
z + geom_text(size=2,aes(label = name), nudge_y = -0.5)

###MA plot
plotMA(dds,ylim=c(-15,10),main='MAplotT1')
dev.copy(png,'deseq2_MAplotT1.png')
dev.off() 

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