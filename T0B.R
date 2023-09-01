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
#head(txi.rsem$counts)
#rownames(txi.rsem$counts)
class <- c( rep("FC",3), rep("FF",3))
condition <- data.frame(condition = factor(c(rep("FC",3), rep("FF",3)))) #UNO; Sampletable
CountData <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
LengthData <- as.data.frame(txi.rsem$length, row.names = rownames(txi.rsem$length))
rownames(condition) = colnames(txi.rsem$counts) #DOS; Sample table rownames
CountData <- CountData[!(CountData$FC0_1 == 0 & CountData$FC0_2 == 0 & CountData$FC0_3 == 0 & CountData$FF0_1 == 0 & CountData$FF0_2 == 0 & CountData$FF0_3 == 0), ]
CountData2 <- CountData
CountData$gene_id <- row.names(CountData)
LengthData <- LengthData[!(LengthData$FC0_1 == 0 | LengthData$FC0_2 == 0 | LengthData$FC0_3 == 0 | LengthData$FF0_1 == 0 | LengthData$FF0_2 == 0 | LengthData$FF0_3 == 0) ,]
LengthData$gene_id <- row.names(LengthData)
DF <- inner_join(CountData,LengthData,"gene_id")
txi.rsem$abundance <- txi.rsem$abundance[DF$gene_id,]
txi.rsem$counts <- txi.rsem$counts[DF$gene_id,]
txi.rsem$length <- txi.rsem$length[DF$gene_id,]
###Crea un DESeqDataSet desde la matriz y las etiquetas
#write.csv(txi.rsem$counts, file = "Mambo5.csv")
dds <- DESeqDataSetFromTximport(txi.rsem,condition, ~ condition)
keep <- rowSums(counts(dds) >= 5) >= 3
###Análisis default de Deseq2 y genera una tabla de resultados
dds <- DESeq(dds[keep,])

#Análisis de expresión
png("qc-dispersionsT0BASD.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot time 0 h Prefilt")
dev.off()
res <- results(dds, c("condition","FF","FC"))
###Acomodar por su p-value y muestra
resOrdered <- res[order(res$padj),]
write.csv(resOrdered,file = "GEDT0B.csv")
res05 <- results(dds, alpha=0.05)
summary(res05)
print("Total")
sum(res05$padj <= 0.05, na.rm=TRUE)
print("Total")
sum(res05$padj <= 0.05 & abs(res05$log2FoldChange) >=2, na.rm=TRUE)
print("UP")
sum(res05$padj <= 0.05 & res05$log2FoldChange >=2, na.rm=TRUE)
print("DN")
sum(res05$padj <= 0.05 & res05$log2FoldChange <=-2, na.rm=TRUE)
####
# crear un MA plot de los resultados 
#####
res05DF <- as.data.frame(res)
res05$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
res05DF$significant <- ifelse(res05$padj <= 0.05, "Significant", NA)
Holakhace <- subset(res05, subset= res05$significant=="Significant")
write.csv(Holakhace, file = "res05DFT0.csv")
Momo <- subset(res05, subset= res05$significant=="Significant" & abs(res05$log2FoldChange)>=2)
write.csv(Momo, file = "res05FCT0.csv", row.names = F)

ggplot(res05DF, aes(baseMean, log2FoldChange, colour=padj)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="Media de conteos normalizados", y="Log Fold Change") + scale_colour_viridis(direction=-1, trans='sqrt') 

################################################################################
#                                     HeatMap
################################################################################
#####
#Tranformaci?n de los datos de conte? con vst (variance satabilizing transform)
#####
res3D <- na.omit(res05,padj)
dds <- DESeq(dds[rownames(res3D[res3D$padj<=0.05,])],)
ddsVST <- vst(dds, nsub = sum( rowMeans( counts(dds, normalized=TRUE)) > 5 ))
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
sampleData_v2 <- data.frame( colnames(CountData2),condition)
colnames(sampleData_v2) <- c("id","condition")
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
down_reg <- row.names(subset(res05, res05$padj <= 0.05 & res$log2FoldChange <=- 0))

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


###PCA
rld<- rlogTransformation(dds, blind=T)
plotPCA(rld, intgroup=c('condition'))

#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=6,aes(label = name), nudge_y = -3)
