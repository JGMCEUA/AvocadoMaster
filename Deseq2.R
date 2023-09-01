library(DESeq2)
###############################################################################
#                             T1                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1P")
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
countDataT1 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
coldata <- data.frame( row.names = colnames(countDataT1), condition)
head(coldata)
countDataT1 <- countDataT1[, rownames(coldata)]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(rownames(coldata) %in% colnames(countDataT1))
all(rownames(coldata)== colnames(countDataT1))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = countDataT1, colData = coldata, design = ~ condition)
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
###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange >= 1))
head(up_reg)
down_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange <=- 1))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj < 0.05))
write.table( DGE_T1, file = "DGE_T1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
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

##Para agregar etiquetas debemos activar la libreria ggplot2
library('ggplot2')
#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=2,aes(label = name), nudge_y = -0.5)

###MA plot
plotMA(dds,ylim=c(-2,2),main='MAplotT1')
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

###############################################################################
#                             T1                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1P2")
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
countDataT1 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
coldata <- data.frame( row.names = colnames(countDataT1), condition)
head(coldata)
countDataT1 <- countDataT1[, rownames(coldata)]
##Checar si todas las ID en colData está también en CountData e indica el orden
all(rownames(coldata) %in% colnames(countDataT1))
all(rownames(coldata)== colnames(countDataT1))
###Crea un DESeqDataSet desde la matriz y las etiquetas
dds <- DESeqDataSetFromMatrix(countData = countDataT1, colData = coldata, design = ~ condition)
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
###extraer el listado de los genes inducidos
up_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange >= 1))
head(up_reg)
down_reg <- row.names(subset(res, res$padj < 0.05 & res$log2FoldChange <=- 1))
head(down_reg)

DGE_T1 <- row.names(subset(res, res$padj < 0.05))
write.table( DGE_T1, file = "DGE_T1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
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

##Para agregar etiquetas debemos activar la libreria ggplot2
library('ggplot2')
#Vamos a crear un respaldo del PCA
z <- plotPCA(rld, intgroup=c('condition'))
##Regraficando pero con etiquetas sin marco y modificando
#el tamano de letra
z + geom_text(size=2,aes(label = name), nudge_y = -0.5)

###MA plot
plotMA(dds,ylim=c(-2,2),main='MAplotT1')
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

