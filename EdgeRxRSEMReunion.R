library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)

################################################################################
#                                 T0                                           #
################################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results","F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results","FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","FH1_3_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results", "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results","F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results","FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results","F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3", "FC1_1","FC1_2","FC1_3" ,"FF1_1","FF1_2","FF1_3", "FC24_1","FC24_2","FC24_3","FF24_1","F24_2","FF24_3", "FC72_1","FC72_2","FC72_3","FF72_1","FF72_2","FF72_3" )
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)

class <- c( rep("FC0",3), rep("FF0",3),rep("FC1",3), rep("FF1",3),rep("FC24",3), rep("FF24",3),rep("FC72",3), rep("FF72",3))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
y <- DGEList(counts = Mat, group = class, genes = Gene)
# Normalizacion TMM
y <- calcNormFactors(y)
#Person <- cor(y$counts,method = "pearson")
#write.csv(Person, file = "Pearson.csv")
##colors <- c("0","1")
##my_palette <- c("Red",colorRampPalette(colors = c("Red", "Blue"))(n = length(colors)), "Blue")
##pheatmap(Person, main = "Coeficiente de Correlaci?n de Pearson", cluster_rows = F, scale = "none" ,cluster_cols = F, color = my_palette, display_numbers = TRUE, number_color = "white", fontsize_number = 8 )

rainbow(9, start = 0, end = 1)
plotMDS(y$counts, top = 5000 ,xlab = "Dim1", ylab = "Dim2"  , main = "MDS conteos",col=c(rep("#FF0000",3),rep("#FFDB00",3),rep("Black",3),rep("#00FF40",3),rep("#8000FF",3),rep("#4900FF",3),rep("#FF00BF",3),rep("#E800FF",3)))

y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
#et <- exactTest(y)
#counts <- getCounts(y)
y$samples
et = exactTest(y, dispersion=y$tagwise.dispersion , pair=c("FC0","FF0")) #Primero COntrol despu?s Tratamiento.
head(et)

y2 <- topTags(et, n=Inf)$table
y2 = rownames(y2)[y2$FDR < 0.05 & abs(y2$logFC) >2]
length(y2)
et2 <- et
summary(de <- decideTestsDGE(et2, adjust.method="BH", FDR=0.05))

out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]

