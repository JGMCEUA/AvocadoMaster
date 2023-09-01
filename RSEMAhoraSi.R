library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)

#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("G:/AgroBio_4TB/ETAPA 2/DEG/AnalisisR")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results",
         "F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results",
         "FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","FH1_3_RHC.genes.results",
         "F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results", 
         "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results",
         "F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results",
         "FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results",
         "F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")

run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results",
         "F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results",
         "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results",
         "F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results",
         "FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results",
         "F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")

files <- file.path(dir,run)
#names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3", "FC1_1","FC1_2","FF1_1","FF1_2","FF1_3", "FC24_1","FC24_2","FC24_3","FF24_1","FF24_2","FF24_3", "FC72_1","FC72_2","FC72_3","FF72_1","FF72_2","FF72_3" )
names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3", "FC24_1","FC24_2","FC24_3","FF24_1","FF24_2","FF24_3", "FC72_1","FC72_2","FC72_3","FF72_1","FF72_2","FF72_3" )
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
#write.csv(txi.rsem$counts,"Conteos.csv")
class <- c( rep("FC0",3), rep("FF0",3),rep("FC24",3), rep("FF24",3),rep("FC72",3), rep("FF72",3))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
#______________________________________________________________________________#
# run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results","F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results","FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","FH1_3_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results", "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results","F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results","FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results","F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")
# files <- file.path(dir,run)
# names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3", "FC1_1","FC1_2", "FC1_3","FF1_1","FF1_2","FF1_3", "FC24_1","FC24_2","FC24_3","FF24_1","FF24_2","FF24_3", "FC72_1","FC72_2","FC72_3","FF72_1","FF72_2","FF72_3" )
# txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
# #write.csv(txi.rsem$counts,"Conteos.csv")
# class <- c( rep("FC0",3), rep("FF0",3),rep("FC1",3), rep("FF1",3),rep("FC24",3), rep("FF24",3),rep("FC72",3), rep("FF72",3))
# Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
# coldata <- data.frame(row.names = colnames(Mat) ,class)
# Gene <- Mat$gene_id
#______________________________________________________________________________#
y <- DGEList(counts = Mat, group = class, genes = Gene)

# Normalizacion TMM
y <- calcNormFactors(y)
Person <- cor(y$counts,method = "pearson")

#write.csv(Person, file = "Pearson.csv")
##colors <- c("0","1")
##my_palette <- c("Red",colorRampPalette(colors = c("Red", "Blue"))(n = length(colors)), "Blue")
##pheatmap(Person, main = "Coeficiente de Correlaci?n de Pearson", cluster_rows = F, scale = "none" ,cluster_cols = F, color = my_palette, display_numbers = TRUE, number_color = "white", fontsize_number = 8 )
#rainbow(8, start = 0, end = 1)
plotMDS(y$counts, xlab = "Dim1", ylab = "Dim2"  , main = "MDS conteos",col=c(rep("#FF0000",3),rep("#FFDB00",3),rep("Black",3),rep("#00FF92",3),rep("#0092FF",3),rep("#4900FF",3),rep("#FF00DB",3),rep("#49FF00",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
#et <- exactTest(y)
#counts <- getCounts(y)
#y$samples


################################################################################
#                             T0
################################################################################
et = exactTest(y, dispersion=y$tagwise.dispersion , pair=c("FC0","FF0")) #Primero COntrol despu?s Tratamiento.
eT0 <- et
yT0 <- topTags(et, n=Inf)$table
#y2 = rownames(y2)[y2$FDR < 0.05 & abs(y2$logFC) >2]
ynames = rownames(yT0)[yT0$FDR <= 0.05]
sum(yT0$logFC>=2 & yT0$FDR<=0.05)
sum(yT0$logFC<=-2 & yT0$FDR<=0.05)
write.csv(yT0, file = "resT0.csv")
write.csv(yT0[yT0$FDR <= 0.05,], file = "res05T0.csv")
summary(de <- decideTestsDGE(eT0, adjust.method="BH", p=0.05))


################################################################################
#                             T1
################################################################################
et = exactTest(y, dispersion=y$tagwise.dispersion , pair=c("FC1","FF1")) #Primero COntrol despu?s Tratamiento.
eT1 <- et
yT1 <- topTags(et, n=Inf)$table
#y2 = rownames(y2)[y2$FDR < 0.05 & abs(y2$logFC) >2]
ynames = rownames(yT1)[yT1$FDR <= 0.05]
sum(yT1$logFC>=2 & yT1$FDR<=0.05)
sum(yT1$logFC<=-2 & yT1$FDR<=0.05)
write.csv(yT1, file = "resT1.csv")
write.csv(yT1[yT1$FDR <= 0.05,], file = "res05T1.csv")
summary(de <- decideTestsDGE(eT1, adjust.method="BH", p=0.05))

################################################################################
#                             T24
################################################################################
et = exactTest(y, dispersion=y$tagwise.dispersion , pair=c("FC24","FF24")) #Primero COntrol despu?s Tratamiento.
eT24 <- et
yT24 <- topTags(et, n=Inf)$table
#y2 = rownames(y2)[y2$FDR < 0.05 & abs(y2$logFC) >2]
ynames = rownames(yT24)[yT24$FDR <= 0.05]
sum(yT24$logFC>=2 & yT24$FDR<=0.05)
sum(yT24$logFC<=-2 & yT24$FDR<=0.05)
write.csv(yT24, file = "resT24.csv")
write.csv(yT24[yT24$FDR <= 0.05,], file = "res05T24.csv")
summary(de <- decideTestsDGE(eT24, adjust.method="BH", p=0.05))

################################################################################
#                             T72
################################################################################
et = exactTest(y, dispersion=y$tagwise.dispersion , pair=c("FC72","FF72")) #Primero COntrol despu?s Tratamiento.
eT72 <- et
yT72 <- topTags(et, n=Inf)$table
#y2 = rownames(y2)[y2$FDR < 0.05 & abs(y2$logFC) >2]
ynames = rownames(yT72)[yT72$FDR <= 0.05]
sum(yT72$logFC>=2 & yT72$FDR<=0.05)
sum(yT72$logFC<=-2 & yT72$FDR<=0.05)
write.csv(yT72, file = "resT72.csv")
write.csv(yT72[yT72$FDR <= 0.05,], file = "res05T72.csv")
summary(de <- decideTestsDGE(eT72, adjust.method="BH", p=0.05))

out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]

#______________________________________________________________________________#

