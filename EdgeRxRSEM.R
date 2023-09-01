library(tximport)
library(edgeR)

################################################################################
#                                 T0                                           #
################################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("F:/ETAPA 2/DEG/AnalisisR")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results","F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC0_1","FC0_2","FC0_3","FF0_1","FF0_2","FF0_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
class <- c( rep("FC",3), rep("FF",3))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
y <- DGEList(counts = Mat, group = class, genes = Gene)
# Normalizaci?n TMM
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")

plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
topTags(et)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))
out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]


write.table(out4, file="DE_genesT0.txt", quote=FALSE, row.names=T, sep="\t")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],1000),file = "T0C1000.csv")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],3000),file = "T0C3000.csv")


T0100 <- as.data.frame(read.csv("T0C1000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0100$X))
y <- DGEList(counts = Mat2, group = class, genes = T0100$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

T0300 <- as.data.frame(read.csv("T0C3000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0300$X))
y <- DGEList(counts = Mat2, group = class, genes = T0300$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

################################################################################
#                                 T1                                           #
################################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","FH1_3_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FC1_3","FF1_1","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
head(txi.rsem$counts)
class <- c( rep("FC",3), rep("FF",3))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
head(Mat)
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
Lol <- as.data.frame(read.csv("Mambo6.csv"))
Gene <- Lol$X
Mat <- Mat[rownames(Mat) %in% Lol$X,]

y <- DGEList(counts = Mat, group = class, genes = Gene)
nrow(y)
# Normalizaci?n TMM
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")

plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
topTags(et)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))
out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]


write.table(out4, file="DE_genesT1.txt", quote=FALSE, row.names=T, sep="\t")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],1000),file = "T1C1000.csv")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],3000),file = "T1C3000.csv")


T0100 <- as.data.frame(read.csv("T1C1000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0100$X))
y <- DGEList(counts = Mat2, group = class, genes = T0100$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

T0300 <- as.data.frame(read.csv("T1C3000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0300$X))
y <- DGEList(counts = Mat2, group = class, genes = T0300$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

################################################################################
#                                 T1-1                                         #
################################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FF1_1","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
head(txi.rsem$counts)
class <- c( rep("FC",2), rep("FF",3))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
head(Mat)
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
y <- DGEList(counts = Mat, group = class, genes = Gene)
nrow(y)
# Normalizaci?n TMM
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")

plotMDS(y, col=c(rep("Red",2),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
topTags(et)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))
out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]


write.table(out4, file="DE_genesT1P.txt", quote=FALSE, row.names=T, sep="\t")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],1000),file = "T1PC1000.csv")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],3000),file = "T1PC3000.csv")


T0100 <- as.data.frame(read.csv("T1PC1000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0100$X))
y <- DGEList(counts = Mat2, group = class, genes = T0100$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",2),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

T0300 <- as.data.frame(read.csv("T1PC3000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0300$X))
y <- DGEList(counts = Mat2, group = class, genes = T0300$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",2),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

################################################################################
#                                 T1-2                                         #
################################################################################
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("E:/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")
files <- file.path(dir,run)
names(files) <- c("FC1_1","FC1_2","FF1_2","FF1_3")
txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
head(txi.rsem$counts)
class <- c( rep("FC",2), rep("FF",2))
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
head(Mat)
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id
y <- DGEList(counts = Mat, group = class, genes = Gene)
nrow(y)
# Normalizaci?n TMM
y <- calcNormFactors(y, col=c(rep("Red",3),rep("blue",3)))
cor(y$counts,method = "pearson")

plotMDS(y,col=c(rep("Red",3),rep("blue",3)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
topTags(et)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))
out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out2 <- cbind(out, counts)
out3 <- out2[as.logical(de),]
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]


write.table(out4, file="DE_genesT1PP.txt", quote=FALSE, row.names=T, sep="\t")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],1000),file = "T1PPC1000.csv")
write.csv(head(out4[order(abs(out4$logFC), decreasing = T), ],3000),file = "T1PPC3000.csv")


T0100 <- as.data.frame(read.csv("T1PPC1000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0100$X))
y <- DGEList(counts = Mat2, group = class, genes = T0100$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",2),rep("blue",2)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))

T0300 <- as.data.frame(read.csv("T1PPC3000.csv"))
Mat2 <- subset(Mat, subset = row.names(Mat) %in% c(T0300$X))
y <- DGEList(counts = Mat2, group = class, genes = T0300$X)
y <- calcNormFactors(y)
cor(y$counts,method = "pearson")
plotMDS(y, col=c(rep("Red",2),rep("blue",2)))
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
counts <- getCounts(y)
summary(de <- decideTestsDGE(et, adjust.method="BH", p=0.05))
