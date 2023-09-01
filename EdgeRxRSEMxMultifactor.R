library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("F:/ETAPA 2/DEG/AnalisisR")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results",
         "F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results",
         "FH1_1_RHC.genes.results","FH1_2_RHC.genes.results",
         "F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results",
         "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results",
         "F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results",
         "FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results",
         "F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")

files <- file.path(dir,run)

names(files) <- c("FC0_1","FC0_2","FC0_3",
                  "FF0_1","FF0_2","FF0_3",
                  "FC1_1","FC1_2",
                  "FF1_1","FF1_2","FF1_3",
                  "FC24_1","FC24_2","FC24_3",
                  "FF24_1","FF24_2","FF24_3",
                  "FC72_1","FC72_2","FC72_3",
                  "FF72_1","FF72_2","FF72_3" )

txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
class <- c( rep("FC0",3), 
            rep("FF0",3),
            rep("FC1",2), 
            rep("FF1",3),
            rep("FC24",3), 
            rep("FF24",3),
            rep("FC72",3), 
            rep("FF72",3))

coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id

Tratamiento <- c( rep("Control",3), 
            rep("Fructanos",3),
            rep("Control",2), 
            rep("Fructanos",3),
            rep("Control",3), 
            rep("Fructanos",3),
            rep("Control",3), 
            rep("Fructanos",3))

Tiempo <- c( rep("0h",3), 
            rep("0h",3),
            rep("1h",2), 
            rep("1h",3),
            rep("24h",3), 
            rep("24h",3),
            rep("72h",3), 
            rep("72h",3))
y <- DGEList(counts = Mat, group = class, genes = Gene)

keep = rowSums(cpm(y)>2) >= 2
y2 = y[keep, , keep.lib.sizes=FALSE]
y3 <- DGEList(counts = y2$counts, group = class, genes = rownames(y2$counts))
y3 <- calcNormFactors(y3)
y3 <- estimateCommonDisp(y3, verbose=TRUE)
y3 <- estimateTagwiseDisp(y3)

targets <- cbind(files,Tratamiento,Tiempo)
targets <- targets[,-1]
targets <- as.data.frame(targets)
Group <- factor(paste(targets$Tratamiento,targets$Tiempo,sep="."))
cbind(targets,Group=Group)

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
fit <- glmQLFit(y3, design, dispersion = 0.1)

my.contrasts <- makeContrasts( FFvsFC.T0 = Fructanos.0h-Control.0h, FFvsFC.T1 = Fructanos.1h-Control.1h, FFvsFC.T24 = Fructanos.24h-Control.24h, FFvsFC.T72 = Fructanos.72h-Control.72h, FF.0vs1 = Fructanos.1h-Fructanos.0h, FF.0vs24 = Fructanos.24h-Fructanos.0h, FF.0vs72 = Fructanos.72h-Fructanos.0h, levels=design)

for(lol in colnames(my.contrasts)){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,lol])
 print(summary(deT0 <- decideTestsDGE(qlf, adjust.method="BH", p=0.05)))
}



