library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
setwd("F:/ETAPA 2/R/DEG/RSEM")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results",
         "F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results")

files <- file.path(dir,run)

names(files) <- c("FC0_1","FC0_2","FC0_3",
                  "FF0_1","FF0_2","FF0_3")

txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
class <- c( rep("FC0",3), 
            rep("FF0",3))

coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id

Tratamiento <- c( rep("Control",3), 
                  rep("Fructanos",3))

Tiempo <- c( rep("0h",3))
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

my.contrasts <- makeContrasts( FFvsFC.T0 = Fructanos.0h-Control.0h, levels=design)
for(lol in colnames(my.contrasts)){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,lol])
  print(summary(deT0 <- decideTestsDGE(qlf, adjust.method="BH", p=0.05)))
}



