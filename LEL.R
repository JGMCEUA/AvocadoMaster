library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
#setwd("/media/alejandra/AgroBio_Tec/ETAPA 2/R/DEG/RSEM")
setwd("F:/ETAPA 2/R/DEG/RSEM")
dir <- getwd()
run <- c("FH1_1_RHC.genes.results","FH1_2_RHC.genes.results",
         "F1_1_RHC.genes.results","F1_2_RHC.genes.results","F1_3_RHC.genes.results")


files <- file.path(dir,run)
names(files) <- c("FC1_1","FC_2",
                  "FF1_1","FF1_2","FF1_3" )

txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
#write.csv(txi.rsem$counts,"Conteos-OK.csv")
class <- c( rep("FC1",2), 
            rep("FF1",3))

Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id



y <- DGEList(counts = Mat, group = class, genes = Gene)

keep = rowSums(cpm(y)>2) >= 2
y2 = y[keep, , keep.lib.sizes=FALSE]
y3 <- DGEList(counts = y2$counts, group = class, genes = rownames(y2$counts))
y3 <- calcNormFactors(y3)
y3 <- estimateCommonDisp(y3, verbose=TRUE)
y3 <- estimateTagwiseDisp(y3)

################################################################################
#                             T0
################################################################################
et = exactTest(y3, dispersion=y3$tagwise.dispersion , pair=c("FC1","FF1")) #Primero COntrol despu?s Tratay3iento.
eT0 <- et
y3T0 <- topTags(et, n=Inf)$table
yUDT0 = rownames(y3T0)[y3T0$FDR < 0.05 & abs(y3T0$logFC) >= 1]
y05T0 = rownames(y3T0)[y3T0$FDR < 0.05]
sum(y3T0$logFC>=1 & y3T0$FDR<0.05)
sum(y3T0$logFC<=-1 & y3T0$FDR<0.05)
sum(abs(y3T0$logFC)>=1 & y3T0$FDR<0.05)
summary(deT0 <- decideTestsDGE(eT0, adjust.method="BH", p=0.05))
sum(y3T0$FDR<0.05)




