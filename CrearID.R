library(readr)
library(tidyverse)
library(dslabs)

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/SampletxtMuestra")

####
DEGT0 <- as.data.frame(read.csv("res05DFT0.csv"))
DEGT1 <- as.data.frame(read.csv("res05DFT1.csv"))
DEGT24 <- as.data.frame(read.csv("res05DFT24.csv"))
DEGT72 <- as.data.frame(read.csv("res05DFT72.csv"))

DT0 <- subset(DEGT0, subset =   DEGT0$lfcSE <= 0.5)
DT1 <- subset(DEGT1, subset =   DEGT1$lfcSE <= 0.5)
DT24 <- subset(DEGT24, subset =  DEGT24$lfcSE <= 0.5)
DT72 <- subset(DEGT72, subset =  DEGT72$lfcSE <= 0.5)

DT01 <- full_join(DT0, DT1, "gene_id")
DT0124 <- full_join(DT01, DT24, "gene_id")
DEGT <- full_join(DT0124, DT72,"gene_id")
write.csv(DEGT,"FinalBoss.csv")


#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('PAH_ittepic.gtf')
gtf_df=as.data.frame(gtf)
gtf_df2 <- subset(gtf_df,subset= gtf_df$type!="exon",)
write.csv(gtf_df2,"Pru.csv")
DF <- read.csv("FinalBoss.csv")

CD <- inner_join(DF, gtf_df2,"gene_id")

write.csv(CD,file = "FinalBoss2.csv")