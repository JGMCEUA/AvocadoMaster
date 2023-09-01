library(readr)
library(tidyverse)
library(dslabs)


setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/SampletxtTiem")

####
DEGT1 <- as.data.frame(read.csv("res05DFT1.csv"))

DEGT24 <- as.data.frame(read.csv("res05DFT24.csv"))

DEGT72 <- as.data.frame(read.csv("res05DFT72.csv"))

DT1 <- subset(DEGT1, subset = abs(DEGT1$log2FoldChange) >= 2 & DEGT1$lfcSE <= 0.5)
DT24 <- subset(DEGT24, subset = abs(DEGT24$log2FoldChange) >= 2 & DEGT24$lfcSE <= 0.5)
DT72 <- subset(DEGT72, subset = abs(DEGT72$log2FoldChange) >= 2 & DEGT72$lfcSE <= 0.5)

DT124 <- full_join(DT1, DT24, "gen_id")
DEGT <- full_join(DT124, DT72,"gen_id")
write.csv(DEGT,"TotaliniC.csv")
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('PAH_ittepic.gtf')
gtf_df=as.data.frame(gtf)
gtf_df2 <- subset(gtf_df,subset= gtf_df$type!="exon",)
head(gtf_df2)
write.csv(gtf_df2,"Pru.csv")
DF <- read.csv("TotaliniC.csv")

CD <- inner_join(DF, gtf_df2,"gene_id")
CD
CD2 <- CD[CD$fill =="hola",]


write.csv(CD2,file = "TotaliniFHDTC2.csv")


Cont <- data.frame(read.csv("TotaliniFHDT.csv"))
ID <- data.frame(read.csv("TotaliniFHDT6.csv"))

Juntos <- inner_join(ID,Cont,"gene_id")
write.csv(Juntos,file = "EsteBueno.csv")
