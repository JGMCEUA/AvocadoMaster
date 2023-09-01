library(readr)
library(tidyverse)
library(dslabs)

################################################################################
#                         RTP0
################################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT0P")
#Cargar resultado de las anotaciones por blastp
Anotaciones <- as.matrix(read.delim("blastx.outfmt6",header = F, sep = "\t"))
colnames(Anotaciones) <- c("seqnames","UniprotID","Identidad%","LongAlineamiento","#NoAlinear",
                           "#GapsAbiertos","InicioConsulta","FinConsulta","InicioMuestra","FinMuestra",
                           "eValue","bitscore")
Anotaciones <- as.data.frame(Anotaciones)
head(Anotaciones)
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('PAH_ittepic.gtf')
gtf_df=as.data.frame(gtf)
gtf_df2 <- subset(gtf_df,subset= gtf_df$type!="exon",)
head(gtf_df2)
#Cargar análisis de expresión diferencial
DEG <- as.data.frame(read.csv("res05DF.csv"))
head(DEG)
#Escribir todos los resultados unidos
ResDEG <- full_join(Anotaciones,gtf_df2, by="seqnames")
head(ResDEG)
ResDEG <- ResDEG[,!(names(ResDEG) %in% c("phase","exon_number"))]
ResDEG <- na.omit(ResDEG)
ResDEG2 <- ResDEG
ResDEG2 <- full_join(ResDEG2, DEG, by="gene_id")
write.csv(ResDEG,"AnotaTotalGTF.csv")
head(ResDEG)
#Escrbir los resultados con Match
borrar <- c("score", "phase", "exon_number", "significant", "source")
ResT <- ResDEG2[ , !(names(ResDEG2) %in% borrar)]
head(ResT)
ResT <- na.omit(ResT)
write.csv(ResT,file = "AnotaRepGTF.csv")
head(ResT)






###Script para eliminar los exones y mantener solo los transcritos




setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT0P")
Gen <- as.data.frame(read.csv("AnotaRep.csv"))
Trans <- subset(Gen, subset = Gen$type != "exon" )
write.csv(Trans, file = "AnotaRepTrans.csv")

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1P")
Gen <- as.data.frame(read.csv("AnotaRep.csv"))
Trans <- subset(Gen, subset = Gen$type != "exon" )
write.csv(Trans, file = "AnotaRepTrans.csv")

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT24P")
Gen <- as.data.frame(read.csv("AnotaRep.csv"))
Trans <- subset(Gen, subset = Gen$type != "exon" )
write.csv(Trans, file = "AnotaRepTrans.csv")

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT72P")
Gen <- as.data.frame(read.csv("AnotaRep.csv"))
Trans <- subset(Gen, subset = Gen$type != "exon" )
write.csv(Trans, file = "AnotaRepTrans.csv")

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1PPP")
Gen <- as.data.frame(read.csv("AnotaRep.csv"))
Trans <- subset(Gen, subset = Gen$type != "exon" )
write.csv(Trans, file = "AnotaRepTrans.csv")

