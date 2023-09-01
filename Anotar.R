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
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('FT0.gtf')
gtf_df=as.data.frame(gtf)
#Cargar análisis de expresión diferencial
DEG <- as.data.frame(read.csv("res05DF.csv"))
DEG
#Escribir todos los resultados unidos
ResDEG <- full_join(Anotaciones,gtf_df, by="seqnames")
ResDEG <- full_join(DEG, ResDEG, by="ref_gene_id")
write.csv(ResDEG,"AnotaTotal.csv")
head(ResDEG)
#Escrbir los resultados con Match
borrar <- c("score", "phase", "exon_number", "significant", "gene_id", "source")
ResT <- ResDEG[ , !(names(ResDEG) %in% borrar)]
head(ResT)
ResT <- na.omit(ResT)
write.csv(ResT,file = "AnotaRep.csv")
head(ResT)
################################################################################
#                         RTP1
################################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1P")
#Cargar resultado de las anotaciones por blastp
Anotaciones <- as.matrix(read.delim("blastx.outfmt6",header = F, sep = "\t"))
colnames(Anotaciones) <- c("seqnames","UniprotID","Identidad%","LongAlineamiento","#NoAlinear",
                           "#GapsAbiertos","InicioConsulta","FinConsulta","InicioMuestra","FinMuestra",
                           "eValue","bitscore")
Anotaciones <- as.data.frame(Anotaciones)
Anotaciones
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('FT1.gtf')
gtf_df=as.data.frame(gtf)
gtf_df
#Cargar análisis de expresión diferencial
DEG <- as.data.frame(read.csv("res05DF.csv"))
DEG
#Escribir todos los resultados unidos
ResDEG <- full_join(Anotaciones,gtf_df, by="seqnames")
ResDEG <- full_join(DEG, ResDEG, by="ref_gene_id")
write.csv(ResDEG,"AnotaTotal.csv")
head(ResDEG)
#Escrbir los resultados con Match
borrar <- c("score", "phase", "exon_number", "significant", "gene_id", "source")
ResT <- ResDEG[ , !(names(ResDEG) %in% borrar)]
head(ResT)
ResT <- na.omit(ResT)
write.csv(ResT,file = "AnotaRep.csv")
head(ResT)



################################################################################
#                         RTP24
################################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT24P")
#Cargar resultado de las anotaciones por blastp
Anotaciones <- as.matrix(read.delim("blastx.outfmt6",header = F, sep = "\t"))
colnames(Anotaciones) <- c("seqnames","UniprotID","Identidad%","LongAlineamiento","#NoAlinear",
                           "#GapsAbiertos","InicioConsulta","FinConsulta","InicioMuestra","FinMuestra",
                           "eValue","bitscore")
Anotaciones <- as.data.frame(Anotaciones)
Anotaciones
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('FT24.gtf')
gtf_df=as.data.frame(gtf)
gtf_df
#Cargar análisis de expresión diferencial
DEG <- as.data.frame(read.csv("res05DF.csv"))
DEG
#Escribir todos los resultados unidos
ResDEG <- full_join(Anotaciones,gtf_df, by="seqnames")
ResDEG <- full_join(DEG, ResDEG, by="ref_gene_id")
write.csv(ResDEG,"AnotaTotal.csv")
head(ResDEG)
#Escrbir los resultados con Match
borrar <- c("score", "phase", "exon_number", "significant", "gene_id", "source")
ResT <- ResDEG[ , !(names(ResDEG) %in% borrar)]
head(ResT)
ResT <- na.omit(ResT)
write.csv(ResT,file = "AnotaRep.csv")
head(ResT)
################################################################################
#                         RTP72
################################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT72P")
#Cargar resultado de las anotaciones por blastp
Anotaciones <- as.matrix(read.delim("blastx.outfmt6",header = F, sep = "\t"))
colnames(Anotaciones) <- c("seqnames","UniprotID","Identidad%","LongAlineamiento","#NoAlinear",
                           "#GapsAbiertos","InicioConsulta","FinConsulta","InicioMuestra","FinMuestra",
                           "eValue","bitscore")
Anotaciones <- as.data.frame(Anotaciones)
Anotaciones
#Cargar archivo GTF con referencia de los contigs ensamblados y los IDS de las lecturas alineadas
gtf <- rtracklayer::import('FT72.gtf')
gtf_df=as.data.frame(gtf)
gtf_df
#Cargar análisis de expresión diferencial
DEG <- as.data.frame(read.csv("res05DF.csv"))
DEG
#Escribir todos los resultados unidos
ResDEG <- full_join(Anotaciones,gtf_df, by="seqnames")
ResDEG <- full_join(DEG, ResDEG, by="ref_gene_id")
write.csv(ResDEG,"AnotaTotal.csv")
head(ResDEG)
#Escrbir los resultados con Match
borrar <- c("score", "phase", "exon_number", "significant", "gene_id", "source")
ResT <- ResDEG[ , !(names(ResDEG) %in% borrar)]
head(ResT)
ResT <- na.omit(ResT)
write.csv(ResT,file = "AnotaRep.csv")
head(ResT)




