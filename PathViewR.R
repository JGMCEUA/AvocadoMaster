library(tidyverse)
library(pathview)

#script para generar los mapas de pathview con los datos de expresión
setwd("F:/ETAPA 2/DEG/AnalisisR")
KODEG <- as.data.frame(read.csv("KOxDEG.csv"))
y <- 3
Name <- colnames(KODEG)
#Mapas <- c("00940","00941","00040","04626","01212","00061","00380","00195","00999","04075","00460", "00710","00196","01200","00630","00780", "04016","00400","00270","00905","04712","00902", "00190", "03040", "04145", "00010", "04144", "04141","00660","03013", "01230", "00620","00770","04146","01110")
# 1: Mapas <- c("00940","00941","00040","04626","01212","00061","00380","00195","00999","04075","00460", "00710","00196","01200","00630","00780", "04016","00400","00270","00905","04712","00902", "00190", "03040", "04145", "00010", "04144", "04141","00660","03013", "01230", "00620","00770","04146","01110")
Mapas <- c("00020","00030")
for(y in c("T0","T24","T72")) {
  if (y == "T0"){
    N<- 3
  }else if (y == "T24"){
    N<- 4
  }else{
    N<- 5
  }
  testVect <- structure(KODEG[KODEG[,N] != "NA",N], .Names =KODEG[KODEG[,N] != "NA",2])
  for(x in Mapas){
    pv.out <- pathview(gene.data = testVect, pathway.id = x, species = "ko", out.suffix = paste("Pathview.Gustavo",y, sep = "."), kegg.native = T)
    str(pv.out)
    head(pv.out$plot.data.gene)
  }}


