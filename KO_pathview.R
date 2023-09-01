#scrip para generar los mapas de pathview con los datos de expresión
setwd("C:/Users/gusta/OneDrive/Escritorio/AngelicaDocs")

library(tidyverse)
library(pathview)

KODEG <- as.data.frame(read.csv("DGE_KO_1.csv"))
colnames(KODEG)
head(KODEG)
testVect <- structure(KODEG[KODEG$T0!="NA",3], .Names =KODEG[KODEG$T0!="NA",2])
Mapas <- c("04626","04016","04075","00999","00270","00905","04712","00902","00941",
           "00940","00190", "03040", "04145", "00010", "04144", "00710", "04141", "00061",
           "00660", "00780","03013", "01230", "00620","01200","00770","04146","01110" )
Mapas <- c("00940","00941","00040","04626","01212","00061","00380","00195","00999","04075","00460", "00710","00196","01200","00630","00780")
for(x in Mapas){
  pv.out <- pathview(gene.data = testVect, pathway.id = x, species = "ko", out.suffix = "Pru.Fin.T1", kegg.native = T)
  str(pv.out)
  head(pv.out$plot.data.gene)
}

