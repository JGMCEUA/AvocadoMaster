library(tidyverse)
library(pathview)

KODEG <- as.data.frame(read.csv("KOxDEG.csv"))
colnames(KODEG)
testVect <- structure(KODEG[KODEG$T72!="NA",5], .Names =KODEG[KODEG$T72!="NA",2])
Mapas <- c("04626","04016","04075","00999","00270","00905","04712","00902","00941",
           "00940","00190", "03040", "04145", "00010", "04144", "00710", "04141", "00061",
           "00660", "00780","03013", "01230", "00620","01200","00770","04146","01110" )

for(x in Mapas){
  pv.out <- pathview(gene.data = testVect, pathway.id = x, species = "ko", out.suffix = "Pru.Fin.T72", kegg.native = T)
  str(pv.out)
  head(pv.out$plot.data.gene)
}
  


