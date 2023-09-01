library(tidyverse)
library(pathview)

setwd("F:/ETAPA 2/R/DEG/RSEM")
T0 <- as.data.frame(read.csv("0613T0KAAS.csv"))
colnames(T0) <- c("genes","ko")
T24 <- as.data.frame(read.csv("0613T24KAAS.csv"))
colnames(T24) <- c("genes","ko")
T72 <- as.data.frame(read.csv("0613T72KAAS.csv"))
colnames(T72) <- c("genes","ko")


ANOT <- as.data.frame(read.csv("ReportexTiempo.csv"))
CT2 <- as.data.frame(read.csv("koTotal.csv"))
RxTxK <- inner_join(ANOT,CT2,"genes")
write.csv(RxTxK,"ReportexTiempoxKO.csv", row.names = F)

Trinot <- as.data.frame(read.csv("ReportexTiempo.csv"))

inner_join(Trinot,CT2,"genes")

DEG <- as.data.frame(read.csv("resDF05TOTALDF.csv"))
DEG <- DEG[, c(1,2)]
head(DEG)

KAAST0 <- inner_join(T72,DEG,"genes")
head(KAAST0)
write.csv(KAAST0,"KT72.csv")

KODEG <- as.data.frame(read.csv("KOxDEG.csv"))

#Tab <- read.csv("KosasT0.csv", header = TRUE, as.is=TRUE,row.names=1,  quote= "")
colnames(KODEG)

testVect <- structure(KODEG[KODEG$T72!="NA",5], .Names =KODEG[KODEG$T72!="NA",2])

#testVect <- structure(KAAST0$T0, .Names = KAAST0$ko)


for(x in Mapas){
  pv.out <- pathview(gene.data = testVect, pathway.id = x, species = "ko", out.suffix = "Pru.Fin.T72", kegg.native = T)
  str(pv.out)
  head(pv.out$plot.data.gene)
}

#scrip para generar los mapas de pathview con los datos de expresión
#setwd("C:/Users/gusta/OneDrive/Escritorio/AngelicaDocs")
setwd("G:/AngelicaDocs/AngelicaDocs")
library(tidyverse)
library(pathview)

KODEG <- as.data.frame(read.csv("DGE_KO_1.csv"))
colnames(KODEG)
head(KODEG)

Mapas <- c("00940","00941","00040","04626","01212","00061","00380","00195","00999","04075","00460", "00710","00196","01200","00630","00780", "04016","00400","00270","00905","04712","00902", "00190", "03040", "04145", "00010", "04144", "04141","00660","03013", "01230", "00620","00770","04146","01110")
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
    pv.out <- pathview(gene.data = testVect, pathway.id = x, species = "ko", out.suffix = paste("Pathview.Angelica",y, sep = "."), kegg.native = T)
    str(pv.out)
    head(pv.out$plot.data.gene)
  }}

  


