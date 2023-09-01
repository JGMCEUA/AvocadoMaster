setwd("F:/ETAPA 2/R/DEG/RSEM")
library(tidyverse)

#Comandos para concatenar las tablas de expresión y las anotaciones
Cont <- as.data.frame(read.csv("resDF05TOTALDF.csv"))
Anot <- as.data.frame(read.csv("Res_TrP_SP_T10.csv"))
Gen <- as.data.frame(read.csv("Genes.csv"))

colnames(Anot)
head(Anot$genes)
head(Cont$genes)
Tot <- inner_join(Cont, Anot, "transcript_id")
head(Tot)
write.csv(Tot,"ReportexTiempo_TrP_SP_T10.csv")

TotGen <- inner_join(Gen, Anot[,c(2,3,7)], "genes")
head(TotGen)
write.csv(TotGen,"GenesxTiempo_TrP_SP_T10.csv")

#______________________________________________________________________________#
Cont <- as.data.frame(read.csv("res05TOTAL.csv"))
Anot <- as.data.frame(read.csv("Res_TrP_SP_T10.csv"))
Tot <- inner_join(Cont, Anot, "transcript_id")
write.csv(Tot,"ReportexTiempoTOT_TrP_SP_T10.csv")

Total <- as.data.frame(read.csv("ID_REP_TOT.csv"))
head(Total)
###Tabla de expresados
colnames(Total)
UPT0 <- Total[ Total$T0 > 0 ,c(1,2,5,9,12)]
DNT0 <- Total[ Total$T0 < 0 ,c(1,2,5,9,12)]

UPT24 <- Total[ Total$T24 > 0 ,c(1,2,5,9,12)]
DNT24 <- Total[ Total$T24 < 0 ,c(1,2,5,9,12)]

UPT72 <- Total[ Total$T72 > 0 ,c(1,2,5,9,12)]
DNT72 <- Total[ Total$T72 < 0 ,c(1,2,5,9,12)]

write.csv(UPT0,"T0UPID.csv",row.names = F)
write.csv(DNT0,"T0DNID.csv",row.names = F)

write.csv(UPT24,"T24UPID.csv",row.names = F)
write.csv(DNT24,"T24DNID.csv",row.names = F)

write.csv(UPT72,"T72UPID.csv",row.names = F)
write.csv(DNT72,"T72DNID.csv",row.names = F)

