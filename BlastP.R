##Filtrar anotaciones de baja caliadad
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT0P")
Blast <- as.data.frame(read.csv("AnotaRepTrans.csv"))
BlastF <- subset(Blast,subset= Blast$eValue<0.0001)
write.csv(BlastF, file = "AnotaRepTransE.csv")
datos <- BlastF[!duplicated(BlastF), ]

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1PPP")
Blast <- as.data.frame(read.csv("AnotaRepTrans.csv"))
BlastF <- subset(Blast,subset= Blast$eValue<0.0001)
write.csv(BlastF, file = "AnotaRepTransE.csv")
datos <- BlastF[!duplicated(BlastF), ]

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT24P")
Blast <- as.data.frame(read.csv("AnotaRepTrans.csv"))
BlastF <- subset(Blast,subset= Blast$eValue<0.0001)
write.csv(BlastF, file = "AnotaRepTransE.csv")

setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT72P")
Blast <- as.data.frame(read.csv("AnotaRepTrans.csv"))
BlastF <- subset(Blast,subset= Blast$eValue<0.0001)
write.csv(BlastF, file = "AnotaRepTransE.csv")
