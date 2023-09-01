################################################################################
# Nombre: José Gustavo Marín Contreras
# Práctica: Enriquecimiento génico e identifcación de rutas metabólicas con cluster pfofile
# Fecha: 01/08/2022
###############################################################################

#Cargar librerías
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)

#Dirección en el disco duro de los archivos
setwd("G:/ETAPA 2/R/DEG/RSEM")
################################################################################
#             ENRIQUECIMIENTO DE LOS PROCESOS BIOLÓGICOS
################################################################################
# Primero se realizará la identifación de las respuestas dada de acuerdo al genoma 
# de arabidopsis con los transcritos ensamblados y anotados con la base de datos de TAIR10
DB <- as.data.frame(read.csv("CoeTAIRF.csv"))
#Se ingresan los resultaddos de únicamnete los transcritos que fueron significativos en su expresión diferencial, al ser solo 3 tiempos se requirieron los 
x <- 9
y <- "ALL"
for (x in colnames(DB[,c(7,8,9)])) {
  DE <- c(DB[DB[,x] <=0.05,"TAIR.locus.name"]) 
  DE <- na.omit(DE)
  for (y in c("BP","MF","CC","ALL")) {
    ego <- enrichGO(DE, OrgDb = "org.At.tair.db", keyType = 'TAIR', ont = y ,readable=TRUE)
    if(y=="ALL"){
      #png(filename = paste(paste("EnrichGO_dotplot",x,y,sep = "_"),"png",sep = "."), width = 1000, height = 700)
      print(dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free"))
      #dev.off()
      
    }else {
      if (y=="BP"){
        ego2 <- simplify(ego)
        geneList <- DB[DB[,x]<=0.05,x-3]
        names(geneList) <- DB[DB[,x]<= 0.05,"TAIR.locus.name"]
        print(cnetplot(ego2, foldChange=geneList))
      }
      #png(filename = paste(paste("EnrichGO_barplot",x,y,sep = "_"),"png",sep = "."), width = 1000, height = 700)
      print(barplot(ego, showCategory=20))
      #dev.off()
      #png(filename = paste(paste("EnrichGO_goplot",x,y,sep = "_"),"png",sep = "."), width = 1000, height = 700)
      print(goplot(ego))
      #dev.off()
    }}}
enrichKEGG()

cnetplot(ego2, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 

ego2 
