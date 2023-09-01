################################################################################
# Nombre: José Gustavo Marín Contreras
# Práctica: Enriquecimiento génico e identifcación de rutas metabólicas comparadas entre los 3 tiempos
# Fecha: 02/08/2022
###############################################################################

library(org.Hs.eg.db)
library(org.At.tair.db)

data(geneList, package="DOSE")
head(geneList)

data(gcSample)
str(gcSample) 
gcSample

ck <- compareCluster(geneCluster = gcSample, fun = enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(formula_res)
dotplot(ck)
dotplot(formula_res)
dotplot(formula_res, x="group") + facet_grid(~othergroup)
cnetplot(ck)

###############################################################################

setwd("G:/ETAPA 2/R/DEG/RSEM")
DB <- as.data.frame(read.csv("CoeTAIRF.csv"))
X0 <- DB[abs(DB$T0exp)>0,3]
X24 <- DB[abs(DB$T24exp)>0,3]
X72 <- DB[abs(DB$T72exp)>0,3]
GOsam <- list(X0,X24,X72)
names(GOsam) <- c("T0","T24","T72")
str(GOsam)

kc <- compareCluster(geneClusters = GOsam, fun = enrichKEGG, organism = "ath")
kc2 <- setReadable(kc, OrgDb = org.At.tair.db, keyType="TAIR")
dotplot(kc2)
cnetplot(kc2)
################################################################################
geneList <- DB[DB$T72sig<=0.05,4]
names(geneList) <- DB[DB$T72sig<= 0.05,"TAIR.locus.name"]

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
#mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
formula_res <- compareCluster(Entrez~group, data=mydf, fun="enrichKEGG", organism = "ath")
head(formula_res)
dotplot(formula_res, x="group")
