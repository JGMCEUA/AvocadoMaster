#Análisis EdgeR
library("baySeq")
library("edgeR")

###############################################################################
#                             T0                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT0")
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
coldata <- data.frame( row.names = colnames(countDataT0), condition)
countDataT0 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
countDataT0 <- countDataT0[, rownames(coldata)]
DgL <- DGEList(counts = countDataT0, group = condition)
DgL <- estimateCommonDisp(DgL)
DgL  #133,237,367
### Se filtran los datos para que solo los que tienen 100 CPM queden en la variable
dim(DgL)
DgL.full <- DgL
head(DgL$counts)
head(cpm(DgL))
apply(DgL$counts, 2, sum)
keep <- rowSums(cpm(DgL)>100) >= 2
DgL <- DgL[keep,]
dim(DgL)
### Comprobar Datos del filtrado
DgL$samples$lib.size <- colSums(DgL$counts)
DgL$samples #48,944,900
DgL <- calcNormFactors(DgL)
DgL
### Gráficar MA
plotMDS(DgL, method = "bcv", col = as.numeric(DgL$samples$group))
legend("bottomfelt", as.character(unique(DgL$samples$group)), col = 1:3, pch = 20)

###Estimando la dispersión [NO SE RECOMIENDA USARLO] Investigar porque es necesario estimar este parámetro y como se puede probar otras distribuciones
d1 <- estimateCommonDisp(DgL, verbose = T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

###Estimar la dispersión por GLM (Generalized Linear Model)
design.mat <- model.matrix(~ 0 + DgL$samples$group)
colnames(design.mat) <- levels(DgL$samples$group)
d2 <- estimateGLMCommonDisp(DgL, design.mat)
d2 <- estimateGLMTrendedDisp(DgL,design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

###Análisis de expresión diferencial
et12 <- exactTest(d1, pair = c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05)
summary(de1)
##Gráfica de análisis de expresión diferencial
de1tags12 <- row.names(d1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2,2), col = "blue")
###GLM test para la expresión diferencial
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast = c(-1, 1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method = "BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags = de2tags12)

###############################################################################
#                             T1                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT1")
countDataT1 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
coldata <- data.frame( row.names = colnames(countDataT1), condition)
countDataT1 <- countDataT1[, rownames(coldata)]
DgL <- DGEList(counts = countDataT1, group = condition)
DgL <- estimateCommonDisp(DgL)
DgL  #133,237,367
### Se filtran los datos para que solo los que tienen 100 CPM queden en la variable
dim(DgL)
DgL.full <- DgL
head(DgL$counts)
head(cpm(DgL))
apply(DgL$counts, 2, sum)
keep <- rowSums(cpm(DgL)>100) >= 2
DgL <- DgL[keep,]
dim(DgL)
### Comprobar Datos del filtrado
DgL$samples$lib.size <- colSums(DgL$counts)
DgL$samples #48,944,900
DgL <- calcNormFactors(DgL)
DgL
### Gráficar MA
plotMDS(DgL, method = "bcv", col = as.numeric(DgL$samples$group))
legend("bottomfelt", as.character(unique(DgL$samples$group)), col = 1:3, pch = 20)

###Estimando la dispersión [NO SE RECOMIENDA USARLO] Investigar porque es necesario estimar este parámetro y como se puede probar otras distribuciones
d1 <- estimateCommonDisp(DgL, verbose = T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

###Estimar la dispersión por GLM (Generalized Linear Model)
design.mat <- model.matrix(~ 0 + DgL$samples$group)
colnames(design.mat) <- levels(DgL$samples$group)
d2 <- estimateGLMCommonDisp(DgL, design.mat)
d2 <- estimateGLMTrendedDisp(DgL,design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

###Análisis de expresión diferencial
et12 <- exactTest(d1, pair = c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05)
summary(de1)
##Gráfica de análisis de expresión diferencial
de1tags12 <- row.names(d1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2,2), col = "blue")
###GLM test para la expresión diferencial
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast = c(-1, 1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method = "BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags = de2tags12)
abline(h= c(-2, 2), col = "blue")

###############################################################################
#                             T24                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT24")
countDataT24 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
coldata <- data.frame( row.names = colnames(countDataT24), condition)
countDataT24 <- countDataT24[, rownames(coldata)]
DgL <- DGEList(counts = countDataT24, group = condition)
DgL <- estimateCommonDisp(DgL)
DgL  #133,237,367
### Se filtran los datos para que solo los que tienen 100 CPM queden en la variable
dim(DgL)
DgL.full <- DgL
head(DgL$counts)
head(cpm(DgL))
apply(DgL$counts, 2, sum)
keep <- rowSums(cpm(DgL)>100) >= 2
DgL <- DgL[keep,]
dim(DgL)
### Comprobar Datos del filtrado
DgL$samples$lib.size <- colSums(DgL$counts)
DgL$samples #48,944,900
DgL <- calcNormFactors(DgL)
DgL
### Gráficar MA
plotMDS(DgL, method = "bcv", col = as.numeric(DgL$samples$group))
legend("bottomfelt", as.character(unique(DgL$samples$group)), col = 1:3, pch = 20)

###Estimando la dispersión [NO SE RECOMIENDA USARLO] Investigar porque es necesario estimar este parámetro y como se puede probar otras distribuciones
d1 <- estimateCommonDisp(DgL, verbose = T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

###Estimar la dispersión por GLM (Generalized Linear Model)
design.mat <- model.matrix(~ 0 + DgL$samples$group)
colnames(design.mat) <- levels(DgL$samples$group)
d2 <- estimateGLMCommonDisp(DgL, design.mat)
d2 <- estimateGLMTrendedDisp(DgL,design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

###Análisis de expresión diferencial
et12 <- exactTest(d1, pair = c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05)
summary(de1)
##Gráfica de análisis de expresión diferencial
de1tags12 <- row.names(d1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2,2), col = "blue")
###GLM test para la expresión diferencial
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast = c(-1, 1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method = "BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags = de2tags12)
abline(h= c(-2, 2), col = "blue")

###############################################################################
#                             T72                                              #
###############################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT72")
countDataT72 <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
coldata <- data.frame( row.names = colnames(countDataT72), condition)
countDataT72 <- countDataT72[, rownames(coldata)]
DgL <- DGEList(counts = countDataT72, group = condition)
DgL <- estimateCommonDisp(DgL)
DgL  #133,237,367
### Se filtran los datos para que solo los que tienen 100 CPM queden en la variable
dim(DgL)
DgL.full <- DgL
head(DgL$counts)
head(cpm(DgL))
apply(DgL$counts, 2, sum)
keep <- rowSums(cpm(DgL)>100) >= 2
DgL <- DgL[keep,]
dim(DgL)
### Comprobar Datos del filtrado
DgL$samples$lib.size <- colSums(DgL$counts)
DgL$samples #48,944,900
DgL <- calcNormFactors(DgL)
DgL
### Gráficar MA
plotMDS(DgL, method = "bcv", col = as.numeric(DgL$samples$group))
legend("bottomfelt", as.character(unique(DgL$samples$group)), col = 1:3, pch = 20)

###Estimando la dispersión [NO SE RECOMIENDA USARLO] Investigar porque es necesario estimar este parámetro y como se puede probar otras distribuciones
d1 <- estimateCommonDisp(DgL, verbose = T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

###Estimar la dispersión por GLM (Generalized Linear Model)
design.mat <- model.matrix(~ 0 + DgL$samples$group)
colnames(design.mat) <- levels(DgL$samples$group)
d2 <- estimateGLMCommonDisp(DgL, design.mat)
d2 <- estimateGLMTrendedDisp(DgL,design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

###Análisis de expresión diferencial
et12 <- exactTest(d1, pair = c(1,2))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05)
summary(de1)
##Gráfica de análisis de expresión diferencial
de1tags12 <- row.names(d1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2,2), col = "blue")
###GLM test para la expresión diferencial
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast = c(-1, 1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method = "BH", p.value = 0.05)
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags = de2tags12)
abline(h= c(-2, 2), col = "blue")

