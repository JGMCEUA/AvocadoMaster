#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
library(limma)
library(tximport)
library(readr)
library(edgeR)
###############################################################################
#                             T0                                              #
###############################################################################

setwd("~/Escritorio/AlinTes2021")
condition <- factor(c(rep("F",3), rep("FH",3)))
condition = relevel(condition, ref = "F")
Dir <- c("F0_1/t_data.ctab","F0_2/t_data.ctab","F0_3/t_data.ctab","FH0_1/t_data.ctab","FH0_2/t_data.ctab","FH0_3/t_data.ctab")
files <- file.path(Dir)
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
y <- DGEList(txi$counts)
### Filtrando usando la información del diseño
SampleTable <- data.frame(condition = condition)
row.names(SampleTable) <- colnames(txi$counts)
desing <- model.matrix(~condition, data = SampleTable)
keep <- filterByExpr(y, condition)
y <- y[keep, ]
###Normalizar y correr Voom transformation
y <- calcNormFactors(y)
v <- voom(y, desing)
