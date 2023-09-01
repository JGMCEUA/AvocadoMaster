library(clusterProfiler)
data(geneList, package="DOSE")
geneList
class(geneList)
de <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

library(enrichplot)
goplot(ego)
ego2 <- simplify(ego)
ego2$ID
cnetplot(ego2, foldChange=geneList)


cnetplot(ego2, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
