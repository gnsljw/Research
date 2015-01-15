library(edgeR)
setwd("D://Work/PTMC/20150112/(3) Hashimoto/")
targets <- readTargets("targets (3-2).txt")
head(targets)
x <- read.delim("(3-2) Thyroiditis(-)_PTMC(24)PTC(339).txt")
dim(x)

y <- DGEList(count=x[,2:364], group=targets$Size, genes=x[,1])
y

library(org.Hs.eg.db)

egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m <- match(y$genes$genes, egSYMBOL$gene_id)
Symbol <- egSYMBOL$symbol[m]
Symbol[which(is.na(Symbol))] <- "?"
head(Symbol)
head(m)
y$genes$symbol <- Symbol
y

rownames(y$counts) <- y$genes$symbol
y

## o <- order(rowSums(y$counts), decreasing=TRUE)
## y <- y[o,]
## d <- duplicated(y$genes$Symbol)
## y <- y[!d,]
## y <- y[is.na==T, ]
nrow(y)

y$samples$lib.size <- colSums(y$counts)
head(y$samples)

keep <- rowSums(cpm(y) > 1) >= 3   # keep genes with at least 1 counts per million (CPM) in at least three samples
y <- y[keep,]
dim(y) # some genes are filtered 

y <- calcNormFactors(y, method="upperquartile")
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)

plotMDS(y, cex=0.5, col=c(rep("red",24),rep("blue",339)))
plotBCV(y)

et <- exactTest(y)
Top <- topTags(et, 20000)
de <- decideTestsDGE(et, adjust.method="BH", p.value=0.01)
summary(de <- decideTestsDGE(et, adjust.method="BH", p.value=0.01))
detags <- rownames(y)[as.logical(de)]

de <- decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(de <- decideTestsDGE(et, adjust.method="BH", p.value=0.05))
detags <- rownames(y)[as.logical(de)]

plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
write.csv(Top, "D://Research//014 MicroPTC//Results/(2-2) EdgeR result.csv")

# heatmap

Top <- topTags(et, 50)
y <- DGEList(count=x[,2:364], group=targets$Size, genes=x[,1])
rownames(y$counts) <- y$genes$Symbol <- Symbol
mat <- y$counts
head(mat[,1:5]); dim(mat)
mat <- mat[as.numeric(rownames(Top)), ]
head(mat[,1:5]); dim(mat)

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
ptcol <- c(rep("red",24), rep("skyblue",339))
hc.complete = function(mat) hclust(mat, method = "complete")
hc.single = function(mat) hclust(mat, method = "single")
hc.avr = function(mat) hclust(mat, method = "average")

heatmap(mat, col=hmcol, cexCol=0.9, cexRow=0.9, ColSideColors = ptcol, hclustfun = hc.complete)
heatmap(mat, col=hmcol, cexCol=0.9, cexRow=0.9, ColSideColors = ptcol, hclustfun = hc.single)
heatmap(mat, col=hmcol, cexCol=0.9, cexRow=0.9, ColSideColors = ptcol, hclustfun = hc.avr)

# GO test
go <- goana(et)
topGO(go,ont="BP",sort="Up",n=30)
