library(edgeR)
setwd("C://Work/PTMC/20150112/(1) All samples/")
targets <- readTargets("targets.txt")
head(targets)
x <- read.delim("(1) PTMC vs PTC gene id.txt")
dim(x)

y <- DGEList(count=x[,2:487], group=targets$Size, genes=x[,1])
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

rownames(y$counts) <- y$genes$genes  # <<- gene names should be converted into gene ID
y

y$samples$lib.size <- colSums(y$counts)
head(y$samples)

keep <- rowSums(cpm(y) > 1) >= 3   # keep genes with at least 1 counts per million (CPM) in at least three samples
y <- y[keep,]
dim(y) # some genes are filtered 

y <- calcNormFactors(y, method="upperquartile")
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)

et <- exactTest(y)

go <- goana(et, FDR = 0.05, species="Hs")

write.csv(go, "G://Research//014 MicroPTC//Results/(2) BRAF/(2-2) GO result.csv")
