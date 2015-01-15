## Data preperation

dat <- read.csv("D://Work/Galectin3/LGALS3 JW.csv")
dat2 <- read.delim("D://Work/Galectin3/galectin related genes.txt")

rownames(dat2) <- dat2[,1]
dat2 <- dat2[,2:135]
dat2.matrix <- as.matrix(dat2)

## t-test

l = c(rep("SNP",6), rep("no SNP",128))
group = factor(l, levels=c("SNP","no SNP"))
t.test(dat2.matrix[1,]~group)

tp = vector()
tf = vector()
for (i in 1:nrow(dat2.matrix)) {
        tmp <- t.test(dat2.matrix[i,] ~ group, paired = FALSE)
        tf[i] <- log2(tmp$estimate[1]/tmp$estimate[2])
        tp[i] <- tmp$p.value
}

tf
tp
id <- which(tp<0.05)
tp[id]
result <- data.frame(rownames(dat2.matrix), tp, tf)

library(dplyr)
result <- tbl_df(result)
result
names(result) <- c("ID","p.value","fold change")
result
signi <- filter(result, p.value < 0.05)

dat2 <- read.delim("D://Work/Galectin3/galectin related genes.txt")
dat3 <- merge(dat2, signi, by=c("ID"))
head(dat3)

dat3 <- dat3[,1:135]
head(dat3)

rownames(dat3) <- dat3[,1]
dat3
dat3 <- dat3[,2:135]
dat3
dat3.matrix <- as.matrix(dat3)
dat3.matrix

## Heatmap 

library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
ptcol <- c(rep("red",6), rep("skyblue",128))
dat3.trans <- t(dat3.matrix)

hc.complete = function(dat3.trans) hclust(dat3.trans, method = "complete")
hc.single = function(dat3.trans) hclust(dat3.trans, method = "single")
hc.avr = function(dat3.trans) hclust(dat3.trans, method = "average")

heatmap(dat3.trans, col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.complete)
heatmap(dat3.trans, col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.single)
heatmap(dat3.trans, col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.avr)

heatmap.2(t(dat3.trans), trace="none", col=hmcol, cexCol=0.9, cexRow=0.8, ColSideColors = ptcol, 
        hclustfun = hc.complete)
heatmap.2(dat3.trans, trace="none", col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.single)
heatmap.2(t(dat3.trans), trace="none", col=hmcol, cexCol=0.9, cexRow=0.8, ColSideColors = ptcol, 
        hclustfun = hc.avr)

dat2.trans <- t(dat2.matrix)
hc.complete2 = function(dat2.trans) hclust(dat2.trans, method = "complete")
hc.single2 = function(dat2.trans) hclust(dat2.trans, method = "single")
hc.avr2 = function(dat2.trans) hclust(dat2.trans, method = "average")
heatmap(dat2.trans, col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.complete2)
heatmap(dat2.trans, col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.single2)
heatmap(t(dat2.trans), col=hmcol, cexCol=0.9, cexRow=0.9, ColSideColors = ptcol, 
        hclustfun = hc.avr2)

heatmap.2(t(dat2.trans), trace="none", col=hmcol, cexCol=0.9, cexRow=0.8, ColSideColors = ptcol, 
        hclustfun = hc.complete2)
heatmap.2(dat2.trans, trace="none", col=hmcol, cexCol=0.9, cexRow=0.9, RowSideColors = ptcol, 
        hclustfun = hc.single2)
heatmap.2(t(dat2.trans), trace="none", col=hmcol, cexCol=0.9, cexRow=0.8, ColSideColors = ptcol, 
          hclustfun = hc.avr2)

## correlation check
LGALS1 <- as.numeric(dat2[3, 2:134])
LGALS3 <- as.numeric(dat2[1, 2:134])
LGALS1.snp <- as.numeric(dat2[3, 2:7])
LGALS1.nosnp <- as.numeric(dat2[3, 8:134])
LGALS3.snp <- as.numeric(dat2[1, 2:7])
LGALS3.nosnp <- as.numeric(dat2[1, 8:134])

cor.test(LGALS1, LGALS3)
cor.test(LGALS1.nosnp, LGALS3.nosnp)
cor.test(LGALS1.snp, LGALS3.snp)

fit1 <- lm(LGALS1 ~ LGALS3)
summary(fit1)
fit2 <- lm(LGALS1.nosnp ~ LGALS3.nosnp)
summary(fit2)
fit3 <- lm(LGALS1.snp ~ LGALS3.snp)
summary(fit3)

par(mfrow=c(1,3))
plot(LGALS3, LGALS1, ylim=c(0,700), xlim=c(0,90), main="All Samples")
abline(fit1, col="red")
plot(LGALS3.nosnp, LGALS1.nosnp, ylim=c(0,700), xlim=c(0,90), 
     xlab="LGALS3", ylab="LGALS1", main="No SNP")
abline(fit2, col="red")
plot(LGALS3.snp, LGALS1.snp, ylim=c(0,700), xlim=c(0,90), 
     xlab="LGALS3", ylab="LGALS1", main="SNP")
abline(fit3, col="red")

## K-means clustering
km.out <- kmeans(dat2.matrix, 4)
sort(km.out$cluster)
par(mfrow=c(1,1))
plot(dat2.matrix, col=km.out$cluster, cex=2, pch=1, lwd=2)

## PCA
pca.out = prcomp(dat2.trans, scale = TRUE)
pca.out
biplot(pca.out, scale = 0, cex = 0.6)

pca.out2 = prcomp(dat3.trans, scale = TRUE)
pca.out2
biplot(pca.out2, scale = 0, cex = 0.6)
