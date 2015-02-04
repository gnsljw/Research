dat <- readRDS("D://Research/014 Kinome, PTMC/Work/kinome clinical info.RDS")
names(dat)
names(dat)[38] <- "cMET"
names(dat)
dat <- dat[-c(111), ]  # removing NA in PTC/PTMC

library(ggplot2); library(gridExtra); library(dplyr); library(reshape2)

dat.mut <- cbind(dat[,5], dat[,17:38])
names(dat.mut)[1] <- "PTC_size"
summary(dat.mut)

dat.long <- melt(dat.mut, id.var = "PTC_size", variable.name = "Genes", value.name="Mutation")
head(dat.long)

dat.long <- tbl_df(dat.long)

fq <- dat.long %>%
        group_by(PTC_size, Genes, Mutation) %>%
        summarise(count=n())
fq

fq.mut <- filter(fq, Mutation == "Yes")
fq.mut <- select(fq.mut, Genes, count)
fq.mut <- tbl_df(fq.mut)
fq.mut <- arrange(fq.mut, desc(count))
fq.mut$Genes <- factor(fq.mut$Genes, 
                       levels=rev(c("BRAF", "cMET", "RAS", "NRAS", "NF1", "APC", "EIF2AK4", "ATM", 
                                "KRAS", "CHEK2", "MST1R", "WEE2", "DAPK1", "IDH2", "HRAS", "KIT", "PIK3CA", 
                                "TP53", "PTEN", "MET")), 
                       labels=rev(c("BRAF", "cMET", "RAS", "NRAS", "NF1", "APC", "EIF2AK4", "ATM", 
                                "KRAS", "CHEK2", "MST1R", "WEE2", "DAPK1", "IDH2", "HRAS", "KIT", "PIK3CA", 
                                "TP53", "PTEN", "MET")))

## draw barplot by ggplot

ggplot(fq.mut, aes(Genes, count, fill=PTC_size)) + geom_bar(stat="identity", position = "dodge") +
        geom_text(stat="identity", color="black", size=3, hjust=0.5,vjust=0.5, aes(y=count, label=count)) +
        coord_flip()


## draw barplot by base plotting system

mut.freq <- dcast(dat.long, Genes + Mutation ~ PTC_size)
mut.freq <- filter(mut.freq, Mutation == "Yes")
mut.freq <- arrange(mut.freq, desc(PTMC), desc(PTC))
mut.freq <- select(mut.freq, Genes, PTMC, PTC)
mut.freq

PTMC <- mut.freq$PTMC
PTC <- mut.freq$PTC
Genes <- mut.freq$Genes

tbl <- as.table(rbind(PTMC, PTC))
colnames(tbl) <- as.character(Genes)
tbl

b <- barplot(tbl, beside=T, las=2, col=c("salmon", "skyblue"), ylim=c(0,80), 
        legend.text = rownames(tbl))
text(b, tbl + 3, as.character(tbl), col="red", cex=0.6)
box()

tbl.pr <- prop.table(tbl)
colnames(tbl.pr)[3] <- "cMET"
tbl.pr
b <- barplot(tbl.pr, beside=T, las=2, col=c("salmon", "skyblue"), ylim=c(0,0.35),
             legend.text = rownames(tbl.pr))
text(b, tbl.pr + 0.01, as.character(round(tbl.pr, digits=2)), col="red", cex=0.5)
