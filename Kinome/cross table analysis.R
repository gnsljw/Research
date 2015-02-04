dat <- readRDS("kinome clinical info.RDS")
names(dat)

dat <- dat[-c(111), ]  # removing NA in PTC/PTMC

library(ggplot2); library(gridExtra); library(dplyr); library(reshape2)
library(gmodels); library(epiR); library(epitools)

dat.mut <- cbind(dat[,5], dat[,17:38])
names(dat.mut)[1] <- "PTC_size"
summary(dat.mut)
attach(dat.mut)

BRAF.tbl <- table(PTC_size, BRAF)
CrossTable(BRAF.tbl, fisher=T, chisq=T, mcnemar=T, expected=T, sresid=T, format="SPSS")
RAS.tbl <- table(PTC_size, RAS)
CrossTable(RAS.tbl, fisher=T, chisq=T, mcnemar=T, expected=T, sresid=T, format="SPSS")
cMET.tbl <- table(PTC_size, cMET.N375S)
CrossTable(cMET.tbl, fisher=T, chisq=T, mcnemar=T, expected=T, sresid=T, format="SPSS")
NF1.tbl <- table(PTC_size, NF1)
CrossTable(NF1.tbl, fisher=T, chisq=T, mcnemar=T, expected=T, sresid=T, format="SPSS")
APC.tbl <- table(PTC_size, APC)
CrossTable(APC.tbl, fisher=T, chisq=T, mcnemar=T, expected=T, sresid=T, format="SPSS")

names(dat)
t.test(dat$Syn ~ dat$PTC_size)
t.test(dat$nonSyn ~ dat$PTC_size)
t.test(dat$Splicing ~ dat$PTC_size)
t.test(dat$Stopgain ~ dat$PTC_size)
t.test(dat$Stoploss ~ dat$PTC_size)
t.test(dat$UTR ~ dat$PTC_size)
t.test(dat$Nonsilent ~ dat$PTC_size)

g1 <- ggplot(dat, aes(PTC_size, nonSyn)) + geom_boxplot() + ggtitle("Non Synonymous") + stat_summary(fun.y=mean, geom="point", shape=5, size=4)
g2 <- ggplot(dat, aes(PTC_size, Splicing)) + geom_boxplot() + ggtitle("Splicing Site")+ stat_summary(fun.y=mean, geom="point", shape=5, size=4)
g3 <- ggplot(dat, aes(PTC_size, Stopgain)) + geom_boxplot() + ggtitle("Stop Gain")+ stat_summary(fun.y=mean, geom="point", shape=5, size=4)
g4 <- ggplot(dat, aes(PTC_size, Stoploss)) + geom_boxplot() + ggtitle("Stop Loss")+ stat_summary(fun.y=mean, geom="point", shape=5, size=4)
g5 <- ggplot(dat, aes(PTC_size, UTR)) + geom_boxplot() + ggtitle("UTR")+ stat_summary(fun.y=mean, geom="point", shape=5, size=4)
g6 <- ggplot(dat, aes(PTC_size, Nonsilent)) + geom_boxplot() + ggtitle("Nonsilent")+ stat_summary(fun.y=mean, geom="point", shape=5, size=4)

grid.arrange(g1, g6, g5, g2, g3, g4, nrow=2)

clinical <- dat[,2:16]
summary(clinical)
clinical.PTMC <- filter(clinical, PTC_size == "PTMC")
summary(clinical.PTMC)
clinical.PTC <- filter(clinical, PTC_size == "PTC")
summary(clinical.PTC)

attach(clinical)
t.test(Age ~ PTC_size)
chisq.test(table(Gender, PTC_size))
chisq.test(table(Variants, PTC_size))
t.test(Max_size ~ PTC_size)
chisq.test(table(Node_meta, PTC_size))
chisq.test(table(ETE, PTC_size))
fisher.test(table(MACIS_level, PTC_size))
fisher.test(table(ATA_risk, PTC_size))
chisq.test(table(Recurrence, PTC_size))
fisher.test(table(T_stage, PTC_size))
fisher.test(table(N_stage, PTC_size))
chisq.test(table(M_stage, PTC_size))
fisher.test(table(AJCC_stage, PTC_size))
