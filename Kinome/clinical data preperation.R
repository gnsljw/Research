clinical <- read.delim("PTC kinome data clinical info.txt")
head(clinical)

ID <- clinical$sampleID
Age <- clinical$oper_age
Gender <- clinical$gender
Variants <- clinical$Variants
PTC_size <- factor(clinical$PTC_size, levels=c(1,2), labels=c("PTMC","PTC"))
Max_size <- clinical$size_max
Node_meta <- factor(clinical$lymphnode, levels=c(0,1), labels=c("No", "Yes"))
ETE <- factor(clinical$extrathyroid_inv, levels=c(0,1), labels=c("No", "Yes"))
Recurrence <- factor(clinical$recur, levels=c(0,1), labels=c("No", "Yes"))
T_stage <- clinical$T
N_stage <- clinical$N
M_stage <- clinical$M
AJCC_stage <- factor(clinical$Stage_TNM, levels=c("I", "II", "III", "Iva", "Ivb"), 
                     labels=c("I","II","III","IVa", "IVb"))
MACIS <- clinical$MACIS
MACIS_level <- factor(clinical$MACIS_3gr, levels=c("low (<6)", "intermediate (6-7)", "High (>7)"), 
                      labels=c("Low","Intermediate","High"))
ATA_risk <- factor(clinical$ATA_risk, levels=c("low","intermediate","High"), labels=c("Low","Intermediate", "High"))
BRAF <- factor(clinical$BRAF, levels=c(0,1), labels=c("No","Yes"))
MST1R <- factor(clinical$MST1R, levels=c(0,1), labels=c("No","Yes"))
NF1 <- factor(clinical$NF1, levels=c(0,1), labels=c("No","Yes"))
TP53 <- factor(clinical$TP53, levels=c(0,1), labels=c("No","Yes"))
WEE2 <- factor(clinical$WEE2, levels=c(0,1), labels=c("No","Yes"))
NRAS <- factor(clinical$NRAS, levels=c(0,1), labels=c("No","Yes"))
HRAS <- factor(clinical$HRAS, levels=c(0,1), labels=c("No","Yes"))
DAPK1 <- factor(clinical$DAPK1, levels=c(0,1), labels=c("No","Yes"))
EIF2AK4 <- factor(clinical$EIF2AK4, levels=c(0,1), labels=c("No","Yes"))
IDH2 <- factor(clinical$IDH2, levels=c(0,1), labels=c("No","Yes"))
MAP3K3 <- factor(clinical$MAP3K3, levels=c(0,1), labels=c("No","Yes"))
KRAS <- factor(clinical$KRAS, levels=c(0,1), labels=c("No","Yes"))
CHEK2 <- factor(clinical$CHEK2, levels=c(0,1), labels=c("No","Yes"))
PTEN <- factor(clinical$PTEN, levels=c(0,1), labels=c("No","Yes"))
ATM <- factor(clinical$ATM, levels=c(0,1), labels=c("No","Yes"))
RB1 <- factor(clinical$RB1, levels=c(0,1), labels=c("No","Yes"))
APC <- factor(clinical$APC, levels=c(0,1), labels=c("No","Yes"))
KIT <- factor(clinical$KIT, levels=c(0,1), labels=c("No","Yes"))
MET <- factor(clinical$MET, levels=c(0,1), labels=c("No","Yes"))
PIK3CA <- factor(clinical$PIK3CA, levels=c(0,1), labels=c("No","Yes"))
cMET.N375S <- factor(clinical$cMET.N375S, levels=c(0,1), labels=c("No","Yes"))
nonSyn <- clinical$nonSyn
Syn <- clinical$syn
Splicing <- clinical$splicing
Stopgain <- clinical$stopgain
Stoploss <- clinical$stoploss
UTR <- clinical$UTR
Nonsilent <- clinical$nonslient

dat <- data.frame(ID, Age, Gender, Variants, PTC_size, Max_size, Node_meta, ETE, Recurrence, 
                  T_stage, N_stage, M_stage, AJCC_stage, MACIS, MACIS_level, ATA_risk, 
                  BRAF, NRAS, KRAS, HRAS, MST1R, NF1, TP53, WEE2, 
                  DAPK1,EIF2AK4, IDH2, MAP3K3, CHEK2,PTEN, 
                  ATM, RB1, APC, KIT, PIK3CA, MET, cMET.N375S,
                  nonSyn, Syn, Splicing, Stopgain, Stoploss, UTR, Nonsilent) 

NRAS <- unclass(dat$NRAS)
HRAS <- unclass(dat$HRAS)
KRAS <- unclass(dat$KRAS)
RAS <- c(NRAS + HRAS + KRAS)
RAS <- factor(RAS, levels=c(3,4), labels=c("No", "Yes"))
summary(RAS)
RAS <- RAS
dat1 <- dat[,c(1:17)]
dat2 <- dat[,c(18:44)]
dat <- cbind(dat1, RAS, dat2)
names(dat)

saveRDS (dat, "kinome clinical info.RDS")
write.csv(dat, "kinome clinical info.csv")
