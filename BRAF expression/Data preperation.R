setwd("D://Research/017 BRAF expression in TCGA/20150226/")
library(dplyr); library(epiR); library(epitools); library(gmodels)
clinical <- read.delim("20150223_TCGA_Cli_Path_BRAF_RAS.txt")
BRAF_tumor_count <- read.delim("BRAF tumor.txt")
BRAF_normal_count <- read.delim("BRAF normal.txt")

dat <- merge(clinical, BRAF_tumor_count, by=c("ID"))
head(dat)
str(dat)

ID <- dat$ID
Gender <- dat$Gender 
Gender <- factor(Gender, levels=c("M","F"))
Vital_Status <- dat$Vital_Status
Recurrence <- dat$Recurrence
FU_days_death <- dat$FW_days
FU_days_recurr <- dat$FW_recurr
T_Stage <- dat$T
summary(T_Stage)
T_Stage[T_Stage == "T1b"] <- "T1a"
summary(T_Stage)
T_Stage[T_Stage == "T4b"] <- "T4a"
summary(T_Stage)
T_Stage <- factor(T_Stage, levels=c("T1a","T2","T3","T4a"), labels=c("T1","T2","T3","T4"))
T_Stage
T_Stage_bin <- T_Stage
T_Stage_bin[T_Stage_bin == "T2"] <- "T1"
T_Stage_bin[T_Stage_bin == "T4"] <- "T3"
T_Stage_bin
T_Stage_bin <- factor(T_Stage_bin, levels=c("T1","T3"), labels=c("T1T2","T3T4"))
T_Stage_bin

N_Stage <- dat$N
N_Stage
summary(N_Stage)
N_Stage[N_Stage == "N1a"] <- "N1"
N_Stage[N_Stage == "N1b"] <- "N1"
N_Stage[N_Stage == "Nx"] <- NA
N_Stage_bin <- factor(N_Stage, levels=c("N0","N1"))
N_Stage_bin

N_Stage <- dat$N
N_Stage[N_Stage == "N1"] <- "N1a"
N_Stage[N_Stage == "Nx"] <- NA
N_Stage <- factor(N_Stage, levels=c("N0","N1a","N1b"))
summary(N_Stage)

ETE <- dat$ETE
ETE[ETE=="Extensive"] <- "Yes"
ETE[ETE=="NO"] <- "No"
ETE <- factor(ETE, levels=c("No","Yes"))

Subtype_TCGA <- dat$Subtype_TCGA
summary(Subtype_TCGA)
Subtype_TCGA <- factor(Subtype_TCGA, levels=c("Classical/Usual (Papillary NOS)", 
                                              "Thyroid Papillary Carcinoma - Follicular (>= 99% follicular patterned)", 
                                              "Thyroid Papillary Carcinoma - Tall Cell (>= 50% tall cell features)"), 
                       labels=c("PTC","FVPTC","TCPTC"))
summary(Subtype_TCGA)
Subtype_TCGA <- unclass(Subtype_TCGA)
Subtype_TCGA[is.na(Subtype_TCGA)] <- "Others"
Subtype_TCGA
Subtype_TCGA <- factor(Subtype_TCGA, levels=c(1,2,3,"Others"), labels=c("PTC","FVPTC","TCPTC","Others"))
Subtype_TCGA

Subtype_path <- dat$Subtype
summary(Subtype_path)
Subtype_path <- factor(Subtype_path, levels=c("Conventional PTC", "FVPTC","TCV"), labels=c("PTC","FVPTC","TCPTC"))
summary(Subtype_path)
Subtype_path <- unclass(Subtype_path)
Subtype_path[is.na(Subtype_path)] <- "Others"
Subtype_path <- factor(Subtype_path, levels=c(1,2,3,"Others"), labels=c("PTC","FVPTC","TCPTC","Others"))
Subtype_path

summary(Subtype_TCGA == Subtype_path)

CLT <- dat$Thyroiditis
CLT[CLT=="Chronic inflam"] <- "No"
CLT[CLT=="Unreported"] <- "No"
CLT[CLT=="Thyroiditis"] <- "No"
summary(CLT)
CLT <- factor(CLT, levels=c("No","CLT"), labels=c("No","Yes"))
summary(CLT)

Thyroiditis <- dat$Thyroiditis
summary(Thyroiditis)
Thyroiditis[Thyroiditis == "Chronic inflam"] <- "Thyroiditis"
Thyroiditis[Thyroiditis == "CLT"] <- "Thyroiditis"
Thyroiditis[Thyroiditis == "Unreported"] <- "No"
summary(Thyroiditis)
Thyroiditis <- factor(Thyroiditis, levels=c("No", "Thyroiditis"), labels=c("No","Yes"))
summary(Thyroiditis)

Meta_nodes <- dat$Node_meta
Size <- dat$Size
Size01 <- factor(ifelse(Size <= 1, "<=1cm",">1cm"))
Size02 <- factor(ifelse(Size <= 2, "<=2cm",">2cm"))
Size04 <- factor(ifelse(Size <= 4, "<=4cm",">4cm"))

LN_ratio <- dat$Node_meta / dat$Node_ret
LN_ratio_01 <- factor(ifelse(LN_raio < 0.33, "Low","High"))

M_Stage <- dat$M_Stage_TCGA 
M_Stage[M_Stage=="MX"] <- "M0"
M_Stage <- factor(M_Stage, levels=c("M0","M1"))
summary(M_Stage)

MACIS_score <- dat$MACIS
MACIS_risk1 <- ifelse(MACIS_score <6, "Low", "Int+High")
MACIS_risk2 <- ifelse(MACIS_score <7, "Low+Int","High")

MACIS_risk3 <- unclass(cut(MACIS_score, breaks=c(0,6,7,8,11)))
MACIS_risk3 <- factor(MACIS_risk3, levels=c(1,2,3,4), labels=c("99%","89%","56%","24%"))
summary(MACIS_risk3)

Stage <- dat$STAGE
Stage[Stage == "IVC"] <- "IVA"
Stage <- factor(Stage, levels=c("I","II","III","IVA"), labels=c("I","II","III","IV") )

Stage_01 <- Stage
Stage_01[Stage == "II"] <- "IV"
Stage_01[Stage == "III"] <- "IV"
Stage_01 <- factor(Stage_01, levels=c("I","IV"), labels=c("I","II,III,IV"))

Stage_02 <- Stage
Stage_02[Stage_02 == "II"] <- "I"
Stage_02[Stage_02 == "III"] <- "I"
Stage_02 <- factor(Stage_02, levels=c("I","IV"), labels=c("I,II,III","IV"))

Stage_03 <- Stage
Stage_03[Stage_03 == "II"] <- "I"
Stage_03[Stage_03 == "III"] <- "IV"
Stage_03 <- factor(Stage_03, levels=c("I","IV"), labels=c("I,II","III,IV"))

summary(Stage);summary(Stage_01);summary(Stage_02);summary(Stage_03)

BRAF <- dat$BRAF
BRAF_Mut <- dat$BRAF_Mutation
BRAF_Mut01 <- factor(BRAF_Mut, levels=c("Wild","V600E"), labels=c("Wild","V600E"))
RAS_Mut <- dat$RAS_Mutation
RAS_Mut01 <- dat$RAS_Mut
RAS_Mut01[RAS_Mut01 == "HRAS"] <- "NRAS"; RAS_Mut01[RAS_Mut01 == "KRAS"] <- "NRAS"
RAS_Mut01 <- factor(RAS_Mut01, levels=c("No","NRAS"), labels=c("No","Yes"))
summary(BRAF_Mut); summary(RAS_Mut)
summary(BRAF_Mut01); summary(RAS_Mut01)

Age <- dat$Age
Age_40 <- factor(ifelse(Age<40, "<40",">=40"))
Age_45 <- factor(ifelse(Age<45, "<45",">=45"))
Age_50 <- factor(ifelse(Age<50, "<50",">=50"))
Age_55 <- factor(ifelse(Age<55, "<55",">=55"))
summary(Age_40);summary(Age_45);summary(Age_50);summary(Age_55)

Node_meta <- dat$Node_meta

dat.all <- data.frame(ID,Gender,Age,Age_40,Age_45,Age_50,Age_55, Subtype_TCGA, Subtype_path, 
                      CLT, Thyroiditis, ETE, Size, Size01, Size02, Size04, Node_meta, LN_ratio_01, 
                      Vital_Status, Recurrence, FU_days_death, FU_days_recurr, 
                      T_Stage, T_Stage_bin, N_Stage, N_Stage_bin, M_Stage, Stage, Stage_01, Stage_02, Stage_03, 
                      MACIS_score, MACIS_risk1, MACIS_risk2, MACIS_risk3, 
                      BRAF_Mut, BRAF_Mut01, RAS_Mut, RAS_Mut01, BRAF)