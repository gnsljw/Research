#=================#
# Library Loading #
#=================#
library(RColorBrewer); library(epiR); library(epitools); library(survival)
library(ggplot2); library(GGally); library(gmodels); library(gridExtra)

#==================#
# Data prepreation #
#==================#

clinical <- read.csv("D://Work/GLUT/clinicalwithras.csv")
summary(clinical)

clinical$Race <- factor(clinical$Race, levels=c("White","Black","Asian","Unknown"))
clinical$Gender <- factor(clinical$Gender, levels=c("Male","Female"))
clinical$Histologic_Dx <- factor(clinical$Histologic_Dx, levels=c("Classic","FV","TCV","OCV","Mix","Other"))
clinical$ETE <- factor(clinical$ETE, levels=c("No","Minimal","Mod/Adv"))
clinical$Stage.binary <- factor(clinical$Stage.binary, levels=c("Early","Advanced"))

Age.binary <- ifelse(clinical$Age_Dx < 45, "<45", "กร45")
clinical$Age.binary <- factor(Age.binary, levels=c("<45","กร45"), labels=c("<45","กร45"))

summary(clinical)
saveRDS(clinical, "clinical.RDS")

clinical <- readRDS("D://Work/GLUT/clinical.RDS")

#==================================#
# Preperation of survival analysis #
#==================================#

Death <- c(rep("FALSE", 486))
Death <- clinical$Vital_Status == "Dead"
summary(Death)

Recurr <- c(rep("FALSE", 486))
Recurr <- clinical$Recurrence == "YES"
summary(Recurr)

FW_Days_Death <- clinical$FW_Days_Death
FW_Days_Recurr <- clinical$FW_Days_Recurr

FW_Months_Death <- round(FW_Days_Death/30, 2)
FW_Months_Recurr <- round(FW_Days_Recurr/30, 2)

summary(FW_Days_Death); summary(FW_Days_Recurr); summary(FW_Months_Death); summary(FW_Months_Recurr)

s1 <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Stage)
p_value1 <- round((1 - pchisq(s1$chisq, length(s1$n) - 1)), 3)
p_value1
s1.fit <- survfit(Surv(FW_Months_Death, Death) ~ clinical$Stage, type="kaplan-meier")
g1 <- ggsurv(s1.fit, xlab="Follow Up Months", main = paste("Survival Curve According to Stage","1.46e-09",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="Stage")
g1

table(Death, Age.binary)
s.age <- survdiff(Surv(FW_Months_Death, Death) ~ Age.binary)
s.age
s.age.fit <- survfit(Surv(FW_Months_Death, Death) ~ Age.binary, type="kaplan-meier")
g2 <- ggsurv(s.age.fit, xlab="Follow Up Months", main = "Survival Curve According to Age\np=3.16e-05") +
        guides(linetype=F) + scale_colour_discrete(name="Age Group")
g2

s.gender <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Gender)
s.gender
s.Race <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Race)
s.Race
s.histol <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Histologic_Dx)
s.histol
s.hashi <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Hashimoto)
s.hashi
s.hashi.fit <- survfit(Surv(FW_Months_Death, Death) ~ clinical$Hashimoto, type="kaplan-meier")
g5 <- ggsurv(s.hashi.fit, xlab="Follow Up Months", main = "Survival Curve According to Hashimoto\np=0.0596") +
        guides(linetype=F) + scale_colour_discrete(name="Hashimoto")
g5

s.braf <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$BRAF_Mutation)
s.braf

s.ras <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$RAS_Mutation)
s.ras

s.ete <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$ETE)
s.ete
s.ete.fit <- survfit(Surv(FW_Months_Death, Death) ~ clinical$ETE, type="kaplan-meier")
g3 <- ggsurv(s.ete.fit, xlab="Follow Up Months", main = "Survival Curve According to ETE\np=1.64e-07") +
        guides(linetype=F) + scale_colour_discrete(name="ETE")
g3

s.t <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$T_Stage)
s.t
s.t.fit <- survfit(Surv(FW_Months_Death, Death) ~ clinical$T_Stage, type="kaplan-meier")
g4 <- ggsurv(s.t.fit, xlab="Follow Up Months", main = "Survival Curve According to T Stage\np=6.61e-07") +
        guides(linetype=F) + scale_colour_discrete(name="T Stage")
g4

s.n <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$N_Stage)
s.n

s.recurr <- survdiff(Surv(FW_Months_Death, Death) ~ clinical$Recurrence)
s.recurr

grid.arrange(g2,g3,g4,g1)


### Recurrence Analysis with clinical variables
r1 <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Stage)
p_value1 <- round((1 - pchisq(r1$chisq, length(r1$n) - 1)), 3)
p_value1
r1
r1.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ clinical$Stage, type="kaplan-meier")
gr1 <- ggsurv(r1.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to Stage","0.17",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="Stage")
gr1

r.age <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ Age.binary)
r.age

r.gender <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Gender)
r.gender
r.Race <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Race)
r.Race
r.histol <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Histologic_Dx)
r.histol
r.hashi <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Hashimoto)
r.hashi
r.hashi.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ clinical$Hashimoto, type="kaplan-meier")
gr5 <- ggsurv(r.hashi.fit, xlab="Follow Up Months", main = "Survival Curve According to Hashimoto") +
        guides(linetype=F) + scale_colour_discrete(name="Hashimoto")
gr5

r.braf <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$BRAF_Mutation)
r.braf

r.ras <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$RAS_Mutation)
r.ras

r.ete <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$ETE)
r.ete

r.t <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$T_Stage)
r.t

r.n <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$N_Stage)
r.n

r.surv <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ clinical$Vital_Status)
r.surv

### GLUT gene analysis
glutGenes <- read.csv("D://Work/GLUT/Target Genes.csv", header=T)
dat1 <- merge(clinical, glutGenes, by=c("ID"))
saveRDS(dat1, "dat1.rds")
names(glutGenes)[-1]


### Univariable Analysis
univariate <- function(i) {
        dat <- readRDS("dat1.rds")
        a <- round(summary(lm(as.formula(dat[,i] ~ Age_Dx), data=dat))$coefficients, 5)
        b <- round(summary(lm(as.formula(dat[,i] ~ Gender), data=dat))$coefficients, 5)
        c <- round(summary(lm(as.formula(dat[,i] ~ Race), data=dat))$coefficients, 5)
        d <- round(summary(lm(as.formula(dat[,i] ~ Histologic_Dx), data=dat))$coefficients, 5)
        e <- round(summary(lm(as.formula(dat[,i] ~ Hashimoto), data=dat))$coefficients, 5)
        f <- round(summary(lm(as.formula(dat[,i] ~ BRAF_Mutation), data=dat))$coefficients, 5)
        g <- round(summary(lm(as.formula(dat[,i] ~ RAS_Mutation), data=dat))$coefficients, 5)
        h <- round(summary(lm(as.formula(dat[,i] ~ Size), data=dat))$coefficients, 5)
        j <- round(summary(lm(as.formula(dat[,i] ~ Metastatic_LN), data=dat))$coefficients, 5) 
        k <- round(summary(lm(as.formula(dat[,i] ~ ETE), data=dat))$coefficients, 5)
        l <- round(summary(lm(as.formula(dat[,i] ~ T_Stage), data=dat))$coefficients, 5)
        m <- round(summary(lm(as.formula(dat[,i] ~ N_Stage), data=dat))$coefficients, 5)
        n <- round(summary(lm(as.formula(dat[,i] ~ M_Stage), data=dat))$coefficients, 5)
        o <- round(summary(lm(as.formula(dat[,i] ~ Stage.binary), data=dat))$coefficients, 5)
        p <- round(summary(lm(as.formula(dat[,i] ~ Recurrence), data=dat))$coefficients, 5)
        q <- round(summary(lm(as.formula(dat[,i] ~ Vital_Status), data=dat))$coefficients, 5)
        
        result <- rbind(a,b,c,d,e,f,g,h,j,k,l,m,n,o,p,q)
        return(result)
}

result.slm <- list()
for (i in 25:53){
        result.slm <- append(result.slm, list(univariate (i)))
}

names(result.slm) <- as.character(names(dat1[25:53]))
result.slm[[1]]
sink("D://Work/GLUT/Univatiate Linear Regression Result.txt")
result.slm
sink()

### Multivariable Analysis
## SLC2A1
SLC2A1.full <- glm(SLC2A1 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A1.full)
SLC2A1.reduced <- glm(SLC2A1 ~  Age_Dx + Histologic_Dx + BRAF_Mutation + 
                              M_Stage + Stage.binary + Vital_Status, data = dat1)
summary(SLC2A1.reduced)
SLC2A1 <- round(summary(SLC2A1.reduced)$coeff, 5)

## SLC2A2
SLC2A2.full <- glm(SLC2A2 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A2.full)
SLC2A2.reduced <- glm(formula = SLC2A2 ~ 1, data = dat1)
summary(SLC2A2.reduced)
SLC2A2 <- round(summary(SLC2A2.reduced)$coeff, 5)

## SLC2A3
SLC2A3.full <- glm(SLC2A3 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A3.full)
SLC2A3.reduced <- glm(formula = SLC2A3 ~ Gender + Race + Histologic_Dx + Hashimoto + 
                              BRAF_Mutation + RAS_Mutation + ETE + N_Stage + M_Stage + 
                              Stage.binary + Vital_Status, data = dat1)
summary(SLC2A3.reduced)
SLC2A3 <- round(summary(SLC2A3.reduced)$coeff, 5)

## SLC2A4
SLC2A4.full <- glm(SLC2A4 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation +Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A4.full)
SLC2A4.reduced <- glm(formula = SLC2A4 ~ Gender + Histologic_Dx + BRAF_Mutation + 
                              RAS_Mutation + N_Stage, data = dat1)
summary(SLC2A4.reduced)
SLC2A4 <- round(summary(SLC2A4.reduced)$coeff, 5)

## SLC2A5
SLC2A5.full <- glm(SLC2A5 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation +Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A5.full)
SLC2A5.reduced <- glm(formula = SLC2A5 ~ Histologic_Dx + Hashimoto + RAS_Mutation + 
                              N_Stage, data = dat1)
summary(SLC2A5.reduced)
SLC2A5 <- round(summary(SLC2A5.reduced)$coeff, 5)

## SLC2A6
SLC2A6.full <- glm(SLC2A6 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A6.full)
SLC2A6.reduced <- glm(formula = SLC2A6 ~ Hashimoto + RAS_Mutation + ETE + M_Stage + 
                              Stage.binary + Vital_Status, data = dat1)
summary(SLC2A6.reduced)
SLC2A6 <- round(summary(SLC2A6.reduced)$coeff, 5)

## SLC2A7
SLC2A7.full <- glm(SLC2A7 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A7.full)
SLC2A7.reduced <- glm(formula = SLC2A7 ~ RAS_Mutation, data = dat1)
summary(SLC2A7.reduced)
SLC2A7 <- round(summary(SLC2A7.reduced)$coeff, 5)

## SLC2A8
SLC2A8.full <- glm(SLC2A8 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A8.full)
SLC2A8.reduced <- glm(formula = SLC2A8 ~ Age_Dx + Gender + Histologic_Dx + BRAF_Mutation + 
                              ETE + N_Stage, data = dat1)
summary(SLC2A8.reduced)
SLC2A8 <- round(summary(SLC2A8.reduced)$coeff, 5)

## SLC2A9
SLC2A9.full <- glm(SLC2A9 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A9.full)
SLC2A9.reduced <- glm(formula = SLC2A9 ~ Age_Dx + Gender + Hashimoto + BRAF_Mutation + 
                              ETE + Stage.binary + Vital_Status, data = dat1)
summary(SLC2A9.reduced)
SLC2A9 <- round(summary(SLC2A9.reduced)$coeff, 5)

## SLC2A10
SLC2A10.full <- glm(SLC2A10 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A10.full)
SLC2A10.reduced <- glm(formula = SLC2A10 ~ Histologic_Dx + RAS_Mutation + Size + 
                               N_Stage + M_Stage, data = dat1)
summary(SLC2A10.reduced)
SLC2A10 <- round(summary(SLC2A10.reduced)$coeff, 5)

## SLC2A11
SLC2A11.full <- glm(SLC2A11 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A11.full)
SLC2A11.reduced <- glm(formula = SLC2A11 ~ Hashimoto + BRAF_Mutation + RAS_Mutation + 
                               Metastatic_LN, data = dat1)
summary(SLC2A11.reduced)
SLC2A11 <- round(summary(SLC2A11.reduced)$coeff, 5)

## SLC2A12
SLC2A12.full <- glm(SLC2A12 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A12.full)
SLC2A12.reduced <- glm(formula = SLC2A12 ~ Age_Dx + Race + Hashimoto + BRAF_Mutation, 
                       data = dat1)
summary(SLC2A12.reduced)
SLC2A12 <- round(summary(SLC2A12.reduced)$coeff, 5)

## SLC2A13
SLC2A13.full <- glm(SLC2A13 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A13.full)
SLC2A13.reduced <- glm(formula = SLC2A13 ~ Gender + Histologic_Dx + Hashimoto + 
                               BRAF_Mutation + N_Stage, data = dat1)
summary(SLC2A13.reduced)
SLC2A13 <- round(summary(SLC2A13.reduced)$coeff, 5)

## SLC2A14
SLC2A14.full <- glm(SLC2A14 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC2A14.full)
SLC2A14.reduced <- glm(formula = SLC2A14 ~ Gender + Race + Histologic_Dx + Hashimoto + 
                               BRAF_Mutation + RAS_Mutation + ETE + N_Stage + M_Stage + 
                               Stage.binary + Vital_Status, data = dat1)
summary(SLC2A14.reduced)
SLC2A14 <- round(summary(SLC2A14.reduced)$coeff, 5)

## SLC25A4
SLC25A4.full <- glm(SLC25A4 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC25A4.full)
SLC25A4.reduced <- glm(formula = SLC25A4 ~ Age_Dx + Histologic_Dx + Hashimoto + 
                               BRAF_Mutation + Metastatic_LN + M_Stage + Stage.binary, data = dat1)
summary(SLC25A4.reduced)
SLC25A4 <- round(summary(SLC25A4.reduced)$coeff, 5)

## SLC26A4
SLC26A4.full <- glm(SLC26A4 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                            BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                            M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC26A4.full)
SLC26A4.reduced <- glm(formula = SLC26A4 ~ Age_Dx + Histologic_Dx + BRAF_Mutation + 
                               Size + ETE + N_Stage + M_Stage, data = dat1)
summary(SLC26A4.reduced)
SLC26A4 <- round(summary(SLC26A4.reduced)$coeff, 5)

## SLC5A5
SLC5A5.full <- glm(SLC5A5 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                           BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                           M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(SLC5A5.full)
SLC5A5.reduced <- glm(formula = SLC5A5 ~ Age_Dx + Hashimoto + BRAF_Mutation, data = dat1)
summary(SLC5A5.reduced)
SLC5A5 <- round(summary(SLC5A5.reduced)$coeff, 5)

## HK1
HK1.full <- glm(HK1 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                        BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                        M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(HK1.full)
HK1.reduced <- glm(formula = HK1 ~ Histologic_Dx + Hashimoto + BRAF_Mutation + 
                           Size + N_Stage + M_Stage + Vital_Status, data = dat1)
summary(HK1.reduced)
HK1 <- round(summary(HK1.reduced)$coeff, 5)

## HK2
HK2.full <- glm(HK2 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                        BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                        M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(HK2.full)
HK2.reduced <- glm(formula = HK2 ~ Race + Hashimoto + BRAF_Mutation + Size + 
                           Metastatic_LN + Recurrence, data = dat1)
summary(HK2.reduced)
HK2 <- round(summary(HK2.reduced)$coeff, 5)

## HK3
HK3.full <- glm(HK3 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                        BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                        M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(HK3.full)
HK3.reduced <- glm(formula = HK3 ~ Hashimoto + BRAF_Mutation + RAS_Mutation + 
                           N_Stage + Vital_Status, data = dat1)
summary(HK3.reduced)
HK3 <- round(summary(HK3.reduced)$coeff, 5)

## HKDC1
HKDC1.full <- glm(HKDC1 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(HKDC1.full)
HKDC1.reduced <- glm(formula = HKDC1 ~ Histologic_Dx + BRAF_Mutation + N_Stage + 
                             M_Stage + Vital_Status, data = dat1)
summary(HKDC1.reduced)
HKDC1 <- round(summary(HKDC1.reduced)$coeff, 5)

## G6PC
G6PC.full <- glm(G6PC ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                         BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                         M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(G6PC.full)
G6PC.reduced <- glm(formula = G6PC ~ Gender + BRAF_Mutation + ETE, data = dat1)
summary(G6PC.reduced)
G6PC <- round(summary(G6PC.reduced)$coeff, 5)

## G6PC2
G6PC2.full <- glm(G6PC2 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(G6PC2.full)
G6PC2.reduced <- glm(formula = G6PC2 ~ BRAF_Mutation, data = dat1)
summary(G6PC2.reduced)
G6PC2 <- round(summary(G6PC2.reduced)$coeff, 5)

## G6PC3
G6PC3.full <- glm(G6PC3 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(G6PC3.full)
G6PC3.reduced <-  glm(formula = G6PC3 ~ Gender + BRAF_Mutation + RAS_Mutation + 
                              Size + ETE + N_Stage + M_Stage, data = dat1)
summary(G6PC3.reduced)
G6PC3 <- round(summary(G6PC3.reduced)$coeff, 5)

## G6PD
G6PD.full <- glm(G6PD ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                         BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                         M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(G6PD.full)
G6PD.reduced <- glm(formula = G6PD ~ Age_Dx + Gender + Histologic_Dx + BRAF_Mutation + 
                            RAS_Mutation + ETE + T_Stage + N_Stage, data = dat1)
summary(G6PD.reduced)
G6PD <- round(summary(G6PD.reduced)$coeff, 5)

## HIF1A
HIF1A.full <- glm(HIF1A ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(HIF1A.full)
HIF1A.reduced <- glm(formula = HIF1A ~ Gender + Race + Histologic_Dx + BRAF_Mutation + 
                             RAS_Mutation + Size + ETE + T_Stage + N_Stage + Recurrence, data = dat1)
summary(HIF1A.reduced)
HIF1A <- round(summary(HIF1A.reduced)$coeff, 5)

## EPAS1
EPAS1.full <- glm(EPAS1 ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(EPAS1.full)
EPAS1.reduced <- glm(formula = EPAS1 ~ Histologic_Dx + BRAF_Mutation + Size + 
                             Metastatic_LN + M_Stage + Stage.binary + Recurrence, data = dat1)
summary(EPAS1.reduced)
EPAS1 <- round(summary(EPAS1.reduced)$coeff, 5)

## EPO
EPO.full <- glm(EPO ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                          BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                          M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(EPO.full)
EPO.reduced <- glm(formula = EPO ~ Gender + BRAF_Mutation + RAS_Mutation + M_Stage, data = dat1)
summary(EPO.reduced)
EPO <- round(summary(EPO.reduced)$coeff, 5)

## EPOR
EPOR.full <- glm(EPOR ~ Age_Dx+Gender+Race+Histologic_Dx+Hashimoto+
                        BRAF_Mutation+RAS_Mutation+Size+Metastatic_LN+ETE+T_Stage+N_Stage+
                        M_Stage+Stage.binary+Recurrence+Vital_Status, data=dat1)
step(EPOR.full)
EPOR.reduced <- glm(formula = EPOR ~ Gender + Histologic_Dx + BRAF_Mutation + 
                            RAS_Mutation + Size + ETE + T_Stage + M_Stage, data = dat1)
summary(EPOR.reduced)
EPOR <- round(summary(EPOR.reduced)$coeff, 5)


## Result
result.mlm <- list()
result.mlm <- list(SLC2A1, SLC2A2, SLC2A3, SLC2A4, SLC2A5, SLC2A6, SLC2A7, SLC2A8, 
                   SLC2A9, SLC2A10, SLC2A11, SLC2A12, SLC2A13, SLC2A14, 
                   SLC25A4, SLC26A4, SLC5A5, HK1, HK2, HK3, HKDC1, G6PC, G6PC2, G6PC3, 
                   G6PD, HIF1A, EPAS1, EPO, EPOR)
names(result.mlm) <- as.character(names(dat1[25:53]))

sink("Multivariate Linear Regression Result.txt")
result.mlm


### Making Heatmap
glut <- read.csv("D://Research/008 Glucose Metabolism/20141229/GLUT genes.csv")

rownames(glut)<-glut[,1]
head(glut)
glut<-glut[,-1]
head(glut)
data_matrix<-t(data.matrix(glut))
head(data_matrix)

library(RColorBrewer)
pal=rev(brewer.pal(7,"Spectral"))

breaks<-c(0, 0.001, 0.01, 0.05, 0.1, 0.5, 1)

#Create layout with 1 row and 2 columns
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(8,1), heights=c(1,1))
par(mar = c(3,15,15,3),oma=c(0.2,0.2,0.2,0.2),mex=0.5) 

image(x=1:nrow(data_matrix),y=1:ncol(data_matrix),
      z=data_matrix,
      xlab="",
      ylab="",
      col=pal[1:(length(breaks)-1)],
      breaks=breaks,
      axes=FALSE)

abline(h=c(1:ncol(data_matrix))+0.5,
       v=c(1:nrow(data_matrix))+0.5,
       col="white",lwd=1,xpd=F)

text(1:nrow(data_matrix), par("usr")[4] + 1, 
     srt = 45, adj = 0, 
     labels = rownames(data_matrix),
     xpd = TRUE, cex=0.85)

axis(side=2,at=1:ncol(data_matrix),
     labels=colnames(data_matrix),
     col="white",las=1, cex.axis=0.85)

text(par("usr")[1]+5, par("usr")[4] + 4,
     "GLUT and Hypoxia Genes", 
     xpd=TRUE,font=2,cex=1.5)

par(mar = c(5,1,4,7)) 
breaks2<-breaks[-length(breaks)]
image(x=1, y=0:length(breaks2),z=t(matrix(breaks2))*1.001,
      col=pal[1:length(breaks)-1],
      axes=FALSE,
      breaks=breaks,
      xlab="", ylab="",
      xaxt="n")

axis(4,at=0:(length(breaks2)-1), labels=breaks2, col="white", las=1)

abline(h=c(1:length(breaks2)),col="white",lwd=2,xpd=F)


### Survival, Cancer Related Genes

mutation <- read.csv("clinicalwith ptc mutation gene.csv")
head(mutation, 2)
mutation$Race <- factor(mutation$Race, levels=c("White","Black","Asian","Unknown"))
mutation$Gender <- factor(mutation$Gender, levels=c("Male","Female"))
mutation$Histologic_Dx <- factor(mutation$Histologic_Dx, levels=c("Classic","FV","TCV","OCV","Mix","Other"))
mutation$ETE <- factor(mutation$ETE, levels=c("No","Minimal","Mod/Adv"))
mutation$Stage.binary <- factor(mutation$Stage.binary, levels=c("Early","Advanced"))

Age.binary <- ifelse(mutation$Age_Dx < 45, "<45", "กร45")
mutation$Age.binary <- factor(Age.binary, levels=c("<45","กร45"), labels=c("<45","กร45"))

Death <- c(rep("FALSE", 486))
Death <- mutation$Vital_Status == "Dead"
summary(Death)

Recurr <- c(rep("FALSE", 486))
Recurr <- mutation$Recurrence == "YES"
summary(Recurr)

FW_Days_Death <- mutation$FW_Days_Death
FW_Days_Recurr <- mutation$FW_Days_Recurr

FW_Months_Death <- round(FW_Days_Death/30, 2)
FW_Months_Recurr <- round(FW_Days_Recurr/30, 2)

summary(FW_Days_Death); summary(FW_Days_Recurr); summary(FW_Months_Death); summary(FW_Months_Recurr)

braf.surv.table <- table(mutation$BRAF_Mutation, mutation$Vital_Status)
braf.surv.table
braf.recurr.table <- table(mutation$BRAF_Mutation, mutation$Recurrence)
braf.recurr.table
CrossTable(braf.surv.table, expected=T, chisq=T, fisher=T, format=c("SPSS"))
CrossTable(braf.recurr.table, expected=T, chisq=T, fisher=T, format=c("SPSS"))

CCDC6.surv.table <- table(mutation$CCDC6, mutation$Vital_Status)
CCDC6.surv.table
CCDC6.recurr.table <- table(mutation$CCDC6, mutation$Recurrence)
CCDC6.recurr.table
CrossTable(CCDC6.surv.table, expected=T, chisq=T, fisher=T, format=c("SPSS"))
CrossTable(CCDC6.recurr.table, expected=T, chisq=T, fisher=T, format=c("SPSS"))

surv.BRAF <- glm(Vital_Status ~ BRAF_Mutation, data=mutation, family=binomial)
summary(surv.BRAF)
recurr.BRAF <- glm(Recurrence ~ BRAF_Mutation, data=mutation, family=binomial)
summary(recurr.BRAF)

Death <- c(rep("FALSE", 486))
Death <- clinical$Vital_Status == "Dead"
summary(Death)

Recurr <- c(rep("FALSE", 486))
Recurr <- clinical$Recurrence == "YES"
summary(Recurr)

FW_Days_Death <- clinical$FW_Days_Death
FW_Days_Recurr <- clinical$FW_Days_Recurr

FW_Months_Death <- round(FW_Days_Death/30, 2)
FW_Months_Recurr <- round(FW_Days_Recurr/30, 2)

s1 <- survdiff(Surv(FW_Months_Death, Death) ~ BRAF_Mutation, data=mutation)
s1
s1.fit <- survfit(Surv(FW_Months_Death, Death) ~ BRAF_Mutation,data=mutation, type="kaplan-meier")
g1 <- ggsurv(s1.fit, xlab="Follow Up Months", main = paste("Survival Curve According to BRAF","p=0.796",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="BRAF Mutation")
g1

r1 <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ BRAF_Mutation, data=mutation)
r1
r1.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ BRAF_Mutation,data=mutation, type="kaplan-meier")
gr1 <- ggsurv(r1.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to BRAF","p=0.225",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="BRAF Mutation")
gr1


#===========================================#
# Survival, ROC Curve of GLUT related Genes #
#===========================================#
### SLC2A1

## Prediction
library(ROCR)
dat1 <- readRDS("dat1.rds")
pred.SLC2A1 <- prediction(dat1$SLC2A1, dat1$Vital_Status)
perf.SLC2A1 <- performance(pred.SLC2A1, "tpr","fpr")
auc.SLC2A1 <- performance(pred.SLC2A1, "auc")@y.values[[1]]
auc.SLC2A1
## ROC curve
roc.DF1 <- data.frame(x=perf.SLC2A1@x.values[[1]], y=perf.SLC2A1@y.values[[1]])
rocr.plot1 <- ggplot(data=roc.DF1, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot1 <- rocr.plot1 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A1,3))), size=4) +
        labs(title="SLC2A1", x="FPR", y="TPR")
rocr.plot1
## Classification, get cutoff value
library(rpart)
SLC2A1.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1$SLC2A1), method="exp", data=dat1)
plot(SLC2A1.rp, branch=0, margin=0.1)
text(SLC2A1.rp, digits=3, use.n=T) # Cutoff : 7.227
## Binary Conversion
SLC2A1.binary <- log(dat1$SLC2A1)
SLC2A1.binary <- ifelse(SLC2A1.binary < 7.227301, 0, 1)
SLC2A1.binary <- factor(SLC2A1.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A1.binary)
## Cross table analysis
SLC2A1.table <- table(dat1$Vital_Status, SLC2A1.binary)
oddsratio(SLC2A1.table)
SLC2A1.OR.xtable <- c(oddsratio(SLC2A1.table)$measure[2,1:3],oddsratio(SLC2A1.table)$p.value[2,2:3])
SLC2A1.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A1.surv <- glm(Vital_Status ~ SLC2A1.binary, data=dat1, family=binomial)
summary(SLC2A1.surv)
summary(SLC2A1.surv)$coeff[,4]
SLC2A1.OR.Logit <- cbind(exp(SLC2A1.surv$coeff), exp(confint(SLC2A1.surv)), summary(SLC2A1.surv)$coeff[,4])
colnames(SLC2A1.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A1.OR.Logit

SLC2A1.result <- list(SLC2A1.OR.xtable, SLC2A1.OR.Logit)
names(SLC2A1.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A1.result
## K-M Curve
SLC2A1.surv <- survdiff(Surv(FW_Months_Death, Death) ~ SLC2A1.binary, data=dat1)
SLC2A1.surv
SLC2A1.fit <- survfit(Surv(FW_Months_Death, Death) ~ SLC2A1.binary,data=dat1, type="kaplan-meier")
g1 <- ggsurv(SLC2A1.fit, xlab="Follow Up Months", main = paste("Survival Curve According to SLC2A1","p=0.105",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A1")
g1

# SLC2A3
## Prediction
pred.SLC2A3 <- prediction(dat1$SLC2A3, dat1$Vital_Status)
perf.SLC2A3 <- performance(pred.SLC2A3, "tpr","fpr")
auc.SLC2A3 <- performance(pred.SLC2A3, "auc")@y.values[[1]]
auc.SLC2A3
## ROC curve
roc.DF2 <- data.frame(x=perf.SLC2A3@x.values[[1]], y=perf.SLC2A3@y.values[[1]])
rocr.plot2 <- ggplot(data=roc.DF2, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot2 <- rocr.plot1 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A3,3))), size=4) +
        labs(title="SLC2A3", x="FPR", y="TPR")
rocr.plot2
## Classification, get cutoff value
SLC2A3.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1$SLC2A3), method="exp", data=dat1)
plot(SLC2A3.rp, branch=0, margin=0.1)
text(SLC2A3.rp, digits=3, use.n=T) # Cutoff : 7.393
## Binary Conversion
SLC2A3.binary <- log(dat1$SLC2A3)
SLC2A3.binary <- ifelse(SLC2A3.binary < 7.39338, 0, 1)
SLC2A3.binary <- factor(SLC2A3.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A3.binary)
## Cross table analysis
SLC2A3.table <- table(dat1$Vital_Status, SLC2A3.binary)
oddsratio(SLC2A3.table)
SLC2A3.OR.xtable <- c(oddsratio(SLC2A3.table)$measure[2,1:3],oddsratio(SLC2A3.table)$p.value[2,2:3])
SLC2A3.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A3.surv <- glm(Vital_Status ~ SLC2A3.binary, data=dat1, family=binomial)
summary(SLC2A3.surv)
SLC2A3.OR.Logit <- cbind(exp(SLC2A3.surv$coeff), exp(confint(SLC2A3.surv)), summary(SLC2A3.surv)$coeff[,4])
colnames(SLC2A3.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A3.OR.Logit

SLC2A3.result <- list(SLC2A3.OR.xtable, SLC2A3.OR.Logit)
names(SLC2A3.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A3.result

## K-M Curve
SLC2A3.surv <- survdiff(Surv(FW_Months_Death, Death) ~ SLC2A3.binary, data=dat1)
SLC2A3.surv
SLC2A3.fit <- survfit(Surv(FW_Months_Death, Death) ~ SLC2A3.binary,data=dat1, type="kaplan-meier")
g2 <- ggsurv(SLC2A3.fit, xlab="Follow Up Months", main = paste("Survival Curve According to SLC2A3","p=0.179",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A3")
g2

## SLC2A14
## Prediction
pred.SLC2A14 <- prediction(dat1$SLC2A14, dat1$Vital_Status)
perf.SLC2A14 <- performance(pred.SLC2A14, "tpr","fpr")
auc.SLC2A14 <- performance(pred.SLC2A14, "auc")@y.values[[1]]
auc.SLC2A14
## ROC curve
roc.DF3 <- data.frame(x=perf.SLC2A14@x.values[[1]], y=perf.SLC2A14@y.values[[1]])
rocr.plot3 <- ggplot(data=roc.DF3, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot3 <- rocr.plot3 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A14,3))), size=4) +
        labs(title="SLC2A14", x="FPR", y="TPR")
rocr.plot3

SLC2A14.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1$SLC2A14), method="exp", data=dat1)
plot(SLC2A14.rp, branch=0, margin=0.1)
text(SLC2A14.rp, digits=3, use.n=T) # Cutoff : 3.174
## Binary Conversion
SLC2A14.binary <- log(dat1$SLC2A14)
SLC2A14.binary <- ifelse(SLC2A14.binary < 3.174006, 0, 1)
SLC2A14.binary <- factor(SLC2A14.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A14.binary)
## Cross table analysis
SLC2A14.table <- table(dat1$Vital_Status, SLC2A14.binary)
oddsratio(SLC2A14.table)
SLC2A14.OR.xtable <- c(oddsratio(SLC2A14.table)$measure[2,1:3],oddsratio(SLC2A14.table)$p.value[2,2:3])
SLC2A14.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A14.surv <- glm(Vital_Status ~ SLC2A14.binary, data=dat1, family=binomial)
summary(SLC2A14.surv)
summary(SLC2A14.surv)$coeff[,4]
SLC2A14.OR.Logit <- cbind(exp(SLC2A14.surv$coeff), exp(confint(SLC2A14.surv)), summary(SLC2A14.surv)$coeff[,4])
colnames(SLC2A14.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A14.OR.Logit

SLC2A14.result <- list(SLC2A14.OR.xtable, SLC2A14.OR.Logit)
names(SLC2A14.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A14.result
## K-M Curve
SLC2A14.surv <- survdiff(Surv(FW_Months_Death, Death) ~ SLC2A14.binary, data=dat1)
SLC2A14.surv
SLC2A14.fit <- survfit(Surv(FW_Months_Death, Death) ~ SLC2A14.binary,data=dat1, type="kaplan-meier")
g3 <- ggsurv(SLC2A14.fit, xlab="Follow Up Months", main = paste("Survival Curve According to SLC2A14","p=0.156",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A14")
g3
grid.arrange(g1,g2,g3)


#=============================================#
# Recurrence, ROC Curve of GLUT related Genes #
#=============================================#
## HIF1A
## Prediction
pred.HIF1A <- prediction(dat1$HIF1A, dat1$Recurrence)
perf.HIF1A <- performance(pred.HIF1A, "tpr","fpr")
auc.HIF1A <- performance(pred.HIF1A, "auc")@y.values[[1]]
auc.HIF1A
## ROC curve
roc.DF4 <- data.frame(x=perf.HIF1A@x.values[[1]], y=perf.HIF1A@y.values[[1]])
rocr.plot4 <- ggplot(data=roc.DF4, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot4 <- rocr.plot4 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.HIF1A,3))), size=4) +
        labs(title="HIF1A", x="FPR", y="TPR")
rocr.plot4

HIF1A.rp <- rpart(Surv(FW_Days_Recurr, Recurrence=="YES") ~ log(dat1$HIF1A), method="exp", data=dat1)
plot(HIF1A.rp, branch=0, margin=0.1)
text(HIF1A.rp, digits=3, use.n=T) # Cutoff : 3.174
## Binary Conversion
HIF1A.binary <- log(dat1$HIF1A)
HIF1A.binary <- ifelse(HIF1A.binary < 8.538076, 0, 1)
HIF1A.binary <- factor(HIF1A.binary, levels = 0:1, labels = c("Low","High"))
summary(HIF1A.binary)
## Cross table analysis
HIF1A.table <- table(dat1$Recurrence, HIF1A.binary)
oddsratio(HIF1A.table)
HIF1A.OR.xtable <- c(oddsratio(HIF1A.table)$measure[2,1:3],oddsratio(HIF1A.table)$p.value[2,2:3])
HIF1A.OR.xtable
## Survival Analysis, Logistic Regression
HIF1A.surv <- glm(Recurrence ~ HIF1A.binary, data=dat1, family=binomial)
summary(HIF1A.surv)
summary(HIF1A.surv)$coeff[,4]
HIF1A.OR.Logit <- cbind(exp(HIF1A.surv$coeff), exp(confint(HIF1A.surv)), summary(HIF1A.surv)$coeff[,4])
colnames(HIF1A.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
HIF1A.OR.Logit

HIF1A.result <- list(HIF1A.OR.xtable, HIF1A.OR.Logit)
names(HIF1A.result) <- c("Cross Table Analysis", "Logistic Regression")
HIF1A.result
## K-M Curve
HIF1A.surv <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ HIF1A.binary, data=dat1)
HIF1A.surv
HIF1A.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ HIF1A.binary,data=dat1, type="kaplan-meier")
g4 <- ggsurv(HIF1A.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to HIF1A","p=0.031",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="HIF1A")
g4


## EPAS1
## Prediction
pred.EPAS1 <- prediction(dat1$EPAS1, dat1$Recurrence)
perf.EPAS1 <- performance(pred.EPAS1, "tpr","fpr")
auc.EPAS1 <- performance(pred.EPAS1, "auc")@y.values[[1]]
auc.EPAS1
## ROC curve
roc.DF5 <- data.frame(x=perf.EPAS1@x.values[[1]], y=perf.EPAS1@y.values[[1]])
rocr.plot5 <- ggplot(data=roc.DF5, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot5 <- rocr.plot5 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.EPAS1,3))), size=4) +
        labs(title="EPAS1", x="FPR", y="TPR")
rocr.plot5

EPAS1.rp <- rpart(Surv(FW_Days_Recurr, Recurrence=="YES") ~ log(dat1$EPAS1), method="exp", data=dat1)
plot(EPAS1.rp, branch=0, margin=0.1)
text(EPAS1.rp, digits=3, use.n=T) # Cutoff : 8.477051
## Binary Conversion
EPAS1.binary <- log(dat1$EPAS1)
EPAS1.binary <- ifelse(EPAS1.binary < 8.477051, 0, 1)
EPAS1.binary <- factor(EPAS1.binary, levels = 0:1, labels = c("Low","High"))
summary(EPAS1.binary)
## Cross table analysis
EPAS1.table <- table(dat1$Recurrence, EPAS1.binary)
oddsratio(EPAS1.table)
EPAS1.OR.xtable <- c(oddsratio(EPAS1.table)$measure[2,1:3],oddsratio(EPAS1.table)$p.value[2,2:3])
EPAS1.OR.xtable
## Survival Analysis, Logistic Regression
EPAS1.surv <- glm(Recurrence ~ EPAS1.binary, data=dat1, family=binomial)
summary(EPAS1.surv)
summary(EPAS1.surv)$coeff[,4]
EPAS1.OR.Logit <- cbind(exp(EPAS1.surv$coeff), exp(confint(EPAS1.surv)), summary(EPAS1.surv)$coeff[,4])
colnames(EPAS1.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
EPAS1.OR.Logit

EPAS1.result <- list(EPAS1.OR.xtable, EPAS1.OR.Logit)
names(EPAS1.result) <- c("Cross Table Analysis", "Logistic Regression")
EPAS1.result
## K-M Curve
EPAS1.surv <- survdiff(Surv(FW_Months_Recurr, Recurr) ~ EPAS1.binary, data=dat1)
EPAS1.surv
EPAS1.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ EPAS1.binary,data=dat1, type="kaplan-meier")
g5 <- ggsurv(EPAS1.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to EPAS1","p=0.040",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="EPAS1")
g5

grid.arrange(g4,g5)

#=========================#
# Effect of BRAF Mutation #
#=========================#
library(dplyr)
dat1.brafmut <- filter(dat1, BRAF_Mutation == "YES")
dat1.brafwild <- filter(dat1, BRAF_Mutation == "NO")

FW_Days_Death <- dat1.brafmut$FW_Days_Death
FW_Days_Recurr <- dat1.brafmut$FW_Days_Recurr

FW_Months_Death <- round(FW_Days_Death/30, 2)
FW_Months_Recurr <- round(FW_Days_Recurr/30, 2)

### SLC2A1 with BRAF mutant group
## Prediction
pred.SLC2A1 <- prediction(dat1.brafmut$SLC2A1, dat1.brafmut$Vital_Status)
perf.SLC2A1 <- performance(pred.SLC2A1, "tpr","fpr")
auc.SLC2A1 <- performance(pred.SLC2A1, "auc")@y.values[[1]]
auc.SLC2A1
## ROC curve
roc.DF1 <- data.frame(x=perf.SLC2A1@x.values[[1]], y=perf.SLC2A1@y.values[[1]])
rocr.plot1 <- ggplot(data=roc.DF1, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot1 <- rocr.plot1 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A1,3))), size=4) +
        labs(title="SLC2A1", x="FPR", y="TPR")
rocr.plot1
## Classification, get cutoff value
library(rpart)
SLC2A1.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1.brafmut$SLC2A1), method="exp", data=dat1.brafmut)
plot(SLC2A1.rp, branch=0, margin=0.1)
text(SLC2A1.rp, digits=3, use.n=T) # Cutoff :6.690549
## Binary Conversion
SLC2A1.binary <- log(dat1.brafmut$SLC2A1)
SLC2A1.binary <- ifelse(SLC2A1.binary < 6.690549, 0, 1)
SLC2A1.binary <- factor(SLC2A1.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A1.binary)
## Cross table analysis
SLC2A1.table <- table(dat1.brafmut$Vital_Status, SLC2A1.binary)
oddsratio(SLC2A1.table)
SLC2A1.OR.xtable <- c(oddsratio(SLC2A1.table)$measure[2,1:3],oddsratio(SLC2A1.table)$p.value[2,2:3])
SLC2A1.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A1.surv <- glm(Vital_Status ~ SLC2A1.binary, data=dat1.brafmut, family=binomial)
summary(SLC2A1.surv)
SLC2A1.OR.Logit <- cbind(exp(SLC2A1.surv$coeff), exp(confint(SLC2A1.surv)), summary(SLC2A1.surv)$coeff[,4])
colnames(SLC2A1.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A1.OR.Logit

SLC2A1.result <- list(SLC2A1.OR.xtable, SLC2A1.OR.Logit)
names(SLC2A1.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A1.result
## K-M Curve
SLC2A1.surv <- survdiff(Surv(FW_Months_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A1.binary, data=dat1.brafmut)
SLC2A1.surv
SLC2A1.fit <- survfit(Surv(FW_Months_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A1.binary,data=dat1, type="kaplan-meier")
g1 <- ggsurv(SLC2A1.fit, xlab="Follow Up Months", main = paste("Survival Curve According to SLC2A1","p=0.105",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A1")
g1

### SLC2A3 with BRAF mutant group
## Prediction
pred.SLC2A3 <- prediction(dat1.brafmut$SLC2A3, dat1.brafmut$Vital_Status)
perf.SLC2A3 <- performance(pred.SLC2A3, "tpr","fpr")
auc.SLC2A3 <- performance(pred.SLC2A3, "auc")@y.values[[1]]
auc.SLC2A3
## ROC curve
roc.DF1 <- data.frame(x=perf.SLC2A3@x.values[[1]], y=perf.SLC2A3@y.values[[1]])
rocr.plot1 <- ggplot(data=roc.DF1, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot1 <- rocr.plot1 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A3,3))), size=4) +
        labs(title="SLC2A3", x="FPR", y="TPR")
rocr.plot1
## Classification, get cutoff value
SLC2A3.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1.brafmut$SLC2A3), method="exp", data=dat1.brafmut)
plot(SLC2A3.rp, branch=0, margin=0.1)
text(SLC2A3.rp, digits=3, use.n=T) # Cutoff :7.716785
## Binary Conversion
SLC2A3.binary <- log(dat1.brafmut$SLC2A3)
SLC2A3.binary <- ifelse(SLC2A3.binary < 7.716785, 0, 1)
SLC2A3.binary <- factor(SLC2A3.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A3.binary)
## Cross table analysis
SLC2A3.table <- table(dat1.brafmut$Vital_Status, SLC2A3.binary)
oddsratio(SLC2A3.table)
SLC2A3.OR.xtable <- c(oddsratio(SLC2A3.table)$measure[2,1:3],oddsratio(SLC2A3.table)$p.value[2,2:3])
SLC2A3.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A3.surv <- glm(Vital_Status ~ SLC2A3.binary, data=dat1.brafmut, family=binomial)
summary(SLC2A3.surv)
SLC2A3.OR.Logit <- cbind(exp(SLC2A3.surv$coeff), exp(confint(SLC2A3.surv)), summary(SLC2A3.surv)$coeff[,4])
colnames(SLC2A3.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A3.OR.Logit

SLC2A3.result <- list(SLC2A3.OR.xtable, SLC2A3.OR.Logit)
names(SLC2A3.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A3.result
## K-M Curve
SLC2A3.surv <- survdiff(Surv(FW_Months_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A3.binary, data=dat1.brafmut)
SLC2A3.surv
SLC2A3.fit <- survfit(Surv(FW_Months_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A3.binary,data=dat1, type="kaplan-meier")
g2 <- ggsurv(SLC2A3.fit, xlab="Follow Up Month", main = paste("Survival Curve According to SLC2A3 with BRAF Mutation","p=0.014",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A3")
g2

## SLC2A14
### SLC2A14 with BRAF mutant group
## Prediction
pred.SLC2A14 <- prediction(dat1.brafmut$SLC2A14, dat1.brafmut$Vital_Status)
perf.SLC2A14 <- performance(pred.SLC2A14, "tpr","fpr")
auc.SLC2A14 <- performance(pred.SLC2A14, "auc")@y.values[[1]]
auc.SLC2A14
## ROC curve
roc.DF1 <- data.frame(x=perf.SLC2A14@x.values[[1]], y=perf.SLC2A14@y.values[[1]])
rocr.plot1 <- ggplot(data=roc.DF1, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot1 <- rocr.plot1 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.SLC2A14,3))), size=4) +
        labs(title="SLC2A14", x="FPR", y="TPR")
rocr.plot1
## Classification, get cutoff value
SLC2A14.rp <- rpart(Surv(FW_Days_Death, Vital_Status=="Dead") ~ log(dat1.brafmut$SLC2A14), method="exp", data=dat1.brafmut)
plot(SLC2A14.rp, branch=0, margin=0.1)
text(SLC2A14.rp, digits=3, use.n=T) # Cutoff :3.854701
## Binary Conversion
SLC2A14.binary <- log(dat1.brafmut$SLC2A14)
SLC2A14.binary <- ifelse(SLC2A14.binary < 3.854701, 0, 1)
SLC2A14.binary <- factor(SLC2A14.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A14.binary)
## Cross table analysis
SLC2A14.table <- table(dat1.brafmut$Vital_Status, SLC2A14.binary)
oddsratio(SLC2A14.table)
SLC2A14.OR.xtable <- c(oddsratio(SLC2A14.table)$measure[2,1:3],oddsratio(SLC2A14.table)$p.value[2,2:3])
SLC2A14.OR.xtable
## Survival Analysis, Logistic Regression
SLC2A14.surv <- glm(Vital_Status ~ SLC2A14.binary, data=dat1.brafmut, family=binomial)
summary(SLC2A14.surv)
SLC2A14.OR.Logit <- cbind(exp(SLC2A14.surv$coeff), exp(confint(SLC2A14.surv)), summary(SLC2A14.surv)$coeff[,4])
colnames(SLC2A14.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
SLC2A14.OR.Logit

SLC2A14.result <- list(SLC2A14.OR.xtable, SLC2A14.OR.Logit)
names(SLC2A14.result) <- c("Cross Table Analysis", "Logistic Regression")
SLC2A14.result
## K-M Curve
SLC2A14.surv <- survdiff(Surv(dat1.brafmut$FW_Days_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A14.binary, data=dat1.brafmut)
SLC2A14.surv
SLC2A14.fit <- survfit(Surv(FW_Months_Death, dat1.brafmut$Vital_Status == "Dead") ~ SLC2A14.binary,data=dat1.brafmut, type="kaplan-meier")
g3 <- ggsurv(SLC2A14.fit, xlab="Follow Up Months", main = paste("Survival Curve According to SLC2A14","p=0.091",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="SLC2A14")
g3


#=============================================#
# Recurrence, ROC Curve of GLUT related Genes #
#=============================================#
## HIF1A
## Prediction
pred.HIF1A <- prediction(dat1.brafmut$HIF1A, dat1.brafmut$Recurrence)
perf.HIF1A <- performance(pred.HIF1A, "tpr","fpr")
auc.HIF1A <- performance(pred.HIF1A, "auc")@y.values[[1]]
auc.HIF1A
## ROC curve
roc.DF4 <- data.frame(x=perf.HIF1A@x.values[[1]], y=perf.HIF1A@y.values[[1]])
rocr.plot4 <- ggplot(data=roc.DF4, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot4 <- rocr.plot4 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.HIF1A,3))), size=4) +
        labs(title="HIF1A", x="FPR", y="TPR")
rocr.plot4

HIF1A.rp <- rpart(Surv(FW_Days_Recurr, Recurrence=="YES") ~ log(dat1.brafmut$HIF1A), method="exp", data=dat1.brafmut)
plot(HIF1A.rp, branch=0, margin=0.1)
text(HIF1A.rp, digits=3, use.n=T) # Cutoff : 8.646985
## Binary Conversion
HIF1A.binary <- log(dat1.brafmut$HIF1A)
HIF1A.binary <- ifelse(HIF1A.binary < 8.646985, 0, 1)
HIF1A.binary <- factor(HIF1A.binary, levels = 0:1, labels = c("Low","High"))
summary(HIF1A.binary)
## Cross table analysis
HIF1A.table <- table(dat1.brafmut$Recurrence, HIF1A.binary)
oddsratio(HIF1A.table)
HIF1A.OR.xtable <- c(oddsratio(HIF1A.table)$measure[2,1:3],oddsratio(HIF1A.table)$p.value[2,2:3])
HIF1A.OR.xtable
## Survival Analysis, Logistic Regression
HIF1A.surv <- glm(Recurrence ~ HIF1A.binary, data=dat1.brafmut, family=binomial)
summary(HIF1A.surv)
summary(HIF1A.surv)$coeff[,4]
HIF1A.OR.Logit <- cbind(exp(HIF1A.surv$coeff), exp(confint(HIF1A.surv)), summary(HIF1A.surv)$coeff[,4])
colnames(HIF1A.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
HIF1A.OR.Logit

HIF1A.result <- list(HIF1A.OR.xtable, HIF1A.OR.Logit)
names(HIF1A.result) <- c("Cross Table Analysis", "Logistic Regression")
HIF1A.result

## K-M Curve
HIF1A.surv <- survdiff(Surv(FW_Months_Recurr, dat1.brafmut$Recurrence=="YES") ~ HIF1A.binary, data=dat1.brafmut)
HIF1A.surv
HIF1A.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ HIF1A.binary,data=dat1, type="kaplan-meier")
g4 <- ggsurv(HIF1A.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to HIF1A","p=0.031",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="HIF1A")
g4


## EPAS1
## Prediction
pred.EPAS1 <- prediction(dat1.brafmut$EPAS1, dat1.brafmut$Recurrence)
perf.EPAS1 <- performance(pred.EPAS1, "tpr","fpr")
auc.EPAS1 <- performance(pred.EPAS1, "auc")@y.values[[1]]
auc.EPAS1
## ROC curve
roc.DF5 <- data.frame(x=perf.EPAS1@x.values[[1]], y=perf.EPAS1@y.values[[1]])
rocr.plot5 <- ggplot(data=roc.DF5, aes(x=x, y=y)) + geom_path(size=1)
rocr.plot5 <- rocr.plot5 + geom_text(aes(x=1, y=0, hjust=1, vjust=0,label=paste(sep="", "AUC = ", round(auc.EPAS1,3))), size=4) +
        labs(title="EPAS1", x="FPR", y="TPR")
rocr.plot5

EPAS1.rp <- rpart(Surv(FW_Days_Recurr, Recurrence=="YES") ~ log(dat1.brafmut$EPAS1), method="exp", data=dat1.brafmut)
plot(EPAS1.rp, branch=0, margin=0.1)
text(EPAS1.rp, digits=3, use.n=T) # Cutoff : 9.011914
## Binary Conversion
EPAS1.binary <- log(dat1.brafmut$EPAS1)
EPAS1.binary <- ifelse(EPAS1.binary < 9.011914, 0, 1)
EPAS1.binary <- factor(EPAS1.binary, levels = 0:1, labels = c("Low","High"))
summary(EPAS1.binary)
## Cross table analysis
EPAS1.table <- table(dat1.brafmut$Recurrence, EPAS1.binary)
oddsratio(EPAS1.table)
EPAS1.OR.xtable <- c(oddsratio(EPAS1.table)$measure[2,1:3],oddsratio(EPAS1.table)$p.value[2,2:3])
EPAS1.OR.xtable
## Survival Analysis, Logistic Regression
EPAS1.surv <- glm(Recurrence ~ EPAS1.binary, data=dat1.brafmut, family=binomial)
summary(EPAS1.surv)
summary(EPAS1.surv)$coeff[,4]
EPAS1.OR.Logit <- cbind(exp(EPAS1.surv$coeff), exp(confint(EPAS1.surv)), summary(EPAS1.surv)$coeff[,4])
colnames(EPAS1.OR.Logit) <- c("OR", "2.5%", "97.5%", "p-value")
EPAS1.OR.Logit

EPAS1.result <- list(EPAS1.OR.xtable, EPAS1.OR.Logit)
names(EPAS1.result) <- c("Cross Table Analysis", "Logistic Regression")
EPAS1.result
## K-M Curve
EPAS1.surv <- survdiff(Surv(FW_Months_Recurr, dat1.brafmut$Recurrence=="YES") ~ EPAS1.binary, data=dat1.brafmut)
EPAS1.surv
EPAS1.fit <- survfit(Surv(FW_Months_Recurr, Recurr) ~ EPAS1.binary,data=dat1, type="kaplan-meier")
g5 <- ggsurv(EPAS1.fit, xlab="Follow Up Months", main = paste("Recurrence Curve According to EPAS1","p=0.040",sep="\n")) +
        guides(linetype=F) + scale_colour_discrete(name="EPAS1")
g5