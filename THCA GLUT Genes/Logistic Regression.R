dat <- readRDS("TCGA GLUT data.rds")

# Survival and Recurrence, Multiple Logistic Regression
surv.all <- glm(Vital_Status ~ SLC2A1+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+
                SLC2A2+SLC2A3+SLC2A4+SLC2A4RG+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+
                SLC25A4+SLC26A4+SLC5A5+HIF1A+HIF1AN+EPAS1+EPO+EPOR+
                HK1+HK2+HK3+HKDC1+G6PC+G6PC2+G6PC3+G6PD+TPO+TSHR+PTEN+PTENP1+BRAF, data=dat, family=binomial)
summary(surv.all)
step(surv.all)
surv.reduced <- glm(Vital_Status ~ SLC2A1 + SLC2A3 + SLC2A4RG + SLC2A5 + HK3 + HKDC1 + PTEN, 
                    family=binomial, data=dat)
summary(surv.reduced)

recurr.all <- glm(Recurrence ~ SLC2A1+SLC2A10+SLC2A11+SLC2A12+SLC2A13+SLC2A14+
                    SLC2A2+SLC2A3+SLC2A4+SLC2A4RG+SLC2A5+SLC2A6+SLC2A7+SLC2A8+SLC2A9+
                    SLC25A4+SLC26A4+SLC5A5+HIF1A+HIF1AN+EPAS1+EPO+EPOR+
                    HK1+HK2+HK3+HKDC1+G6PC+G6PC2+G6PC3+G6PD+TPO+TSHR+PTEN+PTENP1+BRAF, data=dat, family=binomial)
summary(recurr.all)
step(recurr.all)
recurr.reduced <- glm(Recurrence ~ SLC2A12 + SLC2A2 + SLC2A4 + SLC2A8 + 
                      SLC25A4 + SLC26A4 + SLC5A5 + EPAS1 + HK1 + HK2 + G6PD, 
                    family=binomial, data=dat)
summary(recurr.reduced)

# Survival, Univariate Logistic Regression
surv.SLC2A1 <- glm(Vital_Status ~ log(SLC2A1), data=dat, family=binomial)
summary(surv.SLC2A1)
exp(surv.SLC2A1$coeff) # odds ratio
exp(confint(surv.SLC2A1))  # Conf. interval
par(mfrow=c(1,2))
plot(log(dat$SLC2A1), surv.SLC2A1$fitted, ylim=c(0,0.5), yaxt="n", pch=20, col="blue", 
     xlab="log(SLC2A1)", ylab="Survival",main = "SLC2A1\n p=0.010, OR=5.342(1.473-19.604)")  
axis(2, at=c(0,0.5), lab=c("Alive","Dead"), las=1, cex=0.5)
abline(v=7.227, col="salmon")

surv.SLC2A3 <- glm(Vital_Status ~ log(SLC2A3), data=dat, family=binomial)
summary(surv.SLC2A3)
exp(surv.SLC2A3$coeff) # odds ratio
exp(confint(surv.SLC2A3))  # Conf. interval
plot(log(dat$SLC2A3), surv.SLC2A3$fitted, ylim=c(0,0.5), yaxt="n", pch=20, col="blue", 
     xlab="log(SLC2A3)", ylab="Survival",main = "SLC2A3 \n p=0.001, OR=4.092(1.826-10.016) ") 
axis(2, at=c(0,0.5), lab=c("Alive","Dead"), las=1, cex=0.5)
abline(v=7.393, col="salmon")

# Classificatin Analysis
library(ROCR)
par(mfrow=c(1,2))
pred.SLC2A1 <- prediction(dat$SLC2A1, dat$Vital_Status)
perf.SLC2A1 <- performance(pred.SLC2A1, "tpr","fpr")
auc.SLC2A1 <- performance(pred.SLC2A1, "auc")@y.values[[1]]
auc.SLC2A1
plot(perf.SLC2A1, col="blue")
title("SLC2A1\nAUC=0.676")

pred.SLC2A3 <- prediction(dat$SLC2A3, dat$Vital_Status)
perf.SLC2A3 <- performance(pred.SLC2A3, "tpr","fpr")
auc.SLC2A3 <- performance(pred.SLC2A3, "auc")@y.values[[1]]
auc.SLC2A3
plot(perf.SLC2A3, col="blue")
title("SLC2A3\nAUC=0.761")

pred.BRAF <- prediction(dat$BRAF, dat$Vital_Status)
perf.BRAF <- performance(pred.BRAF, "tpr","fpr")
plot(perf.BRAF)
perf.BRAF <- performance(pred.BRAF, "auc")@y.values[[1]]
perf.BRAF

library(rpart)
library(survival)
par(mfrow=c(1,1))
SLC2A1.rp <- rpart(Surv(FW_Days, Vital_Status=="Dead") ~ log(SLC2A1), method="exp", data=dat)
SLC2A1.rp
plot(SLC2A1.rp, branch=0, margin=0.1)
text(SLC2A1.rp, digits=3, use.n=T)

SLC2A3.rp <- rpart(Surv(FW_Days, Vital_Status=="Dead") ~ log(SLC2A3), method="exp", data=dat)
SLC2A3.rp
plot(SLC2A3.rp, branch=0, margin=0.1)
text(SLC2A3.rp, digits=3, use.n=T)

library(epiR)
library(epitools)
library(gmodels)

SLC2A1.binary <- dat$SLC2A1
SLC2A1.binary <- ifelse(log(SLC2A1.binary) < 7.227301, 0, 1)
SLC2A1.binary <- factor(SLC2A1.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A1.binary)

SLC2A3.binary <- dat$SLC2A3
SLC2A3.binary <- ifelse(log(SLC2A3.binary) < 7.39338, 0, 1)
SLC2A3.binary <- factor(SLC2A3.binary, levels = 0:1, labels = c("Low","High"))
summary(SLC2A3.binary)

GLUT1.surv <- table(dat$Vital_Status, SLC2A1.binary)
oddsratio(GLUT1.surv)

GLUT3.surv <- table(dat$Vital_Status, SLC2A3.binary)
oddsratio(GLUT3.surv)

BRAFMUT.surv <- table(dat$Vital_Status, dat$BRAF_V600E)
oddsratio(BRAFMUT.surv)

# Survival Curve, Kaplan Meier
library(survival)
library(GGally)
library(RColorBrewer)

col = brewer.pal(9,"Set1")

death <- c(rep("FALSE", 486))
death <- dat$Vital_Status == "Dead"
FW_Month <- round(dat$FW_Days/30)

s1 <- survdiff(Surv(FW_Month, death) ~ SLC2A1.binary)
p_value1 <- round(1 - pchisq(s1$chisq, length(s1$n) - 1), 3)
p_value1
s1.fit <- survfit(Surv(FW_Month, death) ~ SLC2A1.binary, type="kaplan-meier")
par(mfrow=c(1,1))
plot(s1.fit, main=paste("SLC2A1","\np=",p_value1, sep=""), xlab="Follow Up Months", 
     col=col[2:1], lty=2:3)
legend("bottomright", lty=1, legend = c("Lower Expression", "Higher Expression"), bty = "n", col=col[2:1])

s3 <- survdiff(Surv(FW_Month, death) ~ SLC2A3.binary)
p_value3 <- round(1 - pchisq(s3$chisq, length(s3$n) - 1), 3)
p_value3
s3.fit <- survfit(Surv(FW_Month, death) ~ SLC2A3.binary, type="kaplan-meier")
plot(s3.fit, main=paste("SLC2A3","\np=",p_value3, sep=""), xlab="Follow Up Months", 
     col=col[2:1], lty=2:3)
legend("bottomright", lty=1, legend = c("Lower Expression", "Higher Expression"),  bty = "n",col=col[2:1])

b3 <- survdiff(Surv(FW_Month, death) ~ dat$BRAF_V600E)
p_value5 <- round(1 - pchisq(b3$chisq, length(b3$n) - 1), 3)
b3.fit <- survfit(Surv(FW_Month, death) ~ dat$BRAF_V600E, type="kaplan-meier")
plot(b3.fit, main=paste("BRAF Mutation","\np=",p_value5, sep=""), xlab="Follow Up Months", 
     col=col[2:1], lty=2:3)
legend("bottomright", lty=1, legend = c("Wild", "Mutant"),bty = "n", col=col[2:1])

# Survival Curve, cox proportional hazard model
surv.SLC2A1 <- coxph(Surv(FW_Month, Vital_Status == "Dead") ~ log(SLC2A1), data=dat)
summary(surv.SLC2A1)
plot(survfit(surv.SLC2A1), main="log(SLC2A1)\nHR=2.530(0.706-9.063), p=0.154", col=col, lty=2:3, 
     xlab="Follow Up Months")

surv.SLC2A3 <- coxph(Surv(FW_Month, Vital_Status == "Dead") ~ log(SLC2A3), data=dat)
summary(surv.SLC2A3)
plot(survfit(surv.SLC2A3), main="log(SLC2A3)\nHR=2.108(0.938-4.739), p=0.071", col=col, lty=2:3, 
     xlab="Follow Up Months")

surv.SLC2A1.3 <- coxph(Surv(FW_Month, Vital_Status == "Dead") ~ log(SLC2A1)*log(SLC2A3), data=dat)
summary(surv.SLC2A1.3)
plot(survfit(surv.SLC2A1.3), main="Interaction, SLC2A1 and SLC2A3\nHR=2.971(0.844-10.451), p=0.090", col=col, lty=2:3, 
     xlab="Follow Up Months")

surv.BRAF <- coxph(Surv(FW_Month, Vital_Status == "Dead") ~ log(BRAF), data=dat)
summary(surv.BRAF)
plot(survfit(surv.BRAF), main="log(BRAF)\nHR=1.07(0.297-3.852), p=0.918", col=col, lty=2:3, 
     xlab="Follow Up Months")

