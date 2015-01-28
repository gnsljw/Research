library(dplyr); library(epiR); library(epitools); library(gmodels)
dat1 <- readRDS("D://Work/BRAF exp/brafdata.rds")
summary(dat1)

dat.mutant <- filter(dat1, BRAF_Mutation == "YES") 
summary(dat.mutant)

dat.mutant.old <- filter(dat.mutant, Age.binary == "กร45")
summary(dat.mutant.old)

dat.wild <- filter(dat1, BRAF_Mutation == "NO")
summary(dat.wild)

dat.wild.old <- filter(dat.mutant, Age.binary == "กร45")

### T-test ###

attach(dat1)
summary(dat1)

t.test(BRAF ~ Age.binary)
t.test(BRAF ~ Gender)
t.test(BRAF ~ Hashimoto)
t.test(BRAF ~ ETE.binary)
t.test(BRAF ~ T.binary)
t.test(BRAF ~ N_Stage)
t.test(BRAF ~ Stage.binary2)
t.test(BRAF ~ Recurrence)
t.test(BRAF ~ Vital_Status)

detach(dat1)
attach(dat.mutant)
t.test(BRAF ~ Age.binary)
t.test(BRAF ~ Gender)
t.test(BRAF ~ Hashimoto)
t.test(BRAF ~ ETE.binary)
t.test(BRAF ~ T.binary)
t.test(BRAF ~ N_Stage)
t.test(BRAF ~ Stage.binary2)
t.test(BRAF ~ Recurrence)
t.test(BRAF ~ Vital_Status)

detach(dat.mutant)
attach(dat.mutant.old)
t.test(BRAF ~ Age.binary)
t.test(BRAF ~ Gender)
t.test(BRAF ~ Hashimoto)
t.test(BRAF ~ ETE.binary)
t.test(BRAF ~ T.binary)
t.test(BRAF ~ N_Stage)
t.test(BRAF ~ Stage.binary2)
t.test(BRAF ~ Recurrence)
t.test(BRAF ~ Vital_Status)

detach(dat.mutant.old)
attach(dat.wild)
t.test(BRAF ~ Age.binary)
t.test(BRAF ~ Gender)
t.test(BRAF ~ Hashimoto)
t.test(BRAF ~ ETE.binary)
t.test(BRAF ~ T.binary)
t.test(BRAF ~ N_Stage)
t.test(BRAF ~ Stage.binary2)
t.test(BRAF ~ Recurrence)
t.test(BRAF ~ Vital_Status)


### Cross table analysis ###

attach(dat1)

age.tbl <- table(BRAF.binary2, Age.binary)
oddsratio(age.tbl)
result1 <- c(oddsratio(age.tbl)[[2]][2,], oddsratio(age.tbl)[[3]][6])
names(result1) <- c("OR","2.5%", "97.5%", "p-value")
result1 <- round(result1, 3)

gender.tbl <- table(BRAF.binary2, Gender)
oddsratio(gender.tbl)
result2 <- c(oddsratio(gender.tbl)[[2]][2,], oddsratio(gender.tbl)[[3]][6])
names(result2) <- c("OR","2.5%", "97.5%", "p-value")
result2 <- round(result2, 3)

Hashimoto.tbl <- table(BRAF.binary2, Hashimoto)
oddsratio(Hashimoto.tbl)
result3 <- c(oddsratio(Hashimoto.tbl)[[2]][2,], oddsratio(Hashimoto.tbl)[[3]][6])
names(result3) <- c("OR","2.5%", "97.5%", "p-value")
result3 <- round(result3, 3)

ETE.tbl <- table(BRAF.binary2, ETE.binary)
result4 <- c(oddsratio(ETE.tbl)[[2]][2,], oddsratio(ETE.tbl)[[3]][6])
names(result4) <- c("OR","2.5%", "97.5%", "p-value")
result4 <- round(result4, 3)

T.tbl <- table(BRAF.binary2, T.binary)
result5 <- c(oddsratio(T.tbl)[[2]][2,], oddsratio(T.tbl)[[3]][6])
names(result5) <- c("OR","2.5%", "97.5%", "p-value")
result5 <- round(result5, 3)

N.tbl <- table(BRAF.binary2, N_Stage)
result6 <- c(oddsratio(N.tbl)[[2]][2,], oddsratio(N.tbl)[[3]][6])
names(result6) <- c("OR","2.5%", "97.5%", "p-value")
result6 <- round(result6, 3)
result6

Stage.tbl <- table(BRAF.binary2, Stage.binary2)
result7 <- c(oddsratio(Stage.tbl)[[2]][2,], oddsratio(Stage.tbl)[[3]][6])
names(result7) <- c("OR","2.5%", "97.5%", "p-value")
result7 <- round(result7, 3)
result7

Recurrence.tbl <- table(BRAF.binary2, Recurrence)
result8 <- c(oddsratio(Recurrence.tbl)[[2]][2,], oddsratio(Recurrence.tbl)[[3]][6])
names(result8) <- c("OR","2.5%", "97.5%", "p-value")
result8 <- round(result8, 3)
result8

Vital_Status.tbl <- table(BRAF.binary2, Vital_Status)
result9 <- c(oddsratio(Vital_Status.tbl)[[2]][2,], oddsratio(Vital_Status.tbl)[[3]][6])
names(result9) <- c("OR","2.5%", "97.5%", "p-value")
result9 <- round(result9, 3)
result9

result <- rbind(result1, result2, result3, result4, result5, result6, result7, result8, result9)
rownames(result) <- c("Age","Gender","Thyroiditis", "ETE","T-stage", "N-stage", "Stage", "Recurrence", "Survival")
result
write.table(result, "C://Work/BRAF exp/result1 using mean braf.txt", sep='\t', quote=F, row.names=T)


### Generalized linear model ###

attach(dat1)
fit1 <- glm(BRAF ~ Age.binary, data=dat1)
summary(fit1)
r1 <- summary(fit1)$coeff

fit2 <- glm(BRAF ~ Gender, data=dat1)
summary(fit2)
r2 <- summary(fit2)$coeff

fit3 <- glm(BRAF ~ Hashimoto, data=dat1)
summary(fit3)
r3 <- summary(fit3)$coeff

fit4 <- glm(BRAF ~ ETE.binary, data=dat1)
summary(fit4)
r4 <- summary(fit4)$coeff

fit5 <- glm(BRAF ~ T.binary)
summary(fit5)
r5 <- summary(fit5)$coeff

fit6 <- glm(BRAF ~ N_Stage)
summary(fit6)
r6 <- summary(fit6)$coeff

fit7 <- glm(BRAF ~ Stage.binary2)
summary(fit7)
r7 <- summary(fit7)$coeff
fit8 <- glm(BRAF ~ Recurrence)
summary(fit8)
r8 <- summary(fit8)$coeff
fit9 <- glm(BRAF ~ Vital_Status)
summary(fit9)
r9 <- summary(fit9)$coeff
result.glm.uni <- rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9)
result.glm.uni <- round(result.glm.uni, 3)[c(2,4,6,8,10,12,14,16,18), ]
write.csv(result.glm.uni, "C://Work/BRAF exp/glm result1.csv")

fit10 <- glm(BRAF ~ Age.binary + Gender + Hashimoto + ETE.binary + T.binary + N_Stage + Stage.binary2 + Recurrence + Vital_Status)
summary(fit10)
step(fit10, direction="both")
fit10.reduced <- glm(BRAF ~ Gender + Hashimoto + T.binary + N_Stage + Recurrence)
summary(fit10.reduced)
result.glm.multi <- summary(fit10.reduced)$coeff
result.glm.multi <- round(result.glm.multi, 3)[c(2,3,4,5,6), ]
write.csv(result.glm.multi, "C://Work/BRAF exp/glm result1 multi.csv")


### Plotting OR ###

library(ggplot2); library(gridExtra)

dat <- read.csv("OR result ETE T.csv")
names(dat) <- c("Group", "Variable", "OR", "LCI", "UCI", "p-value")

ggplot(dat, aes(x = Variable, y = OR, ymin = LCI, ymax = UCI)) + 
  geom_pointrange(aes(col = factor(Group)), position=position_dodge(width=0.30)) + 
  ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + scale_color_discrete(name = "Group") + xlab("")

ggplot(dat, aes(x = Variable, y = OR, ymin = LCI, ymax = UCI)) + geom_linerange(aes(col = factor(Group)), position=position_dodge(width=0.30)) +
  geom_point(aes(shape = factor(Group)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")

### Boxplot ###

g1 <- ggplot(dat.mutant, aes(ETE.binary, BRAF, fill=ETE.binary)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  labs(title = "BRAF V600E Mutant\nAll Age") + xlab("ETE") + ylab("BRAF RPKM") + guides(fill=FALSE)

g2 <- ggplot(dat.mutant, aes(T.binary, BRAF, fill=T.binary)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAll Age") + xlab("T Stage") + ylab("BRAF RPKM") + guides(fill=FALSE)

g3 <- ggplot(dat.mutant.old, aes(ETE.binary, BRAF, fill=ETE.binary)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAge Over 45") + xlab("ETE") + ylab("BRAF RPKM") + guides(fill=FALSE)
g4 <- ggplot(dat.mutant.old, aes(T.binary, BRAF, fill=T.binary)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAge Over 45") + xlab("T Stage") + ylab("BRAF RPKM") + guides(fill=FALSE)

grid.arrange(g1, g2, g3, g4)

### Violin Plot ###
g1 <- ggplot(dat.mutant, aes(ETE.binary, BRAF, fill=ETE.binary)) + geom_violin(alpha=0.5, color="gray") + stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  labs(title = "BRAF V600E Mutant\nAll Age") + xlab("ETE") + ylab("BRAF RPKM") + guides(fill=FALSE)

g2 <- ggplot(dat.mutant, aes(T.binary, BRAF, fill=T.binary)) + geom_violin(alpha=0.5, color="gray") + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAll Age") + xlab("T Stage") + ylab("BRAF RPKM") + guides(fill=FALSE)

g3 <- ggplot(dat.mutant.old, aes(ETE.binary, BRAF, fill=ETE.binary)) + geom_violin(alpha=0.5, color="gray") + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAge Over 45") + xlab("ETE") + ylab("BRAF RPKM") + guides(fill=FALSE)
g4 <- ggplot(dat.mutant.old, aes(T.binary, BRAF, fill=T.binary)) + geom_violin(alpha=0.5, color="gray") + stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  labs(title = "BRAF V600E Mutant\nAge Over 45") + xlab("T Stage") + ylab("BRAF RPKM") + guides(fill=FALSE)

grid.arrange(g1, g2, g3, g4)
