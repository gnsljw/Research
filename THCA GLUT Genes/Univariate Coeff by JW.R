univariate <- function(i) {
     dat <- readRDS("TCGA GLUT data.rds")
     a <- round(summary(lm(as.formula(dat[,i] ~ Gender), data=dat))$coefficients, 4)
     b <- round(summary(lm(as.formula(dat[,i] ~ Race), data=dat))$coefficients, 4)
     c <- round(summary(lm(as.formula(dat[,i] ~ Age_Dx), data=dat))$coefficients, 4)
     d <- round(summary(lm(as.formula(dat[,i] ~ Recurrence), data=dat))$coefficients, 4)
     e <- round(summary(lm(as.formula(dat[,i] ~ Vital_Status), data=dat))$coefficients, 4)
     f <- round(summary(lm(as.formula(dat[,i] ~ Hashimoto), data=dat))$coefficients, 4)
     g <- round(summary(lm(as.formula(dat[,i] ~ Histologic_Dx), data=dat))$coefficients, 4)
     h <- round(summary(lm(as.formula(dat[,i] ~ Max_Size), data=dat))$coefficients, 4) 
     j <- round(summary(lm(as.formula(dat[,i] ~ Metastatic_LN), data=dat))$coefficients, 4)
     k <- round(summary(lm(as.formula(dat[,i] ~ ETE), data=dat))$coefficients, 4)
     l <- round(summary(lm(as.formula(dat[,i] ~ T_stage), data=dat))$coefficients, 4)
     m <- round(summary(lm(as.formula(dat[,i] ~ N_stage), data=dat))$coefficients, 4)
     n <- round(summary(lm(as.formula(dat[,i] ~ Stage), data=dat))$coefficients, 4)
     o <- round(summary(lm(as.formula(dat[,i] ~ BRAF_V600E), data=dat))$coefficients, 4)
     
     result <- rbind(a,b,c,d,e,f,g,h,j,k,l,m,n,o)
     return(result)
}

result <- list()
for (i in 28:63){
     result <- append(result, list(univariate (i)))
     }

names(result) <- as.character(names(dat[28:63]))
summary(result)

sink("Univariate Result.txt")
result
sink()


