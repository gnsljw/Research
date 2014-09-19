uni.LR <- function(i) {
     dat <- readRDS("TCGA GLUT data.rds")
     a <- round(summary(glm(as.formula(Vital_Status ~ dat[,i]), data=dat, family="binomial"))$coefficients, 4)
     result <- a
     return(result)
}

result <- list()
for (i in 28:63){
     result <- append(result, list(uni.LR (i)))
}

names(result) <- as.character(names(dat[28:63]))

result

sink("Univatiate Logistic Regression Result.txt")
result
sink()
