setwd("C://Work/Graves")
graves <- read.csv("141226_Graves_comparison_full.csv")
dim(graves)
head(graves)

treatment <- graves$Group
Sex <- graves$Sex
Age <- graves$Age
BMI <- graves$BMI
Indication <- graves$Indication
Dissection <- graves$dissection
ThyroidWt <- graves$ThyWt
FU <- graves$FU
data = cbind(Age,Sex,BMI,Indication,Dissection,ThyroidWt,FU, treatment)
dim(data)

## nearest neighbor matching (1:1)
require(MatchIt)
data1 = data[,c("Age","Sex","BMI","Indication","Dissection",
                "ThyroidWt", "FU", "treatment")]
data1 = as.data.frame(na.omit(data1)) 
m.out = matchit(treatment~Age+Sex+BMI+Indication+Dissection+ThyroidWt+FU,
                method="nearest", data=data1, ratio = 1)
m.out
final_data = match.data(m.out)
write.csv(final_data, file = "matchNN.csv")
plot(m.out)
par(mfrow=c(1,1))
plot(m.out, type = "jitter")
plot(m.out, type = "hist")

# Full Matching
require(optmatch)
data1 = data[,c("Age","Sex","BMI","Indication","Dissection",
                "ThyroidWt", "FU", "treatment")]
data1 = as.data.frame(na.omit(data1)) 
m.out = matchit(treatment~Age+Sex+BMI+Indication+Dissection+ThyroidWt+FU,
                method="full", data=data1)
m.out 
final_data = match.data(m.out) 
write.csv(final_data, file = "matchF.csv")
par(mfrow=c(2,2))
plot(m.out)
plot(m.out, type = "jitter")
plot(m.out, type = "hist")
