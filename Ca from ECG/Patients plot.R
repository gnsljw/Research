library(RColorBrewer)
colorset = brewer.pal(12, "Paired")

patients <- read.csv("D://Work/QTcalcium/20141214/QT complete obs.csv")

### Corr.test
result <- vector()
for (i in 1:nrow(patients)){
        x <- as.numeric(patients[i, 5:7])
        y <- as.numeric(patients[i, 8:10])
        result[i] <- cor.test(x,y)[3]
}
result

### Fig 1

x <- as.numeric(patients[1,5:7])
y <- as.numeric(patients[1,8:10])

plot(x,y, type="o", xlab="QTc Interval", ylab="Calcium", xlim=c(340,480), ylim=c(6.5, 10), 
     col=colorset[1], lwd=2, main="Individual Correlation Plots\nBetween Calcium and QTc")

for (i in 2:12){
        x <- as.numeric(patients[i, 5:7])
        y <- as.numeric(patients[i, 8:10])
        lines(x,y, type="o", col = colorset[i], lwd=2)
}

legend("topright", bty="n", legend=c("Pt.1, p=0.655","Pt.2, p=0.111","Pt.3, p=0.233","Pt.4, p=0.259",
                            "Pt.5, p=0.448","Pt.6, p=0.042","Pt.7, p=0.855","Pt.8, p=0.833",
                            "Pt.9, p=0.574","Pt.10, p=0.846 ","Pt.11, p=0.184 ","Pt.12, p=0.103 "), 
       fill=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
              "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"))

### Fig 2
x <- as.numeric(patients[13,5:7])
y <- as.numeric(patients[13,8:10])

plot(x,y, type="o", xlab="QTc Interval", ylab="Calcium", xlim=c(340,480), ylim=c(6.5, 10), 
     col=colorset[1], lwd=2, main="Individual Correlation Plots\nBetween Calcium and QTc")

for (i in 14:24){
        x <- as.numeric(patients[i, 5:7])
        y <- as.numeric(patients[i, 8:10])
        lines(x,y, type="o", col = colorset[i-12], lwd=2)
}

x <- as.numeric(patients[25,5:7])
y <- as.numeric(patients[25,8:10])
lines(x,y, type="o", col = "black", lwd=2)

legend("topright", legend=c("Pt.13, p=0.163 ","Pt.14, p=0.440 ","Pt.15, p=0.987 ",
                            "Pt.16, p=0.588 ","Pt.17, p=0.212 ","Pt.18, p=0.669 ",
                            "Pt.19, p=0.119 ","Pt.20, p=0.902 ","Pt.21, p=0.030 ",
                            "Pt.22, p=0.891 ","Pt.23, p=0.890 ","Pt.24, p=0.049 ", "Pt.25, p=0.309 "), 
       fill=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
              "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928", "black"), bty="n")

