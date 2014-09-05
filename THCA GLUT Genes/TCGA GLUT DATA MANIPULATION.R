dat <- readRDS("GLUT_TCGA.rds")

dat$N_stage
dat$N_stage[dat$N_stage == "N1a"] <- "N1"
dat$N_stage[dat$N_stage == "N1b"] <- "N1"
dat$N_stage <- factor(dat$N_stage, levels = c("N0","N1","NX"), labels = c("N0","N1","NX"))
summary(dat$N_stage)

dat$T_stage
dat$T_stage[dat$T_stage == "T4a"] <- "T4"
dat$T_stage[dat$T_stage == "T4b"] <- "T4"
dat$T_stage[dat$T_stage == "T1b"] <- "T1a"
dat$T_stage <- factor(dat$T_stage, levels = c("T1a","T2","T3","T4"), labels = c("T1","T2","T3","T4"))
dat$T_stage

dat$M_stage
dat$M_stage[dat$M_stage == "Mx"] <- "MX"
dat$M_stage <- factor(dat$M_stage, levels = c("M0","M1","MX"), labels = c("M0","M1","MX"))

dat$Laterality[dat$Laterality == "Lt lobe"] <- "Left lobe"
dat$Laterality <- factor(dat$Laterality, levels = c("Bilateral","Left lobe","Right lobe","Isthmus"), 
                         labels = c("Bilateral","Left","Right","Isthmus"))

dat$Operation[dat$Operation == "Lobectomy Lt"] <- "Lobectomy Rt"
dat$Operation[dat$Operation == "Completion Lt Thyroidectomy"] <- "Completion Rt Thyroidectomy"
dat$Operation[dat$Operation == "Partial thyroidectomy"] <- "Subtotal thyroidectomy"

dat$Operation <- factor(dat$Operation, levels = c("Total Thyroidectomy","Lobectomy Rt","Completion Rt Thyroidectomy","Subtotal thyroidectomy"), 
                        labels = c("Total", "Lobectomy", "Completion","Other"))
summary(dat$Operation)

summary(dat[,1:27])

saveRDS (dat, "TCGA GLUT data.rds")

