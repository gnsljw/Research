# Email to ehsan at alumni.ubc.ca to report error/suggestions

################################################################
###### Synthetic Data Generation
################################################################

n = 3500
x1 = rnorm(n, 1, 100)
x2 = rnorm(n, 10, 100)
x3 = rnorm(n, 5, 1)
x4 = rnorm(n, 10, 10)
x5 = rbinom(n, 1, .4)
x6 = rnorm(n, 30, 5)
treatment = rbinom(n, 1, .15)
data = cbind(x1,x2,x3,x4,x5,x6, treatment)
dim(data)

################################################################
###### Propensity score matching
###### nearest neighbor matching (1:1)
################################################################

install.packages("MatchIt")
require(MatchIt)
# data1 is the subset of data with only the selected variables mentioned below
data1 = data[,c("x1","x2","x3","x4","x5","x6", "treatment")]
# getting rid of missing values (below)
data1 = as.data.frame(na.omit(data1)) 
# matching is performed below using propensity scores given the covariates mentioned below
m.out = matchit(treatment~x1+x2+x3+x4+x5+x6,method="nearest", data=data1, ratio = 1)
# check the sample sizes (below)
m.out 
# Final matched data saved as final_data
final_data = match.data(m.out) 
# (here distance = propensity score)
# After performing command below check your computer for NN.csv file that can be imported in SPSS
write.csv(final_data, file = "matchNN.csv")
# check balance (below)
plot(m.out) # covariate balance
plot(m.out, type = "jitter") # propensity score locations
plot(m.out, type = "hist") #check matched treated vs matched control


################################################################
###### Propensity score matching
###### caliper matching
################################################################

require(MatchIt)
# data1 is the subset of data with only the selected variables mentioned below
data1 = data[,c("x1","x2","x3","x4","x5","x6", "treatment")]
# getting rid of missing values (below)
data1 = as.data.frame(na.omit(data1)) 
m.out.test = matchit(treatment~x1+x2+x3+x4+x5+x6,method="nearest", data=data1, ratio = 1)
test_data = match.data(m.out.test) 
ps.sd = sd(test_data$distance)
# matching is performed below using propensity scores given the covariates mentioned below
# caliper = 0.25 times sd of propensity scores (optimal)
m.out = matchit(treatment~x1+x2+x3+x4+x5+x6,method="nearest", data=data1, caliper = 0.25*ps.sd)
# check the sample sizes (below)
m.out 
# Final matched data saved as final_data
final_data = match.data(m.out) 
# (here distance = propensity score)
# After performing command below check your computer for matchC.csv file that can be imported in SPSS
write.csv(final_data, file = "matchC.csv")
# check balance (below)
plot(m.out) # covariate balance
plot(m.out, type = "jitter") # propensity score locations
plot(m.out, type = "hist") #check matched treated vs matched control



################################################################
###### Propensity score matching
###### nearest neighbor matching (1:4)
################################################################

require(MatchIt)
# data1 is the subset of data with only the selected variables mentioned below
data1 = data[,c("x1","x2","x3","x4","x5","x6", "treatment")]
# getting rid of missing values (below)
data1 = as.data.frame(na.omit(data1)) 
# matching is performed below using propensity scores given the covariates mentioned below
m.out = matchit(treatment~x1+x2+x3+x4+x5+x6,method="nearest", data=data1, ratio = 4)
# check the sample sizes (below)
m.out 
# Final matched data saved as final_data
final_data = match.data(m.out) 
# (here distance = propensity score)
# After performing command below check your computer for matchNN4.csv file that can be imported in SPSS
write.csv(final_data, file = "matchNN4.csv")
# check balance (below)
plot(m.out) # covariate balance
plot(m.out, type = "jitter") # propensity score locations
plot(m.out, type = "hist") #check matched treated vs matched control

################################################################
###### Propensity score matching 
###### subclassification by Propensity score
################################################################

require(MatchIt)
# data1 is the subset of data with only the selected variables mentioned below
data1 = data[,c("x1","x2","x3","x4","x5","x6", "treatment")]
# getting rid of missing values (below)
data1 = as.data.frame(na.omit(data1)) 
# matching is performed below using propensity scores given the covariates mentioned below
m.out = matchit(treatment~x1+x2+x3+x4+x5+x6, method = "subclass", subclass = 4, data=data1)
# check the sample sizes (below)
m.out 
# Final matched data saved as final_data 
final_data = match.data(m.out) 
# (here distance = propensity score and subclass = the subclass variable)
# ****put this subclass in the final outcome model****
# After performing command below check your computer for matchSC.csv file that can be imported in SPSS
write.csv(final_data, file = "matchSC.csv")
# check balance (below)
plot(m.out) # covariate balance
plot(m.out, type = "jitter") # propensity score locations
plot(m.out, type = "hist") #check matched treated vs matched control

################################################################
###### Propensity score matching
###### full matching
################################################################
install.packages("optmatch")
require(optmatch)
require(MatchIt)
# data1 is the subset of data with only the selected variables mentioned below
data1 = data[,c("x1","x2","x3","x4","x5","x6", "treatment")]
# getting rid of missing values (below)
data1 = as.data.frame(na.omit(data1)) 
# matching is performed below using propensity scores given the covariates mentioned below
m.out = matchit(treatment~x1+x2+x3+x4+x5+x6,method="full", data=data1)
# check the sample sizes (below)
m.out 
# Final matched data saved as final_data
final_data = match.data(m.out) 
# (here distance = propensity score)
# After performing command below check your computer for matchF.csv file that can be imported in SPSS
write.csv(final_data, file = "matchF.csv")
# check balance (below)
plot(m.out) # covariate balance
plot(m.out, type = "jitter") # propensity score locations
plot(m.out, type = "hist") #check matched treated vs matched control
