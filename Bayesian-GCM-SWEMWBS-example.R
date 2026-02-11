
#Author notes: 
#Codes developed based on: Oravecz, Z., & Muth, C. (2018). Fitting growth curve models in the Bayesian framework. Psychonomic Bulletin & Review, 25(1), 235-255.
#If not interesed in fitting the model mannually, please check out the 'brms' package 

# Loading relevant packages
library(logr)
library(foreign)
library(rjags)  

setwd("insert the directory here")

## ----Step 1: load data and create a data vector for JAGS-----------------
# Step 1
# Reading in the data (in Stata format)
Mydata<-read.dta("insert the dataset")

# Visually checking the data
head(Mydata) 

# Checking the number of rows 
N <- dim(Mydata)[1] 
N

# then transform these values into numeric entries of a matrix
data1 <- as.matrix(Mydata[,5:6])

# Creating a variable that indicates the number of time points
nrT <- 2

# Creating a time vector 
time <- c(0, 1) 

# Creating a list of all the variables that we created above
jagsData1 <- list("Y"=data1,"N"=N,"nrT"=nrT,"time"=time)

## ----Step 2: Specify the growth curve model with linear trend, highlight=FALSE----
# Step 2
# Specifying the growth curve model
LinearGrowthCurve = cat("
model {

# Loop over participants
for (i in 1:N) { 
 # Loop over measurement occasions
 for (t in 1:nrT) {  
  # Specifying likelihood function, corresponding to Equation 1:
  Y[i,t] ~ dnorm(betas[i,1]+betas[i,2]*time[t], 1/pow(sdLevel1Error,2))
 }  

 # Specifying the level-2 bivariate distribution of intercepts and slopes (Eq. 6)
 betas[i,1:2] ~ dmnorm(Level2MeanVector[i,1:2], interpersonPrecisionMatrix[1:2,1:2])

 # The mean of the intercept is modeled as a function of positivity group membership
 Level2MeanVector[i,1] <- MedPInt 
 Level2MeanVector[i,2] <- MedPSlope 
}

# Specifying prior distributions
MedPInt ~ dnorm(0,0.01)
MedPSlope ~ dnorm(0,0.01)
sdLevel1Error ~ dunif(0,100)
sdIntercept  ~ dunif(0,100)
sdSlope  ~ dunif(0,100)
corrIntSlope ~ dunif(-1,1)

# Transforming model parameters
# 1. Defining the elements of the level-2 covariance matrix
interpersonCovMatrix[1,1] <- sdIntercept * sdIntercept
interpersonCovMatrix[2,2] <- sdSlope * sdSlope
interpersonCovMatrix[1,2] <- corrIntSlope * sdIntercept* sdSlope
interpersonCovMatrix[2,1] <- interpersonCovMatrix[1,2]
# 2. Taking the inverse of the covariance to get the precision matrix
interpersonPrecisionMatrix <- inverse(interpersonCovMatrix)
# 3. time2 outcome
MedPTime2=MedPInt + MedPSlope
}
",file = "GCM.txt")

## ----Step 3: Collect model parameters, set sampler, tidy = TRUE----------
# Step 3
# Monitoring parameters
parameters  <- c("MedPSlope","MedPInt",
                 "sdIntercept", "sdSlope", 
                 "corrIntSlope", "sdLevel1Error","betas", "MedPTime2")
# Specifying sampler settings
adaptation  <- 2000 
chains  <- 6    
burnin  <- 1000 
thinning    <- 5
# Defining the number of samples drawn from the posterior in each chain
postSamples <- 30000
# Computing the number of posterior samples needed per chain for JAGS  
nrOfIter    <- ceiling((postSamples*thinning)/chains)

## ----echo=FALSE,  eval=TRUE, warning=FALSE, message=FALSE, error=FALSE, include=FALSE----
#, results=hide,>>=
# fixing the random seed for reproducibility
fixedinits<- list(list(.RNG.seed=5,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=6,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=7,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=8,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=9,.RNG.name="base::Mersenne-Twister"),list(.RNG.seed=10,.RNG.name="base::Mersenne-Twister"))


## ----Step 4: Call JAGS and fit the model, tidy = TRUE, eval=TRUE---------
# Step 4

# Creating JAGS model object
jagsModel1<-jags.model("GCM.txt",data=jagsData1,n.chains=chains,n.adapt=adaptation,inits=fixedinits)
# Running burn-in iterations
update(jagsModel1,n.iter=burnin)
# Drawing posterior samples
codaSamples1<-coda.samples(jagsModel1,variable.names=parameters,n.iter=nrOfIter, thin = thinning,seed=15)

## ----Step 5: Explore results, tidy = TRUE, eval = TRUE-------------------
# Step 5 
source("posteriorSummaryStats.R")
# Part 1: Check convergence
resulttable <- summarizePost(codaSamples1)
saveNonConverged <- resulttable[resulttable$RHAT>1.1,]
if (nrow(saveNonConverged) == 0){
  print("Convergence criterion was met for every parameter.")
}else{ 
  print("Not converged parameter(s):")
  show(saveNonConverged)
}

sink("log/SWEMWBS.txt")
show(summarizePost(codaSamples1, filters =  c("^Med", "^sd","^corr"))) 
sink()

# export the results table in a cvs file 
write.csv(resulttable, "SWEMWBS_resulttable.csv")
