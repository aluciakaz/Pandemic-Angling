# Pandemic-Angling
Code for models and analyses presented in Kaz et al. (in prep).

####Load Packages & Data####
#setwd("C:/Users/aluci/Dropbox/COVID Project/Manuscripts/UrbanRural Analysis/Analysis")

# Load packages
library(tidyverse)
library(VGAM)
library(Hmisc)
library(reshape2)
library(plyr)
library(maps)
library(fabricatr)
library(lme4)
library(rcompanion)
library(dplyr)
library(MuMIn)
library(emmeans)
library(patchwork)
library(jagsUI)
library(effsize)

# Load data and workspace
allDat <- read.csv("allDat.csv",header=TRUE)
#load("workspace_08.RData")

# Add log-transformed variables to data
allDat$PD.log <- log(allDat$POPdensity)
allDat$CD.log <- log(allDat$COVIDdensity)

# Get Differences in trips
allDat$diff <- allDat$Q16_1 - allDat$Q15_1

# COVID density scales to integers
allDat$COVIDdensity_1 <- allDat$COVIDdensity*100

#Low correlation: 0.18 untransformed, 0.37 log-transformed
ggplot(NULL) + 
  theme_classic(base_size = 25) +
  theme(legend.position = "none")+
  labs(x = "Population density", y = "COVID-19 cases per capita") +
  geom_point(aes(x=PD.log, y=CD.log, color = as.factor(STATEFP)), allDat) 

cor.test(allDat$PD.log, allDat$CD.log,  method = "pearson")
cor.test(allDat$POPdensity, allDat$COVIDdensity,  method = "pearson")

c=c("ARKANSAS", "CONNECTICUT","FLORIDA","IOWA","MISSOURI","NORTH CAROLINA","SOUTH CAROLINA", "TEXAS", "UTAH", "WYOMING" )
namevec <- map(database = "state", col = "white",fill = F, namesonly=TRUE)
map(database = "state",col = c("gray90", "skyblue4")[1+(namevec %in% tolower(c) )],fill=T)




####Q1: Were license sales motivated by pop density / covid case density####

#remove non-responses
Dat=allDat[which(allDat$Q9!=""),]

#recode character responses as ordered factors
Dat$Q9_Score <- revalue(Dat$Q9,
                        c("I decided not to purchase a license in this time frame (but I typically would)"="5",
                          "The pandemic delayed my purchased"="4",
                          "No change"="3",
                          "The pandemic hastened my purchase"="2",
                          "I decided to purchase a license (when I wasn’t planning to purchase one this year)"="1"))
Dat$Q9_Score <- as.factor(Dat$Q9_Score)
Dat <- within(Dat, {Q9_Score <- factor(Q9_Score, ordered = TRUE)})

#Code ordinal model
vglm.parallel <- vglm(Q9_Score ~  POPdensity + COVIDdensity_1, 
                      Dat,
                      family = cumulative(link = "logitlink", 
                                          parallel = TRUE, 
                                          reverse = TRUE))

#Code null model and models for likelihood ratio testing
vglm.parallel0 <- vglm(Q9_Score ~ + 1, 
                       Dat,
                       family = cumulative(link = "logitlink", 
                                           parallel = TRUE, 
                                           reverse = TRUE))
vglm.parallel1 <- vglm(Q9_Score ~  POPdensity, 
                      Dat,
                      family = cumulative(link = "logitlink", 
                                          parallel = TRUE, 
                                          reverse = TRUE))
vglm.parallel2 <- vglm(Q9_Score ~  COVIDdensity_1, 
                      Dat,
                      family = cumulative(link = "logitlink", 
                                          parallel = TRUE, 
                                          reverse = TRUE))


#Likelihood ratio test: Population density is not an important predictor of license purchase behavior (lrt test stat: 0.3154)
lrtest(vglm.parallel, vglm.parallel0)
lrtest(vglm.parallel1, vglm.parallel0)
lrtest(vglm.parallel2, vglm.parallel0)

#Pull odds ratios and associated confidence intervals
exp(cbind("Odds ratio" = coef(vglm.parallel), confint.default(vglm.parallel, level = 0.95)))

#Examine summary output and calculated R-squared values
summary(vglm.parallel)
nagelkerke(vglm.parallel) 

#Create new dataset of predicted probabilities for a figure
newdat <- data.frame(POPdensity = rep(seq(from = min(Dat$POPdensity), to = max(Dat$POPdensity), length.out = 100)),
                     COVIDdensity_1 = rep(seq(from = min(Dat$COVIDdensity_1), to = max(Dat$COVIDdensity_1), length.out = 100), 3))
newdat <- cbind(newdat, predict(vglm.parallel, newdat, type = "response"))
lnewdat <- melt(newdat, id.vars = c("COVIDdensity_1", "POPdensity"),
                variable.name = "Level", value.name="Probability")

#Add 95% confidence intervals for each group
lnewdat1 <- filter(lnewdat,Level==1)
lnewdat1$lb <- lnewdat1$Probability - (3.588/100)
lnewdat1$ub <- lnewdat1$Probability + (3.588/100)

lnewdat2 <- filter(lnewdat,Level==2)
lnewdat2$lb <- lnewdat2$Probability - (1.405/100)
lnewdat2$ub <- lnewdat2$Probability + (1.405/100)

lnewdat3 <- filter(lnewdat,Level==3)
lnewdat3$lb <- lnewdat3$Probability - (0.005/100)
lnewdat3$ub <- lnewdat3$Probability + (0.005/100)

lnewdat4 <- filter(lnewdat,Level==4)
lnewdat4$lb <- lnewdat4$Probability - (0.002/100)
lnewdat4$ub <- lnewdat4$Probability + (0.002/100)

lnewdat5 <- filter(lnewdat,Level==5)
lnewdat5$lb <- lnewdat5$Probability - (0.045/100)
lnewdat5$ub <- lnewdat5$Probability + (0.045/100)

lnewdat <- rbind(lnewdat1, lnewdat2, lnewdat3, lnewdat4, lnewdat5)

#Create figure of OLR probability curves
(p1 <- ggplot(data = lnewdat) + 
    labs(x = "COVID-19 cases per 100 people", y = "Probability of license purchase", size = 15) +
    geom_line(mapping = aes(x = COVIDdensity_1, y = Probability, color = Level), size = 1) +
    geom_ribbon(aes(ymin=lb, ymax=ub, x = COVIDdensity_1, fill = Level), alpha=0.2) +
    scale_color_manual(values = c("skyblue4", "powderblue", "gray35", "indianred", "indianred4")) +
    scale_fill_manual(values = c("skyblue4", "powderblue", "gray35", "indianred", "indianred4")) +
    theme_classic(base_size = 15) +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none"))


#Get averaged stats from 1st and 4th quantile to compare purchase probabilities at different COVID caseloads
dat <- lnewdat[which(lnewdat$Level==2),]
dat2 <- dat %>%  mutate(quantile = ntile(COVIDdensity_1, 4))
dat3 <- dat2[which(dat2$quantile==1),]
first <- mean(dat3$Probability)
dat3 <- dat2[which(dat2$quantile==4),]
fourth <- mean(dat3$Probability)
(hasten <- fourth/first)

dat <- lnewdat[which(lnewdat$Level==1),]
dat2 <- dat %>%  mutate(quantile = ntile(COVIDdensity_1, 4))
dat3 <- dat2[which(dat2$quantile==1),]
first <- mean(dat3$Probability)
dat3 <- dat2[which(dat2$quantile==4),]
fourth <- mean(dat3$Probability)
(new <- fourth/first)



####Q2 - Overall: Differenced fishing trips & demographic info#### 

#Write Bayesian linear model for overall trend
sink("lm1.txt")
cat("
model {

# Priors
 alpha ~ dnorm(0,0.001)
 beta1 ~ dnorm(0,0.001)
 beta2 ~ dnorm(0,0.001)
 sigma ~ dunif(0, 100)
 
# Likelihood
 for (i in 1:n) {
    diff[i] ~ dnorm(mu[i], tau) 
    mu[i] <- alpha + beta1*CD.log[i] + beta2*PD.log[i]}
    
# Derived quantities
 tau <- 1/ (sigma * sigma)}
",fill=TRUE)
sink()

# Bundle data
data1 <- list(diff = as.numeric(allDat$diff), 
              CD.log = allDat$CD.log, 
              PD.log = allDat$PD.log,
              n = nrow(allDat))

# Inits function
inits <- function(){ list(alpha=rnorm(1), beta1=rnorm(1), beta2 = rnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha","beta1", "beta2", "sigma")

# MCMC settings
nc = 3  ;  ni=1200  ;  nb=200  ;  nt=1

# Start Gibbs sampler
out1 <- autojags(data = data1, inits = inits, parameters = params, model = "lm1.txt", 
                 n.thin = nt, n.chains = nc)
print(out1, dig = 3)

#Visualize
overall <-  ggplot() +
  labs(x= "log(COVID-19 cases per capita)", y = "Differenced fishing trips", size = 15) +
  geom_point(aes(x=CD.log, y=diff), allDat, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.276, intercept=1.670, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5))




####Q2 - Gender: Differenced fishing trips & demographic info####

#Recode character responses as factors
Dat <- allDat
Dat$gender <- recode(Dat$Q23, 
                  " " = "NA",  
                  'Prefer not to answer' = "NA",
                  'Female' = "1",
                  'Male' = "2")

#Create new df with necessary columns and omit NAs
Dat <- Dat[,c("CD.log","diff","gender")]
Dat$gender <- as.numeric(Dat$gender)
Dat <- na.omit(Dat)

# Bundle data
str(bdata <- list(diff = Dat$diff, gender = Dat$gender, 
                  CD.log = Dat$CD.log, n.group = 2, n = nrow(Dat)))

#Write Bayesian interaction model for trend by gender
cat(file = "lm.txt", "
model {

# Priors
 for (i in 1:n.group){		
    alpha[i] ~ dnorm(0, 0.001)		# Intercepts
    beta[i] ~ dnorm(0, 0.001)		# Slopes
 }
 sigma ~ dunif(0, 100)			# Residual standard deviation
 tau <- 1 / ( sigma * sigma)

# Likelihood
 for (i in 1:n) {
    diff[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[gender[i]] + beta[gender[i]]* CD.log[i]
 }

}
")

# Inits function
inits <- function(){ list(alpha = rnorm(2, 0, 2), 
                          beta = rnorm(2, 1, 1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta", "sigma")

# MCMC settings
na <- 1000  ;  nc <- 3  ;  ni <- 1200  ;  nb <- 200  ;  nt <- 2

# Call JAGS, check convergence and summarize posteriors
out <- autojags(bdata, inits, params, "lm.txt", n.adapt = na, n.thin = nt, n.chains = nc, parallel = TRUE)
print(out, dig = 3)

#Visualize
female <-  ggplot() +
  labs(x= "log(per capita COVID-19 cases)", y = "Differenced \nfishing trips", size = 15) +
  geom_point(aes(x=CD.log, y=diff), allDat[which(allDat$Q23=="Female"),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.371, intercept=2.423, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("Female") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
male <-  ggplot() +
  labs(x= "log(per capita COVID-19 cases)", size = 15) +
  geom_point(aes(x=CD.log, y=diff), allDat[which(allDat$Q23=="Male"),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.230, intercept=1.136, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("Male") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))




####Q2 - Income: Differenced fishing trips & demographic info#### 

#Recode character responses as factors
Dat <- allDat
Dat$income <- recode(allDat$Q24, 
                     " " = "NA",  
                     'Less than $10,000' = "1",
                     '$10,000 - $19,999' = "1",
                     '$20,000 - $29,999' = "1",
                     '$30,000 - $39,999' = "2",
                     '$40,000 - $49,999' = "2",
                     '$50,000 - $59,999' = "3",
                     '$60,000 - $69,999' = "3",
                     '$70,000 - $79,999' = "3",
                     '$80,000 - $89,999' = "4",
                     '$90,000 - $99,999' = "4",
                     '$100,000 - $149,999' = "4",
                     'More than $150,000' = "4")

#Create new df with necessary columns and omit NAs
Dat <- Dat[,c("CD.log","diff","income")]
Dat$income <- as.numeric(Dat$income)
Dat <- na.omit(Dat)

#lmer model for comparison
lme1 <- lmer(diff ~ CD.log + (CD.log | income), data = Dat)
summary(lme1)

#Write Bayesian random effects model for trend by income group
sink("lme.model1.txt")
cat("
model {

# Likelihood
 for (i in 1:n) {
    diff[i] ~ dnorm(mu[i], tau)		# The 'residual' random variable
    mu[i] <- alpha[income[i]] + beta[income[i]]* CD.log[i]  # Expectation
 }
 
# Priors
 for (i in 1:ngroups){
    alpha[i] <- B[i,1]
    beta[i] <- B[i,2]
    B[i,1:2] ~ dmnorm(B.hat[i,], Tau.B[,])
    B.hat[i,1] <- mu.int
    B.hat[i,2] <- mu.slope}

 mu.int ~ dnorm(0, 0.001)		# Hyperpriors for random intercepts
 mu.slope ~ dnorm(0, 0.001)		# Hyperpriors for random slopes

 Tau.B[1:2,1:2] <- inverse(Sigma.B[,])
 Sigma.B[1,1] <- pow(sigma.int,2)
 sigma.int ~ dunif(0, 100)		# SD of intercepts
 Sigma.B[2,2] <- pow(sigma.slope,2)
 sigma.slope ~ dunif(0, 100)		# SD of slopes
 Sigma.B[1,2] <- rho*sigma.int*sigma.slope
 Sigma.B[2,1] <- Sigma.B[1,2]
 rho ~ dunif(-1,1)
 covariance <- Sigma.B[1,2]

 tau <- 1 / ( sigma * sigma)		# Residual
 sigma ~ dunif(0, 100)}			# Residual standard deviation
",fill=TRUE)
sink()

# Bundle data
data3 <- list(diff = as.numeric(Dat$diff), 
              income = as.numeric(Dat$income), 
              CD.log = Dat$CD.log, 
              ngroups = 4, 
              n = nrow(Dat))

# Inits function
inits <- function(){ list(mu.int = rnorm(1, 0, 1), 
                          sigma.int = rlnorm(1), 
                          mu.slope = rnorm(1, 0, 1), 
                          sigma.slope = rlnorm(1), 
                          rho = runif(1, -1, 1), 
                          sigma = rlnorm(1))}

# Parameters to estimate
parameters <- c("alpha", "beta", "mu.int", "sigma.int", "mu.slope", 
                "sigma.slope", "rho", "covariance", "sigma")

# MCMC settings
ni <- 30000
nb <- 5000
nt <- 2
nc <- 3

# Start Gibbs sampler
out3 <- autojags(data3, 
                 inits, 
                 parameters, 
                 "lme.model1.txt", 
                 n.thin=nt, 
                 n.chains=nc, 
                 #n.burnin=nb, 
                 #n.iter=ni,
                 parallel = TRUE, 
                 n.cores = 2
)

#Look at Bayesian and Frequentist analyses
print(out3, dig = 2)
lme1

#Visualize
income1 <-  ggplot() +
  labs(y = "Differenced \nfishing trips", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$income==1),], color = "grey65", pch=21, alpha = .7, size = 1.5, fill = "grey90") +
  geom_abline(slope=-0.06, intercept=-0.8, color = "grey55", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("< $30K") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
income2 <-   ggplot() +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$income==2),], color = "grey65", pch=21, alpha = .7, size = 1.5, fill = "grey90") +
  geom_abline(slope=0.2, intercept=0.3, color = "grey55", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("$30K - 49K") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
income3 <-  ggplot() +
  labs(y = "Differenced \nfishing trips", x= "log(per capita COVID-19 cases)", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$income==3),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.33, intercept=1.62, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("$50K - 79K") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
income4 <-  ggplot() +
  labs(x= "log(per capita COVID-19 cases)", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$income==4),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.22, intercept=1.43, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("> $79K") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))



####Q2 - Age: Differenced fishing trips & demographic info####

#Recode character responses as factors
Dat <- allDat
Dat$age <- recode(Dat$Q22, 
                  " " = "NA",  
                  '17 or younger' = "1",
                  '18–20' = "2",
                  '21–29' = "2",
                  '30–39' = "2",
                  '40–49' = "3",
                  '50–59' = "3",
                  '60 or older' = "4")

#Create new df with necessary columns and omit NAs
Dat <- Dat[,c("CD.log","diff","age")]
Dat$age <- as.numeric(Dat$age)
Dat <- na.omit(Dat)

#lmer model for comparison
lme2 <- lmer(diff ~ CD.log + (CD.log | age), data = Dat)

#Write Bayesian random effects model for trend by age group
sink("lme.model2.txt")
cat("
model {

# Likelihood
 for (i in 1:n) {
    diff[i] ~ dnorm(mu[i], tau)		# The 'residual' random variable
    mu[i] <- alpha[age[i]] + beta[age[i]]* CD.log[i]  # Expectation
 }
 
# Priors
 for (i in 1:ngroups){
    alpha[i] <- B[i,1]
    beta[i] <- B[i,2]
    B[i,1:2] ~ dmnorm(B.hat[i,], Tau.B[,])
    B.hat[i,1] <- mu.int
    B.hat[i,2] <- mu.slope}

 mu.int ~ dnorm(0, 0.001)		# Hyperpriors for random intercepts
 mu.slope ~ dnorm(0, 0.001)		# Hyperpriors for random slopes

 Tau.B[1:2,1:2] <- inverse(Sigma.B[,])
 Sigma.B[1,1] <- pow(sigma.int,2)
 sigma.int ~ dunif(0, 100)		# SD of intercepts
 Sigma.B[2,2] <- pow(sigma.slope,2)
 sigma.slope ~ dunif(0, 100)		# SD of slopes
 Sigma.B[1,2] <- rho*sigma.int*sigma.slope
 Sigma.B[2,1] <- Sigma.B[1,2]
 rho ~ dunif(-1,1)
 covariance <- Sigma.B[1,2]

 tau <- 1 / ( sigma * sigma)		# Residual
 sigma ~ dunif(0, 100)}			# Residual standard deviation
",fill=TRUE)
sink()

# Bundle data
data4 <- list(diff = as.numeric(Dat$diff), 
              age = as.numeric(Dat$age), 
              CD.log = Dat$CD.log, # try log, or standardize
              ngroups = 4, 
              n = nrow(Dat))

# Inits function
inits <- function(){ list(mu.int = rnorm(1, 0, 1), 
                          sigma.int = rlnorm(1), 
                          mu.slope = rnorm(1, 0, 1), 
                          sigma.slope = rlnorm(1), 
                          rho = runif(1, -1, 1), 
                          sigma = rlnorm(1))}

# Parameters to estimate
parameters <- c("alpha", "beta", "mu.int", "sigma.int", "mu.slope", 
                "sigma.slope", "rho", "covariance", "sigma")

# MCMC settings
ni <- 30000
nb <- 5000
nt <- 2
nc <- 3

# Start Gibbs sampler
out4 <- autojags(data4, 
                 inits, 
                 parameters, 
                 "lme.model2.txt", 
                 n.thin=nt, 
                 n.chains=nc, 
                 #n.burnin=nb, 
                 #n.iter=ni,
                 parallel = TRUE, 
                 n.cores = 2)

#Look at Bayesian and Frequentist analyses
print(out4, dig = 2)			
lme2	

#Visualize
age1 <-  ggplot() +
  labs(y = "Differenced \nfishing trips", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$age==1),], color = "grey65", pch=21, alpha = .7, size = 1.5, fill = "grey90") +
  geom_abline(slope=0.89, intercept=7.10, color = "grey55", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20) ) +
  theme_classic(base_size = 15) +
  ggtitle("< 18") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
age2 <-  ggplot() +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$age==2),], color = "grey65", pch=21, alpha = .7, size = 1.5, fill = "grey90") +
  geom_abline(slope=0.13, intercept=1.19, color = "grey55", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("18 - 39") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
age3 <-  ggplot() +
  labs(y = "Differenced \nfishing trips", x= "log(COVID-19 cases \n per capita)", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$age==3),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.25, intercept=1.32, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("40 - 59") +
  theme(plot.title = element_text(hjust = 0.5))
age4 <-  ggplot() +
  labs(x= "log(COVID-19 cases \n per capita)", size = 15) +
  geom_point(aes(x=CD.log, y=diff), Dat[which(Dat$age==4),], color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  geom_abline(slope=0.4, intercept=1.32, color = "skyblue4", size = 2)+
  scale_x_continuous(limits=c(-9, -3) ) +
  scale_y_continuous(limits=c( -20, 20 ) ) +
  theme_classic(base_size = 15) +
  ggtitle("> 60") +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))



####Q3: Fishing motivations before and during pandemic####

#During pandemic: look at fishing motivation data
table(allDat$Q14_1) #Food
table(allDat$Q14_2) #Sport/thrill
table(allDat$Q14_3) #Nature
table(allDat$Q14_4) #Social/Family
table(allDat$Q14_5) #Stress
table(allDat$Q14_6) #Competition
table(allDat$Q14_7) #Free time
table(allDat$Q14_8) #Get Away
table(allDat$Q14_10) #Don't normally fish

#During pandemic: Code ranked responses as binary: reason or not reason for fishing
allDat$Q14_1b <- ifelse(allDat$Q14_1 > 0, 1, 0)
allDat$Q14_2b <- ifelse(allDat$Q14_2 > 0, 1, 0)
allDat$Q14_3b <- ifelse(allDat$Q14_3 > 0, 1, 0)
allDat$Q14_4b <- ifelse(allDat$Q14_4 > 0, 1, 0)
allDat$Q14_5b <- ifelse(allDat$Q14_5 > 0, 1, 0)
allDat$Q14_6b <- ifelse(allDat$Q14_6 > 0, 1, 0)
allDat$Q14_7b <- ifelse(allDat$Q14_7 > 0, 1, 0)
allDat$Q14_8b <- ifelse(allDat$Q14_8 > 0, 1, 0)
allDat$Q14_10b <- ifelse(allDat$Q14_10 > 0, 1, 0)

#During pandemic: Code non-responses as zeros
allDat$Q14_1b <- allDat$Q14_1b %>% replace(is.na(.), 0)
allDat$Q14_2b <- allDat$Q14_2b %>% replace(is.na(.), 0)
allDat$Q14_3b <- allDat$Q14_3b %>% replace(is.na(.), 0)
allDat$Q14_4b <- allDat$Q14_4b %>% replace(is.na(.), 0)
allDat$Q14_5b <- allDat$Q14_5b %>% replace(is.na(.), 0)
allDat$Q14_6b <- allDat$Q14_6b %>% replace(is.na(.), 0)
allDat$Q14_7b <- allDat$Q14_7b %>% replace(is.na(.), 0)
allDat$Q14_8b <- allDat$Q14_8b %>% replace(is.na(.), 0)
allDat$Q14_10b <- allDat$Q14_10b %>% replace(is.na(.), 0)

#During pandemic: Code logistic regression models - log transformed
d1 <- glm(Q14_1b ~ CD.log + PD.log, data = allDat, family = "binomial") #Food 
d2 <- glm(Q14_2b ~ CD.log + PD.log, data = allDat, family = "binomial") #Sport/thrill 
d3 <- glm(Q14_3b ~ CD.log + PD.log, data = allDat, family = "binomial") #Nature 
d4 <- glm(Q14_4b ~ CD.log + PD.log, data = allDat, family = "binomial") #Social/Family
d5 <- glm(Q14_5b ~ CD.log + PD.log, data = allDat, family = "binomial") #Stress 
d6 <- glm(Q14_6b ~ CD.log + PD.log, data = allDat, family = "binomial") #Competition 
d7 <- glm(Q14_7b ~ CD.log + PD.log, data = allDat, family = "binomial") #Free time 
d8 <- glm(Q14_8b ~ CD.log + PD.log, data = allDat, family = "binomial") #Get Away 
d9 <- glm(Q14_10b ~ CD.log + PD.log, data = allDat, family = "binomial") #Don't normally fish

#Look at output
summary(d1)
summary(d2)
summary(d3)
summary(d4)
summary(d5)
summary(d6)
summary(d7)
summary(d8)
summary(d9)

#Before pandemic: look at fishing motivation data
table(allDat$Q12_1) #Food
table(allDat$Q12_2) #Sport/thrill
table(allDat$Q12_3) #Nature
table(allDat$Q12_4) #Social/Family
table(allDat$Q12_5) #Stress
table(allDat$Q12_6) #Competition
table(allDat$Q12_7) #Free time
table(allDat$Q12_8) #Get Away
table(allDat$Q12_10) #Don't normally fish

#Before pandemic: Code ranked responses as binary: reason or not reason for fishing
allDat$Q12_1b <- ifelse(allDat$Q12_1 > 0, 1, 0)
allDat$Q12_2b <- ifelse(allDat$Q12_2 > 0, 1, 0)
allDat$Q12_3b <- ifelse(allDat$Q12_3 > 0, 1, 0)
allDat$Q12_4b <- ifelse(allDat$Q12_4 > 0, 1, 0)
allDat$Q12_5b <- ifelse(allDat$Q12_5 > 0, 1, 0)
allDat$Q12_6b <- ifelse(allDat$Q12_6 > 0, 1, 0)
allDat$Q12_7b <- ifelse(allDat$Q12_7 > 0, 1, 0)
allDat$Q12_8b <- ifelse(allDat$Q12_8 > 0, 1, 0)
allDat$Q12_10b <- ifelse(allDat$Q12_10 > 0, 1, 0)

#Before pandemic: Code non-responses as zeros
allDat$Q12_1b <- allDat$Q12_1b %>% replace(is.na(.), 0)
allDat$Q12_2b <- allDat$Q12_2b %>% replace(is.na(.), 0)
allDat$Q12_3b <- allDat$Q12_3b %>% replace(is.na(.), 0)
allDat$Q12_4b <- allDat$Q12_4b %>% replace(is.na(.), 0)
allDat$Q12_5b <- allDat$Q12_5b %>% replace(is.na(.), 0)
allDat$Q12_6b <- allDat$Q12_6b %>% replace(is.na(.), 0)
allDat$Q12_7b <- allDat$Q12_7b %>% replace(is.na(.), 0)
allDat$Q12_8b <- allDat$Q12_8b %>% replace(is.na(.), 0)
allDat$Q12_10b <- allDat$Q12_10b %>% replace(is.na(.), 0)

#Before pandemic: Code logistic regression models - log transformed
b1 <- glm(Q12_1b ~ CD.log + PD.log, data = allDat, family = "binomial") #Food 
b2 <- glm(Q12_2b ~ CD.log + PD.log, data = allDat, family = "binomial") #Sport/thrill
b3 <- glm(Q12_3b ~ CD.log + PD.log, data = allDat, family = "binomial") #Nature
b4 <- glm(Q12_4b ~ CD.log + PD.log, data = allDat, family = "binomial") #Social/Family
b5 <- glm(Q12_5b ~ CD.log + PD.log, data = allDat, family = "binomial") #Stress
b6 <- glm(Q12_6b ~ CD.log + PD.log, data = allDat, family = "binomial") #Competition
b7 <- glm(Q12_7b ~ CD.log + PD.log, data = allDat, family = "binomial") #Free time 
b8 <- glm(Q12_8b ~ CD.log + PD.log, data = allDat, family = "binomial") #Get Away
b9 <- glm(Q12_10b ~ CD.log + PD.log, data = allDat, family = "binomial") #Don't normally fish

#Look at output
summary(b1)
summary(b2)
summary(b3)
summary(b4)
summary(b5)
summary(b6)
summary(b7)
summary(b8)
summary(b9)

#Calculate odds ratios
exp(cbind("Odds ratio" = coef(b1), confint.default(b1, level = 0.95)))
exp(cbind("Odds ratio" = coef(b2), confint.default(b2, level = 0.95)))
exp(cbind("Odds ratio" = coef(b3), confint.default(b3, level = 0.95)))
exp(cbind("Odds ratio" = coef(b4), confint.default(b4, level = 0.95)))
exp(cbind("Odds ratio" = coef(b5), confint.default(b5, level = 0.95)))
exp(cbind("Odds ratio" = coef(b6), confint.default(b6, level = 0.95)))
exp(cbind("Odds ratio" = coef(b7), confint.default(b7, level = 0.95)))
exp(cbind("Odds ratio" = coef(b8), confint.default(b8, level = 0.95)))
exp(cbind("Odds ratio" = coef(b9), confint.default(b9, level = 0.95)))

exp(cbind("Odds ratio" = coef(d1), confint.default(d1, level = 0.95)))
exp(cbind("Odds ratio" = coef(d2), confint.default(d2, level = 0.95)))
exp(cbind("Odds ratio" = coef(d3), confint.default(d3, level = 0.95)))
exp(cbind("Odds ratio" = coef(d4), confint.default(d4, level = 0.95)))
exp(cbind("Odds ratio" = coef(d5), confint.default(d5, level = 0.95)))
exp(cbind("Odds ratio" = coef(d6), confint.default(d6, level = 0.95)))
exp(cbind("Odds ratio" = coef(d7), confint.default(d7, level = 0.95)))
exp(cbind("Odds ratio" = coef(d8), confint.default(d8, level = 0.95)))
exp(cbind("Odds ratio" = coef(d9), confint.default(d9, level = 0.95)))


#Visualize
p1 <- ggplot(data = allDat) +
  labs(y = "Sport", title = "Pre-pandemic", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_2b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_2b), method = "glm", se = TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p2 <- ggplot(data = allDat) +
  labs(y = "Sport", title = "During pandemic", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_2b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_2b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.y = element_text(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p3 <- ggplot(data = allDat) +
  labs(y = "Stress relief", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_5b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_5b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p4 <- ggplot(data = allDat) +
  labs(y = "Stress relief", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_5b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_5b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.y = element_text(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p5 <- ggplot(data = allDat) +
  labs(y = "Competition", x = "log(Population Density)", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_6b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_6b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))
p6 <- ggplot(data = allDat) +
  labs(y = "Competition", x = "log(Population Density)", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_6b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_6b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.y = element_text(),
        plot.title = element_text(hjust = 0.5))

#Visualize
food <- ggplot(data = allDat) +
  labs(y = "Food") +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_1b), height = .15, color = "indianred4", pch=21, alpha = .3, size = 1.5, fill = "indianred1", alpha = .5) +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_1b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "indianred4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
sport <- ggplot(data = allDat) +
  labs(y = "Sport", title = "During pandemic") +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_2b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_2b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
freetime <-  ggplot(data = allDat) +
  labs(y = "Free time", x = "log(COVID-19 cases \n per capita)") +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_7b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_7b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))




####Additional Manuscript Figures####

female | male
(income1 + income2) / (income3 + income4)
(age1 + age2) / (age3 + age4)
(female | male) / (income3 | income4) / (age3 | age4)
(p1 + p2 + sport) / (p3 + p4 + food) / (p5 + p6 + freetime)

#Figure of all fishing motivations
pb1 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Pre-pandemic", title = "Food", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_1b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_1b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
pb2 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Sport", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_2b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_2b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb3 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Nature", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_3b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_3b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb4 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Social", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_4b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_4b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb5 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Stress relief", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_5b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_5b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb6 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Competition", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_6b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_6b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb7 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Free time", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_7b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_7b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb8 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Get away", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_8b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_8b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
pb9 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", title = "Don't normally fish", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q12_10b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q12_10b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

####

pd1 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "During pandemic", title = "During pandemic", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_1b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_1b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(size = 15),
        plot.title = element_blank())
pd2 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Sport", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_2b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_2b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd3 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Nature", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_3b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_3b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd4 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Social", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_4b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_4b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd5 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Stress relief", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_5b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_5b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd6 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Competition", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_6b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_6b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd7 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Free time", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_7b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_7b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd8 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Get away", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_8b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_8b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
pd9 <- ggplot(data = allDat) +
  labs(x = "log(Population \n density)", y = "Don't normally fish", size = 15) +
  geom_jitter(mapping = aes(x = PD.log, y = Q14_10b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = PD.log, y = Q14_10b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())

####

c1 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "During pandemic", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_1b), height = .15, color = "indianred4", pch=21, alpha = .3, size = 1.5, fill = "indianred1", alpha = .5) +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_1b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "indianred4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(size = 15),
        plot.title = element_blank())
c2 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Sport", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_2b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_2b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c3 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Nature", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_3b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_3b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c4 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Social", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_4b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_4b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c5 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Stress relief", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_5b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_5b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c6 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Competition", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_6b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_6b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c7 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Free time", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_7b), height = .15, color = "skyblue3", pch=21, alpha = .3, size = 1.5, fill = "powderblue") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_7b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "skyblue4", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c8 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Get away", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_8b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_8b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())
c9 <- ggplot(data = allDat) +
  labs(x = "log(COVID- 19 cases \n per capita)", y = "Don't normally fish", size = 15) +
  geom_jitter(mapping = aes(x = CD.log, y = Q14_10b), height = .15, color = "grey65", pch=21, alpha = .3, size = 1.5, fill = "grey90") +
  stat_smooth(mapping = aes(x = CD.log, y = Q14_10b), method = "glm", se=TRUE, method.args = list(family=binomial), color = "grey40", size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_text(),
        axis.title.y = element_blank(),
        plot.title = element_blank())

(pb1 | pb2 | pb3 | pb4 | pb5 | pb6 | pb7 | pb8 | pb9) / (pd1 | pd2 | pd3 | pd4 | pd5 | pd6 | pd7 | pd8 | pd9) / (c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8 | c9)
