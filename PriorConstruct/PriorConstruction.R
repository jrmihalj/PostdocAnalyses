# Constructing prior distributions from field experiment.

# I'm going to construct priors for:
# - nu.bar
# - CV
# I will base these priors on three of the isolates:
# - KLIP, TMB-1, NM
# I will exclude LOVE, because it seems like such an outlier virus

##################################################################
#-----------------------------------------------------------------
##################################################################

######################
## PRIORS OF NU.BAR ##
######################
library(ggmcmc)

klip.nu <- ggs(out.null, family="nu") #best model = NULL
tmb.nu <- ggs(out.noCV, family="nu") #best model = NO CV
nm.nu <- ggs(out.full, family="nu.bar") #best model = FULL (CV + MU)

nu.combined <- c(klip.nu$value, tmb.nu$value, nm.nu$value)

# IMPORTANT:
# Need to re-scale these estimates to # 4th instars per meter^2, not 1st per cm^2
# So, divide nu.bar by ratio (~0.01) and then divide by 10^4 to convert to m^2


hist(nu.combined[nu.combined <= 25], breaks=40)
# Distribution looks like either gamma or lognormal

# Fit a gamma model to the data to estimate scale and shape
setwd("~/Documents/PostdocResearch/2015Summer/FieldDataAnalyses/PriorConstruct")
cat("model {
  #### PRIORS #####

  # FOR GAMMA #
  scale ~ dgamma(0.001, 0.001)
  shape ~ dgamma(0.001, 0.001)
  rate <- 1/scale

  # FOR LNORM #
  mu ~ dnorm(0,.01)
  sd ~ dunif(0,100)
  tau <- 1 / (sd^2)
  
  ### LIKELIHOOD ###
  for(i in 1:length(nu.combined1)){
      nu.combined1[i] ~ dgamma( rate , shape )
      nu.combined2[i] ~ dlnorm( mu, tau )
  }

}", fill=T, file="PriorFit.txt")

library(rjags)
prior.dat <- list(nu.combined1 = sample(nu.combined, 1500),
                  nu.combined2 = sample(nu.combined, 1500))
prior.mod <- jags.model(file="PriorFit.txt", data=prior.dat, n.chains = 3, n.adapt=1000)
out.prior <- coda.samples(prior.mod, variable.names=c("scale", "shape", "mu", "sd"), n.iter=10000, thin=10)

quartz()
plot(out.prior)
summary(out.prior)

# Check model fit
shape.fit <- median(ggs(out.prior, family="shape")$value)
scale.fit <- median(ggs(out.prior, family="scale")$value)

G <- Gammad(scale= scale.fit, shape= shape.fit)
fitted.gam <- r(G)(1500)

mu.fit <- median(ggs(out.prior, family="mu")$value)
sd.fit <- median(ggs(out.prior, family="sd")$value)

L <- Lnorm(meanlog=mu.fit, sdlog = sd.fit)
fitted.lnorm<- r(L)(1500)

fit.df <- data.frame(prior.dat[[1]], fitted.gam, fitted.lnorm)
colnames(fit.df) <- c("raw", "fitted.g", "fitted.l")

# Plot the fitted gam and the data
plot.fit.nu <- ggplot(fit.df, aes(x=raw))+
  geom_histogram(aes(y=..density..), binwidth=.5, fill="gray")+
  #geom_density(aes(x=fitted.g), color="red")+
  geom_density(aes(x=fitted.l), color="blue")+
  theme_classic()

quartz()
plot.fit.nu

# Log norm fits well, while gamma is very poor fit.

# So the priors for nu are:
# Lognormal distribution with mu = 1.12, sd = 0.84

##################################################################
#-----------------------------------------------------------------
##################################################################


#########################
## PRIORS OF Coef.Var. ##
#########################

library(ggmcmc)

# Only NM had a significant CV, but I'll use the estimates from the full model for all three isolates
CV.all <- ggs(out.full, family="C") 

CV.combined <- subset(CV.all, Parameter != "C[2]")$value
hist(CV.combined, breaks=40)
# Distribution looks like a half-normal, but I need k, which is 1/CV^2
K <- CV.combined^-2
hist(K[K<100], breaks=40)
# maybe follows case 2 in the C code, where
# exp(N(mu,sd))

# Fit a gamma model to the data to estimate scale and shape
setwd("~/Documents/PostdocResearch/2015Summer/FieldDataAnalyses/PriorConstruct")
cat("model {
  #### PRIORS #####

  # FOR STUDENT-t #
  mu ~ dnorm(0,.01)
  sd ~ dunif(0,100)
  tau <- 1 / (sd^2)

  ### LIKELIHOOD ###
  for(i in 1:length(K)){
      K[i] ~ exp(dnorm( mu, tau))
  }

}", fill=T, file="CVPriorFit.txt")

library(rjags)
prior.dat.CV <- list(K = sample(K, 1000))
prior.mod.CV <- jags.model(file="CVPriorFit.txt", data=prior.dat.CV, n.chains = 3, n.adapt=1000)
out.prior.CV <- coda.samples(prior.mod.CV, variable.names=c("mu", "sd"), n.iter=10000, thin=10)

quartz()
plot(out.prior.CV)
summary(out.prior.CV)


# Check model fit

mu.fit.CV <- mean(ggs(out.prior.CV, family="mu")$value)
sd.fit.CV <- mean(ggs(out.prior.CV, family="sd")$value)

fitted.CV <- abs(rnorm(length(CV.combined), mu.fit.CV, sd.fit.CV))

hist(CV.combined, col="gray", breaks=40)
hist(fitted.CV, main="Fit Compare", xlab="CV", breaks=40, add =T)

# The distribution is very well modeled by a half-normal distribution:
# mu: 1.89
# sd: 1.25
# Problem...the model uses k, which is 1/CV^2

# Solution: 
# Just fit to a normal-log: exp(mean + N(0,sd))
K.rand <- exp(1e-4 + rnorm(1000, 0, 1))
hist(K.rand, breaks=40)
median(K.rand)
hist(K[K<50], breaks=40)

# Try a log-normal
K.rand <- rlnorm(1000,meanlog=0,sdlog=1)
hist(K.rand, breaks=40)
median(K.rand)
hist(K[K<50], breaks=40)



