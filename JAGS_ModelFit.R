# JAGS model fitting of Field Experiment 2015 data
# Author: JR Mihaljevic
# Date: July 2015

##################################################################
#-----------------------------------------------------------------
##################################################################

setwd("~/Documents/PostdocResearch/2015Summer")


# Need to format "Prop.Inf" dataframe using the following source:
# 'AutopsyDataFormat.R'

head(Prop.Inf)

# Sort the data appropriately:
Prop.Inf <- Prop.Inf[order(Prop.Inf$P.Virus,Prop.Inf$P.Decay,Prop.Inf$Cad.Den.Adj),]

# make the cadaver density per cm^2
Prop.Inf$Cad.Den.cm <- Prop.Inf$Cad.Den.Adj / 1e4
# Create necessary model variables:

n.inf <- Prop.Inf$N.Inf
n.recov <- Prop.Inf$N.Recov
Virus <- as.factor(as.numeric(Prop.Inf$P.Virus))
P <- Prop.Inf$Cad.Den.cm
decay <- ifelse(Prop.Inf$P.Decay == "0D", 0, 3)
Cad.Den <-as.numeric(factor(as.character(Prop.Inf$Cad.Den)))

N.obs <- length(P)

# I want to show the error in n.inf and P
# Need a bunch of NA's

P <- c(P, rep("NA", N.obs))

n.inf <- rep(n.inf,2)
n.recov <- rep(n.recov,2)
Virus <- rep(Virus,2)
decay <- rep(decay,2)
Cad.Den <-rep(Cad.Den,2)

N.obs <- length(P)

n.inf.na <- rep("NA", N.obs)



###################
## RUN THE MODEL ##
###################

library(rjags)

j.data <- list(N.obs = N.obs,
			n.inf = n.inf,
			n.recov = n.recov,
			Virus = Virus,
			P = P,
			decay = decay,
			Cad.Den = Cad.Den,
			n.inf.na = n.inf.na)

params <- c("C", "nu", "mu", "P", "n.inf.na")

mod.inits <- list(C = c(1,1,1,1),
		  	nu = c(.005,.005,.005,.005),
			mu = c(0,0,0,0))


mod <- jags.model(file="InfModel.txt", inits=mod.inits,
			data=j.data, n.chains = 3, n.adapt=2000)

out <- coda.samples(mod, variable.names=params, n.iter=10000, n.thin=10)

#windows(record=TRUE)
plot(out)
summary(out)

######################
# SUMMARY #
######################

library(ggmcmc)
source("HDI.R")

C.data <- ggs(out, family="C")
nu.data <- ggs(out, family="nu")
mu.data <- ggs(out, family="mu")

param.sum <- mat.or.vec(nr=4,nc=9)

colnames(param.sum) <- c("C.med","C.lo","C.hi","nu.med","nu.lo","nu.hi","mu.med","mu.lo","mu.hi")

for(i in 1:4){
	param.sum[i,1] <- median(C.data[which(C.data$Parameter == paste("C[",i,"]",sep="")),]$value)
	param.sum[i,2:3] <- HDI(C.data[which(C.data$Parameter == paste("C[",i,"]",sep="")),]$value)

	param.sum[i,4] <- median(nu.data[which(nu.data$Parameter == paste("nu[",i,"]",sep="")),]$value)
	param.sum[i,5:6] <- HDI(nu.data[which(nu.data$Parameter == paste("nu[",i,"]",sep="")),]$value)

	param.sum[i,7] <- median(mu.data[which(mu.data$Parameter == paste("mu[",i,"]",sep="")),]$value)
	param.sum[i,8:9] <- HDI(mu.data[which(mu.data$Parameter == paste("mu[",i,"]",sep="")),]$value)
}

param.sum

