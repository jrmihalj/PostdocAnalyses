# JAGS model fitting of Field Experiment 2015 data
# Author: JR Mihaljevic
# Date: July 2015

##################################################################
#-----------------------------------------------------------------
##################################################################

setwd("~/Documents/PostdocResearch/2015Summer/FieldDataAnalyses")


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

##################################################################
#-----------------------------------------------------------------
# I'M GOING TO RUN 4 MODELS AND COMPARE:
# 1. FULL MODEL W/ CV AND MU
# 2. REDUCED MODEL W/ CV, NO MU
# 3. REDUCED MODEL W/ MU, NO CV
# 4. NULL MODEL W/ NEITHER CV NOR MU
#-----------------------------------------------------------------
##################################################################

#####################################
## RUN THE FULL MODEL W/ CV AND MU ##
#####################################

library(rjags)

j.data <- list(N.obs = N.obs,
			n.inf = n.inf,
			n.recov = n.recov,
			Virus = Virus,
			P = P,
			decay = decay,
			Cad.Den = Cad.Den)

params.full <- c("C", "nu.bar", "mu", "prob")

mod.inits.full <- list(C = c(1,1,1,1),
		  	nu.bar = c(.005,.005,.005,.005),
			mu = c(0,0,0,0))


mod.full <- jags.model(file="InfModel_CVmu.txt", inits=mod.inits.full,
			data=j.data, n.chains = 3, n.adapt=10000)

out.full <- coda.samples(mod.full, variable.names=params.full, n.iter=50000, thin=50)

plot(out.full[,"mu[4]"])
summary(out.full)

######################
# SUMMARY #
######################

library(ggmcmc)
source("HDI.R")

C.data <- ggs(out.full, family="C")
nu.bar.data <- ggs(out.full, family="nu.bar")
mu.full.data <- ggs(out.full, family="mu")

param.sum.full <- mat.or.vec(nr=4,nc=9)

colnames(param.sum.full) <- c("C.med","C.lo","C.hi","nu.med","nu.lo","nu.hi","mu.med","mu.lo","mu.hi")

for(i in 1:4){
  param.sum.full[i,1] <- median(C.data[which(C.data$Parameter == paste("C[",i,"]",sep="")),]$value)
  param.sum.full[i,2:3] <- HDI(C.data[which(C.data$Parameter == paste("C[",i,"]",sep="")),]$value)

  param.sum.full[i,4] <- median(nu.bar.data[which(nu.bar.data$Parameter == paste("nu.bar[",i,"]",sep="")),]$value)
  param.sum.full[i,5:6] <- HDI(nu.bar.data[which(nu.bar.data$Parameter == paste("nu.bar[",i,"]",sep="")),]$value)

  param.sum.full[i,7] <- median(mu.full.data[which(mu.full.data$Parameter == paste("mu[",i,"]",sep="")),]$value)
  param.sum.full[i,8:9] <- HDI(mu.full.data[which(mu.full.data$Parameter == paste("mu[",i,"]",sep="")),]$value)
}

param.sum.full

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

##########################################
## RUN THE NULL MODEL WITH CV BUT NO MU ##
##########################################
params.noMu <- c("C", "nu.bar", "prob")

mod.inits.noMu <- list(C = c(1,1,1,1),
                       nu.bar = c(.005,.005,.005,.005))


mod.noMu <- jags.model(file="InfModel_CV.txt", inits=mod.inits.noMu,
                       data=j.data, n.chains = 3, n.adapt=10000)

out.noMu <- coda.samples(mod.noMu, variable.names=params.noMu, n.iter=50000, thin=50)

plot(out.noMu)
summary(out.noMu)

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

##########################################
## RUN THE NULL MODEL WITH MU BUT NO CV ##
##########################################
params.noCV <- c("mu", "nu", "prob")

mod.inits.noCV <- list(nu= c(.005,.005,.005,.005),
                       mu = c(0,0,0,0))


mod.noCV <- jags.model(file="InfModel_mu.txt", inits=mod.inits.noCV,
                       data=j.data, n.chains = 3, n.adapt=10000)

out.noCV <- coda.samples(mod.noCV, variable.names=params.noCV, n.iter=50000, thin=50)

plot(out.noCV[])
summary(out.noCV)

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

#########################################
## RUN THE NULL MODEL WITHOUT CV OR MU ##
#########################################

params.null <- c("nu", "prob")

mod.inits.null<- list(nu = c(1,1,.1,1))


mod.null <- jags.model(file="InfModel_Null.txt", inits=mod.inits.null,
                     data=j.data, n.chains = 3, n.adapt=10000)

out.null <- coda.samples(mod.null, variable.names=params.null, n.iter=50000, thin=50)

#windows(record=TRUE)
plot(out.null)
summary(out.null)

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

####################
#  CALCULATE WAIC  #
####################
source("waic_inf.R")
WAIC_tab <- array(dim=c(4,4))

WAIC_tab[1,] <- waic_inf(out.full, j.data)$WAIC
WAIC_tab[2,] <- waic_inf(out.noMu, j.data)$WAIC
WAIC_tab[3,] <- waic_inf(out.noCV, j.data)$WAIC
WAIC_tab[4,] <- waic_inf(out.null, j.data)$WAIC

rownames(WAIC_tab) <- c("Full", "CV,NoMu", "MU,NoCV", "Null")
colnames(WAIC_tab) <- c("KLIP", "LOVE", "NM", "TMB-1")

WAIC_tab
#             KLIP      LOVE       NM    TMB-1
# Full    147.2398  92.22460 135.9228 130.1110
# CV,NoMu 142.7814  90.49347 145.1884 130.0893
# MU,NoCV 141.4899 102.56395 142.8584 129.8757
# Null    138.0584 100.79667 162.4113 132.2143
