#########################
# A FEW DIFFERENT PLOTS #
#########################

##########################################
# 1. All data points with median model fit
##########################################

library(ggmcmc)
source("HDI.R")

plots.fitted1 <- list()

for(i in 1:4){
	C.vals <- NULL; nu.vals <- NULL; mu.vals <- NULL;

	C.vals <- subset(C.data, Parameter==paste("C[", i, "]", sep=""))$value
	nu.vals <- subset(nu.bar.data, Parameter==paste("nu.bar[", i, "]", sep=""))$value
	mu.vals <- subset(mu.cv.data, Parameter==paste("mu[", i, "]", sep=""))$value
	
	C.mode <- NULL; nu.mode <- NULL; mu.mode <- NULL;
	C.mode <- median(C.vals)
	nu.mode <- median(nu.vals)
	mu.mode <- median(mu.vals)

	fake.P <- seq(0, 0.05, length.out=300)

	fit.0D <- NULL; fit.3D <- NULL;
	fit.0D <- 1 - (1 + C.mode^2 * nu.mode * fake.P * 7)^(-1/(C.mode^2))
	fit.3D <- 1 - (1 + C.mode^2 * nu.mode * exp(-(mu.mode) * 3) * fake.P * 7)^(-1/(C.mode^2))

	data.fitted <- NULL;
	data.fitted <- data.frame(fit.0D, fit.3D, fake.P)
	data.fitted <- data.fitted[order(data.fitted$fake.P), ]
	
	sub <- NULL
	sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i])

	fit.plot <- NULL
	fit.plot <- ggplot(data.fitted, aes(x=fake.P,y=fit.0D))+
		geom_line(linetype=2, size=1)+
		geom_line(aes(y=fit.3D), linetype=1, size=1)+
		geom_point(data=sub, aes(x=Cad.Den.cm, y=P.Inf, shape=P.Decay))+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
		scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
		scale_x_continuous(limits=c(0,0.05))+
		#theme(legend.position="none")+
		ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


plots.fitted1[[i]] <- fit.plot


}

library(gridExtra)

names(plots.fitted1)=paste("plot", 1:4)
arg.list <- c(plots.fitted1, list(nrow=2))
do.call(grid.arrange, arg.list)


################################################
# 2. Only 3 data points with X and Y error bars
################################################

library(ggplot2)
library(ggmcmc)
library(modeest)
source("HDI.R")

C.data <- ggs(out, family="C")
nu.bar.data <- ggs(out, family="nu")
mu.cv.data <- ggs(out, family="mu")

plots.fitted2 <- list()


for(i in 1:4){
	C.vals <- NULL; nu.vals <- NULL; mu.vals <- NULL;

	C.vals <- subset(C.data, Parameter==paste("C[", i, "]", sep=""))$value
	nu.vals <- subset(nu.bar.data, Parameter==paste("nu.bar[", i, "]", sep=""))$value
	mu.vals <- subset(mu.cv.data, Parameter==paste("mu[", i, "]", sep=""))$value

	C.mode <- NULL; nu.mode <- NULL; mu.mode <- NULL;
	C.mode <- median(C.vals)
	nu.mode <- median(nu.vals)
	mu.mode <- median(mu.vals)

	fake.P <- seq(0, 0.05, length.out=300)
	fit.0D.mean <- 1 - (1 + C.mode^2 * nu.mode * fake.P * 7)^(-1/(C.mode^2))
	fit.3D.mean <- 1 - (1 + C.mode^2 * nu.mode * exp(-(mu.mode) * 3) * fake.P * 7)^(-1/(C.mode^2))


	# Bootstrap storage:
	n.iter=1000
	fit.0D <- mat.or.vec(nr=length(fake.P), nc=n.iter)
	fit.3D <- mat.or.vec(nr=length(fake.P), nc=n.iter)

	for(j in 1:n.iter){
		c.temp <- sample(C.vals,1); nu.temp <- sample(nu.vals,1); mu.temp <- sample(mu.vals,1)
		fit.0D[,j] <- 1 - (1 + c.temp^2 * nu.temp * fake.P * 7)^(-1/(c.temp^2))
		fit.3D[,j] <- 1 - (1 + c.temp^2 * nu.temp * exp(-(mu.temp) * 3) * fake.P * 7)^(-1/(c.temp^2))
	}

	fit.0D.range <- mat.or.vec(nr=length(fake.P), nc=3)
	fit.3D.range <- mat.or.vec(nr=length(fake.P), nc=3)

	for(j in 1:length(fake.P)){
		fit.0D.range[j,1] <- fit.0D.mean[j] ; fit.3D.range[j,1] <- fit.3D.mean[j]
		fit.0D.range[j,2:3] <- as.vector(quantile(fit.0D[j,], prob=c(0.025,0.975))[1:2])
		fit.3D.range[j,2:3] <- as.vector(quantile(fit.3D[j,], prob=c(0.025,0.975))[1:2])

	}

	data.fitted <- NULL;
	data.fitted <- data.frame(fit.0D.range, fit.3D.range, fake.P)
	colnames(data.fitted) <- c("fit.0Dmean","fit.0Dlow","fit.0Dhi","fit.3Dmean","fit.3Dlow","fit.3Dhi","fake.P")
	data.fitted <- data.fitted[order(data.fitted$fake.P), ]
	
	sub <- NULL
	sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i])

	# Find Y error:
	se <- function(x) sqrt(var(x)/length(x))
	se.y <- NULL; mean.y <- NULL
	se.y[1] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 0)$P.Inf); mean.y[1] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 0)$P.Inf)
	se.y[2] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 10)$P.Inf); mean.y[2] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 10)$P.Inf)
	se.y[3] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 40)$P.Inf); mean.y[3] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 40)$P.Inf)
	se.y[4] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 0)$P.Inf); mean.y[4] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 0)$P.Inf)
	se.y[5] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 10)$P.Inf); mean.y[5] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 10)$P.Inf)
	se.y[6] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 40)$P.Inf); mean.y[6] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 40)$P.Inf)   
	

	# Find X error:
	se.x <- NULL; mean.x <- NULL
	se.x[1] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 0)$Cad.Den.cm); mean.x[1] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 0)$Cad.Den.cm)
	se.x[2] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 10)$Cad.Den.cm); mean.x[2] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 10)$Cad.Den.cm)
	se.x[3] <- se(subset(sub, P.Decay=="0D" & Cad.Den == 40)$Cad.Den.cm); mean.x[3] <- mean(subset(sub, P.Decay=="0D" & Cad.Den == 40)$Cad.Den.cm)
	se.x[4] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 0)$Cad.Den.cm); mean.x[4] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 0)$Cad.Den.cm)
	se.x[5] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 10)$Cad.Den.cm); mean.x[5] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 10)$Cad.Den.cm)
	se.x[6] <- se(subset(sub, P.Decay=="3D" & Cad.Den == 40)$Cad.Den.cm); mean.x[6] <- mean(subset(sub, P.Decay=="3D" & Cad.Den == 40)$Cad.Den.cm)  
	
	errs <- data.frame(se.y,mean.y,se.x,mean.x)

	fit.plot <- NULL
	fit.plot <- ggplot(data.fitted, aes(x=fake.P))+
		geom_line(aes(y=fit.0Dmean),linetype=5, size=1)+
		#geom_line(aes(y=fit.0Dlow), linetype=5, size=.85, color="gray")+
		#geom_line(aes(y=fit.0Dhi), linetype=5, size=.85, color="gray")+
		geom_line(aes(y=fit.3Dmean), linetype=1, size=1)+
		#geom_line(aes(y=fit.3Dlow), linetype=1, size=.85, color="gray")+
		#geom_line(aes(y=fit.3Dhi), linetype=1, size=.85, color="gray")+
		geom_point(data = errs[1:3,], aes(x = mean.x, y = mean.y), shape=1, size=3)+
		geom_point(data = errs[4:6,], aes(x = mean.x, y = mean.y), shape=17, size=3)+
		geom_segment(data=errs, aes(y=mean.y-se.y, yend=mean.y+se.y, x=mean.x, xend=mean.x),width=.5)+
		geom_segment(data=errs, aes(y=mean.y, yend=mean.y, x=mean.x-se.x, xend=mean.x+se.x),width=.5)+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
		scale_x_continuous(limits=c(0,0.035))+
		#theme(legend.position="none")+
		ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


	plots.fitted2[[i]] <- fit.plot

}

library(gridExtra)

names(plots.fitted2)=paste("plot", 1:4)
arg.list <- c(plots.fitted2, list(nrow=2))
#windows(record=T)
do.call(grid.arrange, arg.list)




#####################
# WITH ERROR BOUNDS #
#####################

library(ggplot2)
library(ggmcmc)
library(modeest)
source("HDI.R")

C.data <- ggs(out, family="C")
nu.bar.data <- ggs(out, family="nu")
mu.cv.data <- ggs(out, family="mu")

plots.fitted3 <- list()


for(i in 1:4){
	C.vals <- NULL; nu.vals <- NULL; mu.vals <- NULL;

	C.vals <- subset(C.data, Parameter==paste("C[", i, "]", sep=""))$value
	nu.vals <- subset(nu.bar.data, Parameter==paste("nu.bar[", i, "]", sep=""))$value
	mu.vals <- subset(mu.cv.data, Parameter==paste("mu[", i, "]", sep=""))$value

	C.mode <- NULL; nu.mode <- NULL; mu.mode <- NULL;
	C.mode <- median(C.vals)
	nu.mode <- median(nu.vals)
	mu.mode <- median(mu.vals)

	fake.P <- seq(0, 0.05, length.out=300)
	fit.0D.mean <- 1 - (1 + C.mode^2 * nu.mode * fake.P * 7)^(-1/(C.mode^2))
	fit.3D.mean <- 1 - (1 + C.mode^2 * nu.mode * exp(-(mu.mode) * 3) * fake.P * 7)^(-1/(C.mode^2))


	# Bootstrap storage:
	n.iter=1000
	fit.0D <- mat.or.vec(nr=length(fake.P), nc=n.iter)
	fit.3D <- mat.or.vec(nr=length(fake.P), nc=n.iter)

	for(j in 1:n.iter){
		c.temp <- sample(C.vals,1); nu.temp <- sample(nu.vals,1); mu.temp <- sample(mu.vals,1)
		fit.0D[,j] <- 1 - (1 + c.temp^2 * nu.temp * fake.P * 7)^(-1/(c.temp^2))
		fit.3D[,j] <- 1 - (1 + c.temp^2 * nu.temp * exp(-(mu.temp) * 3) * fake.P * 7)^(-1/(c.temp^2))
	}

	fit.0D.range <- mat.or.vec(nr=length(fake.P), nc=3)
	fit.3D.range <- mat.or.vec(nr=length(fake.P), nc=3)

	for(j in 1:length(fake.P)){
		fit.0D.range[j,1] <- fit.0D.mean[j] ; fit.3D.range[j,1] <- fit.3D.mean[j]
		fit.0D.range[j,2:3] <- as.vector(quantile(fit.0D[j,], prob=c(0.025,0.975))[1:2])
		fit.3D.range[j,2:3] <- as.vector(quantile(fit.3D[j,], prob=c(0.025,0.975))[1:2])

	}

	data.fitted <- NULL;
	data.fitted <- data.frame(fit.0D.range, fit.3D.range, fake.P)
	colnames(data.fitted) <- c("fit.0Dmean","fit.0Dlow","fit.0Dhi","fit.3Dmean","fit.3Dlow","fit.3Dhi","fake.P")
	data.fitted <- data.fitted[order(data.fitted$fake.P), ]
	
	sub <- NULL
	sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i])

	fit.plot <- NULL
	fit.plot <- ggplot(data.fitted, aes(x=fake.P))+
		geom_line(aes(y=fit.0Dmean),linetype=5, size=1)+
		geom_line(aes(y=fit.0Dlow), linetype=5, size=.85, color="gray")+
		geom_line(aes(y=fit.0Dhi), linetype=5, size=.85, color="gray")+
		geom_line(aes(y=fit.3Dmean), linetype=1, size=1)+
		geom_line(aes(y=fit.3Dlow), linetype=1, size=.85, color="gray")+
		geom_line(aes(y=fit.3Dhi), linetype=1, size=.85, color="gray")+
		geom_point(data=sub, aes(x=Cad.Den.cm, y=P.Inf, shape=P.Decay))+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
		scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
		scale_x_continuous(limits=c(0,0.05))+
		#theme(legend.position="none")+
		ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))

plots.fitted3[[i]] <- fit.plot

}

library(gridExtra)

names(plots.fitted3)=paste("plot", 1:4)
arg.list <- c(plots.fitted3, list(nrow=2))
#windows(record=T)
do.call(grid.arrange, arg.list)





###################
# Most Complicated:
###################

library(ggmcmc)
library(modeest)
source("HDI.R")

C.data <- ggs(out, family="C")
nu.bar.data <- ggs(out, family="nu")
mu.cv.data <- ggs(out, family="mu")

P.error <- ggs(out, family="P")
Y.error <- ggs(out, family="n.inf.na")

plots.fitted4 <- list()

for(i in 1:4){
	# Fitted values:
	C.vals <- NULL; nu.vals <- NULL; mu.vals <- NULL;

	C.vals <- subset(C.data, Parameter==paste("C[", i, "]", sep=""))$value
	nu.vals <- subset(nu.bar.data, Parameter==paste("nu[", i, "]", sep=""))$value
	mu.vals <- subset(mu.cv.data, Parameter==paste("mu[", i, "]", sep=""))$value
	
	C.mode <- NULL; nu.mode <- NULL; mu.mode <- NULL;
	C.mode <- median(C.vals)
	nu.mode <- median(nu.vals)
	mu.mode <- median(mu.vals)

	fit.0D <- NULL; fit.3D <- NULL;
	fake.P <- seq(0, 0.05, length.out=300)
	fit.0D <- 1 - (1 + C.mode^2 * nu.mode * fake.P * 7)^(-1/(C.mode^2))
	fit.3D <- 1 - (1 + C.mode^2 * nu.mode * exp(-(mu.mode) * 3) * fake.P * 7)^(-1/(C.mode^2))

	data.fitted <- NULL;
	data.fitted <- data.frame(fit.0D, fit.3D, fake.P)
	colnames(data.fitted) <- c("fit.0D","fit.3D","P")
	data.fitted <- data.fitted[order(data.fitted$P), ]
	
	# Subset by virus:
	sub <- NULL
	sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i])

	# Error in Prop Inf
	which.y <- paste("n.inf.na[", which(Prop.Inf$P.Virus == levels(Prop.Inf$P.Virus)[i]), "]", sep="")
	hdi.y <- as.data.frame(mat.or.vec(nr=length(which.y), nc=2))
	
	for(j in 1:length(which.y)){
		hdi.y[j,1:2] <- HDI(
					(Y.error[which(Y.error$Parameter == which.y[j]),]$value / n.recov[which(Prop.Inf$P.Virus == levels(Prop.Inf$P.Virus)[i])[j]]), 
					.95)
	}
	
	y.error.df <- data.frame(sub$Cad.Den.cm, hdi.y[,1], hdi.y[,2])
	colnames(y.error.df) <- c("Den", "low", "hi")

	# Error in Cad Density	
	which.x <- paste("P[", which(Prop.Inf$P.Virus == levels(Prop.Inf$P.Virus)[i])+(N.obs/2), "]", sep="")
	hdi.x <- as.data.frame(mat.or.vec(nr=length(which.x), nc=2))

	for(j in 1:length(which.x)){
		hdi.x[j,1:2] <- HDI(P.error[which(P.error$Parameter == which.x[j]),]$value, .95)
	}
	
	x.error.df <- data.frame(sub$P.Inf, hdi.x[,1], hdi.x[,2])
	colnames(x.error.df) <- c("P.Inf", "low", "hi")


	fit.plot <- NULL
	fit.plot <- ggplot(data.fitted, aes(x=P,y=fit.0D))+
		geom_line(linetype=2, size=1)+
		geom_line(aes(y=fit.3D), linetype=1, size=1)+
		geom_point(data=sub, aes(x=Cad.Den.cm, y=P.Inf, shape=P.Decay))+
		geom_segment(data=y.error.df, aes(y=low, yend=hi, x=Den, xend=Den))+
		geom_segment(data=x.error.df, aes(x=low, xend=hi, y=P.Inf, yend=P.Inf))+
		scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
		scale_x_continuous(limits=c(0,0.05))+
		#theme(legend.position="none")+
		ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))

plots.fitted4[[i]] <- fit.plot


}

library(gridExtra)

names(plots.fitted4)=paste("plot", 1:4)
arg.list <- c(plots.fitted4, list(nrow=2))
do.call(grid.arrange, arg.list)


