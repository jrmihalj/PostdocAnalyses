# Figures for each virus isolate and it's appropriate best-fit model

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

#####################
# NM # FULL MODEL
#####################


#####################
# FIGURE 1
#####################

library(ggplot2)
library(ggmcmc)
source("HDI.R")

C.nm <- ggs(out.full, family="C")
nu.nm <- ggs(out.full, family="nu.bar")
mu.nm <- ggs(out.full, family="mu")

C.vals <- subset(C.nm, Parameter=="C[3]")$value
nu.vals <- subset(nu.nm, Parameter=="nu.bar[3]")$value
mu.vals <- subset(mu.nm, Parameter=="mu[3]")$value

C.mode <- NULL; nu.mode <- NULL; mu.mode <- NULL;
C.mode <- mean(C.vals)
nu.mode <- mean(nu.vals)
mu.mode <- mean(mu.vals)

fake.P <- seq(0, 0.05, length.out=500)
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
sub <- subset(Prop.Inf, P.Virus=="NM")

nm.plot1a <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=fit.0Dmean),linetype=5, size=1)+
  geom_line(aes(y=fit.0Dlow), linetype=5, size=.85, color="gray")+
  geom_line(aes(y=fit.0Dhi), linetype=5, size=.85, color="gray")+
  geom_line(aes(y=fit.3Dmean), linetype=1, size=1)+
  geom_line(aes(y=fit.3Dlow), linetype=1, size=.85, color="gray")+
  geom_line(aes(y=fit.3Dhi), linetype=1, size=.85, color="gray")+
  geom_point(data=sub, aes(x=Cad.Den.cm, y=P.Inf, shape=P.Decay))+
  theme_classic()+
  ggtitle("NM")+
  scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
  scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
  scale_x_continuous(limits=c(0,0.05), breaks=c(0,.02,.04))+
  #theme(legend.position="none")+
  ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))

nm.plot1a

nm.plot1b <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=-log(1-fit.0Dmean)),linetype=5, size=1, color="gray")+
  #geom_line(aes(y=-log(1-fit.0Dlow)), linetype=5, size=.85, color="gray")+
  #geom_line(aes(y=-log(1-fit.0Dhi)), linetype=5, size=.85, color="gray")+
  geom_line(aes(y=-log(1-fit.3Dmean)), linetype=1, size=1)+
  #geom_line(aes(y=-log(1-fit.3Dlow)), linetype=1, size=.85, color="gray")+
  #geom_line(aes(y=-log(1-fit.3Dhi)), linetype=1, size=.85, color="gray")+
  geom_point(data=sub, aes(x=Cad.Den.cm, y=-log(1-P.Inf), shape=P.Decay))+
  theme_classic()+
  ggtitle("NM")+
  scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
  scale_y_continuous(limits=c(0,2),breaks=c(0,1,2))+
  scale_x_continuous(limits=c(0,0.05), breaks=c(0,.02,.04))+
  theme(legend.position="none")+
  ylab("-ln(1 - Prop.Inf.)") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))

nm.plot1b

#####################
# FIGURE 2
#####################

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

nm.plot2a <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=fit.0Dmean),linetype=1, size=1, color="gray")+
  #geom_line(aes(y=fit.0Dlow), linetype=5, size=.5, color="red")+
  #geom_line(aes(y=fit.0Dhi), linetype=5, size=.5, color="red")+
  geom_line(aes(y=fit.3Dmean), linetype=1, size=1)+
  #geom_line(aes(y=fit.3Dlow), linetype=3, size=.5)+
  #geom_line(aes(y=fit.3Dhi), linetype=3, size=.5)+
  geom_segment(data=errs, aes(y=mean.y-se.y, yend=mean.y+se.y, x=mean.x, xend=mean.x),width=.5)+
  geom_segment(data=errs, aes(y=mean.y, yend=mean.y, x=mean.x-se.x, xend=mean.x+se.x),width=.5)+
  geom_point(data = errs[1:3,], aes(x = mean.x, y = mean.y), shape=19, size=3, color="gray")+
  geom_point(data = errs[4:6,], aes(x = mean.x, y = mean.y), shape=17, size=3)+
  theme_classic()+
  ggtitle("NM")+
  scale_y_continuous(limits=c(min(sub$P.Inf),.6),breaks=c(0,.3,.6))+
  scale_x_continuous(limits=c(0,0.035), breaks=c(0,0.01,0.03))+
  ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


quartz()
plot.new()
nm.plot2a
legend(x=c("topleft"), legend=c("0 Day", "3 Day"), pch=c(19,17), pt.cex=c(1.25,1.25), 
       lty=c(1,1), lwd=c(2,2), col=c("gray","black"), bty="n")


### For plot2b, need to re-calculate the std.errs. ###

# Find Y error:
se.y <- NULL; mean.y <- NULL
se.y[1] <- se(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 0)$P.Inf)); mean.y[1] <- mean(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 0)$P.Inf))
se.y[2] <- se(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 10)$P.Inf)); mean.y[2] <- mean(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 10)$P.Inf))
se.y[3] <- se(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 40)$P.Inf)); mean.y[3] <- mean(-log(1-subset(sub, P.Decay=="0D" & Cad.Den == 40)$P.Inf))
se.y[4] <- se(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 0)$P.Inf)); mean.y[4] <- mean(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 0)$P.Inf))
se.y[5] <- se(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 10)$P.Inf)); mean.y[5] <- mean(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 10)$P.Inf))
se.y[6] <- se(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 40)$P.Inf)); mean.y[6] <- mean(-log(1-subset(sub, P.Decay=="3D" & Cad.Den == 40)$P.Inf))   

errs <- data.frame(se.y,mean.y,se.x,mean.x)

nm.plot2b <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=-log(1-fit.0Dmean)),linetype=1, size=1, color="gray")+
  #geom_line(aes(y=fit.0Dlow), linetype=5, size=.5, color="red")+
  #geom_line(aes(y=fit.0Dhi), linetype=5, size=.5, color="red")+
  geom_line(aes(y=-log(1-fit.3Dmean)), linetype=1, size=1)+
  #geom_line(aes(y=fit.3Dlow), linetype=3, size=.5)+
  #geom_line(aes(y=fit.3Dhi), linetype=3, size=.5)+
  geom_segment(data=errs, aes(y=mean.y-se.y, yend=mean.y+se.y, x=mean.x, xend=mean.x),width=.5)+
  geom_segment(data=errs, aes(y=mean.y, yend=mean.y, x=mean.x-se.x, xend=mean.x+se.x),width=.5)+
  geom_point(data = errs[1:3,], aes(x = mean.x, y = mean.y), shape=19, size=3, color="gray")+
  geom_point(data = errs[4:6,], aes(x = mean.x, y = mean.y), shape=17, size=3)+
  theme_classic()+
  ggtitle("NM")+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
  scale_x_continuous(limits=c(0,0.035), breaks=c(0,0.01,0.03))+
  #theme(legend.position="none")+
  ylab("-ln(1 - Prop.Inf.)") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


quartz()
plot.new()
nm.plot2b
legend(x=c("topleft"), legend=c("0 Day", "3 Day"), pch=c(19,17), pt.cex=c(1.25,1.25), 
       lty=c(1,1), lwd=c(2,2), col=c("gray","black"), bty="n")
