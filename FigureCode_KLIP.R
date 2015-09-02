# Figures for each virus isolate and it's appropriate best-fit model

##################################################################
#-----------------------------------------------------------------
#-----------------------------------------------------------------
##################################################################

#####################
# KLIP # NULL MODEL
#####################


#####################
# FIGURE 1
#####################

library(ggplot2)
library(ggmcmc)
source("HDI.R")

nu.klip <- ggs(out.null, family="nu")
nu.klip <- subset(nu.klip, Parameter=="nu[1]")$value
nu.klip.mode <- median(nu.klip)

fake.P <- seq(0, 0.05, length.out=500)
fit.mean <- 1 - exp( -nu.klip.mode * fake.P * 7 )
  
# Bootstrap storage:
n.iter <- 1000
fit.klip <- mat.or.vec(nr=length(fake.P), nc=n.iter)
  
for(j in 1:n.iter){
    nu.temp <- sample(nu.klip,1)
    fit.klip[,j] <- 1 - exp( -nu.temp * fake.P * 7 )
}
  
fit.range <- mat.or.vec(nr=length(fake.P), nc=3)
  
for(j in 1:length(fake.P)){
    fit.range[j,1] <- fit.mean[j]
    fit.range[j,2:3] <- as.vector(quantile(fit.klip[j,], prob=c(0.025,0.975))[1:2])
    
}
  
data.fitted <- NULL;
data.fitted <- data.frame(fit.range, fake.P)
colnames(data.fitted) <- c("fit.mean","fit.low","fit.hi", "fake.P")
data.fitted <- data.fitted[order(data.fitted$fake.P), ]

sub <- NULL
sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[1]) #Subset out Klip

klip.plot1a <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=fit.mean), linetype=1, size=1)+
  geom_line(aes(y=fit.low), linetype=1, size=.85, color="gray")+
  geom_line(aes(y=fit.hi), linetype=1, size=.85, color="gray")+
  geom_point(data=sub, aes(x=Cad.Den.cm, y=P.Inf), shape=1)+
  theme_classic()+
  ggtitle("KLIP")+
  scale_y_continuous(limits=c(min(sub$P.Inf),1),breaks=c(0,.5,1))+
  scale_x_continuous(limits=c(0,0.05), breaks=c(0,.02,.04))+
  #theme(legend.position="none")+
  ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))
  
klip.plot1a

klip.plot1b <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=-log(1-fit.mean)), linetype=1, size=1)+
  geom_line(aes(y=-log(1-fit.hi)), linetype=1, size=.85, color="gray")+
  geom_line(aes(y=-log(1-fit.low)), linetype=1, size=.85, color="gray")+
  geom_point(data=sub, aes(x=Cad.Den.cm, y=-log(1-P.Inf)), shape=1)+
  theme_classic()+
  ggtitle("KLIP")+
  scale_y_continuous(limits=c(0,2),breaks=c(0,1,2))+
  scale_x_continuous(limits=c(0,0.05), breaks=c(0,.02,.04))+
  #theme(legend.position="none")+
  ylab("-ln(1 - Prop.Inf.)") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))

klip.plot1b

#####################
# FIGURE 2
#####################

# Find Y error:
se <- function(x) sqrt(var(x)/length(x))
se.y <- NULL; mean.y <- NULL
se.y[1] <- se(subset(sub, Cad.Den == 0)$P.Inf); mean.y[1] <- mean(subset(sub, Cad.Den == 0)$P.Inf)
se.y[2] <- se(subset(sub, Cad.Den == 10)$P.Inf); mean.y[2] <- mean(subset(sub, Cad.Den == 10)$P.Inf)
se.y[3] <- se(subset(sub, Cad.Den == 40)$P.Inf); mean.y[3] <- mean(subset(sub, Cad.Den == 40)$P.Inf)
se.y[4] <- se(-log(1-subset(sub, Cad.Den == 0)$P.Inf)); mean.y[4] <- mean(-log(1-subset(sub, Cad.Den == 0)$P.Inf))
se.y[5] <- se(-log(1-subset(sub, Cad.Den == 10)$P.Inf)); mean.y[5] <- mean(-log(1-subset(sub, Cad.Den == 10)$P.Inf))
se.y[6] <- se(-log(1-subset(sub, Cad.Den == 40)$P.Inf)); mean.y[6] <- mean(-log(1-subset(sub, Cad.Den == 40)$P.Inf))


# Find X error:
se.x <- NULL; mean.x <- NULL
se.x[1] <- se(subset(sub, Cad.Den == 0)$Cad.Den.cm); mean.x[1] <- mean(subset(sub, Cad.Den == 0)$Cad.Den.cm)
se.x[2] <- se(subset(sub, Cad.Den == 10)$Cad.Den.cm); mean.x[2] <- mean(subset(sub, Cad.Den == 10)$Cad.Den.cm)
se.x[3] <- se(subset(sub, Cad.Den == 40)$Cad.Den.cm); mean.x[3] <- mean(subset(sub, Cad.Den == 40)$Cad.Den.cm)
se.x[4] <- se(subset(sub, Cad.Den == 0)$Cad.Den.cm); mean.x[4] <- mean(subset(sub, Cad.Den == 0)$Cad.Den.cm)
se.x[5] <- se(subset(sub, Cad.Den == 10)$Cad.Den.cm); mean.x[5] <- mean(subset(sub, Cad.Den == 10)$Cad.Den.cm)
se.x[6] <- se(subset(sub, Cad.Den == 40)$Cad.Den.cm); mean.x[6] <- mean(subset(sub, Cad.Den == 40)$Cad.Den.cm)

errs <- data.frame(se.y,mean.y,se.x,mean.x)

klip.plot2a <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=fit.mean), linetype=1, size=1)+
  geom_line(aes(y=fit.low), linetype=1, size=.85, color="gray")+
  geom_line(aes(y=fit.hi), linetype=1, size=.85, color="gray")+
  geom_point(data = errs[1:3,], aes(x = mean.x, y = mean.y), shape=19, size=3)+
  geom_segment(data=errs[1:3,], aes(y=mean.y-se.y, yend=mean.y+se.y, x=mean.x, xend=mean.x),width=.5)+
  geom_segment(data=errs[1:3,], aes(y=mean.y, yend=mean.y, x=mean.x-se.x, xend=mean.x+se.x),width=.5)+
  theme_classic()+
  ggtitle("KLIP")+
  scale_y_continuous(limits=c(min(sub$P.Inf),.6),breaks=c(0,.3,.6))+
  scale_x_continuous(limits=c(0,0.035), breaks=c(0,0.01,0.03))+
  #theme(legend.position="none")+
  ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


klip.plot2a

klip.plot2b <- ggplot(data.fitted, aes(x=fake.P))+
  geom_line(aes(y=-log(1-fit.mean)), linetype=1, size=1)+
  geom_line(aes(y=-log(1-fit.hi)), linetype=5, size=.65, color="gray")+
  geom_line(aes(y=-log(1-fit.low)), linetype=5, size=.65, color="gray")+
  geom_point(data = errs[4:6,], aes(x = mean.x, y = mean.y), shape=19, size=3)+
  geom_segment(data=errs[4:6,], aes(y=mean.y-se.y, yend=mean.y+se.y, x=mean.x, xend=mean.x),width=.5)+
  geom_segment(data=errs[4:6,], aes(y=mean.y, yend=mean.y, x=mean.x-se.x, xend=mean.x+se.x),width=.5)+
  theme_classic()+
  ggtitle("KLIP")+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
  scale_x_continuous(limits=c(0,0.035), breaks=c(0,0.01,0.03))+
  #theme(legend.position="none")+
  ylab("-ln(1 - Prop.Inf.)") + xlab(expression(paste("Cad. Density (cadavers*", cm^-2, ")")))


klip.plot2b
