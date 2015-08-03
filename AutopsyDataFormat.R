###################
#IMPORT AND FORMAT#
###################

recov <- read.csv(file.choose(), header=T)
head(recov)

recov$Decay[1:7] <- "0D"
recov$Dose[1:7] <- "HI"

vir <- NULL
Virus <- NULL
dec <- NULL
Decay <- NULL
dos <- NULL
Dose <- NULL
rep <- NULL
Rep <- NULL
ind <- NULL
Indiv <- NULL

n.bags <-119

for(i in 1:n.bags){
	
	n.indiv <- recov[i,11]
	if(n.indiv > 0){
		vir <- NULL
		dec <- NULL
		dos <- NULL
		rep <- NULL
		ind <- NULL

		vir <- paste(rep(recov[i,3], times=n.indiv))
		dec <- paste(rep(recov[i,4], times=n.indiv))
		dos <- paste(rep(recov[i,5], times=n.indiv))
		rep <- paste(rep(recov[i,6], times=n.indiv))
		ind <- 1:n.indiv

		Virus <- c(Virus, vir)
		Decay <- c(Decay, dec)
		Dose <- c(Dose, dos)
		Rep <- c(Rep, rep)
		Indiv <- c(Indiv, ind)
	}
}

Autopsy.Data <- data.frame(Virus,Decay,Dose,Rep,Indiv)

write.table(Autopsy.Data, file="Autopsy.Data.csv", row.names=F, sep=",") 
getwd()

#########################################
# CREATE FIGURES OF PROPORTION INFECTED #
#########################################

data <- read.csv(file.choose(), header=T)
head(data)

# Fill in 0's for NA's in Virus.Pres
data$Virus.Pres[is.na(data$Virus.Pres)] <- 0

# Create a new table, which has the proportion infected for each treatment
Treat <- levels(data$Virus)
Dec <- levels(data$Decay)
Dos <- levels(data$Dose)

# STORAGE
P.Virus <- NULL
P.Decay <- NULL
P.Dose <- NULL
P.Rep <- NULL
P.Inf <- NULL
N.Inf <- NULL
N.Recov <- NULL
p.count <- 1



for(t in 1:length(Treat)){
	if(t == 1){
		sub <- NULL
		sub <- subset(data, Virus == Treat[1])
		n.rep1 <- length(levels(factor(sub$Rep)))
		for(r1 in 1:n.rep1){
			pres <- NULL
			pres <- subset(data, Virus == Treat[1] & Rep == r1)$Virus.Pres
			P.Inf[p.count] <- length(which(pres == 1)) / length(pres)
			P.Virus[p.count] <- Treat[1]
			P.Decay[p.count] <- NA
			P.Dose[p.count] <- NA
			P.Rep[p.count] <- r1
			N.Inf[p.count] <- length(which(pres == 1))
			N.Recov[p.count] <- length(pres)
			p.count <- p.count + 1
		}
	}else{
		for(de in 1:length(Dec)){
			for(do in 1:length(Dos)){
				n.reps <- max(unique(subset(data, Virus == Treat[t] & Decay == Dec[de] & Dose == Dos[do])$Rep)) 
				for(r in 1:n.reps){
					pres <- NULL
					pres <- subset(data, Virus == Treat[t] & Decay == Dec[de] & Dose == Dos[do] & Rep == r)$Virus.Pres
					P.Inf[p.count] <- length(which(pres == 1)) / length(pres)
					P.Virus[p.count] <- Treat[t]
					P.Decay[p.count] <- Dec[de]
					P.Dose[p.count] <- Dos[do]
					P.Rep[p.count] <- r
					N.Inf[p.count] <- length(which(pres == 1))
					N.Recov[p.count] <- length(pres)
					p.count <- p.count + 1
				}
			}
		}


	}
}

Prop.Inf <- data.frame(P.Virus, P.Decay, P.Dose, P.Rep, P.Inf, N.Inf, N.Recov)
Prop.Inf

#Remove bags that had zero individuals recovered
Prop.Inf <- Prop.Inf[-which(is.na(Prop.Inf$P.Inf)), ]
Prop.Inf

#Need to make the Controls = to Zero Cadavers for each virus
cont <- subset(Prop.Inf, P.Virus == "CONTROL")
cont$P.Decay <- "0D"
cont$P.Dose <- 0

cont <- rbind(cont,cont)
cont$P.Decay[8:14] <- "3D"

fix.cont <- mat.or.vec(nr=1, nc=ncol(cont))
colnames(fix.cont) <- colnames(cont)

# Need to repeat this "cont" for each virus and relabel P.Virus
for(v in 2:length(levels(Prop.Inf$P.Virus))){
	temp <- cont
	temp$P.Virus <- levels(Prop.Inf$P.Virus)[v]
	fix.cont <-rbind(fix.cont, temp)
}

fix.cont <- fix.cont[-1, ]

# Add to Prop.Inf

Prop.Inf <- rbind(Prop.Inf, fix.cont)
Prop.Inf <- Prop.Inf[-c(1:7),]

###Format to plot

#1. Reformat P.Dose

Prop.Inf$Cad.Den <- 0
for(i in 1:nrow(Prop.Inf)){
	if(Prop.Inf$P.Dose[i] == "LOW" | is.na(Prop.Inf$P.Dose[i])){
		Prop.Inf$Cad.Den[i] <- 10} else {Prop.Inf$Cad.Den[i] <- 40}	
}

#2.Rename rows
rownames(Prop.Inf) <- 1:nrow(Prop.Inf)

#3. Fix Lowest Density (Controls)
Prop.Inf$Cad.Den[106:nrow(Prop.Inf)] <- 0

#4 Re-factor P.Virus
Prop.Inf$P.Virus <- factor(Prop.Inf$P.Virus)

# GENERATE MEANS AND SE FOR PROP INF
errors <- summarySE(Prop.Inf, measurevar="P.Inf", groupvars=c("P.Virus","P.Decay","Cad.Den"))
errors

################################
# GRAPH THE PROP INF PER VIRUS 
# ASSUMING CAD. DENSITY IS FIXED (0, 10, 40)
################################

library(ggplot2)
library(scales)

plots <- list()


for(i in 1:4){

sub.main <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i])
sub.error <- subset(errors, P.Virus == levels(Prop.Inf$P.Virus)[i])

PI.plot <- ggplot(sub.main, aes(x=Cad.Den, y=P.Inf, shape=P.Decay))+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		#geom_jitter(aes(alpha=0.4), position = position_jitter(width = .5, height=0))+
		geom_errorbar(data=sub.error, aes(ymin=P.Inf-se, ymax=P.Inf+se), width=1)+
		geom_point(data=sub.error, aes(size=2))+
		scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
		geom_line(data=sub.error)+
		scale_y_continuous(limits=c(min(P.Inf),1),breaks=c(0,.5,1))+
		guides(alpha=FALSE, size=FALSE)+
		#theme(legend.position="none")+
		ylab("Prop. Infected") + xlab("Cadaver Density")

plots[[i]] <- PI.plot
}


#NOW PLOT USING GRID

library(gridExtra)

names(plots)=paste("plot", 1:4)
arg.list <- c(plots, list(nrow=2))
windows(record=T)
do.call(grid.arrange, arg.list)

####################################
# GRAPH THE PROP INF PER VIRUS 
# ASSUMING CAD. DENSITY IS CONTINUOUS
####################################

# First, import branch area data:
b.area <- read.csv(file.choose(), header=T)
head(b.area)

# Need to put in the branch area data for each
# bag in the Prop.Inf table

# Storage:
Prop.Inf$Picture.ID <- 0
Prop.Inf$Area.in2 <- 0

#Deal with controls after...
for (i in 1:105){
	sub <- subset(b.area, Virus==paste(Prop.Inf$P.Virus[i]) & Decay == paste(Prop.Inf$P.Decay[i])
					& Dose == paste(Prop.Inf$P.Dose[i]) & Rep == paste(Prop.Inf$P.Rep[i]))

	Prop.Inf$Picture.ID[i] <- sub$Picture.ID
	Prop.Inf$Fan.Cyl[i] <- paste(sub$Fan.Cyl)
	Prop.Inf$Area.in2[i] <- sub$Area.in2
}

#Control Bags
cont.area <- b.area[which(b.area$Virus=="CONTROL"),c(5,7:9)]
cont.area <- cont.area[order(cont.area$Rep),]

Prop.Inf[106:nrow(Prop.Inf), 9:10] <- rep(cont.area[,c(2,4)], times=8)
Prop.Inf[106:nrow(Prop.Inf), 11] <- rep(paste(cont.area[,3]),times=8)

# Convert in2 to m2
Prop.Inf$Area.m2 <- Prop.Inf$Area.in2 * 0.00064516 

# CREATE A NEW CAD.DEN COLUMN
Prop.Inf$Cad.Den.Adj <- Prop.Inf$Cad.Den / Prop.Inf$Area.m2

# Center and scale?
#Prop.Inf$Cad.Den.Adj <- (Prop.Inf$Cad.Den.Adj - mean(Prop.Inf$Cad.Den.Adj, na.rm=T)) / sd(Prop.Inf$Cad.Den.Adj, na.rm=T)

##################################
# PLOT DATA WITH LOESS SMOOTHER? #
##################################

library(ggplot2)

plots.2 <- list()

for(i in 1:4){

	sub <- subset(Prop.Inf, P.Virus==levels(Prop.Inf$P.Virus)[i] & P.Inf < 1)

	PI.plot2 <- ggplot(sub, aes(x=Cad.Den.Adj, y=P.Inf, shape=P.Decay))+
		theme_classic()+
		ggtitle(paste(levels(Prop.Inf$P.Virus)[i]))+
		geom_point()+
		scale_shape_manual(name="Decay", labels=c("0 Days","3 Days"),values=c(1,17))+
		scale_y_continuous(limits=c(min(P.Inf),1),breaks=c(0,.5,1))+
		scale_x_continuous(limits=c(0,450))+
		guides(alpha=FALSE, size=FALSE)+
		scale_linetype_manual(name="Decay", labels=c("0 Days","3 Days"), values=c(1,2))+
		#theme(legend.position="none")+
		geom_smooth(method="loess", se=F, aes(linetype=P.Decay), color="black")+
		ylab("Prop. Infected") + xlab(expression(paste("Cad. Density (cadavers*", m^-2, ")")))

plots.2[[i]] <- PI.plot2


}

#NOW PLOT USING GRID

library(gridExtra)

names(plots.2)=paste("plot", 1:4)
arg.list.2 <- c(plots.2, list(nrow=2))
do.call(grid.arrange, arg.list.2)


#############
# FUNCTIONS #
#############

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


