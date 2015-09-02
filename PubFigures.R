# Publication quality figures

#####################
# Figure 1, version 1
#####################

library(gridExtra)

quartz(height=8, width=8, dpi=300)
grid.arrange(klip.plot1b, tmb.plot1b, love.plot1b, nm.plot1b, nrow=2)

quartz(height=8, width=8, dpi=300)
plot.new()
legend(locator(1), legend=c("0 Day", "3 Day"), pch=c(1,17), pt.cex=c(1.15,1.15), 
       lty=c(1,1), col=c("gray","black"), lwd=c(2,2), bty="n")

graphics.off()
#####################
# Figure 1, version 2
#####################

library(gridExtra)

quartz(height=8, width=8, dpi=300)
grid.arrange(klip.plot2b, tmb.plot2b, love.plot2b, nm.plot2b, nrow=2)

quartz(height=8, width=8, dpi=300)
plot.new()
legend(locator(1), legend=c("0 Day", "3 Day"), pch=c(19,17), pt.cex=c(1.15,1.15), 
       lty=c(1,1), col=c("gray","black"), lwd=c(2,2), bty="n")

graphics.off()
