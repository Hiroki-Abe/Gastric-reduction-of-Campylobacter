setwd("data/")
d1 <- read.csv("yakitori.csv")
d2 <- read.csv("yakitori temp.csv")

#作図
par(mar=c(5,5,2,5.5))
par(mgp=c(2.6,1,0))
plot(d1$TIME[1:8], d1$MEAN[1:8], type="b", xlab="Heating time (min)", ylab="Survival cell count (log CFU/g)",
     lwd=3, cex=2.5, cex.lab=2, cex.axis=1.6, ylim=c(0, 6), xlim=c(0,4), yaxt="n")
axis(2,at=seq(0,6,1),las=1,cex.lab=2,cex.axis=1.6)
arrows(d1$TIME[1:8],d1$MEAN[1:8],d1$TIME[1:8],d1$MEAN[1:8]+d1$STDEV[1:8],angle=90,length=0.1, lwd=3)
arrows(d1$TIME[1:8],d1$MEAN[1:8],d1$TIME[1:8],d1$MEAN[1:8]-d1$STDEV[1:8],angle=90,length=0.1, lwd=3)
points(d1$TIME[9], d1$MEAN[9], cex=4.5, pch="*", lwd=3)


par(new=T)
plot(d2$TIME, d2$TEMP, xlim=c(0,4), ylim=c(0,80), xlab = "",ylab = "", 
     type="l", lty=2 , col=1,lwd=3,axes = FALSE)
axis(4,at=seq(0,80,20), ylim=c(0,80), las=1,cex.lab=2,cex.axis=1.6) 
mtext("Internal temperature (°C)",side=4,line=3.7, cex=2)
