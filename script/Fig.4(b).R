library(MASS)

min <- 0.025
max <- 0.975
#実測値のデータ
ds <- read.csv("simulation_liver_data.csv")

N_mcmc <- length(ms$lp__)
alld_mcmc <- NULL
#pHデータの取得
dpH <- read.csv("pH_liver_data.csv")

#0分を入力するための操作
y <- -0
min_y <- quantile(y, probs=min)
min_y <- as.numeric(min_y)
max_y <- quantile(y, probs=max)
max_y <- as.numeric(max_y)
med_y <- median(y)
d_mcmc <- data.frame(time=0, min_y, max_y, med_y)
alld_mcmc <- rbind(alld_mcmc, d_mcmc)

r_abef <- matrix(0,length(ms$lp__),4)
for(n in 1:length(ms$lp__))
{
  r_abef[n,] <- mvrnorm(1,mu=ms$abef0[n,], Sigma=ms$cov[n,,])
}
a <- r_abef[,1]
b <- r_abef[,2]
e <- r_abef[,3]
f <- r_abef[,4]

for(i in seq(1,180,1)){
  #pHデータから各時間のpHを推定
  pH1 <- dpH$pH[floor(i/10)+1]-(dpH$pH[floor(i/10)+1]-dpH$pH[ceiling(i/10)+1])*((i-10*floor(i/10))/10)
  pH2 <- dpH$pH[floor((i-1)/10)+1]-(dpH$pH[floor((i-1)/10)+1]-dpH$pH[ceiling((i-1)/10)+1])*((i-10*floor((i-1)/10))/10)
  pH <- (pH1+pH2)/2
  d <- exp(a*pH+b)
  p <- exp(e*pH+f)
  #前回の生残率から今回のpHの場合何分に相当する減少率か計算
  medt <- d*((-alld_mcmc[i/1, 4])^(1/p))
  y <- -((medt+1)/exp(d))^p + (medt/exp(d))^p 
  rowmin_y <- quantile(y, probs=min)
  min_y <- as.numeric(rowmin_y) + alld_mcmc[i/1, 2]
  rowmax_y <- quantile(y, probs=max)
  max_y <- as.numeric(rowmax_y) + alld_mcmc[i/1, 3]
  med_y <- median(y) +  alld_mcmc[i/1, 4]
  d_mcmc <- data.frame(time=i, min_y, max_y, med_y)
  alld_mcmc <- rbind(alld_mcmc, d_mcmc)
}


#作図
par(mar=c(5,5,2,5.5))
par(mgp=c(3.4,1,0))
plot(alld_mcmc[,1], alld_mcmc[,4], type="l", xlab="Model digestion duration (min)", ylab="Survival ratio (log N/N0)",
     lwd=3, cex.lab=2, cex.axis=1.6, ylim=c(-6.2, 0), xlim=c(0,180), xaxp=c(0,180,3), yaxt="n")
axis(2,at=seq(-6,0,1),las=1,cex.lab=2,cex.axis=1.6)
points(alld_mcmc[,1], alld_mcmc[,2], type="l", lty=2, lwd=3)
points(alld_mcmc[,1], alld_mcmc[,3], type="l", lty=2, lwd=3)
points(ds$TIME[1:7], ds$MEAN[1:7], cex=2.5, lwd=3)
arrows(ds$TIME[1:7],ds$MEAN[1:7],ds$TIME[1:7],ds$MEAN[1:7]+ds$STDEV[1:7],angle=90,length=0.1, lwd=3)
arrows(ds$TIME[1:7],ds$MEAN[1:7],ds$TIME[1:7],ds$MEAN[1:7]-ds$STDEV[1:7],angle=90,length=0.1, lwd=3)

dpH[1,] <- c(0,1.5,0)
par(new=T)
plot(dpH$TIME, dpH$pH, xlim=c(0,180), ylim=c(1.5,6.2), xlab = "",ylab = "", 
     type="l", lty=3 , col=1, lwd=3,axes = FALSE)
axis(4,at=seq(2,6,1), ylim=c(1.5,6.2), las=1, cex.axis=1.6) 
mtext("pH",side=4,line=3, cex=2)


#RMSE calculations
(mean((ds$MEAN-alld_mcmc[ds$TIME+1,4])^2))^0.5

