source("SMLE.r") 

#****************DATA*******************#
data_firmcost <- read.table("data_RiskSurvey.txt",h=F)
attach(data_firmcost)

#***Covariates***#
#----------------#
SIZELOG <- as.numeric(V4);
INDCOST <- as.numeric(V5);
#----------------#

#**************************************************#
#******Beta regression with varying precision*****#
#**************************************************#
y <- as.numeric(V1/100) #response variable
X <- matrix(c(rep(1,73), INDCOST, SIZELOG), ncol=3, byrow=F) #regressor matrix for the mean submodel
Z <- X #regressor matrix for the precision submodel

#******************************Subsets without outliers*******************************************#
y15 <- y[-c(15)]; X15 <- X[-c(15),];Z15 <- as.matrix(Z[-c(15),])
y16 <- y[-c(16)]; X16 <- X[-c(16),];Z16 <- as.matrix(Z[-c(16),])
y72 <- y[-c(72)]; X72 <- X[-c(72),];Z72 <- as.matrix(Z[-c(72),])
y15_16 <- y[-c(15,16)]; X15_16 <- X[-c(15,16),];Z15_16 <- as.matrix(Z[-c(15,16),])
y15_72 <- y[-c(15,72)]; X15_72 <- X[-c(15,72),];Z15_72 <- as.matrix(Z[-c(15,72),])
y16_72 <- y[-c(16,72)]; X16_72 <- X[-c(16,72),];Z16_72 <- as.matrix(Z[-c(16,72),])
y15_16_72 <- y[-c(15,16,72)]; X15_16_72 <- X[-c(15,16,72),];Z15_16_72 <- as.matrix(Z[-c(15,16,72),])
#**************************************************************************************************#

#********************************************MLE FITS*****************************************************#
# MLE FIT FULL DATA
fitMLE <- SMLE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE

# MLE FIT WITHOUT 15
fitMLE_15 <- SMLE_BETA(y=y15, X=X15, Z=Z15, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",  
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15

# MLE FIT WITHOUT 16
fitMLE_16 <- SMLE_BETA(y=y16, X=X16, Z=Z16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_16

# MLE FIT WITHOUT 72
fitMLE_72 <- SMLE_BETA(y=y72, X=X72, Z=Z72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_72

# MLE FIT WITHOUT 15 AND 16
fitMLE_15_16 <- SMLE_BETA(y=y15_16, X=X15_16, Z=Z15_16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_16

# MLE FIT WITHOUT 15 AND 72
fitMLE_15_72 <- SMLE_BETA(y=y15_72, X=X15_72, Z=Z15_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_72

# MLE FIT WITHOUT 16 AND 72
fitMLE_16_72 <- SMLE_BETA(y=y16_72, X=X16_72, Z=Z16_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",  
startV="CP", linkmu="logit", linkphi="log")
fitMLE_16_72

# MLE FIT WITHOUT 15, 16 AND 72
fitMLE_15_16_72 <- SMLE_BETA(y=y15_16_72, X=X15_16_72, Z=Z15_16_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,
method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_16_72

#************************SMLE*********************#

# SMLE FIT FULL DATA
fitRob <- SMLE_BETA(y=y, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob

# SMLE FIT WITHOUT 15
fitRob_15 <- SMLE_BETA(y=y15, X=X15, Z=Z15, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15

# SMLE FIT WITHOUT 16
fitRob_16 <- SMLE_BETA(y=y16, X=X16, Z=Z16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_16

# SMLE FIT WITHOUT 72
fitRob_72 <- SMLE_BETA(y=y72, X=X72, Z=Z72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_72

# SMLE FIT WITHOUT 15 AND 16
fitRob_15_16 <- SMLE_BETA(y=y15_16, X=X15_16, Z=Z15_16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_16

# SMLE FIT WITHOUT 15 AND 72
fitRob_15_72 <- SMLE_BETA(y=y15_72, X=X15_72, Z=Z15_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_72

# SMLE FIT WITHOUT 16 AND 72
fitRob_16_72 <- SMLE_BETA(y=y16_72, X=X16_72, Z=Z16_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_16_72

# SMLE FIT WITHOUT 15, 16 AND 72
fitRob_15_16_72 <- SMLE_BETA(y=y15_16_72, X=X15_16_72, Z=Z15_16_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,
method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_16_72

#*************************SCATTER PLOTS*******************************#
value_indcost <- seq(0.09,1.22,length.out=1000)
value_sizelog <- seq(5.27,10.6,length.out=1000)

#*************INDCOST*************#

#*****MLE*****#

#full data
beta_mle <- fitMLE$beta; gama_mle <- fitMLE$gama
eta_mle <- beta_mle[1] + beta_mle[2]*value_indcost + beta_mle[3]*median(SIZELOG)
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
#phi_mle <- exp(gama_mle)

#without 15
beta_mle_w15 <- fitMLE_15$beta; gama_mle_w15 <- fitMLE_15$gama
eta_mle_w15 <- beta_mle_w15[1] + beta_mle_w15[2]*value_indcost + beta_mle_w15[3]*median(SIZELOG)
mu_mle_w15 <- exp(eta_mle_w15)/(1+exp(eta_mle_w15))
#phi_mle_w15 <- exp(gama_mle_w15)

#without 16
beta_mle_w16 <- fitMLE_16$beta; gama_mle_w16 <- fitMLE_16$gama
eta_mle_w16 <- beta_mle_w16[1] +  beta_mle_w16[2]*value_indcost + beta_mle_w16[3]*median(SIZELOG)
mu_mle_w16 <- exp(eta_mle_w16)/(1+exp(eta_mle_w16))
#phi_mle_w16 <- exp(gama_mle_w16)

#without 72
beta_mle_w72 <- fitMLE_72$beta; gama_mle_w72 <- fitMLE_72$gama
eta_mle_w72 <- beta_mle_w72[1] + beta_mle_w72[2]*value_indcost + beta_mle_w72[3]*median(SIZELOG)
mu_mle_w72 <- exp(eta_mle_w72)/(1+exp(eta_mle_w72))
#phi_mle_w72 <- exp(gama_mle_w72)

#without 15 and 16
beta_mle_w1516 <- fitMLE_15_16$beta; gama_mle_w1516 <- fitMLE_15_16$gama
eta_mle_w1516 <- beta_mle_w1516[1] + beta_mle_w1516[2]*value_indcost + beta_mle_w1516[3]*median(SIZELOG)
mu_mle_w1516 <- exp(eta_mle_w1516)/(1+exp(eta_mle_w1516))
#phi_mle_w1516 <- exp(gama_mle_w1516)

#without 15 and 72
beta_mle_w1572 <- fitMLE_15_72$beta; gama_mle_w1572 <- fitMLE_15_72$gama
eta_mle_w1572 <- beta_mle_w1572[1] + beta_mle_w1572[2]*value_indcost + beta_mle_w1572[3]*median(SIZELOG)
mu_mle_w1572 <- exp(eta_mle_w1572)/(1+exp(eta_mle_w1572))
#phi_mle_w1572 <- exp(gama_mle_w1572)

#without 16 and 72
beta_mle_w1672 <- fitMLE_16_72$beta; gama_mle_w1672 <- fitMLE_16_72$gama
eta_mle_w1672 <- beta_mle_w1672[1] + beta_mle_w1672[2]*value_indcost + beta_mle_w1672[3]*median(SIZELOG)
mu_mle_w1672 <- exp(eta_mle_w1672)/(1+exp(eta_mle_w1672))
#phi_mle_w1672 <- exp(gama_mle_w1672)

#without 15, 16 and 72
beta_mle_w151672 <- fitMLE_15_16_72$beta; gama_mle_w151672 <- fitMLE_15_16_72$gama
eta_mle_w151672 <- beta_mle_w151672[1] + beta_mle_w151672[2]*value_indcost + beta_mle_w151672[3]*median(SIZELOG)
mu_mle_w151672 <- exp(eta_mle_w151672)/(1+exp(eta_mle_w151672))
phi_mle_w151672 <- exp(gama_mle_w151672)


#****SMLE***#

#full data
beta_smle <- fitRob$beta; gama_smle <- fitRob$gama
eta_smle <- beta_smle[1] +  beta_smle[2]*value_indcost + beta_smle[3]*median(SIZELOG)
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
#phi_smle <- exp(gama_smle)

#without 15
beta_smle_w15 <- fitRob_15$beta; gama_smle_w15 <- fitRob_15$gama
eta_smle_w15 <- beta_smle_w15[1] + beta_smle_w15[2]*value_indcost + beta_smle_w15[3]*median(SIZELOG)
mu_smle_w15 <- exp(eta_smle_w15)/(1+exp(eta_smle_w15))
#phi_smle_w15 <- exp(gama_smle_w15)

#without 16
beta_smle_w16 <- fitRob_16$beta; gama_smle_w16 <- fitRob_16$gama
eta_smle_w16 <- beta_smle_w16[1] +  beta_smle_w16[2]*value_indcost + beta_smle_w16[3]*median(SIZELOG)
mu_smle_w16 <- exp(eta_smle_w16)/(1+exp(eta_smle_w16))
#phi_smle_w16 <- exp(gama_smle_w16)

#without 72
beta_smle_w72 <- fitRob_72$beta; gama_smle_w72 <- fitRob_72$gama
eta_smle_w72 <- beta_smle_w72[1] + beta_smle_w72[2]*value_indcost + beta_smle_w72[3]*median(SIZELOG)
mu_smle_w72 <- exp(eta_smle_w72)/(1+exp(eta_smle_w72))
#phi_smle_w72 <- exp(gama_smle_w72)

#without 15 and 16
beta_smle_w1516 <- fitRob_15_16$beta; gama_smle_w1516 <- fitRob_15_16$gama
eta_smle_w1516 <- beta_smle_w1516[1] + beta_smle_w1516[2]*value_indcost + beta_smle_w1516[3]*median(SIZELOG)
mu_smle_w1516 <- exp(eta_smle_w1516)/(1+exp(eta_smle_w1516))
#phi_smle_w1516 <- exp(gama_smle_w1516)

#without 15 and 72
beta_smle_w1572 <- fitRob_15_72$beta; gama_smle_w1572 <- fitRob_15_72$gama
eta_smle_w1572 <- beta_smle_w1572[1] + beta_smle_w1572[2]*value_indcost + beta_smle_w1572[3]*median(SIZELOG)
mu_smle_w1572 <- exp(eta_smle_w1572)/(1+exp(eta_smle_w1572))
#phi_smle_w1572 <- exp(gama_smle_w1572)

#without 16 and 72
beta_smle_w1672 <- fitRob_16_72$beta; gama_smle_w1672 <- fitRob_16_72$gama
eta_smle_w1672 <- beta_smle_w1672[1] + beta_smle_w1672[2]*value_indcost + beta_smle_w1672[3]*median(SIZELOG)
mu_smle_w1672 <- exp(eta_smle_w1672)/(1+exp(eta_smle_w1672))
#phi_smle_w1672 <- exp(gama_smle_w1672)

#without 15, 16 and 72
beta_smle_w151672 <- fitRob_15_16_72$beta; gama_smle_w151672 <- fitRob_15_16_72$gama
eta_smle_w151672 <- beta_smle_w151672[1] + beta_smle_w151672[2]*value_indcost + beta_smle_w151672[3]*median(SIZELOG)
mu_smle_w151672 <- exp(eta_smle_w151672)/(1+exp(eta_smle_w151672))
#phi_smle_w151672 <- exp(gama_smle_w151672)

#******SMLE****#

par(mar=c(5.0,5.0,4.0,2.0))
plot(INDCOST,y,pch=16,xlab="Ind cost",ylab="Firm cost",ylim=c(0,1), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(value_indcost,mu_smle,         lwd=2, col=1, lty=1)
lines(value_indcost,mu_smle_w15,     lwd=4, col=2, lty=2)
lines(value_indcost,mu_smle_w16,     lwd=2, col=6, lty=6)
lines(value_indcost,mu_smle_w72,     lwd=2, col=4, lty=4)
lines(value_indcost,mu_smle_w1516,   lwd=2, col=5, lty=5)
lines(value_indcost,mu_smle_w1572,   lwd=3, col=3, lty=1)
lines(value_indcost,mu_smle_w1672,   lwd=2, col="orangered", lty=3)
lines(value_indcost,mu_smle_w151672, lwd=3, col=8, lty=4)
identify(INDCOST,y,cex=1.3, n=3)
legend(0.4, 0.85, c("SMLE","SMLE w/o 15","SMLE w/o 16","SMLE w/o 72", "SMLE w/o 15, 16",
"SMLE w/o 15, 72", "SMLE w/o 16, 72", "SMLE w/o 15, 16, 72"),
col=c(1,2,6,4,5,3,"orangered",8), lty=c(1,2,6,4,5,1,3,4),cex=1.0,lwd=c(2,4,2,2,2,3,2,3))


#******MLE****#

par(mar=c(5.0,5.0,4.0,2.0))
plot(INDCOST,y,pch=16,xlab="Ind cost",ylab="Firm cost",ylim=c(0,1), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(value_indcost,mu_mle,         lwd=2, col=1, lty=1)
lines(value_indcost,mu_mle_w15,     lwd=4, col=2, lty=2)
lines(value_indcost,mu_mle_w16,     lwd=2, col=6, lty=6)
lines(value_indcost,mu_mle_w72,     lwd=2, col=4, lty=4)
lines(value_indcost,mu_mle_w1516,   lwd=2, col=5, lty=5)
lines(value_indcost,mu_mle_w1572,   lwd=3, col=3, lty=1)
lines(value_indcost,mu_mle_w1672,   lwd=2, col="orangered", lty=3)
lines(value_indcost,mu_mle_w151672, lwd=3, col=8, lty=4)
identify(INDCOST,y,cex=1.3, n=3)
legend(0.4, 0.85,c("MLE","MLE w/o 15","MLE w/o 16","MLE w/o 72", "MLE w/o 15, 16",
"MLE w/o 15, 72", "MLE w/o 16, 72", "MLE w/o 15, 16, 72"),
col=c(1,2,6,4,5,3,"orangered",8), lty=c(1,2,6,4,5,1,3,4),cex=1.0,lwd=c(2,4,2,2,2,3,2,3))

#***********scatter plot in introduction section*******#
par(mar=c(5.0,5.0,4.0,2.0))
plot(INDCOST,y,pch=16,xlab="Ind cost",ylab="Firm cost",ylim=c(0,1), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(value_indcost,mu_mle,    lwd=2,col=1,lty=1)
lines(value_indcost,mu_mle_w15,lwd=4,col=2,lty=2)
lines(value_indcost,mu_mle_w16,lwd=2,col=6,lty=6)
lines(value_indcost,mu_mle_w72,lwd=2,col=4,lty=4)
identify(INDCOST,y,cex=1.3)
legend(0.45,0.85,c("MLE","MLE w/o 15","MLE w/o 16","MLE w/o 72"),
col=c(1,2,6,4),lty=c(1,2,6,4),cex=1.0,lwd=c(2,4,2,2))


#******************SIZELOG****************#

#*****MLE*****#

#full data
beta_mle <- fitMLE$beta; gama_mle <- fitMLE$gama
eta_mle <- beta_mle[1] + beta_mle[2]*median(INDCOST) + beta_mle[3]*value_sizelog
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
#phi_mle <- exp(gama_mle)

#without 15
beta_mle_w15 <- fitMLE_15$beta; gama_mle_w15 <- fitMLE_15$gama
eta_mle_w15 <- beta_mle_w15[1] + beta_mle_w15[2]*median(INDCOST) + beta_mle_w15[3]*value_sizelog
mu_mle_w15 <- exp(eta_mle_w15)/(1+exp(eta_mle_w15))
#phi_mle_w15 <- exp(gama_mle_w15)

#without 16
beta_mle_w16 <- fitMLE_16$beta; gama_mle_w16 <- fitMLE_16$gama
eta_mle_w16 <- beta_mle_w16[1] +  beta_mle_w16[2]*median(INDCOST) + beta_mle_w16[3]*value_sizelog
mu_mle_w16 <- exp(eta_mle_w16)/(1+exp(eta_mle_w16))
#phi_mle_w16 <- exp(gama_mle_w16)

#without 72
beta_mle_w72 <- fitMLE_72$beta; gama_mle_w72 <- fitMLE_72$gama
eta_mle_w72 <- beta_mle_w72[1] + beta_mle_w72[2]*median(INDCOST) + beta_mle_w72[3]*value_sizelog
mu_mle_w72 <- exp(eta_mle_w72)/(1+exp(eta_mle_w72))
#phi_mle_w72 <- exp(gama_mle_w72)

#without 15 and 16
beta_mle_w1516 <- fitMLE_15_16$beta; gama_mle_w1516 <- fitMLE_15_16$gama
eta_mle_w1516 <- beta_mle_w1516[1] + beta_mle_w1516[2]*median(INDCOST) + beta_mle_w1516[3]*value_sizelog
mu_mle_w1516 <- exp(eta_mle_w1516)/(1+exp(eta_mle_w1516))
#phi_mle_w1516 <- exp(gama_mle_w1516)

#without 15 and 72
beta_mle_w1572 <- fitMLE_15_72$beta; gama_mle_w1572 <- fitMLE_15_72$gama
eta_mle_w1572 <- beta_mle_w1572[1] + beta_mle_w1572[2]*median(INDCOST) + beta_mle_w1572[3]*value_sizelog
mu_mle_w1572 <- exp(eta_mle_w1572)/(1+exp(eta_mle_w1572))
#phi_mle_w1572 <- exp(gama_mle_w1572)

#without 16 and 72
beta_mle_w1672 <- fitMLE_16_72$beta; gama_mle_w1672 <- fitMLE_16_72$gama
eta_mle_w1672 <- beta_mle_w1672[1] + beta_mle_w1672[2]*median(INDCOST) + beta_mle_w1672[3]*value_sizelog
mu_mle_w1672 <- exp(eta_mle_w1672)/(1+exp(eta_mle_w1672))
#phi_mle_w1672 <- exp(gama_mle_w1672)

#without 15, 16 and 72
beta_mle_w151672 <- fitMLE_15_16_72$beta; gama_mle_w151672 <- fitMLE_15_16_72$gama
eta_mle_w151672 <- beta_mle_w151672[1] + beta_mle_w151672[2]*median(INDCOST) + beta_mle_w151672[3]*value_sizelog
mu_mle_w151672 <- exp(eta_mle_w151672)/(1+exp(eta_mle_w151672))
#phi_mle_w151672 <- exp(gama_mle_w151672)

#*****SMLE*****#

#full data
beta_smle <- fitRob$beta; gama_smle <- fitRob$gama
eta_smle <- beta_smle[1] + beta_smle[2]*median(INDCOST) + beta_smle[3]*value_sizelog
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
#phi_smle <- exp(gama_smle)

#without 15
beta_smle_w15 <- fitRob_15$beta; gama_smle_w15 <- fitRob_15$gama
eta_smle_w15 <- beta_smle_w15[1] + beta_smle_w15[2]*median(INDCOST) + beta_smle_w15[3]*value_sizelog
mu_smle_w15 <- exp(eta_smle_w15)/(1+exp(eta_smle_w15))
#phi_smle_w15 <- exp(gama_smle_w15)

#without 16
beta_smle_w16 <- fitRob_16$beta; gama_smle_w16 <- fitRob_16$gama
eta_smle_w16 <- beta_smle_w16[1] +  beta_smle_w16[2]*median(INDCOST) + beta_smle_w16[3]*value_sizelog
mu_smle_w16 <- exp(eta_smle_w16)/(1+exp(eta_smle_w16))
#phi_smle_w16 <- exp(gama_smle_w16)

#without 72
beta_smle_w72 <- fitRob_72$beta; gama_smle_w72 <- fitRob_72$gama
eta_smle_w72 <- beta_smle_w72[1] + beta_smle_w72[2]*median(INDCOST) + beta_smle_w72[3]*value_sizelog
mu_smle_w72 <- exp(eta_smle_w72)/(1+exp(eta_smle_w72))
#phi_smle_w72 <- exp(gama_smle_w72)

#without 15 and 16
beta_smle_w1516 <- fitRob_15_16$beta; gama_smle_w1516 <- fitRob_15_16$gama
eta_smle_w1516 <- beta_smle_w1516[1] + beta_smle_w1516[2]*median(INDCOST) + beta_smle_w1516[3]*value_sizelog
mu_smle_w1516 <- exp(eta_smle_w1516)/(1+exp(eta_smle_w1516))
#phi_smle_w1516 <- exp(gama_smle_w1516)

#without 15 and 72
beta_smle_w1572 <- fitRob_15_72$beta; gama_smle_w1572 <- fitRob_15_72$gama
eta_smle_w1572 <- beta_smle_w1572[1] + beta_smle_w1572[2]*median(INDCOST) + beta_smle_w1572[3]*value_sizelog
mu_smle_w1572 <- exp(eta_smle_w1572)/(1+exp(eta_smle_w1572))
#phi_smle_w1572 <- exp(gama_smle_w1572)

#without 16 and 72
beta_smle_w1672 <- fitRob_16_72$beta; gama_smle_w1672 <- fitRob_16_72$gama
eta_smle_w1672 <- beta_smle_w1672[1] + beta_smle_w1672[2]*median(INDCOST) + beta_smle_w1672[3]*value_sizelog
mu_smle_w1672 <- exp(eta_smle_w1672)/(1+exp(eta_smle_w1672))
#phi_smle_w1672 <- exp(gama_smle_w1672)

#without 15, 16 and 72
beta_smle_w151672 <- fitRob_15_16_72$beta; gama_smle_w151672 <- fitRob_15_16_72$gama
eta_smle_w151672 <- beta_smle_w151672[1] + beta_smle_w151672[2]*median(INDCOST) + beta_smle_w151672[3]*value_sizelog
mu_smle_w151672 <- exp(eta_smle_w151672)/(1+exp(eta_smle_w151672))
#phi_smle_w151672 <- exp(gama_smle_w151672)


#******SMLE****#

par(mar=c(5.0,5.0,4.0,2.0))
plot(SIZELOG,y,pch=16,xlab="Sizelog",ylab="Firm cost",ylim=c(0,1), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(value_sizelog,mu_smle,         lwd=2, col=1, lty=1)
lines(value_sizelog,mu_smle_w15,     lwd=4, col=2, lty=2)
lines(value_sizelog,mu_smle_w16,     lwd=2, col=6, lty=6)
lines(value_sizelog,mu_smle_w72,     lwd=2, col=4, lty=4)
lines(value_sizelog,mu_smle_w1516,   lwd=2, col=5, lty=5)
lines(value_sizelog,mu_smle_w1572,   lwd=3, col=3, lty=1)
lines(value_sizelog,mu_smle_w1672,   lwd=2, col="orangered", lty=3)
lines(value_sizelog,mu_smle_w151672, lwd=3, col=8, lty=4)
identify(SIZELOG,y,cex=1.3)
legend(8, 0.8, c("SMLE","SMLE w/o 15","SMLE w/o 16","SMLE w/o 72", "SMLE w/o 15, 16",
"SMLE w/o 15, 72", "SMLE w/o 16, 72", "SMLE w/o 15, 16, 72"),
col=c(1,2,6,4,5,3,"orangered",8), lty=c(1,2,6,4,5,1,3,4),cex=1.0,lwd=c(2,4,2,2,2,3,2,3))

#*****MLE****#

par(mar=c(5.0,5.0,4.0,2.0))
plot(SIZELOG,y,pch=16,xlab="Sizelog",ylab="Firm cost",ylim=c(0,1), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(value_sizelog,mu_mle,         lwd=2, col=1, lty=1)
lines(value_sizelog,mu_mle_w15,     lwd=4, col=2, lty=2)
lines(value_sizelog,mu_mle_w16,     lwd=2, col=6, lty=6)
lines(value_sizelog,mu_mle_w72,     lwd=2, col=4, lty=4)
lines(value_sizelog,mu_mle_w1516,   lwd=2, col=5, lty=5)
lines(value_sizelog,mu_mle_w1572,   lwd=3, col=3, lty=1)
lines(value_sizelog,mu_mle_w1672,   lwd=2, col="orangered", lty=3)
lines(value_sizelog,mu_mle_w151672, lwd=3, col=8, lty=4)
identify(SIZELOG,y,cex=1.3)
legend(8, 0.8,c("MLE","MLE w/o 15","MLE w/o 16","MLE w/o 72", "MLE w/o 15, 16",
"MLE w/o 15, 72", "MLE w/o 16, 72", "MLE w/o 15, 16, 72"),
col=c(1,2,6,4,5,3,"orangered",8), lty=c(1,2,6,4,5,1,3,4),cex=1.0,lwd=c(2,4,2,2,2,3,2,3))


#*********************************BOOTSTRAP PVALUE*************************#
source("pvalue_boot.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=X,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y,X=X,Z=X,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15***#
pvalue_boot(y=y15,X=X15,Z=Z15,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y15,X=X15,Z=Z15,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16***#
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 72***#
pvalue_boot(y=y72,X=X72,Z=Z72,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y72,X=X72,Z=Z72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15 AND 16***#
pvalue_boot(y=y15_16,X=X15_16,Z=Z15_16,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y15_16,X=X15_16,Z=Z15_16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15 AND 72***#
pvalue_boot(y=y15_72,X=X15_72,Z=Z15_72,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y15_72,X=X15_72,Z=Z15_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16 AND 72***#
pvalue_boot(y=y16_72,X=X16_72,Z=Z16_72,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y16_72,X=X16_72,Z=Z16_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15, 16 AND 72***#
pvalue_boot(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE


#*************************Normal probability plots of residuals with simulated envelope**********************#
source("envelope_function.r")

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T,main.title = "SMLE full data",
 faixa.fixed = c(-5,8), labels.fixed = c(1:73))

#***DATA WITHOUT 15***#
envelope_SMLE(y=y15,X=X15,Z=Z15,theta=c(fitRob_15$beta,fitRob_15$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 15",
 faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:73))

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 16", 
faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:73))

#***DATA WITHOUT 72***#
envelope_SMLE(y=y72,X=X72,Z=Z72,theta=c(fitRob_72$beta,fitRob_72$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 72",
 faixa.fixed = c(-5,8), labels.fixed = c(1:71,73))

#***DATA WITHOUT 15 AND 16***#
envelope_SMLE(y=y15_16,X=X15_16,Z=Z15_16,theta=c(fitRob_15_16$beta,fitRob_15_16$gama), linkmu="logit", linkphi="log", SMLE=T, 
main.title = "SMLE without observations 15 \nand 16", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:73))

#***DATA WITHOUT 15 AND 72***#
envelope_SMLE(y=y15_72,X=X15_72,Z=Z15_72,theta=c(fitRob_15_72$beta,fitRob_15_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 15 \nand 72", faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:71,73))

#***DATA WITHOUT 16 AND 72***#
envelope_SMLE(y=y16_72,X=X16_72,Z=Z16_72,theta=c(fitRob_16_72$beta,fitRob_16_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 16 \nand 72", faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:71,73))

#***DATA WITHOUT 15, 16 AND 72***#
envelope_SMLE(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,theta=c(fitRob_15_16_72$beta,fitRob_15_16_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 15, \n16 and 72", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:71,73))


#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F, main.title = "MLE full data",
 faixa.fixed = c(-5,8), labels.fixed = c(1:73))

#***DATA WITHOUT 15***#
envelope_SMLE(y=y15,X=X15,Z=Z15,theta=c(fitMLE_15$beta,fitMLE_15$gama), linkmu="logit", linkphi="log",
 SMLE=F, main.title = "MLE without observation 15", faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:73))

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitMLE_16$beta,fitMLE_16$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observation 16", faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:73))

#***DATA WITHOUT 72***#
envelope_SMLE(y=y72,X=X72,Z=Z72,theta=c(fitMLE_72$beta,fitMLE_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observation 72", faixa.fixed = c(-5,8), labels.fixed = c(1:71,73))

#***DATA WITHOUT 15 AND 16***#
envelope_SMLE(y=y15_16,X=X15_16,Z=Z15_16,theta=c(fitMLE_15_16$beta,fitMLE_15_16$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 15 \nand 16", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:73))

#***DATA WITHOUT 15 AND 72***#
envelope_SMLE(y=y15_72,X=X15_72,Z=Z15_72,theta=c(fitMLE_15_72$beta,fitMLE_15_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 15 \nand 72",
 faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:71,73))

#***DATA WITHOUT 16 AND 72***#
envelope_SMLE(y=y16_72,X=X16_72,Z=Z16_72,theta=c(fitMLE_16_72$beta,fitMLE_16_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 16 \nand 72", faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:71,73))

#***DATA WITHOUT 15, 16 AND 72***#
envelope_SMLE(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,theta=c(fitMLE_15_16_72$beta,fitMLE_15_16_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 15, \n16 and 72", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:71,73))


#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")
#***FULL DATA***#
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
RP2_mle <- residuals_beta(y,X,Z,c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.2))
#identify(RP2_smle,weights,cex=1.3)
identify(RP2_smle[c(15,16,72)],weights[c(15,16,72)],cex=1.3,labels=c("15","16","72"))

#***DATA WITHOUT 15***#
weights_15 <- fitRob_15$weights
RP2_smle_w15 <- residuals_beta(y15,X15,Z15,c(fitRob_15$beta,fitRob_15$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15,weights_15,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 15"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15,weights_15,cex=1.3)

#***DATA WITHOUT 16***#
weights_16 <- fitRob_16$weights
RP2_smle_w16 <- residuals_beta(y16,X16,Z16,c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16,weights_16,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 16"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16,weights_16,cex=1.3)

#***DATA WITHOUT 16***#
weights_72 <- fitRob_72$weights
RP2_smle_w72 <- residuals_beta(y72,X72,Z72,c(fitRob_72$beta,fitRob_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w72,weights_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w72,weights_72,cex=1.3)

#***DATA WITHOUT 15 AND 16***#
weights_15_16 <- fitRob_15_16$weights
RP2_smle_w15_16 <- residuals_beta(y15_16,X15_16,Z15_16,c(fitRob_15_16$beta,fitRob_15_16$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_16,weights_15_16,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 15 \nand 16"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15_16,weights_15_16,cex=1.3)

#***DATA WITHOUT 15 AND 72***#
weights_15_72 <- fitRob_15_72$weights
RP2_smle_w15_72 <- residuals_beta(y15_72,X15_72,Z15_72,c(fitRob_15_72$beta,fitRob_15_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_72,weights_15_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 15 \nand 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15_72,weights_15_72,cex=1.3)

#***DATA WITHOUT 16 AND 72***#
weights_16_72 <- fitRob_16_72$weights
RP2_smle_w16_72 <- residuals_beta(y16_72,X16_72,Z16_72,c(fitRob_16_72$beta,fitRob_16_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16_72,weights_16_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 16 \nand 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16_72,weights_16_72,cex=1.3)

#***DATA WITHOUT 15, 16 AND 72***#
weights_15_16_72 <- fitRob_15_16_72$weights
RP2_smle_w15_16_72 <- residuals_beta(y15_16_72,X15_16_72,Z15_16_72,c(fitRob_15_16_72$beta,fitRob_15_16_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_16_72,weights_15_16_72,pch=16,xlab="Residuals", ylab="Weights",
main="Without observations 15, \n16 and 72", cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16_72,weights_16_72,cex=1.3)


#********ZOOMED VERSIONS*********#

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T,main.title = "SMLE full data",
 faixa.fixed = c(-4.5,4.5), labels.fixed = c(1:73))

#***DATA WITHOUT 15***#
envelope_SMLE(y=y15,X=X15,Z=Z15,theta=c(fitRob_15$beta,fitRob_15$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 15",
 faixa.fixed = c(-4,4), labels.fixed = c(1:14,16:73))

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 16", 
faixa.fixed = c(-4,4), labels.fixed = c(1:15,17:73))

#***DATA WITHOUT 72***#
envelope_SMLE(y=y72,X=X72,Z=Z72,theta=c(fitRob_72$beta,fitRob_72$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 72",
 faixa.fixed = c(-4,4), labels.fixed = c(1:71,73))

#***DATA WITHOUT 15 AND 16***#
envelope_SMLE(y=y15_16,X=X15_16,Z=Z15_16,theta=c(fitRob_15_16$beta,fitRob_15_16$gama), linkmu="logit", linkphi="log", SMLE=T, 
main.title = "SMLE without observations 15 \nand 16", faixa.fixed = c(-4,4), labels.fixed = c(1:14,17:73))

#***DATA WITHOUT 15 AND 72***#
envelope_SMLE(y=y15_72,X=X15_72,Z=Z15_72,theta=c(fitRob_15_72$beta,fitRob_15_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 15 \nand 72", faixa.fixed = c(-4,4), labels.fixed = c(1:14,16:71,73))

#***DATA WITHOUT 16 AND 72***#
envelope_SMLE(y=y16_72,X=X16_72,Z=Z16_72,theta=c(fitRob_16_72$beta,fitRob_16_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 16 \nand 72", faixa.fixed = c(-4,4), labels.fixed = c(1:15,17:71,73))

#***DATA WITHOUT 15, 16 AND 72***#
envelope_SMLE(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,theta=c(fitRob_15_16_72$beta,fitRob_15_16_72$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 15, \n16 and 72", faixa.fixed = c(-4,4), labels.fixed = c(1:14,17:71,73))


#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F, main.title = "MLE full data",
 faixa.fixed = c(-4.5,4.5), labels.fixed = c(1:73))

#***DATA WITHOUT 15***#
envelope_SMLE(y=y15,X=X15,Z=Z15,theta=c(fitMLE_15$beta,fitMLE_15$gama), linkmu="logit", linkphi="log",
 SMLE=F, main.title = "MLE without observation 15", faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:73))

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitMLE_16$beta,fitMLE_16$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observation 16", faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:73))

#***DATA WITHOUT 72***#
envelope_SMLE(y=y72,X=X72,Z=Z72,theta=c(fitMLE_72$beta,fitMLE_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observation 72", faixa.fixed = c(-5,8), labels.fixed = c(1:71,73))

#***DATA WITHOUT 15 AND 16***#
envelope_SMLE(y=y15_16,X=X15_16,Z=Z15_16,theta=c(fitMLE_15_16$beta,fitMLE_15_16$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 15 \nand 16", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:73))

#***DATA WITHOUT 15 AND 72***#
envelope_SMLE(y=y15_72,X=X15_72,Z=Z15_72,theta=c(fitMLE_15_72$beta,fitMLE_15_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 15 \nand 72",
 faixa.fixed = c(-5,8), labels.fixed = c(1:14,16:71,73))

#***DATA WITHOUT 16 AND 72***#
envelope_SMLE(y=y16_72,X=X16_72,Z=Z16_72,theta=c(fitMLE_16_72$beta,fitMLE_16_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 16 \nand 72", faixa.fixed = c(-5,8), labels.fixed = c(1:15,17:71,73))

#***DATA WITHOUT 15, 16 AND 72***#
envelope_SMLE(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,theta=c(fitMLE_15_16_72$beta,fitMLE_15_16_72$gama), linkmu="logit", linkphi="log", SMLE=F,
 main.title = "MLE without observations 15, \n16 and 72", faixa.fixed = c(-5,8), labels.fixed = c(1:14,17:71,73))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")
#***FULL DATA***#
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
RP2_mle <- residuals_beta(y,X,Z,c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4.5,4.5), ylim=c(0.8,1.2))
#identify(RP2_smle,weights,cex=1.3)
identify(RP2_smle[c(15,16,72)],weights[c(15,16,72)],cex=1.3,labels=c("15","16","72"))

#***DATA WITHOUT 15***#
weights_15 <- fitRob_15$weights
RP2_smle_w15 <- residuals_beta(y15,X15,Z15,c(fitRob_15$beta,fitRob_15$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15,weights_15,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 15"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15,weights_15,cex=1.3)

#***DATA WITHOUT 16***#
weights_16 <- fitRob_16$weights
RP2_smle_w16 <- residuals_beta(y16,X16,Z16,c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16,weights_16,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 16"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16,weights_16,cex=1.3)

#***DATA WITHOUT 16***#
weights_72 <- fitRob_72$weights
RP2_smle_w72 <- residuals_beta(y72,X72,Z72,c(fitRob_72$beta,fitRob_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w72,weights_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w72,weights_72,cex=1.3)

#***DATA WITHOUT 15 AND 16***#
weights_15_16 <- fitRob_15_16$weights
RP2_smle_w15_16 <- residuals_beta(y15_16,X15_16,Z15_16,c(fitRob_15_16$beta,fitRob_15_16$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_16,weights_15_16,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 15 \nand 16"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15_16,weights_15_16,cex=1.3)

#***DATA WITHOUT 15 AND 72***#
weights_15_72 <- fitRob_15_72$weights
RP2_smle_w15_72 <- residuals_beta(y15_72,X15_72,Z15_72,c(fitRob_15_72$beta,fitRob_15_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_72,weights_15_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 15 \nand 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w15_72,weights_15_72,cex=1.3)

#***DATA WITHOUT 16 AND 72***#
weights_16_72 <- fitRob_16_72$weights
RP2_smle_w16_72 <- residuals_beta(y16_72,X16_72,Z16_72,c(fitRob_16_72$beta,fitRob_16_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16_72,weights_16_72,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 16 \nand 72"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16_72,weights_16_72,cex=1.3)

#***DATA WITHOUT 15, 16 AND 72***#
weights_15_16_72 <- fitRob_15_16_72$weights
RP2_smle_w15_16_72 <- residuals_beta(y15_16_72,X15_16_72,Z15_16_72,c(fitRob_15_16_72$beta,fitRob_15_16_72$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w15_16_72,weights_15_16_72,pch=16,xlab="Residuals", ylab="Weights",
main="Without observations 15, \n16 and 72", cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-5,8), ylim=c(0,1.4))
identify(RP2_smle_w16_72,weights_16_72,cex=1.3)













