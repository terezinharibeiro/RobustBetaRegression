source("SMLE.r") 

#*************packages and data required************#
require(sn); data(ais); attach(ais)
aisrow <- ais[ which(sport=='Row'), ]
attach(aisrow)
BFP <- Bfat/100
#***************************************************#

#******Beta regression with constant precision******#
y <- BFP #response variable
X <- matrix(c(rep(1,37),LBM),ncol=2,byrow=F); #regressor matrix for the mean submodel
Z <- matrix(c(rep(1,37)),ncol=1,byrow=F); #regressor matrix for the precision submodel
#*********************Subsets without outliers**********************************#
y16 <- y[-c(16)]; X16<- X[-c(16),];Z16 <- as.matrix(Z[-c(16),])
y30<- y[-c(30)]; X30<- X[-c(30),];Z30 <- as.matrix(Z[-c(30),])
y16_30 <- y[-c(16,30)]; X16_30<- X[-c(16,30),];Z16_30 <- as.matrix(Z[-c(16,30),])
#********************************************************************************#

#********************************************MLE FITS*********************************************#
# MLE FIT FULL DATA
fitMLE<- SMLE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE

# MLE FIT WITHOUT 16 
fitMLE_16 <- SMLE_BETA(y=y16, X=X16, Z=Z16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_16

# MLE FIT WITHOUT 30
fitMLE_30<- SMLE_BETA(y=y30, X=X30, Z=Z30, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_30

# MLE FIT WITHOUT 16 AND 30
fitMLE_16_30<- SMLE_BETA(y=y16_30, X=X16_30, Z=Z16_30, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_16_30

#*******************************************SMLE FITS******************************************************#
# SMLE FIT FULL DATA
fitRob <- SMLE_BETA(y=y, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit", linkphi="log", weights=T)
fitRob

# SMLE FIT WITHOUT 16
fitRob_16 <- SMLE_BETA(y=y16, X=X16, Z=Z16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_16

# SMLE FIT WITHOUT 30
fitRob_30 <- SMLE_BETA(y=y30, X=X30, Z=Z30, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_30

# SMLE FIT WITHOUT 16 AND 30
fitRob_16_30 <- SMLE_BETA(y=y16_30, X=X16_30, Z=Z16_30, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_16_30


#*********************SCATTER PLOTS***********************#
values <- seq(30,130,length.out=1000)

#*****MLE fit full data*****#
beta_mle <- fitMLE$beta; gama_mle <- fitMLE$gama
eta_mle <- beta_mle[1] +  beta_mle[2]*values 
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
phi_mle <- exp(gama_mle)

#***MLE fit without 16***#
beta_mle_w16 <- fitMLE_16$beta; gama_mle_w16 <- fitMLE_16$gama
eta_mle_w16 <- beta_mle_w16[1] +  beta_mle_w16[2]*values 
mu_mle_w16 <- exp(eta_mle_w16)/(1+exp(eta_mle_w16))
phi_mle_w16 <- exp(gama_mle_w16)

#***MLE fit without 30***#
beta_mle_w30 <- fitMLE_30$beta; gama_mle_w30 <- fitMLE_30$gama
eta_mle_w30 <- beta_mle_w30[1] +  beta_mle_w30[2]*values 
mu_mle_w30 <- exp(eta_mle_w30)/(1+exp(eta_mle_w30))
phi_mle_w30 <- exp(gama_mle_w30)

#***MLE fit without 16 and 30***#
beta_mle_w16_30 <- fitMLE_16_30$beta; gama_mle_w16_30 <- fitMLE_16_30$gama
eta_mle_w16_30 <- beta_mle_w16_30[1] +  beta_mle_w16_30[2]*values 
mu_mle_w16_30 <- exp(eta_mle_w16_30)/(1+exp(eta_mle_w16_30))
phi_mle_w16_30 <- exp(gama_mle_w16_30)

#*****SMLE fit full data*****#
beta_smle <- fitRob$beta; gama_smle <- fitRob$gama
eta_smle <- beta_smle[1] +  beta_smle[2]*values 
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
phi_smle <- exp(gama_smle)

#***SMLE fit without 16***#
beta_smle_w16 <- fitRob_16$beta; gama_smle_w16 <- fitRob_16$gama
eta_smle_w16 <- beta_smle_w16[1] +  beta_smle_w16[2]*values 
mu_smle_w16 <- exp(eta_smle_w16)/(1+exp(eta_smle_w16))
phi_smle_w16 <- exp(gama_smle_w16)

#***SMLE fit without 30***#
beta_smle_w30 <- fitRob_30$beta; gama_smle_w30 <- fitRob_30$gama
eta_smle_w30 <- beta_smle_w30[1] +  beta_smle_w30[2]*values 
mu_smle_w30 <- exp(eta_smle_w30)/(1+exp(eta_smle_w30))
phi_smle_w30 <- exp(gama_smle_w30)

#***SMLE fit without 16 and 30***#
beta_smle_w16_30 <- fitRob_16_30$beta; gama_smle_w16_30 <- fitRob_16_30$gama
eta_smle_w16_30 <- beta_smle_w16_30[1] +  beta_smle_w16_30[2]*values 
mu_smle_w16_30 <- exp(eta_smle_w16_30)/(1+exp(eta_smle_w16_30))
phi_smle_w16_30 <- exp(gama_smle_w16_30)

#***SCATTER PLOT FOR MLE***#
par(mar=c(5.0,5.0,4.0,2.0))
plot(LBM,y,pch=16,xlab="LBM",ylab="BFP",ylim=c(0,0.5),xlim=c(30,90), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(values, mu_mle,       lwd=2,col=2,lty=2)
lines(values,mu_mle_w16,    lwd=2.5,col=3,lty=1)
lines(values,mu_mle_w30,    lwd=4,col=6,lty=3)
lines(values,mu_mle_w16_30, lwd=2,col=5,lty=4)
identify(LBM,y, cex=1.3, n=2)
legend(60,0.5,c("MLE","MLE w/o 16","MLE w/o 30","MLE w/o 16, 30"),
col=c(2,3,6,5),lty=c(2,1,3,4),cex=1.0,lwd=c(2,2.5,4,2))

#***SCATTER PLOT FOR SMLE****#
par(mar=c(5.0,5.0,4.0,2.0))
plot(LBM,y,pch=16,xlab="LBM",ylab="BFP",ylim=c(0,0.5),xlim=c(30,90), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(values, mu_smle,       lwd=5,col=2,lty=2)
lines(values,mu_smle_w16,    lwd=2.5,col=3,lty=1)
lines(values,mu_smle_w30,    lwd=4,col=6,lty=3)
lines(values,mu_smle_w16_30, lwd=2,col=5,lty=4)
identify(LBM,y, cex=1.3, n=2)
legend(60,0.5,c("SMLE","SMLE w/o 16","SMLE w/o 30","SMLE w/o 16, 30"),
col=c(2,3,6,5),lty=c(2,1,3,4),cex=1.0,lwd=c(2,2.5,4,2))


#********************BOOTSTRAP PVALUE*************************#
source("pvalue_boot.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16***#
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 30***#
pvalue_boot(y=y30,X=X30,Z=Z30,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y30,X=X30,Z=Z30,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16 AND 30***#
pvalue_boot(y=y16_30,X=X16_30,Z=Z16_30,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y16_30,X=X16_30,Z=Z16_30,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE
#*********************************************************************************#

#*************************Normal probability plots of residuals with simulated envelope**********************#
source("envelope_function.r")

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z, theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE full data",
faixa.fixed = c(-8,4), labels.fixed = 1:37)

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log", SMLE=T,
main.title = "SMLE without observation 16", faixa.fixed = c(-8,4), labels.fixed = c(1:15,17:37))

#***DATA WITHOUT 30***#
envelope_SMLE(y=y30,X=X30,Z=Z30,theta=c(fitRob_30$beta,fitRob_30$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE without observation 30",
faixa.fixed = c(-8,4), labels.fixed = c(1:29,31:37))

#***DATA WITHOUT 16 AND 30***#
envelope_SMLE(y=y16_30,X=X16_30,Z=Z16_30,theta=c(fitRob_16_30$beta,fitRob_16_30$gama), linkmu="logit", linkphi="log", SMLE=T,
main.title = "SMLE without observations 16 \nand 30", faixa.fixed = c(-8,4), labels.fixed = c(1:15,17:29, 31:37))

#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", 
SMLE=F,main.title = "MLE full data", faixa.fixed = c(-8,4),labels.fixed = 1:37)

#***DATA WITHOUT 16***#
envelope_SMLE(y=y16,X=X16,Z=Z16,theta=c(fitMLE_16$beta,fitMLE_16$gama), linkmu="logit", linkphi="log", SMLE=F, main.title = "MLE without observation 16",
faixa.fixed = c(-8,4), labels.fixed = c(1:15,17:37))

#***DATA WITHOUT 30***#
envelope_SMLE(y=y30,X=X30,Z=Z30,theta=c(fitMLE_30$beta,fitMLE_30$gama), linkmu="logit", linkphi="log", SMLE=F, main.title = "MLE without observation 30",
faixa.fixed = c(-8,4), labels.fixed = c(1:29,31:37))

#***DATA WITHOUT 16 AND 30***#
envelope_SMLE(y=y16_30,X=X16_30,Z=Z16_30,theta=c(fitMLE_16_30$beta,fitMLE_16_30$gama), linkmu="logit", linkphi="log", SMLE=F,
main.title = "MLE without observations 16 \nand 30", faixa.fixed = c(-8,4), labels.fixed = c(1:15,17:29, 31:37))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")

#***FULL DATA***#
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-8,4), ylim=c(0,1.9))
identify(RP2_smle,weights,cex=1.3, n=2)

#***DATA WITHOUT 16***#
weights_16 <- fitRob_16$weights
RP2_smle_w16 <- residuals_beta(y16,X16,Z16,c(fitRob_16$beta,fitRob_16$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16,weights_16,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 16"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-8,4), ylim=c(0,1.9))
identify(RP2_smle_w16,weights_16,c(1:15,17:37),cex=1.3)

#***DATA WITHOUT 30***#
weights_30 <- fitRob_30$weights
RP2_smle_w30 <- residuals_beta(y30,X30,Z30,c(fitRob_30$beta,fitRob_30$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w30,weights_30,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 30"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-8,4), ylim=c(0,1.9))
identify(RP2_smle_w30,weights_30,cex=1.3)

#***DATA WITHOUT 16 AND 30***#
weights_16_30 <- fitRob_16_30$weights
RP2_smle_w16_30 <- residuals_beta(y16_30,X16_30,Z16_30,c(fitRob_16_30$beta,fitRob_16_30$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w16_30,weights_16_30,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 16 and 30"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-8,4), ylim=c(0,1.9))
identify(RP2_smle_w16_30,weights_16_30,cex=1.3)
