source("SMLE.r") 

data_tuna <- read.table("data_2000_Indian.txt",h=T) #original data
attach(data_tuna)
tropicalP <-  trop/100 
index0<- 1:length(tropicalP)
index0[tropicalP==0] #no observation equal to zero
index1<- 1:length(tropicalP)
index1[tropicalP==1] #one observation equal to one: #46
tropicalP[c(46)] <- 0.999 #replacement value
Southern_Indian <- tropicalP[lat<0] #response variable
SST_Southern_Indian <- SST[lat<0] #covariate

#**************************************************#
#******Beta regression with constant precision*****#
#**************************************************#
y <- Southern_Indian 
X <- matrix(c(rep(1,77), SST_Southern_Indian), ncol=2, byrow=F);
Z <- matrix(c(rep(1,77)), ncol=1, byrow=F);
#*************************Subsets without outliers******************************#
y46 <- y[-c(46)]; X46<- X[-c(46),];Z46 <- as.matrix(Z[-c(46),])
y25 <- y[-c(25)]; X25<- X[-c(25),];Z25 <- as.matrix(Z[-c(25),])
y25_46 <- y[-c(25,46)]; X25_46<- X[-c(25,46),];Z25_46 <- as.matrix(Z[-c(25,46),])
#********************************************************************************#

#********************************************MLE FITS*****************************************************#
# MLE FIT FULL DATA
fitMLE<- SMLE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02,method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE

# MLE FIT WITHOUT 46 
fitMLE_46 <- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_46

# MLE FIT WITHOUT 25 
fitMLE_25 <- SMLE_BETA(y=y25, X=X25, Z=Z25, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25

# MLE FIT WITHOUT 25 and 46
fitMLE_25_46 <- SMLE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25_46

#********************************************SMLE FITS*****************************************************#
# SMLE FIT FULL DATA
fitRob <- SMLE_BETA(y=y, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob

# SMLE FIT WITHOUT 46
fitRob_46 <- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_46

# SMLE FIT WITHOUT 25
fitRob_25 <- SMLE_BETA(y=y25, X=X25, Z=Z25, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5,spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25

# SMLE FIT WITHOUT 25 and 46
fitRob_25_46 <- SMLE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25_46

#*********************SCATTER PLOTS***********************#
values <- seq(15,30,length.out=1000)
#*****SMLE fit full data*****#
beta_smle <- fitRob$beta; gama_smle <- fitRob$gama
eta_smle <- beta_smle[1] +  beta_smle[2]*values 
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
phi_smle <- exp(gama_smle)

#***MLE fit full data***#
beta_mle <- fitMLE$beta; gama_mle <- fitMLE$gama
eta_mle <- beta_mle[1] +  beta_mle[2]*values 
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
phi_mle <- exp(gama_mle)

#***SMLE fit without 46***#
beta_smle_w46 <- fitRob_46$beta; gama_smle_w46 <- fitRob_46$gama
eta_smle_w46 <- beta_smle_w46[1] +  beta_smle_w46[2]*values 
mu_smle_w46 <- exp(eta_smle_w46)/(1+exp(eta_smle_w46))
phi_smle_w46 <- exp(gama_smle_w46)

#***MLE fit without 46***#
beta_mle_w46 <- fitMLE_46$beta; gama_mle_w46 <- fitMLE_46$gama
eta_mle_w46 <- beta_mle_w46[1] +  beta_mle_w46[2]*values 
mu_mle_w46 <- exp(eta_mle_w46)/(1+exp(eta_mle_w46))
phi_mle_w46 <- exp(gama_mle_w46)

#***SCATTER PLOT****#
par(mar=c(5.0,5.0,4.0,2.0))
yoriginal <- y; yoriginal[46] <- 1
plot(SST_Southern_Indian, yoriginal, pch=16, xlab="SST", ylab="TTP", ylim=c(0,1), xlim=c(15,30), cex=1.5, cex.lab=2.0, cex.axis=1.5, 
cex.main=2.0)
lines(values, mu_mle,       lwd=2,col=2,lty=1)
lines(values,mu_mle_w46,    lwd=2,col=6,lty=4)
lines(values,mu_smle,       lwd=2,col=4,lty=3)
identify(SST_Southern_Indian, yoriginal, cex=1.3)
legend(16,0.8,c("MLE","MLE, SMLE w/o 46", "SMLE"),
col=c(2,6,4),lty=c(1,4,3), cex=1.5, lwd=c(2,2,2))

#************************************BOOTSTRAP PVALUE*************************#
source("pvalue_boot.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y, X=X, Z=Z, B=B, optimal=T, linkmu="logit", linkphi="log") #SMLE
pvalue_boot(y=y, X=X, Z=Z, B=B, optimal=F, linkmu="logit", linkphi="log") #MLE

#***DATA WITHOUT 46***#
pvalue_boot(y=y46, X=X46, Z=Z46, B=B, optimal=T, linkmu="logit", linkphi="log") #SMLE
pvalue_boot(y=y46, X=X46, Z=Z46, B=B, optimal=F, linkmu="logit", linkphi="log") #MLE
#*********************************************************************************#

#*************************Normal probability plots of residuals with simulated envelope**********************#
source("envelope_function.r")

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y, X=X, Z=Z, theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE full data",
faixa.fixed = c(-6,18), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log", SMLE=T,
main.title = "SMLE without observation 46", faixa.fixed = c(-6,18), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 25", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 25 \nand 46", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:45,47:77))

#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y, X=X, Z=Z, theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F,
main.title = "MLE full data", faixa.fixed = c(-6,18), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitMLE_46$beta,fitMLE_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 46", faixa.fixed = c(-6,18), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitMLE_25$beta,fitMLE_25$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 25", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitMLE_25_46$beta,fitMLE_25_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 25 \nand 46", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:45,47:77))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")
#***FULL DATA***
weights <- fitRob$weights
RP2_smle <- residuals_beta(y, X, Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle,weights,cex=1.3)

#***DATA WITHOUT 46***#
weights_46 <- fitRob_46$weights
RP2_smle_w46 <- residuals_beta(y46,X46,Z46,c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w46,weights_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w46,weights_46,c(1:45,47:77),cex=1.3)

#***DATA WITHOUT 25***#
weights_25 <- fitRob_25$weights
RP2_smle_w25 <- residuals_beta(y25,X25,Z25,c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25,weights_25,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 25"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w25,weights_25,c(1:24,26:77),cex=1.3)

#***DATA WITHOUT 25 AND 46***#
weights_25_46 <- fitRob_25_46$weights
RP2_smle_w25_46 <- residuals_beta(y25_46,X25_46,Z25_46,c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25_46,weights_25_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 25 and 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w25_46,weights_25_46,c(1:24,26:45,47:77),cex=1.3)


#***ZOOMED VERSIONS***#

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T,main.title = "SMLE full data",
 faixa.fixed = c(-4,4), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 46", faixa.fixed = c(-4,4), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 25", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 25 \nand 46", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:45,47:77))

#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y, X=X, Z=Z, theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F,
main.title = "MLE full data", faixa.fixed = c(-4,4), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitMLE_46$beta,fitMLE_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 46", faixa.fixed = c(-4,4), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitMLE_25$beta,fitMLE_25$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 25", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitMLE_25_46$beta,fitMLE_25_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 25 \nand 46", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:45,47:77))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")
#***FULL DATA***
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0.8,1.2))
identify(RP2_smle,weights,cex=1.3)

#***DATA WITHOUT 46***#
weights_46 <- fitRob_46$weights
RP2_smle_w46 <- residuals_beta(y46,X46,Z46,c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w46,weights_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w46,weights_46,c(1:45,47:77),cex=1.3)

#***DATA WITHOUT 25***#
weights_25 <- fitRob_25$weights
RP2_smle_w25 <- residuals_beta(y25,X25,Z25,c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25,weights_25,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 25"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w25,weights_25,c(1:24,26:77),cex=1.3)

#***DATA WITHOUT 25 AND 46***#
weights_25_46 <- fitRob_25_46$weights
RP2_smle_w25_46 <- residuals_beta(y25_46,X25_46,Z25_46,c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25_46,weights_25_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 25 and 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w25_46,weights_25_46,c(1:24,26:45,47:77),cex=1.3)


#**************************************************#
#******Beta regression with varying precision******#
#**************************************************#
y <- Southern_Indian 
X <- matrix(c(rep(1,77), SST_Southern_Indian), ncol=2, byrow=F);
Z <- X
#*************************Subsets without outliers******************************#
y46 <- y[-c(46)]; X46<- X[-c(46),];Z46 <- as.matrix(Z[-c(46),])
y25 <- y[-c(25)]; X25<- X[-c(25),];Z25 <- as.matrix(Z[-c(25),])
y25_46 <- y[-c(25,46)]; X25_46<- X[-c(25,46),];Z25_46 <- as.matrix(Z[-c(25,46),])
#********************************************************************************#

#********************************************MLE FITS*****************************************************#
# MLE FIT FULL DATA
fitMLE<- SMLE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE

# MLE FIT WITHOUT 46 
fitMLE_46 <- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS",
startV="CP", linkmu="logit",linkphi="log")
fitMLE_46

# MLE FIT WITHOUT 25 
fitMLE_25 <- SMLE_BETA(y=y25, X=X25, Z=Z25, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25

# MLE FIT WITHOUT 25 and 46
fitMLE_25_46 <- SMLE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5,spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25_46

#********************************************SMLE FITS*****************************************************#
# SMLE FIT FULL DATA
fitRob <- SMLE_BETA(y=y, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob

# SMLE FIT WITHOUT 46
fitRob_46 <- SMLE_BETA(y=y46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_46

# SMLE FIT WITHOUT 25
fitRob_25 <- SMLE_BETA(y=y25, X=X25, Z=Z25, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25

# SMLE FIT WITHOUT 25 and 46
fitRob_25_46 <- SMLE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25_46

#*********************SCATTER PLOTS***********************#
values <- seq(15,30,length.out=1000)
#*****SMLE fit full data*****#
beta_smle <- fitRob$beta; gama_smle <- fitRob$gama
eta_smle <- beta_smle[1] +  beta_smle[2]*values 
mu_smle <- exp(eta_smle)/(1+exp(eta_smle))
phi_smle <- exp(gama_smle)

#*****MLE fit full data*****#
beta_mle <- fitMLE$beta; gama_mle <- fitMLE$gama
eta_mle <- beta_mle[1] +  beta_mle[2]*values 
mu_mle <- exp(eta_mle)/(1+exp(eta_mle))
phi_mle <- exp(gama_mle)

#***SMLE fit without 46***#
beta_smle_w46 <- fitRob_46$beta; gama_smle_w46 <- fitRob_46$gama
eta_smle_w46 <- beta_smle_w46[1] +  beta_smle_w46[2]*values 
mu_smle_w46 <- exp(eta_smle_w46)/(1+exp(eta_smle_w46))
phi_smle_w46 <- exp(gama_smle_w46)

#***MLE fit without 46***#
beta_mle_w46 <- fitMLE_46$beta; gama_mle_w46 <- fitMLE_46$gama
eta_mle_w46 <- beta_mle_w46[1] +  beta_mle_w46[2]*values 
mu_mle_w46 <- exp(eta_mle_w46)/(1+exp(eta_mle_w46))
phi_mle_w46 <- exp(gama_mle_w46)

#***SCATTER PLOT****#
par(mar=c(5.0,5.0,4.0,2.0))
plot(SST_Southern_Indian,y,pch=16,xlab="SST",ylab="TTP",ylim=c(0,1),xlim=c(15,30), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(values, mu_mle,       lwd=2,col=2,lty=1)
lines(values,mu_mle_w46,    lwd=2,col=6,lty=6)
lines(values,mu_smle,       lwd=2,col=4,lty=4)
identify(SST_Southern_Indian,y, cex=1.3)
legend(16,0.8,c("MLE","MLE, SMLE w/o 46", "SMLE"),
col=c(2,6,4),lty=c(1,6,4),cex=1.5,lwd=c(2,2,2))

#************************************BOOTSTRAP PVALUE*************************#
source("pvalue_boot.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 46***#
pvalue_boot(y=y46,X=X46,Z=Z46,B=B, optimal=T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y46,X=X46,Z=Z46,B=B, optimal=F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 25***#
pvalue_boot(y=y25,X=X25,Z=Z25,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y25,X=X25,Z=Z25,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 25 AND 46***#
pvalue_boot(y=y25_46,X=X25_46,Z=Z25_46,B=B, optimal =T, linkmu="logit",linkphi="log") #SMLE
pvalue_boot(y=y25_46,X=X25_46,Z=Z25_46,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE
#*********************************************************************************#

#*************************Normal probability plots of residuals with simulated envelope**********************#
source("envelope_function.r")

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T, main.title = "SMLE full data",
 faixa.fixed = c(-6,18), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 46", faixa.fixed = c(-6,18), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 25", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 25 \nand 46", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:45,47:77))

#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F,
main.title = "MLE full data", faixa.fixed = c(-6,18), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitMLE_46$beta,fitMLE_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 46", faixa.fixed = c(-6,18), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitMLE_25$beta,fitMLE_25$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 25", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitMLE_25_46$beta,fitMLE_25_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 25 \nand 46", faixa.fixed = c(-6,18), labels.fixed = c(1:24,26:45,47:77))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")
#***FULL DATA***#
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle,weights,cex=1.3)

#***DATA WITHOUT 46***#
weights_46 <- fitRob_46$weights
RP2_smle_w46 <- residuals_beta(y46,X46,Z46,c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w46,weights_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w46,weights_46,c(1:45,47:77),cex=1.3)

#***DATA WITHOUT 25***#
weights_25 <- fitRob_25$weights
RP2_smle_w25 <- residuals_beta(y25,X25,Z25,c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25,weights_25,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 25"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w25,weights_25,c(1:24,26:77),cex=1.3)

#***DATA WITHOUT 25 AND 46***#
weights_25_46 <- fitRob_25_46$weights
RP2_smle_w25_46 <- residuals_beta(y25_46,X25_46,Z25_46,c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25_46,weights_25_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 25 \nand 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-6,18), ylim=c(0,1.3))
identify(RP2_smle_w25_46,weights_25_46,c(1:24,26:45,47:77),cex=1.3)

#***ZOOMED VERSIONS***#

#************SMLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log", SMLE=T,main.title = "SMLE full data",
 faixa.fixed = c(-4,4), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 46", faixa.fixed = c(-4,4), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observation 25", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log", SMLE=T,
 main.title = "SMLE without observations 25 \nand 46", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:45,47:77))

#************MLE***********#

#***FULL DATA***#
envelope_SMLE(y=y,X=X,Z=Z,theta=c(fitMLE$beta,fitMLE$gama), linkmu="logit", linkphi="log", SMLE=F,
main.title = "MLE full data", faixa.fixed = c(-4,4), labels.fixed =1:77)

#***DATA WITHOUT 46***#
envelope_SMLE(y=y46,X=X46,Z=Z46,theta=c(fitMLE_46$beta,fitMLE_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 46", faixa.fixed = c(-4,4), labels.fixed = c(1:45,47:77))

#***DATA WITHOUT 25***#
envelope_SMLE(y=y25,X=X25,Z=Z25,theta=c(fitMLE_25$beta,fitMLE_25$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observation 25", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:77))

#***DATA WITHOUT 25 AND 46***#
envelope_SMLE(y=y25_46,X=X25_46,Z=Z25_46,theta=c(fitMLE_25_46$beta,fitMLE_25_46$gama), linkmu="logit", linkphi="log", SMLE=F, 
main.title = "MLE without observations 25 \nand 46", faixa.fixed = c(-4,4), labels.fixed = c(1:24,26:45,47:77))

#*****************RESIDUALS AGAINST WEIGHTS******************#
source("Resfunction.r")

#***FULL DATA***#
weights <- fitRob$weights
RP2_smle <- residuals_beta(y,X,Z,c(fitRob$beta,fitRob$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle,weights,pch=16,xlab="Residuals", ylab="Weights",main="Full data"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0.8,1.2))
identify(RP2_smle,weights,cex=1.3)

#***DATA WITHOUT 46***#
weights_46 <- fitRob_46$weights
RP2_smle_w46 <- residuals_beta(y46,X46,Z46,c(fitRob_46$beta,fitRob_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w46,weights_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w46,weights_46,c(1:45,47:77),cex=1.3)

#***DATA WITHOUT 25***#
weights_25 <- fitRob_25$weights
RP2_smle_w25 <- residuals_beta(y25,X25,Z25,c(fitRob_25$beta,fitRob_25$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25,weights_25,pch=16,xlab="Residuals", ylab="Weights",main="Without observation 25"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w25,weights_25,c(1:24,26:77),cex=1.3)

#***DATA WITHOUT 25 AND 46***#
weights_25_46 <- fitRob_25_46$weights
RP2_smle_w25_46 <- residuals_beta(y25_46,X25_46,Z25_46,c(fitRob_25_46$beta,fitRob_25_46$gama), linkmu="logit", linkphi="log") 
par(mar=c(5.0,5.0,4.0,2.0))
plot(RP2_smle_w25_46,weights_25_46,pch=16,xlab="Residuals", ylab="Weights",main="Without observations 25 \nand 46"
, cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, xlim=c(-4,4), ylim=c(0,1.3))
identify(RP2_smle_w25_46,weights_25_46,c(1:24,26:45,47:77),cex=1.3)


#**************************SCATTER PLOTS FOR DIFFERENT REPLACEMENT*******************#

#*****************Varying Precision***********#

y1 <- y; y2 <- y; y3 <- y
y1[c(46)] <- 0.999 #observation equal to one
y2[c(46)] <- 0.99 #observation equal to one
y3[c(46)] <- 0.95 #observation equal to one

X <- matrix(c(rep(1,77), SST_Southern_Indian),ncol=2,byrow=F);
Z <-X;
X46<- X[-c(46),];Z46 <- as.matrix(Z[-c(46),])
y1_46 <- y1[-c(46)]; y2_46 <- y2[-c(46)]; y3_46 <- y3[-c(46)];

#***************MLE FIT**********************#

# y[46] = 0.999
fitMLE1<- SMLE_BETA(y=y1, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5,spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE1

# y[46] = 0.99
fitMLE2<- SMLE_BETA(y=y2, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",linkmu="logit",linkphi="log")
fitMLE2
pvalue_boot(y=y2, X=X, Z=Z, B=B, optimal =F, linkmu="logit", linkphi="log") 

# y[46] = 0.95
fitMLE3<- SMLE_BETA(y=y3, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE3
pvalue_boot(y=y3, X=X, Z=Z, B=B, optimal=F, linkmu="logit",linkphi="log") 

# MLE FIT WITHOUT 46 

# y[46] = 0.999
fitMLE1_46 <- SMLE_BETA(y=y1_46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE1_46

# y[46] = 0.99
fitMLE2_46 <- SMLE_BETA(y=y2_46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE2_46

# y[46] = 0.95
fitMLE3_46 <- SMLE_BETA(y=y3_46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE3_46

#************************SMLE FIT*********************#

# y[46] = 0.999
fitRob1 <- SMLE_BETA(y=y1, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob1

# y[46] = 0.99
fitRob2 <- SMLE_BETA(y=y2, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob2

# y[46] = 0.95
fitRob3 <- SMLE_BETA(y=y3, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob3

# SMLE FIT WITHOUT 46

# y[46] = 0.999
fitRob1_46 <- SMLE_BETA(y=y1_46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob1_46

# y[46] = 0.99
fitRob2_46 <- SMLE_BETA(y=y2_46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob2_46

# y[46] = 0.95
fitRob3_46 <- SMLE_BETA(y=y3_46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob3_46


#********SCATTER PLOT*********#
values <- seq(15,30,length.out=1000)

#***FULL DATA****#

#**y[46] = 0.999**#

beta_smle1 <- fitRob1$beta; gama_smle1 <- fitRob1$gama
eta_smle1 <- beta_smle1[1] +  beta_smle1[2]*values 
mu_smle1 <- exp(eta_smle1)/(1+exp(eta_smle1))
phi_smle1 <- exp(gama_smle1)

beta_mle1 <- fitMLE1$beta; gama_mle1 <- fitMLE1$gama
eta_mle1 <- beta_mle1[1] +  beta_mle1[2]*values 
mu_mle1 <- exp(eta_mle1)/(1+exp(eta_mle1))
phi_mle1 <- exp(gama_mle1)

#**y[46] = 0.99**#

beta_smle2 <- fitRob2$beta; gama_smle2 <- fitRob2$gama
eta_smle2 <- beta_smle2[1] +  beta_smle2[2]*values 
mu_smle2 <- exp(eta_smle2)/(1+exp(eta_smle2))
phi_smle2 <- exp(gama_smle2)

beta_mle2 <- fitMLE2$beta; gama_mle2 <- fitMLE2$gama
eta_mle2 <- beta_mle2[1] +  beta_mle2[2]*values 
mu_mle2 <- exp(eta_mle2)/(1+exp(eta_mle2))
phi_mle2 <- exp(gama_mle2)

#**y[46] = 0.95**#

beta_smle3 <- fitRob3$beta; gama_smle3 <- fitRob3$gama
eta_smle3 <- beta_smle3[1] +  beta_smle3[2]*values 
mu_smle3 <- exp(eta_smle3)/(1+exp(eta_smle3))
phi_smle3 <- exp(gama_smle3)

beta_mle3 <- fitMLE3$beta; gama_mle3 <- fitMLE3$gama
eta_mle3 <- beta_mle3[1] +  beta_mle3[2]*values 
mu_mle3 <- exp(eta_mle3)/(1+exp(eta_mle3))
phi_mle3 <- exp(gama_mle3)

#without 46

#**y[46] = 0.999**#

beta_smle1_w46 <- fitRob1_46$beta; gama_smle1_w46 <- fitRob1_46$gama
eta_smle1_w46 <- beta_smle1_w46[1] +  beta_smle1_w46[2]*values 
mu_smle1_w46 <- exp(eta_smle1_w46)/(1+exp(eta_smle1_w46))
phi_smle1_w46 <- exp(gama_smle1_w46)

beta_mle1_w46 <- fitMLE1_46$beta; gama_mle1_w46 <- fitMLE1_46$gama
eta_mle1_w46 <- beta_mle1_w46[1] +  beta_mle1_w46[2]*values 
mu_mle1_w46 <- exp(eta_mle1_w46)/(1+exp(eta_mle1_w46))
phi_mle1_w46 <- exp(gama_mle1_w46)


#**y[46] = 0.99**#

beta_smle2_w46 <- fitRob2_46$beta; gama_smle2_w46 <- fitRob2_46$gama
eta_smle2_w46 <- beta_smle2_w46[1] +  beta_smle2_w46[2]*values 
mu_smle2_w46 <- exp(eta_smle2_w46)/(1+exp(eta_smle2_w46))
phi_smle2_w46 <- exp(gama_smle2_w46)

beta_mle2_w46 <- fitMLE2_46$beta; gama_mle2_w46 <- fitMLE2_46$gama
eta_mle2_w46 <- beta_mle2_w46[1] +  beta_mle2_w46[2]*values 
mu_mle2_w46 <- exp(eta_mle2_w46)/(1+exp(eta_mle2_w46))
phi_mle2_w46 <- exp(gama_mle2_w46)

#**y[46] = 0.95**#

beta_smle3_w46 <- fitRob3_46$beta; gama_smle3_w46 <- fitRob3_46$gama
eta_smle3_w46 <- beta_smle3_w46[1] +  beta_smle3_w46[2]*values 
mu_smle3_w46 <- exp(eta_smle3_w46)/(1+exp(eta_smle3_w46))
phi_smle3_w46 <- exp(gama_smle3_w46)

beta_mle3_w46 <- fitMLE3_46$beta; gama_mle3_w46 <- fitMLE3_46$gama
eta_mle3_w46 <- beta_mle3_w46[1] +  beta_mle3_w46[2]*values 
mu_mle3_w46 <- exp(eta_mle3_w46)/(1+exp(eta_mle3_w46))
phi_mle3_w46 <- exp(gama_mle3_w46)

#*** MLE***#
par(mar=c(5.0,5.0,4.0,2.0))
yoriginal <- y; yoriginal[46] <- 1
plot(SST_Southern_Indian,yoriginal,pch=16,xlab="SST",ylab="TTP",ylim=c(0,1),xlim=c(15,30), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(values, mu_mle1,       lwd=2,col=2,lty=1)
lines(values, mu_mle2,       lwd=2.5,col=3,lty=2)
lines(values, mu_mle3,       lwd=2.5,col=5,lty=5)
identify(SST_Southern_Indian,yoriginal, cex=1.3, n=1)
legend(16,0.8,c("MLE - 0.999","MLE - 0.99", "MLE - 0.95"),
col=c(2,3,5),lty=c(1,2,5),cex=1.5,lwd=c(2,2.5,2.5))

#*** SMLE***#
par(mar=c(5.0,5.0,4.0,2.0))
yoriginal <- y; yoriginal[46] <- 1
plot(SST_Southern_Indian,yoriginal,pch=16,xlab="SST",ylab="TTP",ylim=c(0,1),xlim=c(15,30), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0)
lines(values, mu_smle1,       lwd=3,col=4,lty=4)
lines(values, mu_smle2,       lwd=2,col="orangered",lty=7)
lines(values, mu_smle3,       lwd=2,col=1,lty=8)
identify(SST_Southern_Indian,yoriginal, cex=1.3, n=1)
legend(16,0.8,c("SMLE - 0.999","SMLE - 0.99", "SMLE - 0.95"),
col=c(4,"orangered",1),lty=c(4,7,8),cex=1.5,lwd=c(2,2,2))


















