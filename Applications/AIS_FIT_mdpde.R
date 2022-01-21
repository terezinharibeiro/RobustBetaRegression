source("MDPDE.r") 

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
fitMLE<- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE

# MLE FIT WITHOUT 16 
fitMLE_16 <- MDPDE_BETA(y=y16, X=X16, Z=Z16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_16

# MLE FIT WITHOUT 30
fitMLE_30<- MDPDE_BETA(y=y30, X=X30, Z=Z30, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_30

# MLE FIT WITHOUT 16 AND 30
fitMLE_16_30<- MDPDE_BETA(y=y16_30, X=X16_30, Z=Z16_30, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log")
fitMLE_16_30

#*******************************************MDPDE FITS******************************************************#
# MDPDE FIT FULL DATA
fitRob <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit", linkphi="log", weights=T)
fitRob

#  MDPDE FIT WITHOUT 16
fitRob_16 <- MDPDE_BETA(y=y16, X=X16, Z=Z16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_16

#  MDPDE FIT WITHOUT 30
fitRob_30 <- MDPDE_BETA(y=y30, X=X30, Z=Z30, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_30

#  MDPDE FIT WITHOUT 16 AND 30
fitRob_16_30 <- MDPDE_BETA(y=y16_30, X=X16_30, Z=Z16_30, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
linkmu="logit",linkphi="log", weights=T)
fitRob_16_30


#********************BOOTSTRAP PVALUE*************************#
source("pvalue_boot_MDPDE.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16***#
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 30***#
pvalue_boot(y=y30,X=X30,Z=Z30,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y30,X=X30,Z=Z30,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16 AND 30***#
pvalue_boot(y=y16_30,X=X16_30,Z=Z16_30,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y16_30,X=X16_30,Z=Z16_30,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE
#*********************************************************************************#