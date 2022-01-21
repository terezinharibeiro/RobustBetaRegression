source("MDPDE.r") 

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
fitMLE<- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE

# MLE FIT WITHOUT 46 
fitMLE_46 <- MDPDE_BETA(y=y46, X=X46, Z=Z46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS",
startV="CP", linkmu="logit",linkphi="log")
fitMLE_46

# MLE FIT WITHOUT 25 
fitMLE_25 <- MDPDE_BETA(y=y25, X=X25, Z=Z25, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25

# MLE FIT WITHOUT 25 and 46
fitMLE_25_46 <- MDPDE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5,spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_25_46

#********************************************MDPDE FITS*****************************************************#
# MDPDE FIT FULL DATA
fitRob <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob

# MDPDE FIT WITHOUT 46
fitRob_46 <- MDPDE_BETA(y=y46, X=X46, Z=Z46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_46

# MDPDE FIT WITHOUT 25
fitRob_25 <- MDPDE_BETA(y=y25, X=X25, Z=Z25, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac =0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25

# MDPDE FIT WITHOUT 25 and 46
fitRob_25_46 <- MDPDE_BETA(y=y25_46, X=X25_46, Z=Z25_46, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_25_46

#************************************BOOTSTRAP PVALUE*************************#
source("pvalue_boot_MDPDE.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y,X=X,Z=Z,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 46***#
pvalue_boot(y=y46,X=X46,Z=Z46,B=B, optimal=T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y46,X=X46,Z=Z46,B=B, optimal=F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 25***#
pvalue_boot(y=y25,X=X25,Z=Z25,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y25,X=X25,Z=Z25,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 25 AND 46***#
pvalue_boot(y=y25_46,X=X25_46,Z=Z25_46,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y25_46,X=X25_46,Z=Z25_46,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE
#*********************************************************************************#