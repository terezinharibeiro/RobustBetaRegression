source("MDPDE.r") 

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
fitMLE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE

# MLE FIT WITHOUT 15
fitMLE_15 <- MDPDE_BETA(y=y15, X=X15, Z=Z15, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",  
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15

# MLE FIT WITHOUT 16
fitMLE_16 <- MDPDE_BETA(y=y16, X=X16, Z=Z16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_16

# MLE FIT WITHOUT 72
fitMLE_72 <- MDPDE_BETA(y=y72, X=X72, Z=Z72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_72

# MLE FIT WITHOUT 15 AND 16
fitMLE_15_16 <- MDPDE_BETA(y=y15_16, X=X15_16, Z=Z15_16, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", 
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_16

# MLE FIT WITHOUT 15 AND 72
fitMLE_15_72 <- MDPDE_BETA(y=y15_72, X=X15_72, Z=Z15_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_72

# MLE FIT WITHOUT 16 AND 72
fitMLE_16_72 <- MDPDE_BETA(y=y16_72, X=X16_72, Z=Z16_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",  
startV="CP", linkmu="logit", linkphi="log")
fitMLE_16_72

# MLE FIT WITHOUT 15, 16 AND 72
fitMLE_15_16_72 <- MDPDE_BETA(y=y15_16_72, X=X15_16_72, Z=Z15_16_72, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,
method="BFGS", startV="CP", linkmu="logit", linkphi="log")
fitMLE_15_16_72

#************************MDPDE*********************#

# MDPDE FIT FULL DATA
fitRob <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob

# MDPDE FIT WITHOUT 15
fitRob_15 <- MDPDE_BETA(y=y15, X=X15, Z=Z15, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15

# MDPDE FIT WITHOUT 16
fitRob_16 <- MDPDE_BETA(y=y16, X=X16, Z=Z16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_16

# MDPDE FIT WITHOUT 72
fitRob_72 <- MDPDE_BETA(y=y72, X=X72, Z=Z72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_72

# MDPDE FIT WITHOUT 15 AND 16
fitRob_15_16 <- MDPDE_BETA(y=y15_16, X=X15_16, Z=Z15_16, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_16

# MDPDE FIT WITHOUT 15 AND 72
fitRob_15_72 <- MDPDE_BETA(y=y15_72, X=X15_72, Z=Z15_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_72

# MDPDE FIT WITHOUT 16 AND 72
fitRob_16_72 <- MDPDE_BETA(y=y16_72, X=X16_72, Z=Z16_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_16_72

# MDPDE FIT WITHOUT 15, 16 AND 72
fitRob_15_16_72 <- MDPDE_BETA(y=y15_16_72, X=X15_16_72, Z=Z15_16_72, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02,
method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=T)
fitRob_15_16_72

#*********************************BOOTSTRAP PVALUE*************************#
source("pvalue_boot_MDPDE.r")
B <- 500 #number of bootstrap replicates
#***FULL DATA***#
pvalue_boot(y=y,X=X,Z=X,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y,X=X,Z=X,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15***#
pvalue_boot(y=y15,X=X15,Z=Z15,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y15,X=X15,Z=Z15,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16***#
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y16,X=X16,Z=Z16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 72***#
pvalue_boot(y=y72,X=X72,Z=Z72,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y72,X=X72,Z=Z72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15 AND 16***#
pvalue_boot(y=y15_16,X=X15_16,Z=Z15_16,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y15_16,X=X15_16,Z=Z15_16,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15 AND 72***#
pvalue_boot(y=y15_72,X=X15_72,Z=Z15_72,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y15_72,X=X15_72,Z=Z15_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 16 AND 72***#
pvalue_boot(y=y16_72,X=X16_72,Z=Z16_72,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y16_72,X=X16_72,Z=Z16_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE

#***DATA WITHOUT 15, 16 AND 72***#
pvalue_boot(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,B=B, optimal =T, linkmu="logit",linkphi="log") #MDPDE
pvalue_boot(y=y15_16_72,X=X15_16_72,Z=Z15_16_72,B=B, optimal =F, linkmu="logit",linkphi="log") #MLE





