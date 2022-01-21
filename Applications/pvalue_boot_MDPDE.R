#*************************************************************************Description*************************************************************************************************#
# Bootstrap p-values for beta regression fit under MLE and MDPDE
#***ARGUMENTS***#
# y - response variable.
# X -  regressor matrix for the mean submodel.
# Z -  regressor matrix for the precision submodel.
# B -  number of boostrap replicates. Default is 500. 
# optimal - logical; if TRUE, the MDPDE is used to calculate the p-values; if FALSE, the MLE is used to calculate the p-values. Default is TRUE.
# linkmu - character specification of the link function in the mean submodel. Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Default is "logit".
# linkphi - character specification of the link function in the precision submodel. Currently, "identity", "log", "sqrt" are supported. Default is "log".
#**************************************************************************************************************************************************************************************#

pvalue_boot <- function(y, X, Z, B=500, optimal = T, linkmu="logit", linkphi="log") { 
source("MDPDE.r")
n <- nrow(X)
kk1 <- ncol(X)
kk2 <- ncol(Z)
#****fit observed sample****#
fit<- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=optimal, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu=linkmu, linkphi=linkphi)
#***observed sample test statistic (includes intercepts)***#
chi2_original_beta <- (fit$beta/fit$se_beta)^2
chi2_original_gama <- (fit$gama/fit$se_gama)^2
chi2_original <- c(chi2_original_beta, chi2_original_gama) 
  if(kk1>1){
    #without intercept#
    chi2_B_beta<- matrix(numeric(0), ncol=kk1-1, nrow=B);
    muhat_fitx_aux <- matrix(numeric(0), ncol=kk1-1, nrow=n)
    phihat_fitx_aux <- matrix(numeric(0), ncol=kk1-1, nrow=n) 

    #*** Estimating the regression models without one covariate at a time ***#
    #*** Comment: This is made to generate the bootstrap response variable y under H0, i.e., under the respective coefficient equal to zero***#
           for(j in 2:kk1){
                #excluding one column at a time in X
                fitx_aux <- suppressWarnings(MDPDE_BETA(y=y, X=as.matrix(X[,-j]), Z=Z, qoptimal= optimal, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",
                startV="CP", linkmu=linkmu, linkphi=linkphi))
                beta_fitx_aux <- fitx_aux$beta
                gama_fitx_aux <- fitx_aux$gama
                eta_fitx_aux <- as.matrix(X[,-j])%*%as.vector(beta_fitx_aux)
                delta_fitx_aux <- Z%*%as.vector(gama_fitx_aux) 
                muhat_fitx_aux[,j-1] <- exp(eta_fitx_aux)/(1+exp(eta_fitx_aux))
                phihat_fitx_aux[,j-1] <- exp(delta_fitx_aux)  
          }
  }
  if(kk2>1){
    chi2_B_gama<- matrix(numeric(0),ncol=kk2-1,nrow=B)
    muhat_fitz_aux <- matrix(numeric(0), ncol=kk2-1, nrow=n) 
    phihat_fitz_aux <- matrix(numeric(0), ncol=kk2-1, nrow=n) 
          for(l in 2:kk2){
               #excluding one column at a time in Z
               fitz_aux <- suppressWarnings(MDPDE_BETA(y=y, X=X, Z=as.matrix(Z[,-l]), qoptimal= optimal, q0=1, m=3, L=0.02, qmin=0.5, spac = 0.02, method="BFGS",
               startV="CP", linkmu=linkmu, linkphi=linkphi))
               beta_fitz_aux <- fitz_aux$beta
               gama_fitz_aux <- fitz_aux$gama
               eta_fitz_aux <- X%*%as.vector(beta_fitz_aux)
               delta_fitz_aux <- as.matrix(Z[,-l])%*%as.vector(gama_fitz_aux) 
               muhat_fitz_aux[,l-1] <- exp(eta_fitz_aux)/(1+exp(eta_fitz_aux))
               phihat_fitz_aux[,l-1] <- exp(delta_fitz_aux) 
          }
  }
set.seed(c(1994,1991), kind="Marsaglia-Multicarry")
   #*****opens bootstrap*****#
  for(i in 1:B){
      if(kk1>1){
            for(j in 2:kk1){
              #generating under sob H0
              ygen_fitx_aux <- rbeta(n, muhat_fitx_aux[,j-1]*phihat_fitx_aux[,j-1], (1.0-muhat_fitx_aux[,j-1])*phihat_fitx_aux[,j-1])
                while(any(round(ygen_fitx_aux,5)==0|round(ygen_fitx_aux,5)==1)){	
                     ygen_fitx_aux <- rbeta(n, muhat_fitx_aux[,j-1]*phihat_fitx_aux[,j-1], (1.0-muhat_fitx_aux[,j-1])*phihat_fitx_aux[,j-1])			
                }
             #estimate with all covariates
             fitx <- suppressWarnings(MDPDE_BETA(y=ygen_fitx_aux, X=X, Z=Z, qoptimal= optimal, q0=1, m=3, L=0.02, qmin=0.5, spac = 0.02, method="BFGS",
             startV="CP", linkmu=linkmu, linkphi=linkphi))
             chi2_B_beta[i,j-1] <- (fitx$beta[j]/fitx$se_beta[j])^2
            }
      } 
      if(kk2>1){
            for(l in 2:kk2){
              #generating under H0
              ygen_fitz_aux <- rbeta(n, muhat_fitz_aux[,l-1]*phihat_fitz_aux[,l-1], (1.0-muhat_fitz_aux[,l-1])*phihat_fitz_aux[,l-1])	
                while(any(round(ygen_fitz_aux,4)==0|round(ygen_fitz_aux,4)==1)){
                     ygen_fitz_aux <- rbeta(n, muhat_fitz_aux[,l-1]*phihat_fitz_aux[,l-1], (1.0-muhat_fitz_aux[,l-1])*phihat_fitz_aux[,l-1])			
                }
             #estimate with all covariates
             fitz <- suppressWarnings(MDPDE_BETA(y=ygen_fitz_aux, X=X, Z=Z, qoptimal= optimal, q0=1, m=3, L=0.02, qmin=0.5, spac = 0.02, method="BFGS",
             startV="CP", linkmu=linkmu, linkphi=linkphi))
             chi2_B_gama[i,l-1] <- (fitz$gama[l]/fitz$se_gama[l])^2
            }
      }
  }#closes bootstrap
  #*******Calculating the bootstrap p-values*******#
  if(kk1>1) {
  pvalue_x_chi2 <- numeric(kk1-1)
       for(l in 2:kk1){
           pvalue_x_chi2[l-1] <- mean(chi2_B_beta[,l-1] > chi2_original_beta[l]);
       }
  }
  if(kk2>1) {
  pvalue_z_chi2 <- numeric(kk2-1)
      for(l in 2:kk2){
          pvalue_z_chi2[l-1] <- mean(chi2_B_gama[,l-1] > chi2_original_gama[l]);
      }
  }
  if(kk1==1){
     pvalue_x_chi2 <- NULL
  }
  if(kk2==1){
     pvalue_z_chi2 <- NULL
  }
results <- list(pvalue_beta = pvalue_x_chi2, pvalue_gama = pvalue_z_chi2)
return(results)
}


