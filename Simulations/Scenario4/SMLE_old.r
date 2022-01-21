#*************************************************************************Description*************************************************************************************************#
# SMLE for beta regression 
# Reference
# Ribeiro, T. K. A.; Ferrari, S. L. P. (2020). Robust estimation in beta regression via maximum Lq-likelihood.
#***ARGUMENTS***#
# y - response variable.
# X -  regressor matrix for the mean submodel.
# Z -  regressor matrix for the precision submodel.
# qoptimal - logical; if TRUE, the tuning constant q is chosen via the data-driven algorithm proposed by Ribeiro and Ferrari (2020)
# if FALSE, a fixed 0 < q0 <= 1 is required. Defaults are TRUE and q0 = 1, respectively.
# m - size of the grids for the data-driven algorithm. Default is m = 3.
# L - threshold for the data-driven algorithm. Default is L = 0.02.
# qmin -  minimum value for the tuning constant q. Default is 0.5.
# spac -  grid spacing for the data-driven algorithm. Default is 0.02.
# method - method passed to optim. Default is BFGS.
# maxit - maximum number of iterations. Default is 200.
# startV - character specification of the starting values. Currently, "CP" (constant precision model) and "VP" (varying precision model) are supported. Default is "CP".
# linkmu - character specification of the link function in the mean submodel. Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Default is "logit".
# linkphi - character specification of the link function in the precision model. Currently, "identity", "log", "sqrt" are supported. Default is "log".
# weights - logical; if TRUE, the weights given to each observation are printed. Default is FALSE.
#**************************************************************************************************************************************************************************************#

#------------------------->FUNCTION STARTS HERE<--------------------------------#
SMLE_BETA <- function(y, X, Z, qoptimal=TRUE, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", maxit=200, startV ="CP", linkmu="logit", linkphi="log", weights=FALSE){
#****Packages required****#
if(!suppressWarnings(require(betareg))){suppressWarnings(install.packages("betareg"));suppressWarnings(library(betareg))}#used for MLE fit
if(!suppressWarnings(require(robustbase))){suppressWarnings(install.packages("robustbase"));suppressWarnings(library(robustbase))}#used for robust starting values
if(!suppressWarnings(require(matlib))){suppressWarnings(install.packages("matlib"));suppressWarnings(library(matlib))}#used for Ginv function
if(!suppressWarnings(require(BBmisc))){suppressWarnings(install.packages("BBmisc"));suppressWarnings(library(BBmisc))} #used for is.error function

#**********link functions**********#
mean_link_functions <- function(beta, y, X, linkmu="logit"){
  kk1 <- ncol(X);  n <- nrow(X);  eta <- as.vector(X%*%beta)	 
 if(linkmu == "logit"){
  mu <- exp(eta)/(1.0+exp(eta))
  T_1 <- diag(mu*(1.0-mu))
  yg1 <- log(y/(1-y))
 }
 if(linkmu == "probit"){
  mu <- pnorm(eta) 
  T_1 <- diag(dnorm(qnorm(mu)))
  yg1 <- qnorm(y)
 }
 if(linkmu == "cloglog"){
  mu <- 1.0 - exp(-exp(eta)) 
  T_1 <- diag((mu - 1)*log(1 - mu))
  yg1 <- log(-log(1-y))
 }
 if(linkmu == "log"){
  mu <- exp(eta) 
  T_1 <- diag(mu)
  yg1 <- log(y)
 }
 if(linkmu == "loglog"){
  mu <- exp(-exp(-eta)) 
  T_1 <- diag(-mu*log(mu))
  yg1 <- -log(-log(y))
 }
 if(linkmu == "cauchit"){
  mu <- (pi^(-1))*atan(eta) + 0.5   
  T_1 <- diag((1/pi)*((cospi(mu-0.5))^2))     
  yg1 <- tan(pi*(y-0.5))                          
 }
  results <- list(mu= mu, T_1 = T_1, yg1=yg1)
  return(results)
}#ends mean link function

precision_link_functions <- function(gama, Z, linkphi="log"){
  kk2 <- ncol(Z);  n <- nrow(Z);  delta <- as.vector(Z%*%gama) 
 if(linkphi == "log"){
  phi <- exp(delta) 
  T_2 <- diag(phi)
 }
 if(linkphi == "identify"){
  phi <- delta 
  T_2 <- diag(rep(1,n))
 }
 if(linkphi == "sqrt"){
  phi <- delta^2 
  T_2 <- diag(2*sqrt(phi))
 }
  results <- list(phi=phi, T_2=T_2)
  return(results)
}#ends precision link function

#*************Function to be maximized**********************#
log_liksurrogate <- function(theta) 
{
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]                                                  
  mu_q <- mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$mu
  phi_q <- precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$phi
  phi <-(1.0/q_const)*(phi_q - 2.0) + 2.0
  mu <- ((1.0/q_const)*(mu_q*phi_q - 1.0) + 1.0)/phi                                  
  a <- mu*phi						 
  b <- (1.0 - mu)*phi		
  log_likS <- sum(dbeta(y, a, b, log = F)^(1.0 - q_const))#function to be maximized  
  return(log_likS)
}

#********************Score-type function*******************#
Score <- function(theta) { 
  avScore = numeric(0) 
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]                                                  
  mu_q <- mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$mu
  phi_q <- precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$phi
  T_1_q <- mean_link_functions(beta = beta, y=y, X=X, linkmu=linkmu)$T_1
  T_2_q <- precision_link_functions(gama = gama, Z=Z, linkphi=linkphi)$T_2
  phi <-(1.0/q_const)*(phi_q - 2.0)+2.0
  mu <- ((1.0/q_const)*(mu_q*phi_q - 1.0) + 1.0)/phi                                            
  a <- mu*phi						 
  b <- (1.0 - mu)*phi
  m_phiq <- diag(phi_q)
  F_q <- diag(dbeta(y, a, b)^(1.0 - q_const))
  ystar <- log(y/(1.0 - y))
  ydagger <- log(1.0 - y)
  mustar <- psigamma(a, 0) - psigamma(b, 0)
  mudagger <-  psigamma(b, 0) - psigamma(a + b, 0)
  #****vector type Score*****#
  avScore[1:kk1] <- (q_const^(-1.0))*t(X)%*%m_phiq%*%F_q%*%T_1_q%*%(ystar - mustar)
  avScore[(kk1+1.0):(kk1+kk2)] <- (q_const^(-1.0))*t(Z)%*%F_q%*%T_2_q%*%(mu_q*(ystar - mustar) + ydagger - mudagger)
  return(avScore)
}#ends score-type function

#********************Covariance matrix****************#
V_matrix <- function(theta) {  
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta_p <- theta[1:kk1]
  gama_p <- theta[(kk1+1.0):(kk1+kk2)]    	   
  muhat_q <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$mu
  phihat_q <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$phi  
  ahat_q <- muhat_q*phihat_q
  bhat_q <- (1.0 - muhat_q)*phihat_q
  T_1_q <- mean_link_functions(beta = beta_p, y=y, X=X, linkmu=linkmu)$T_1
  T_2_q <- precision_link_functions(gama = gama_p, Z=Z, linkphi=linkphi)$T_2
  phihat_n <-   (1.0/q_const)*(phihat_q - 2.0) + 2.0
  muhat_n <- 	((1.0/q_const)*(muhat_q*phihat_q  - 1.0) + 1.0)/phihat_n
  ahat_n <- muhat_n*phihat_n
  bhat_n <- (1.0 - muhat_n)*phihat_n
  phihat_2_q <- (2.0 - q_const)*(phihat_n - 2.0) + 2.0;	#expression of phi_(2-q)
  muhat_2_q <- ((2.0 - q_const)*(ahat_n - 1.0) + 1.0)/phihat_2_q; #expression of mu_(2-q)
  a2_qhat <- muhat_2_q*phihat_2_q
  b2_qhat <- (1.0 - muhat_2_q)*phihat_2_q
  mustarhat_n <- psigamma(ahat_n, 0) - psigamma(bhat_n, 0)
  mudaggerhat_n <-  psigamma(bhat_n, 0) - psigamma(phihat_n, 0)		
  mustarhat_2_q <- psigamma(a2_qhat, 0) - psigamma(b2_qhat, 0)
  mudaggerhat_2_q <-  psigamma(b2_qhat, 0) - psigamma(phihat_2_q, 0)	
  muhat_d_2_q <- muhat_q*(mustarhat_2_q - mustarhat_n) + mudaggerhat_2_q - mudaggerhat_n
  m_phiq <- diag(phihat_q)
  psi1_n <- psigamma(ahat_n, 1.0) 
  psi2_n <- psigamma(bhat_n, 1.0) 
  psi3_n <- psigamma(phihat_n, 1.0) 
  V_n <- diag(psi1_n + psi2_n)	
  B1 <- diag(exp(q_const*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - (lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))
  B2 <- diag(exp(lgamma(a2_qhat) + lgamma(b2_qhat) - lgamma(phihat_2_q) - (2.0*(1.0-q_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n))  +
  lgamma(ahat_q) + lgamma(bhat_q) - lgamma(phihat_q))))	
  C_q_0 <- diag(phihat_q*(muhat_q*psi1_n - (1.0 - muhat_q)*psi2_n))	
  D_q_0 <- diag((muhat_q^2.0)*psi1_n + ((1.0 - muhat_q)^2.0)*psi2_n - psi3_n) 
  psi1_2_q <- psigamma(a2_qhat, 1.0)	
  psi2_2_q <- psigamma(b2_qhat, 1.0)
  psi3_2_q <- psigamma(phihat_2_q, 1.0)		
  V_2_q <- diag(psi1_2_q + psi2_2_q)
  C_q_2_q <- diag(phihat_q*(muhat_q*psi1_2_q - (1.0 - muhat_q)*psi2_2_q))
  D_q_2_q <- diag((muhat_q^2.0)*psi1_2_q + ((1.0 - muhat_q)^2.0)*psi2_2_q - psi3_2_q)	
  M1 <- diag(mustarhat_2_q - mustarhat_n)
  M2 <- diag(muhat_d_2_q)
  M3 <- diag((mustarhat_2_q - mustarhat_n)*muhat_d_2_q)									
  Jq_betabeta <- as.matrix(t(X)%*%B1%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%V_n%*%X)
  Jq_betagamma <- as.matrix(t(X)%*%B1%*%T_1_q%*%T_2_q%*%C_q_0%*%Z)
  Jq_gammagamma <- as.matrix(t(Z)%*%B1%*%(T_2_q^2.0)%*%D_q_0%*%Z)  
  Jq <- matrix(numeric(0), kk1+kk2, kk1+kk2) 
  Jq[1:kk1,1:kk1] <- Jq_betabeta
  Jq[1:kk1,(kk1+1):(kk1+kk2)] <- Jq_betagamma
  Jq[(kk1+1):(kk1+kk2),1:kk1] <- t(Jq_betagamma)
  Jq[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Jq_gammagamma)
  Jq <- -(q_const^(-1.0))*Jq
  Kq_betabeta <- as.matrix(t(X)%*%B2%*%(T_1_q^2.0)%*%(m_phiq^2.0)%*%(V_2_q+M1^2.0)%*%X) 
  Kq_betagamma <- as.matrix(t(X)%*%B2%*%T_1_q%*%T_2_q%*%m_phiq%*%((C_q_2_q*(1/phihat_q)) + M3)%*%Z)	 
  Kq_gammagamma <- as.matrix(t(Z)%*%B2%*%(T_2_q^2.0)%*%(D_q_2_q+M2^2.0)%*%Z)
  Kq <- matrix(numeric(0),kk1+kk2,kk1+kk2)
  Kq[1:kk1,1:kk1] <- Kq_betabeta
  Kq[1:kk1,(kk1+1):(kk1+kk2)] <- Kq_betagamma
  Kq[(kk1+1):(kk1+kk2),1:kk1] <- t(Kq_betagamma)
  Kq[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Kq_gammagamma)
  Kq <- (q_const^(-2.0))*Kq	
  Vq <- tryCatch( solve(Jq)%*%Kq%*%t(solve(Jq)), error=function(e) {e})  #asymptotic covariance matrix
   if(is.error(Vq)){
     Vq <- Ginv(Jq)%*%Kq%*%t(Ginv(Jq)) 
   }
return(Vq)
}#ends Covariance matrix function

#**********************Starting Values function*****************#
Starting_values <- function(y, X, Z, startV="CP", linkmu, linkphi){
  kk1 <- ncol(X);kk2 <- ncol(Z);n <- nrow(X)
#********initial values based on constant precision*******#
 if(startV=="CP"){# opens startV "CP"
  #*******IID CASE******#
  if(kk1==1 && kk2==1){#opens 1
   fit_mle_start <- betareg(y~1|1,link = linkmu, link.phi= linkphi)   #MLE FIT
   initials_mle <- as.numeric(c(coef(fit_mle_start)[1:kk1],coef(fit_mle_start)[kk1+1])) #starting values iid case
   fit_mle <- fit_mle_start 
   thetaStart <- initials_mle 
   #***opens robust startV***#
   beta_mle_start <- initials_mle[1:kk1]
   gama_mle_start <- initials_mle[(kk1+1.0):(kk1+kk2)]    	   
   muhat_mle_start <- mean_link_functions(beta = beta_mle_start,y=y, X=X, linkmu=linkmu)$mu
   phihat_mle_start <- precision_link_functions(gama = gama_mle_start, Z=Z, linkphi=linkphi)$phi
   yg <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$yg1
      if(any(muhat_mle_start*phihat_mle_start<1) | any((1-muhat_mle_start)*phihat_mle_start<1)){### evaluating condition muphi<1 and/or (1-mu)phi<1
         fit_lm_Rob <- suppressWarnings(lmrob(yg~1)) #robust regression to obtain a starting value for beta
         initials_beta_rob <-  as.numeric(coef(fit_lm_Rob)) #beta estimates from robustbase package 
         muhat_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$mu  #mu estimate from robust method
         if(linkmu == "logit") muhatg <- log(muhat_Rob/(1-muhat_Rob))
         if(linkmu == "probit") muhatg <- qnorm(muhat_Rob)
         if(linkmu == "cloglog") muhatg <- log(-log(1-muhat_Rob))
         if(linkmu == "log") muhatg <- log(muhat_Rob)
         if(linkmu == "loglog") muhatg <- -log(-log(muhat_Rob))
         if(linkmu == "cauchit") muhatg <- tan(pi*(muhat_Rob-0.5))       
         T_1_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$T_1
         sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*diag(T_1_Rob^2)) #estimate of variance of y based on robustbase package
         phihat_Rob <- mean((muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob) #estimate of phi from robust method
         if(linkphi == "log") gama1hat_rob <- log(phihat_Rob)
         if(linkphi == "identify") gama1hat_rob <- phihat_Rob
         if(linkphi == "sqrt") gama1hat_rob <- sqrt(phihat_Rob)
         initials_gama_rob <-  as.numeric(gama1hat_rob) #gamma estimates from robustbase package
         thetaStart <- c(initials_beta_rob, initials_gama_rob) #robust starting values for beta and gamma
      }#closes if for the evaluation of the condition
  }#closes 1
  #***constant precision case***#
  if(kk1>1&&kk2==1){ #opens 2
   fit_mle_start <- betareg(y~X[,2:kk1]|1,link = linkmu, link.phi= linkphi)# MLE FIT
   initials_mle <- as.numeric(c(coef(fit_mle_start)[1:kk1],coef(fit_mle_start)[kk1+1])) #starting values for constant precision
   fit_mle <- fit_mle_start
   thetaStart <- initials_mle
   #***starts robust startV***#
   beta_mle_start <- initials_mle[1:kk1]
   gama_mle_start <- initials_mle[(kk1+1.0):(kk1+kk2)]    	   
   muhat_mle_start <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$mu
   phihat_mle_start <- precision_link_functions(gama = gama_mle_start, Z=Z, linkphi=linkphi)$phi
   yg <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$yg1
      if(any(muhat_mle_start*phihat_mle_start<1) | any((1-muhat_mle_start)*phihat_mle_start<1)  ){### evaluating condition muphi<1 and (1-mu)phi<1
         fit_lm_Rob <- suppressWarnings(lmrob(yg~X[,2:kk1])) #robust regression to obtain a starting value for beta
         initials_beta_rob <-  as.numeric(coef(fit_lm_Rob)) #beta estimates from robustbase package 
         muhat_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$mu  #mu estimate from robust method
         if(linkmu == "logit") muhatg <- log(muhat_Rob/(1-muhat_Rob))
         if(linkmu == "probit") muhatg <- qnorm(muhat_Rob)
         if(linkmu == "cloglog") muhatg <- log(-log(1-muhat_Rob))
         if(linkmu == "log") muhatg <- log(muhat_Rob)
         if(linkmu == "loglog") muhatg <- -log(-log(muhat_Rob))
         if(linkmu == "cauchit") muhatg <- tan(pi*(muhat_Rob-0.5))                          
         T_1_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$T_1
         sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*diag(T_1_Rob^2)) #estimate of variance of y based on robustbase package
         phihat_Rob <- mean((muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob) #phi estimate from robust method
         if(linkphi == "log") gama1hat_rob <- log(phihat_Rob)
         if(linkphi == "identify") gama1hat_rob <- phihat_Rob
         if(linkphi == "sqrt") gama1hat_rob <- sqrt(phihat_Rob)
         initials_gama_rob <-  as.numeric(gama1hat_rob) #gamma estimates from robustbase package 
         thetaStart <- c(initials_beta_rob, initials_gama_rob) #robust starting values for beta and gamma
      }#closes if for the evaluation of the condition
  }# closes 2
  if(kk1==1&&kk2>1){ #opens 3
   fit_mle <- betareg(y~1|Z[,2:kk2],link = linkmu, link.phi= linkphi)
   fit_mle_start <- betareg(y~1|1,link = linkmu, link.phi= linkphi)
   initials_mle <- as.numeric(c(coef(fit_mle_start)[1:kk1], coef(fit_mle_start)[kk1+1], rep(0,kk2-1))) #starting values for varying precision based on constant precision
   #***starts robust startV***#
   beta_mle_start <- initials_mle[1:kk1]
   gama_mle_start <- initials_mle[(kk1+1.0)]
   Xaux <- matrix(c(rep(1,n)),ncol=1,byrow=F);
   Zaux <- matrix(c(rep(1,n)),ncol=1,byrow=F);
   thetaStart <- c(beta_mle_start, gama_mle_start, rep(0,kk2-1))	   
   muhat_mle_start <- mean_link_functions(beta = beta_mle_start, y=y, X=Xaux, linkmu=linkmu)$mu
   phihat_mle_start <- precision_link_functions(gama = gama_mle_start, Z=Zaux, linkphi=linkphi)$phi
   yg <- mean_link_functions(beta = beta_mle_start, y=y, X=Xaux, linkmu=linkmu)$yg1
      if(any(muhat_mle_start*phihat_mle_start<1) | any((1-muhat_mle_start)*phihat_mle_start<1)  ){### evaluating condition muphi<1 and (1-mu)phi<1
         fit_lm_Rob <- suppressWarnings(lmrob(yg~1)) #robust regression to obtain a starting value for beta
         initials_beta_rob <-  as.numeric(coef(fit_lm_Rob)) #beta estimates from robustbase package 
         muhat_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=Xaux, linkmu=linkmu)$mu  #mu estimate from robust method
         if(linkmu == "logit") muhatg <- log(muhat_Rob/(1-muhat_Rob))
         if(linkmu == "probit") muhatg <- qnorm(muhat_Rob)
         if(linkmu == "cloglog") muhatg <- log(-log(1-muhat_Rob))
         if(linkmu == "log") muhatg <- log(muhat_Rob)
         if(linkmu == "loglog") muhatg <- -log(-log(muhat_Rob))
         if(linkmu == "cauchit") muhatg <- tan(pi*(muhat_Rob-0.5))        
         T_1_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=Xaux, linkmu=linkmu)$T_1
         sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*diag(T_1_Rob^2)) #estimate of variance of y based on robust package
         phihat_Rob <- mean((muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob) #phi estimate from robust method
         if(linkphi == "log") gama1hat_rob <- log(phihat_Rob) 
         if(linkphi == "identify") gama1hat_rob <- phihat_Rob
         if(linkphi == "sqrt") gama1hat_rob <- sqrt(phihat_Rob)
         initials_gama_rob <-  c(as.numeric(gama1hat_rob), rep(0,kk2-1)) #gamma estimates from robustbase package 
         thetaStart <- c(initials_beta_rob, initials_gama_rob) #robust starting values for beta and gamma
      }#closes if for the evaluation of the condition
  }# closes 3
  #****varying precision****#
  if(kk1>1&&kk2>1){#opens 4
   fit_mle <- betareg(y~X[,2:kk1]|Z[,2:kk2],link = linkmu, link.phi= linkphi)
   fit_mle_start <- betareg(y~X[,2:kk1]|1,link = linkmu, link.phi= linkphi)
   initials_mle <- as.numeric(c(coef(fit_mle_start)[1:kk1], coef(fit_mle_start)[kk1+1], rep(0,kk2-1))) #starting values for varying precision based on constant precision
   #***starts robust startV***#
   beta_mle_start <- initials_mle[1:kk1]
   gama_mle_start <- initials_mle[(kk1+1.0)]
   thetaStart <- c(beta_mle_start, gama_mle_start, rep(0,kk2-1))	   
   muhat_mle_start <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$mu
   Zaux <- matrix(c(rep(1,n)),ncol=1,byrow=F);
   phihat_mle_start <- precision_link_functions(gama = gama_mle_start, Z=Zaux, linkphi=linkphi)$phi
   yg <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$yg1
      if(any(muhat_mle_start*phihat_mle_start<1) | any((1-muhat_mle_start)*phihat_mle_start<1)  ){### evaluating condition muphi<1 and (1-mu)phi<1
         fit_lm_Rob <- suppressWarnings(lmrob(yg~X[,2:kk1])) #robust regression to obtain a starting value for beta
         initials_beta_rob <-  as.numeric(coef(fit_lm_Rob)) #beta estimates from robustbase package 
         muhat_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$mu  #mu estimate from robust method
         if(linkmu == "logit") muhatg <- log(muhat_Rob/(1-muhat_Rob))
         if(linkmu == "probit") muhatg <- qnorm(muhat_Rob)
         if(linkmu == "cloglog") muhatg <- log(-log(1-muhat_Rob))
         if(linkmu == "log") muhatg <- log(muhat_Rob)
         if(linkmu == "loglog") muhatg <- -log(-log(muhat_Rob))
         if(linkmu == "cauchit") muhatg <- tan(pi*(muhat_Rob-0.5))        
         T_1_Rob <- mean_link_functions(beta = initials_beta_rob, y=y, X=X, linkmu=linkmu)$T_1
         sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*diag(T_1_Rob^2)) #estimate of variance of y based on robustbase package
         phihat_Rob <- mean((muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob) #phi estimate from robust method
         if(linkphi == "log") gama1hat_rob <- log(phihat_Rob) 
         if(linkphi == "identify") gama1hat_rob <- phihat_Rob
         if(linkphi == "sqrt") gama1hat_rob <- sqrt(phihat_Rob)
         initials_gama_rob <-  c(as.numeric(gama1hat_rob), rep(0,kk2-1)) #gamma estimates from robustbase package
         thetaStart <- c(initials_beta_rob, initials_gama_rob) #robust starting values for beta and gama
      }#closes if for the evaluation of the condition
  }#closes 4
 }#closes if "CP"
 if(startV=="VP"){
   fit_mle <- betareg(y~X[,2:kk1]|Z[,2:kk2],link = linkmu, link.phi= linkphi)
   fit_mle_start <- fit_mle # MLE FIT
   initials_mle <- as.numeric(coef(fit_mle_start)) #starting values for varying precision based on varying precision ("VP")
   thetaStart <- initials_mle
   #***starts robust startV***#
   beta_mle_start <- initials_mle[1:kk1]
   gama_mle_start <- initials_mle[(kk1+1.0):(kk1+kk2)]    	   
   muhat_mle_start <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$mu
   phihat_mle_start <- precision_link_functions(gama = gama_mle_start, Z=Z, linkphi=linkphi)$phi
   yg <- mean_link_functions(beta = beta_mle_start, y=y, X=X, linkmu=linkmu)$yg1
      if(any(muhat_mle_start*phihat_mle_start<1) | any((1-muhat_mle_start)*phihat_mle_start<1)  ){### evaluating condition muphi<1 and (1-mu)phi<1
         fit_lm_Rob <- suppressWarnings(lmrob(yg~X[,2:kk1])) #robust regression to obtain a starting value for beta
         initials_beta_rob <-  as.numeric(coef(fit_lm_Rob)) #beta estimates from robustbase package 
         muhat_Rob <- mean_link_functions(beta = initials_beta_rob,y=y, X=X, linkmu=linkmu)$mu  #mu estimate from robust method
         if(linkmu == "logit") muhatg <- log(muhat_Rob/(1-muhat_Rob))
         if(linkmu == "probit") muhatg <- qnorm(muhat_Rob)
         if(linkmu == "cloglog") muhatg <- log(-log(1-muhat_Rob))
         if(linkmu == "log") muhatg <- log(muhat_Rob)
         if(linkmu == "loglog") muhatg <- -log(-log(muhat_Rob))
         if(linkmu == "cauchit") muhatg <- tan(pi*(muhat_Rob-0.5))       
         T_1_Rob <- mean_link_functions(beta = initials_beta_rob, y=y,X=X, linkmu=linkmu)$T_1
         sigma2hat_Rob <- ((fit_lm_Rob$scale^2)*diag(T_1_Rob^2)) #estimate of variance of y based on robustbase package
         phihat_Rob <- (muhat_Rob*(1-muhat_Rob))/sigma2hat_Rob #phi estimate from robust method
         if(linkphi == "log") phihatg <- log(phihat_Rob)
         if(linkphi == "identify") phihatg <- phihat_Rob
         if(linkphi == "sqrt") phihatg <- sqrt(phihat_Rob)
         fit2_lm_Rob <- suppressWarnings(lmrob(phihatg~Z[,2:kk1])) #robust regression to obtain a starting value for gama
         initials_gama_rob <-  as.numeric(coef(fit2_lm_Rob)) #gamma estimates from robustbase package 
         thetaStart <- c(initials_beta_rob, initials_gama_rob) #robust starting values for beta and gama
      }#closes if for the evaluation of the condition
 }#closes if "VP"
  results <- list(thetaStart=thetaStart, MLEfit= fit_mle)
  return(results)
}#ends function of starting values

#----------------------------------------------------------------------------------------------------------------#
#------------------------------------->Estimation procedure starts here<-----------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
   
  #***Evaluating MLE and starting values for SMLE***#
  kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X)
  SV <- Starting_values(y=y, X=X, Z=Z, startV=startV, linkmu=linkmu, linkphi=linkphi)
  thetaStart <- SV$thetaStart 
  fit_mle <- SV$MLEfit
  theta_mle <- as.numeric(coef(fit_mle))
  se_mle <- as.numeric(sqrt(diag(fit_mle$vcov)))
  trace_mle <-  sum(se_mle^2)	
  z_mle <- as.matrix(theta_mle/se_mle)
  pvalue_mle <- apply(z_mle,1,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
  etahat_mle <- X%*%theta_mle[1:kk1]
  muhat_mle <- mean_link_functions(beta = theta_mle[1:kk1], y=y, X=X, linkmu=linkmu)$mu
  phihat_mle <- precision_link_functions(gama = theta_mle[(kk1+1):(kk1+kk2)], Z=Z, linkphi=linkphi)$phi
  ygm <- mean_link_functions(beta = theta_mle[1:kk1], y=y, X=X, linkmu=linkmu)$yg1
  R2_mle <- cor(ygm, etahat_mle)^2.0 #pseudo R2
  #********RESULTS for q0=1********#
  results_mle_w <- list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
  se_gama <- se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
  pvalue_gama <- pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle), weights= rep(1,n)) #results with weights
  results_mle <- list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
  se_gama = se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
  pvalue_gama = pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle)) #results without weights

  #------------------------------->If q is fixed<----------------------------------#
     if(qoptimal==FALSE){#starts option "qoptimal==FALSE"
          if(q0==1){#******q0 equal to 1 (MLE)*******#
                     if(weights==T){
                          results <- results_mle_w
                     }
                     else{
                        results <- results_mle
                     }#closes else of weights
              return(results) #returns results based on MLE
                   }#closes if for q0 = 1 
          else{ #opens if for q0 < 1 
              q_const <- q0
              #-------->maximization with q fixed<------------#
              res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate,hessian = F, control=list(fnscale=-1,maxit=maxit),
              gr = Score, method=method), error=function(e) {e}))#maximization
                   if(is.error(res)){#starts else in case of error in the estimation procedure
                       if(weights==T){
                           results <- results_mle_w
                       }
                       else{
                           results <- results_mle
                       }#ends else of weights
                     return(results)
                   }   
                   else{
                        if(res$convergence == 0||res$convergence == 1){ 
                               estimates <- res$par
                               log.lik=res$value
                               Vq <- V_matrix(estimates)	
                               se_smle <-  suppressWarnings(t(sqrt(diag(Vq))))	 
                               trace_smle <-  sum(diag(Vq))	
                               z_smle <- as.matrix(estimates/se_smle)
                               pvalue_smle <- apply(z_smle,2,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
                               etahat_smle <- X%*%estimates[1:kk1]
                               deltahat_smle <- Z%*%estimates[(kk1+1):(kk1+kk2)] 
                               muhat_smle <- mean_link_functions(beta = estimates[1:kk1], y=y, X=X, linkmu=linkmu)$mu
                               phihat_smle <- precision_link_functions(gama = estimates[(kk1+1):(kk1+kk2)], Z=Z, linkphi=linkphi)$phi
                               ygm <- mean_link_functions(beta = estimates[1:kk1], y=y, X=X, linkmu=linkmu)$yg1
                               R2_smle <- cor(ygm, etahat_smle)^2.0 #pseudo R2
                                   if(weights==T){
                                      w <- (dbeta(y,muhat_smle*phihat_smle,  (1-muhat_smle)*phihat_smle, log = F))^(1.0 - q0)
                                      results <- list(beta = estimates[1:kk1], gama = estimates[(kk1+1):(kk1+kk2)], se_beta = se_smle[1:kk1],
                                      se_gama <- se_smle[(kk1+1):(kk1+kk2)], z_beta = z_smle[1:kk1], z_gama = z_smle[(kk1+1):(kk1+kk2)],
                                      pvalue_beta = pvalue_smle[1:kk1],
                                      pvalue_gama <- pvalue_smle[(kk1+1):(kk1+kk2)], q_value = q0, R2 = as.numeric(R2_smle), weights = as.numeric(t(w)),
                                      Initial_Values = thetaStart)
                                   }
                                   else{
                                   results <- list(beta = estimates[1:kk1], gama = estimates[(kk1+1):(kk1+kk2)], se_beta = se_smle[1:kk1],
                                   se_gama = se_smle[(kk1+1):(kk1+kk2)], z_beta = z_smle[1:kk1], z_gama = z_smle[(kk1+1):(kk1+kk2)],
                                   pvalue_beta = pvalue_smle[1:kk1],
                                   pvalue_gama = pvalue_smle[(kk1+1):(kk1+kk2)], q_value = q0, R2 = as.numeric(R2_smle),Initial_Values = thetaStart)
                                   }
                          return(results)
                        }
                        else {# in case of no convergence, take the MLE
                            if(weights==T){
                                results <- results_mle_w
                            }
                            else{
                                results <- results_mle
                            }#closes else of weights 
                                return(results)
                        }
                   }
          }#closes else of q0 < 1 

  #------------------------------->If q is not fixed<----------------------------------#

     }else {#starts option "qoptimal==TRUE"
         #***grid of values for q***#
         q_values <- seq(from=qmin,to=1.0,by=spac)
         q_values <- sort(q_values, decreasing = TRUE)
         nq <- length(q_values)
         thetahat_n <- matrix(numeric(0),nrow=kk1+kk2,ncol= nq)
         se_smle <- matrix(numeric(0),nrow=kk1+kk2,ncol= nq)
         trace_smle <- matrix(numeric(0),nrow=1,ncol= nq)		  
         thetahat_n[,1] =  theta_mle
         se_smle[,1] =  se_mle
         trace_smle[,1] = trace_mle	  	 
         counter <- 1
         grid = m
         f1 <- 0.0; fc <- 0; f1e <- 0.0; cfailure <- 0.0
         q_const <- 1.0 
         #***************Searching for the optimal q starts here****************#    
               while(grid > 0 && q_const >= qmin){
                    q_const <- q_values[1.0 + counter]
                    res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=-1, maxit=maxit), 
                    gr = Score, method=method), error=function(e) {e}))#maximization
                          if(is.error(res)&& q_const >= qmin){
                                 counter <- counter + 1.0
                          }
                          else{# if estimation procedure converged
                                 estimates <- res$par
                                 log.lik=res$value
                                 thetahat_n[,counter + 1.0] <- estimates
                                 Vq <- V_matrix(thetahat_n[,counter + 1.0])	
                                 se_smle[,counter + 1.0] =  suppressWarnings(t(sqrt(diag(Vq))))	  
                                 trace_smle[,counter + 1.0] =  sum(diag(Vq))	 
                                 #****checking stability of the estimates****#
                                 SQV <- round((1.0/(kk1+kk2))*sqrt(sum( (thetahat_n[,counter]/(sqrt(n)*se_smle[,counter]) - 
                                 thetahat_n[,counter+1.0]/(sqrt(n)*se_smle[,counter + 1.0]))^2)),5) 

                                        #***In case of no convergence or NaN in SQV***#
                                       if(is.nan(SQV)|| res$convergence == 10 || res$convergence == 51|| res$convergence == 52){
                                                grid = m # if no convergence or SQV=NaN, take one more q for the grid
                                                f1= f1 + 1.0  # sum of failures
                                                      if(f1 == nq - 1.0) {#if the number of failures is equal to the grid size, take SMLE = MLE  
                                                              grid = 0;	  #grid = 0 means that the search for the optimal value of q is over
                                                              q_optimal =	  1.0	 #q_optimal = 1 means that SMLE = MLE 
                                                      }
	                                 }	else{  	
                                       if(SQV < L){# if stability condition is satisfied
                                                grid = grid - 1.0 # subtract 1 in the grid	
                                                      if(grid==0) q_optimal = q_values[counter-m+1.0] #if grid = 0 take the maximum m (i.e, take the q-max of the grid)
	                                 }  else{  # if  stability condition is not satisfied	   
                                                      grid = m
                                                      f1e = f1e + 1.0
                                                                 if(f1e == nq - 1.0) {	# if for all grids, the stability condition is not satisfied
                                                                           grid = 0 #grid = 0 means that the search for the optimal value of q is over	  
                                                                           q_optimal =	  1.0	 #q_optimal = 1 means that SMLE = MLE 
     	                                                           }	
	                                    }
	                                    }

                                                     #*Checking convergence of SMLE*#
                                       if(res$convergence == 10 || res$convergence == 51|| res$convergence == 52){
                                                fc = fc + 1.0
                                                            if(fc==nq-1.0) {
                                                                      cfailure =  cfailure + 1.0

                                                            }
                                       }
                                 counter = counter + 1.0
                 
                          }#closes else of "the estimation procedure converged"

         if(is.error(res)&&q_const==qmin){#if qmin is reached and there is no convergence for some q
                    if(weights==T){
                          results <- list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
                          se_gama = se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
                          pvalue_gama = pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle),
                          Warning = "Lack of stability and non-convergence for some q values", weights = rep(1,n), 
                          Initial_Values = thetaStart)
                    }
                    else{
                          results <- list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
                          se_gama = se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
                          pvalue_gama = pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle), 
                          Warning = "Lack of stability and non-convergence for some q values", Initial_Values = thetaStart)
                    }
             return(results)
         }
        if((grid>0)&&q_const==qmin){#if qmin is reached and there is no stability of the estimates
                    if(weights==T){
                          results <- list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
                          se_gama = se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
                          pvalue_gama = pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle), Warning = "Lack of stability", 
                          weights = rep(1,n), Initial_Values = thetaStart)
                    }
                    else{
                          results <-list(beta = theta_mle[1:kk1], gama = theta_mle[(kk1+1):(kk1+kk2)], se_beta = se_mle[1:kk1],
                          se_gama = se_mle[(kk1+1):(kk1+kk2)], z_beta = z_mle[1:kk1], z_gama = z_mle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_mle[1:kk1],
                          pvalue_gama = pvalue_mle[(kk1+1):(kk1+kk2)], q_value = 1.0, R2 = as.numeric(R2_mle), Warning= "Lack of stability",
                          Initial_Values = thetaStart)
                    }
              return(results)
        }
               }
           #***************Searching for the optimal q ends here****************# 

         #****Selecting estimates corresponding to the optimal q****#
         seq <- 1:nq
         index_op <- seq[q_values==q_optimal]
         thehat_n_optimal <- thetahat_n[,index_op]
         se_smle_otimo <- se_smle[,index_op]
         trace_smle_optimal <- trace_smle[,index_op]
         z_smle <- as.matrix(thetahat_n[,index_op]/se_smle[,index_op])
         pvalue_smle <- apply(z_smle,1,function(x)(2.0*(1.0-pnorm(mean(abs(x))))))
         etahat_optimal <- X%*%thehat_n_optimal[1:kk1]
         deltahat_optimal <- Z%*%thehat_n_optimal[(kk1+1):(kk1+kk2)] 
         muhat_optimal <- mean_link_functions(beta = thehat_n_optimal[1:kk1], y=y, X=X, linkmu=linkmu)$mu
         phihat_optimal <- precision_link_functions(gama = thehat_n_optimal[(kk1+1):(kk1+kk2)], Z=Z, linkphi=linkphi)$phi
         ygm <- mean_link_functions(beta = thehat_n_optimal[1:kk1], y=y, X=X, linkmu=linkmu)$yg1
         R2 <- cor(ygm, etahat_optimal)^2.0 #pseudo R2
         if(weights==T){
                w <- (dbeta(y,muhat_optimal*phihat_optimal,  (1-muhat_optimal)*phihat_optimal, log = F))^(1.0 - q_optimal)
                results <- list(beta = thehat_n_optimal[1:kk1], gama = thehat_n_optimal[(kk1+1):(kk1+kk2)], se_beta = se_smle_otimo[1:kk1],
                se_gama = se_smle_otimo[(kk1+1):(kk1+kk2)], z_beta = z_smle[1:kk1], z_gama = z_smle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_smle[1:kk1],
                pvalue_gama = pvalue_smle[(kk1+1):(kk1+kk2)],q_value = q_optimal, R2 = as.numeric(R2), weights=as.numeric(t(w)), Initial_Values = thetaStart)
         }
         else{
                results <- list(beta = thehat_n_optimal[1:kk1], gama = thehat_n_optimal[(kk1+1):(kk1+kk2)], se_beta = se_smle_otimo[1:kk1],
                se_gama = se_smle_otimo[(kk1+1):(kk1+kk2)], z_beta = z_smle[1:kk1], z_gama = z_smle[(kk1+1):(kk1+kk2)], pvalue_beta = pvalue_smle[1:kk1],
                pvalue_gama = pvalue_smle[(kk1+1):(kk1+kk2)], q_value = q_optimal, R2 = as.numeric(R2), Initial_Values = thetaStart)
         }
           return(results)
      }# ends option "qoptimal==TRUE"
}# ends function
#------------------------->FUNCTION ENDS HERE<--------------------------------#

