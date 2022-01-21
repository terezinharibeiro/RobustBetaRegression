source("SMLE.r")
source("MDPDE.r")
source("SAMPLE.r")

sink("Results_Scenario5_5p.txt")

names_results1 <- c("Estimates_Scenario5_n40", "Estimates_Scenario5_n80", "Estimates_Scenario5_n160", "Estimates_Scenario5_n320")
names_results2 <- c("SE_Estimates_Scenario5_n40", "SE_Estimates_Scenario5_n80", "SE_Estimates_Scenario5_n160", "SE_Estimates_Scenario5_n320")
names_results3 <- c("Qoptimal_Scenario5_n40", "Qoptimal_Scenario5_n80", "Qoptimal_Scenario5_n160", "Qoptimal_Scenario5_n320")

########################################################################
#######################Global Variables#################################
########################################################################
 
VN <- c(40, 80, 160, 320) #Sample sizes
VBETA <- c(-1.2, -2.5, -2.5) #true beta values
VGAMA <- c(3.8, 0.7, 0.7) #true gama values
VTHETA <- c(VBETA, VGAMA) # true theta values
NREP <- 1000 #number of Monte Carlo replicates
kk1 <- length(VBETA); kk2 <- length(VGAMA)

#*************generated covariates******************#
se1 = 2 ; se2 = 3 #random seed
set.seed(c(se1,se2), kind="Marsaglia-Multicarry") #To ensure repeatability of the experiment
x1 <- runif(VN[1]); x2 <- runif(VN[1]); z1 <- x1; z2 <- x2 
X <- matrix(c(rep(1,VN[1]), x1, x2), ncol=3, byrow=F); #regressor matrix for the mean submodel
Z <- matrix(c(rep(1,VN[1]), z1, z2), ncol=3, byrow=F); #regressor matrix for the precision submodel

for(l in 1:length(VN))#Loop for the sample size
{
N = VN[l]

#***Initializing vectors to store the estimates (MLE, SMLE and MDPDE)***#
thetahat <- matrix(NA, NREP, kk1 + kk2)
thetahat_C <- matrix(NA, NREP, kk1 + kk2)
thetahat_smle <- matrix(NA, NREP, kk1 + kk2)
thetahat_smle_C <- matrix(NA, NREP, kk1 + kk2)
thetahat_mdpde <- matrix(NA, NREP, kk1 + kk2)
thetahat_mdpde_C <- matrix(NA, NREP, kk1 + kk2)

se_thetahat <- matrix(NA, NREP, kk1 + kk2)
se_thetahat_C <- matrix(NA, NREP, kk1 + kk2)
se_thetahat_smle <- matrix(NA, NREP, kk1 + kk2)
se_thetahat_smle_C <- matrix(NA, NREP, kk1 + kk2)
se_thetahat_mdpde <- matrix(NA, NREP, kk1 + kk2)
se_thetahat_mdpde_C <- matrix(NA, NREP, kk1 + kk2)

z_thetahat <- matrix(NA, NREP, kk1 + kk2)
z_thetahat_C <- matrix(NA, NREP, kk1 + kk2)
z_thetahat_smle <- matrix(NA, NREP, kk1 + kk2)
z_thetahat_smle_C <- matrix(NA, NREP, kk1 + kk2)
z_thetahat_mdpde <- matrix(NA, NREP, kk1 + kk2)
z_thetahat_mdpde_C <- matrix(NA, NREP, kk1 + kk2)

trace_thetahat <- matrix(NA, NREP, 1.0)
trace_thetahat_C <- matrix(NA, NREP, 1.0)
trace_thetahat_smle <- matrix(NA, NREP, 1.0)
trace_thetahat_smle_C <- matrix(NA, NREP, 1.0)
trace_thetahat_mdpde <- matrix(NA, NREP, 1.0)
trace_thetahat_mdpde_C <- matrix(NA, NREP, 1.0)

qopt_thetahat_smle <- matrix(NA, NREP, 1.0)
qopt_thetahat_smle_C <- matrix(NA, NREP, 1.0)
qopt_thetahat_mdpde <- matrix(NA, NREP, 1.0)
qopt_thetahat_mdpde_C <- matrix(NA, NREP, 1.0)

cont <- 0 # Counter for Monte Carlo replicates 
while(cont < NREP)
{
cont <- cont + 1 
perc <- cont/NREP
#***To print the progress of the simulation***#
if(perc == 0.25 || perc == 0.5 || perc == 0.75 || perc ==1) cat("Perc. Replic. MC =",perc*100,"%","\n")

 SampleG <- FunSample(n = N, mXini= X, mZini = Z, size = VN[1], theta = c(VBETA, VGAMA), linkmu="logit", linkphi="log")
 y <- SampleG$y; mu <- SampleG$mu; phi <- SampleG$phi
 #********contamination process************#
 PC <- 0.05*N #percentage of contamination
 ind_C <- sort(mu, index.return=TRUE)$ix[1:PC] 
 mu_pc <- (1.0 + mu[ind_C])/2.0 
 yaux_c <- rbeta(PC, mu_pc*phi[ind_C], (1.0- mu_pc)*phi[ind_C])  
 y_c <- y; y_c[ind_C] <- yaux_c 

 # with contamination
 fitMLE   <- SMLE_BETA(y=SampleG$y, X=SampleG$X, Z=SampleG$Z, qoptimal= F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",startV="CP",linkmu="logit",linkphi="log", weights=F)
 fitSMLE  <- SMLE_BETA(y=SampleG$y, X=SampleG$X, Z=SampleG$Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",startV="CP",linkmu="logit",linkphi="log", weights=F)
 fitMDPDE <-MDPDE_BETA(y=SampleG$y, X=SampleG$X, Z=SampleG$Z, qoptimal= T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS",startV="CP",linkmu="logit",linkphi="log", weights=F)

 #without contamination
 fitMLE_C   <- SMLE_BETA(y=y_c, X=SampleG$X, Z=SampleG$Z, qoptimal=F, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=F)
 fitSMLE_C  <- SMLE_BETA(y=y_c, X=SampleG$X, Z=SampleG$Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=F)
 fitMDPDE_C<- MDPDE_BETA(y=y_c, X=SampleG$X, Z=SampleG$Z, qoptimal=T, q0=1, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP", linkmu="logit", linkphi="log", weights=F)

 #***************parameter estimates*************************#
 thetahat[cont,] <- c(fitMLE$beta, fitMLE$gama)
 thetahat_C[cont,] <- c(fitMLE_C$beta, fitMLE_C$gama)
 thetahat_smle[cont,] <- c(fitSMLE$beta, fitSMLE$gama)
 thetahat_smle_C[cont,] <- c(fitSMLE_C$beta, fitSMLE_C$gama)
 thetahat_mdpde[cont,] <- c(fitMDPDE$beta, fitMDPDE$gama)
 thetahat_mdpde_C[cont,] <- c(fitMDPDE_C$beta, fitMDPDE_C$gama)

  #***************asymptotic standard error estimates*************************#
 se_thetahat[cont,] <- c(fitMLE$se_beta, fitMLE$se_gama)
 se_thetahat_C[cont,] <- c(fitMLE_C$se_beta, fitMLE_C$se_gama)
 se_thetahat_smle[cont,] <- c(fitSMLE$se_beta, fitSMLE$se_gama)
 se_thetahat_smle_C[cont,] <- c(fitSMLE_C$se_beta, fitSMLE_C$se_gama)
 se_thetahat_mdpde[cont,] <- c(fitMDPDE$se_beta, fitMDPDE$se_gama)
 se_thetahat_mdpde_C[cont,] <- c(fitMDPDE_C$se_beta, fitMDPDE_C$se_gama)
 
  #***************z-statistics*************************#
 z_thetahat[cont, ] <- (thetahat[cont,] - VTHETA)/se_thetahat[cont,]
 z_thetahat_C[cont, ] <- (thetahat_C[cont,] - VTHETA)/se_thetahat_C[cont,]
 z_thetahat_smle[cont, ] <- (thetahat_smle[cont,] - VTHETA)/se_thetahat_smle[cont,]
 z_thetahat_smle_C[cont, ] <- (thetahat_smle_C[cont,] - VTHETA)/se_thetahat_smle_C[cont,]
 z_thetahat_mdpde[cont, ] <- (thetahat_mdpde[cont,] - VTHETA)/se_thetahat_mdpde[cont,]
 z_thetahat_mdpde_C[cont, ] <- (thetahat_mdpde_C[cont,] - VTHETA)/se_thetahat_mdpde_C[cont,]

  #***************trace of covariance matrix (sum of asymptotic variances)*************************#
 trace_thetahat[cont, ] <-  sum(se_thetahat[cont,]^2)
 trace_thetahat_C[cont, ] <-  sum(se_thetahat_C[cont,]^2)
 trace_thetahat_smle[cont, ] <-  sum(se_thetahat_smle[cont,]^2)
 trace_thetahat_smle_C[cont, ] <-  sum(se_thetahat_smle_C[cont,]^2) 
 trace_thetahat_mdpde[cont, ] <-  sum(se_thetahat_mdpde[cont,]^2)
 trace_thetahat_mdpde_C[cont, ] <-  sum(se_thetahat_mdpde_C[cont,]^2) 

#***************optimal q values*************************#
 qopt_thetahat_smle[cont, ] <- fitSMLE$q_value
 qopt_thetahat_smle_C[cont, ] <- fitSMLE_C$q_value
 qopt_thetahat_mdpde[cont, ] <- fitMDPDE$q_value
 qopt_thetahat_mdpde_C[cont, ] <- fitMDPDE_C$q_value


}#closes MC replicates

#******************results***************************#

 #**************mean of estimates*************#
M_thetahat <- apply(thetahat,2,mean)
M_thetahat_C <- apply(thetahat_C,2,mean)
M_thetahat_smle <- apply(thetahat_smle,2,mean)
M_thetahat_smle_C <- apply(thetahat_smle_C,2,mean)
M_thetahat_mdpde <- apply(thetahat_mdpde,2,mean)
M_thetahat_mdpde_C <- apply(thetahat_mdpde_C,2,mean)

#**************median of estimates*************#
Med_thetahat <- apply(thetahat,2,median)
Med_thetahat_C <- apply(thetahat_C,2,median)
Med_thetahat_smle <- apply(thetahat_smle,2,median)
Med_thetahat_smle_C <- apply(thetahat_smle_C,2,median)
Med_thetahat_mdpde <- apply(thetahat_mdpde,2,median)
Med_thetahat_mdpde_C <- apply(thetahat_mdpde_C,2,median)

#**************relative bias of estimates*************#
RB_thetahat <- (M_thetahat - VTHETA)/VTHETA
RB_thetahat_C <- (M_thetahat_C - VTHETA)/VTHETA
RB_thetahat_smle <- (M_thetahat_smle - VTHETA)/VTHETA
RB_thetahat_smle_C <- (M_thetahat_smle_C - VTHETA)/VTHETA
RB_thetahat_mdpde <- (M_thetahat_mdpde - VTHETA)/VTHETA
RB_thetahat_mdpde_C <- (M_thetahat_mdpde_C - VTHETA)/VTHETA

#**************relative median bias of estimates*************#
RMB_thetahat <- (Med_thetahat - VTHETA)/VTHETA
RMB_thetahat_C <- (Med_thetahat_C - VTHETA)/VTHETA
RMB_thetahat_smle <- (Med_thetahat_smle - VTHETA)/VTHETA
RMB_thetahat_smle_C <- (Med_thetahat_smle_C - VTHETA)/VTHETA
RMB_thetahat_mdpde <- (Med_thetahat_mdpde - VTHETA)/VTHETA
RMB_thetahat_mdpde_C <- (Med_thetahat_mdpde_C - VTHETA)/VTHETA

 #**************mean of s.e. estimates*************#
M_se_thetahat <- apply(se_thetahat,2,mean)
M_se_thetahat_C <- apply(se_thetahat_C,2,mean)
M_se_thetahat_smle <- apply(se_thetahat_smle,2,mean)
M_se_thetahat_smle_C <- apply(se_thetahat_smle_C,2,mean)
M_se_thetahat_mdpde <- apply(se_thetahat_mdpde,2,mean)
M_se_thetahat_mdpde_C <- apply(se_thetahat_mdpde_C,2,mean)

#**************median of s.e. estimates*************#
Med_se_thetahat <- apply(se_thetahat,2,median)
Med_se_thetahat_C <- apply(se_thetahat_C,2,median)
Med_se_thetahat_smle <- apply(se_thetahat_smle,2,median)
Med_se_thetahat_smle_C <- apply(se_thetahat_smle_C,2,median)
Med_se_thetahat_mdpde <- apply(se_thetahat_mdpde,2,median)
Med_se_thetahat_mdpde_C <- apply(se_thetahat_mdpde_C,2,median)

#**************ratio of variances of estimates*************#
RV_MLE_SMLE <- mean(trace_thetahat/trace_thetahat_smle)
RV_MLE_SMLE_C <- mean(trace_thetahat_C/trace_thetahat_smle_C)
RV_MLE_MDPDE <- mean(trace_thetahat/trace_thetahat_mdpde)
RV_MLE_MDPDE_C <- mean(trace_thetahat_C/trace_thetahat_mdpde_C)
RV_SMLE_MDPDE <- mean(trace_thetahat_smle/trace_thetahat_mdpde)
RV_SMLE_MDPDE_C <- mean(trace_thetahat_smle_C/trace_thetahat_mdpde_C)

#**************ratio of MSE of estimates*************#
RMSE_MLE_SMLE <-     sum(apply((t(thetahat)        -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_smle)    -  VTHETA)^2.0, 1,mean)) 
RMSE_MLE_SMLE_C <-   sum(apply((t(thetahat_C)      -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_smle_C)  -  VTHETA)^2.0, 1,mean))
RMSE_MLE_MDPDE <-    sum(apply((t(thetahat)        -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_mdpde)   -  VTHETA)^2.0, 1,mean))
RMSE_MLE_MDPDE_C <-  sum(apply((t(thetahat_C)      -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_mdpde_C) -  VTHETA)^2.0, 1,mean))
RMSE_SMLE_MDPDE <-   sum(apply((t(thetahat_smle)   -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_mdpde)   -  VTHETA)^2.0, 1,mean))
RMSE_SMLE_MDPDE_C <- sum(apply((t(thetahat_smle_C) -  VTHETA)^2.0, 1,mean))/sum(apply((t(thetahat_mdpde_C) -  VTHETA)^2.0, 1,mean))

#**************empirical null levels of Wald-type tests*************#
LEVEL_MLE <- apply(abs(z_thetahat)>1.959964, 2, mean)
LEVEL_MLE_C <- apply(abs(z_thetahat_C)>1.959964, 2, mean)
LEVEL_SMLE <- apply(abs(z_thetahat_smle)>1.959964, 2, mean)
LEVEL_SMLE_C <- apply(abs(z_thetahat_smle_C)>1.959964, 2, mean)		 
LEVEL_MDPDE <- apply(abs(z_thetahat_mdpde)>1.959964, 2, mean)
LEVEL_MDPDE_C <- apply(abs(z_thetahat_mdpde_C)>1.959964, 2, mean)		

#**********************mean of optimal q****************#
QOPTIMAL_SMLE <- mean(qopt_thetahat_smle)
QOPTIMAL_SMLE_C <- mean(qopt_thetahat_smle_C)
QOPTIMAL_MDPDE <- mean(qopt_thetahat_mdpde)
QOPTIMAL_MDPDE_C <- mean(qopt_thetahat_mdpde_C)

###############################################Outputs################################################

out_estimates_MLE <- cbind(M_thetahat, RB_thetahat, M_se_thetahat, LEVEL_MLE, Med_thetahat, RMB_thetahat, Med_se_thetahat)
out_estimates_MLE_C <- cbind(M_thetahat_C, RB_thetahat_C, M_se_thetahat_C, LEVEL_MLE_C, Med_thetahat_C, RMB_thetahat_C, Med_se_thetahat_C)
out_estimates_SMLE <- cbind(M_thetahat_smle, RB_thetahat_smle, M_se_thetahat_smle, LEVEL_SMLE, Med_thetahat_smle, RMB_thetahat_smle, Med_se_thetahat_smle)
out_estimates_SMLE_C <- cbind(M_thetahat_smle_C, RB_thetahat_smle_C, M_se_thetahat_smle_C, LEVEL_SMLE_C, Med_thetahat_smle_C, RMB_thetahat_smle_C, Med_se_thetahat_smle_C)
out_estimates_MDPDE <- cbind(M_thetahat_mdpde, RB_thetahat_mdpde, M_se_thetahat_mdpde, LEVEL_MDPDE, Med_thetahat_mdpde, RMB_thetahat_mdpde, Med_se_thetahat_mdpde)
out_estimates_MDPDE_C <- cbind(M_thetahat_mdpde_C, RB_thetahat_mdpde_C, M_se_thetahat_mdpde_C, LEVEL_MDPDE_C, Med_thetahat_mdpde_C, RMB_thetahat_mdpde_C, Med_se_thetahat_mdpde_C)
out_ratio_RV <- rbind(cbind(RV_MLE_SMLE, RV_MLE_MDPDE, RV_SMLE_MDPDE), cbind(RV_MLE_SMLE_C, RV_MLE_MDPDE_C, RV_SMLE_MDPDE_C))
out_ratio_MSE <- rbind(cbind(RMSE_MLE_SMLE, RMSE_MLE_MDPDE, RMSE_SMLE_MDPDE), cbind(RMSE_MLE_SMLE_C, RMSE_MLE_MDPDE_C, RMSE_SMLE_MDPDE_C))
out_PMLE <- rbind(cbind(mean(qopt_thetahat_smle==1.0), mean(qopt_thetahat_mdpde==1.0)), cbind(mean(qopt_thetahat_smle_C==1.0), mean(qopt_thetahat_mdpde_C==1.0)))

cat(" ========================================================================================== \n")
cat(" TRUE VALUES ","\n")
cat(" Sample size ",N ,"\n")
out_VTHETA <- rbind(VTHETA);colnames(out_VTHETA)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");rownames(out_VTHETA)=c("");print(out_VTHETA)
Summ_mu <- round(c(min(mu), max(mu), median(mu), min(phi), max(phi), median(phi), max(phi)/min(phi)),3)
out_Smu <- rbind(Summ_mu);colnames(out_Smu)=c("Mu_min", "Mu_max", "Mu_median", "Phi_min", "Phi_max", "Phi_median", "Het_Intensity");rownames(out_Smu)=c("");print(out_Smu)
cat(" ========================================================================================== \n")
cat(" ================================ MLE estimates without contamination ================================ \n")
output_mle = rbind(out_estimates_MLE);colnames(output_mle)=c("Mean estimates","Bias","Mean SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_mle)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_mle)
cat(" ================================ MLE estimates with contamination ================================ \n")
output_mle_c = rbind(out_estimates_MLE_C);colnames(output_mle_c)=c("Estimates","Bias","SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_mle_c)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_mle_c)
cat(" ================================ SMLE estimates without contamination ================================\n")
cat(" Mean of q-optimal is ",QOPTIMAL_SMLE ,"\n")
output_smle = rbind(out_estimates_SMLE);colnames(output_smle)=c("Estimates","Bias","SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_smle)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_smle)
cat(" ================================ SMLE estimates with contamination ================================ \n")
cat(" Mean of q-optimal is ",QOPTIMAL_SMLE_C ,"\n")
output_smle_c = rbind(out_estimates_SMLE_C);colnames(output_smle_c)=c("Estimates","Bias","SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_smle_c)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_smle_c)
cat(" ================================ MDPDE estimates without contamination ================================ \n")
cat(" Mean of q-optimal is ",QOPTIMAL_MDPDE ,"\n")
output_mdpde = rbind(out_estimates_MDPDE);colnames(output_mdpde)=c("Estimates","Bias","SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_mdpde)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_mdpde)
cat(" ================================ MDPDE estimates with contamination ================================ \n")
cat(" Mean of q-optimal is ",QOPTIMAL_MDPDE_C ,"\n")
output_mdpde_c = rbind(out_estimates_MDPDE_C);colnames(output_mdpde_c)=c("Estimates","Bias","SE", "Null level", "Median estimates","Median Bias","Median SE");rownames(output_mdpde_c)=c("beta1", "beta2", "beta3", "gama1", "gama2", "gama3");print(output_mdpde_c)
cat(" ================================ Ratio of MSE ================================ \n")
output_mse = rbind(out_ratio_MSE);colnames(output_mse)=c("MLE/SMLE","MLE/MDPDE","SMLE/MDPDE");rownames(output_mse)=c("Without contamination", "With contamination");print(output_mse)
cat(" ================================ Ratio of sum of variances ================================ \n")
output_rv = rbind(out_ratio_RV);colnames(output_rv)=c("MLE/SMLE","MLE/MDPDE","SMLE/MDPDE");rownames(output_rv)=c("Without contamination", "With contamination");print(output_rv)
cat(" ========================================================================================== \n")
cat(" ================================ Percentage of times that MLE is choosen ================================ \n")
output_Pmle = rbind(out_PMLE);colnames(output_Pmle)=c("SMLE","MDPDE");rownames(output_Pmle)=c("Without contamination", "With contamination");print(output_Pmle)
cat(" ========================================================================================== \n")

#******Print replicates*****#

write.table(cbind(thetahat, thetahat_C, thetahat_smle, thetahat_smle_C, thetahat_mdpde, thetahat_mdpde_C), 
file= paste(names_results1[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
col.names = c("MLE_beta1", "MLE_beta2", "MLE_beta3", "MLE_gama1", "MLE_gama2", "MLE_gama3", "MLE_C_beta1", "MLE_C_beta2", "MLE_C_beta3", 
"MLE_C_gama1", "MLE_C_gama2", "MLE_C_gama3","SMLE_beta1", "SMLE_beta2", "SMLE_beta3", "SMLE_gama1", "SMLE_gama2", "SMLE_gama3",
"SMLE_C_beta1", "SMLE_C_beta2", "SMLE_C_beta3", "SMLE_C_gama1", "SMLE_C_gama2", "SMLE_C_gama3", "MDPDE_beta1", "MDPDE_beta2", "MDPDE_beta3",
"MDPDE_gama1", "MDPDE_gama2", "MDPDE_gama3", "MDPDE_C_beta1", "MDPDE_C_beta2", "MDPDE_C_beta3", "MDPDE_C_gama1", "MDPDE_C_gama2", "MDPDE_C_gama3"))

write.table(cbind(se_thetahat, se_thetahat_C, se_thetahat_smle, se_thetahat_smle_C, se_thetahat_mdpde, se_thetahat_mdpde_C), 
file= paste(names_results2[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
col.names = c("se_MLE_beta1", "se_MLE_beta2", "se_MLE_beta3", "se_MLE_gama1", "se_MLE_gama2", "se_MLE_gama3", "se_MLE_C_beta1", "se_MLE_C_beta2", "se_MLE_C_beta3", 
"se_MLE_C_gama1", "se_MLE_C_gama2", "se_MLE_C_gama3","se_SMLE_beta1", "se_SMLE_beta2", "se_SMLE_beta3", "se_SMLE_gama1", "se_SMLE_gama2", "se_SMLE_gama3",
"se_SMLE_C_beta1", "se_SMLE_C_beta2", "se_SMLE_C_beta3", "se_SMLE_C_gama1", "se_SMLE_C_gama2", "se_SMLE_C_gama3", "se_MDPDE_beta1", "se_MDPDE_beta2", "se_MDPDE_beta3",
"se_MDPDE_gama1", "se_MDPDE_gama2", "se_MDPDE_gama3", "se_MDPDE_C_beta1", "se_MDPDE_C_beta2", "se_MDPDE_C_beta3", "se_MDPDE_C_gama1", "se_MDPDE_C_gama2", "se_MDPDE_C_gama3"))

write.table(cbind(qopt_thetahat_smle, qopt_thetahat_smle_C, qopt_thetahat_mdpde, qopt_thetahat_mdpde_C), 
file= paste(names_results3[l], ".txt", sep = ""), append = FALSE, sep = " ", dec = ".", row.names = FALSE,
col.names = c("qoptimal_SMLE", "qoptimal_SMLE_C", "qoptimal_MDPDE", "qoptimal_MDPDE_C"))

}

sink()
