# R codes and datasets for the paper "Robust estimation in beta regression via maximum Lq-likelihood"

This repository contains R codes and datasets used in the applications and simulations of the paper "Robust estimation in beta regression via maximum Lq-likelihood" by Ribeiro and Ferrari (2020).

### Applications

The directory Applications contains the three datasets and the R scripts to replicate the results presented in the Real data applications section of the paper.

- AIS_FIT.R: R script to replicate the robust inference results for the AIS dataset. 
- FIRMCOST_FIT.R: R script for replicate the robust inference results for the Firm cost dataset.
- data_RiskSurvey.txt: File with Firm cost dataset.
- TUNA_FIT.R: R script for replicate the robust inference results for the Tuna dataset.
- data_2000_Indian.txt: File with Tuna dataset.
- envelope_function.R: R function that creates normal probability plots of residuals with simulated envelope under SMLE and MLE.
- Resfunction.R: R function that calculates the residuals under MLE and SMLE.
- pvalue_boot.R: R function that calculates the bootstrap p-values under MLE and SMLE.
- SMLE.R: R function that fits beta regression models via SMLE and MLE.

### Simulations

The directory Simulations contains five folders corresponding to the simulation scenarios presented in the Monte Carlo simulation results section of the paper and the supplementary material. For example, folder Scenario1 contains the following files: 

- Simulation_Scenario1.R: R script to generate the results and files.
- Estimates_Scenario1_n40.txt: File with 1000 estimates for each parameter with n = 40 (analogously with n= 80, 160, 320).
- SE_Estimates_Scenario1_n40.txt: File with 1000 standard-error estimates for each parameter with n = 40 (analogously with n= 80, 160, 320).
- Qoptimal_Scenario1_n40.txt: File with 1000 q-optimal values for the 1000 samples with n = 40 (analogously with n= 80, 160, 320).
- SMLE.R: R function that fits beta regression models via SMLE and MLE.
- MDPDE.R: R function that fits beta regression models via MDPDE and MLE.
- SAMPLE.R: R function that replicates X and Z of the smallest sample and generate the response variable y.
- BOXPLOTS_ESTIMATES.R: R function to obtain the boxplots of the parameter estimates.
- BOXPLOTS_TUNING.R: R function to obtain the boxplots of the optimal values of q.
  
The general results are displayed in the R console.
