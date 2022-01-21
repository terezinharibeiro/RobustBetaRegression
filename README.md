# R codes and datasets for the paper "Robust estimation in beta regression via maximum Lq-likelihood"

This repository contains R codes and datasets used in the applications and simulations of the paper "Robust estimation in beta regression via maximum Lq-likelihood" by Ribeiro and Ferrari (2022).

### Applications

The directory Applications contains the three datasets and the R scripts to replicate the results presented in the Real data applications section of the paper.

- AIS_FIT.R: R script to replicate the robust inference results for the AIS dataset using the SMLE. 
- AIS_FIT_mdpde.R: R script to replicate the robust inference results for the AIS dataset using the MDPDE. 
- FIRMCOST_FIT.R: R script to replicate the robust inference results for the Firm cost dataset using the SMLE.
- FIRMCOST_FIT_mdpde.R: R script to replicate the robust inference results for the Firm cost dataset using the MDPDE.
- MDPDE.R: R function that fits beta regression models via MDPDE and MLE.
- Resfunction.R: R function that computes the residuals under MLE and SMLE.
- SMLE.R: R function that fits beta regression models via SMLE and MLE.
- SMLE_old.R: Old version of the R function that fits beta regression models via SMLE and MLE.
- TUNA_FIT.R: R script to replicate the robust inference results for the Tuna dataset using the SMLE.
- TUNA_FIT_mdpde.R: R script to replicate the robust inference results for the Tuna dataset using the MDPDE.
- data_2000_Indian.txt: Tuna dataset.
- data_RiskSurvey.txt: Firm cost dataset.
- envelope_function.R: R function that creates normal probability plots of residuals with a simulated envelope under SMLE and MLE.
- pvalue_boot.R: R function that computes the bootstrap p-values under MLE and SMLE.
- pvalue_boot_MDPDE.R: R function that computes the bootstrap p-values under MLE and MDPDE.

### Simulations

The directory Simulations contains five folders corresponding to the simulation scenarios presented in the paper's Monte Carlo simulation results section and the supplementary material. For example, folder Scenario1 contains the following files: 

- BOXPLOTS_ESTIMATES_10p.R: R function to obtain the boxplots of the parameter estimates in the absence of contamination and the presence of 10% of contamination.
- BOXPLOTS_ESTIMATES_2_5p.R: R function to obtain the boxplots of the parameter estimates in the absence of contamination and the presence of 2.5% of contamination.
- BOXPLOTS_ESTIMATES_5p.R: R function to obtain the boxplots of the parameter estimates in the absence of contamination and the presence of 5% of contamination.
- BOXPLOTS_ESTIMATES_old.R: Old version of the R function to obtain the boxplots of the parameter estimates in the absence of contamination and the presence of 5% of contamination.
- BOXPLOTS_TUNING.R: R function to obtain the boxplots of the optimal values of q.
- MDPDE.R: R function that fits beta regression models via MDPDE and MLE.
- MDPDE_old.R: Old version of the R function that fits beta regression models via MDPDE and MLE.
- Results_Scenario1_10p.txt: File with general simulation results in the absence of contamination and the presence of 10% of contamination.
- Results_Scenario1_2_5p.txt: File with general simulation results in the absence of contamination and the presence of 2.5% of contamination.
- Results_Scenario1_5p.txt: File with general simulation results in the absence of contamination and the presence of 5% of contamination.
- Results_Scenario1_old.txt: Old version of the file with general simulation results in the absence of contamination and the presence of 5% of contamination.
- SAMPLE.R: R function that replicates X and Z of the smallest sample and generates the response variable y.
- SMLE.R: R function that fits beta regression models via SMLE and MLE.
- SMLE_old.R: Old version of the R function that fits beta regression models via SMLE and MLE.
- Simulation_Scenario1_10p.R: R script to generate the results and files in the absence of contamination and the presence of 10% of contamination.
- Simulation_Scenario1_2_5p.R: R script to generate the results and files in the absence of contamination and the presence of 2.5% of contamination.
- Simulation_Scenario1_5p.R: R script to generate the results and files in the absence of contamination and the presence of 5% of contamination.
