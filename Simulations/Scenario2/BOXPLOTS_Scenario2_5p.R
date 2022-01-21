rm(list=ls())

data40 <- read.table("Estimates_Scenario2_n40.txt", h=T)
data80 <- read.table("Estimates_Scenario2_n80.txt", h=T)
data160 <- read.table("Estimates_Scenario2_n160.txt", h=T)
data320 <- read.table("Estimates_Scenario2_n320.txt", h=T)

#**************************MLE estimates******************************#
MLE_beta1_n40 <- data40$MLE_beta1; MLEC_beta1_n40 <- data40$MLE_C_beta1
MLE_beta2_n40 <- data40$MLE_beta2; MLEC_beta2_n40 <- data40$MLE_C_beta2
MLE_beta3_n40 <- data40$MLE_beta3; MLEC_beta3_n40 <- data40$MLE_C_beta3
MLE_gama1_n40 <- data40$MLE_gama1; MLEC_gama1_n40 <- data40$MLE_C_gama1
MLE_gama2_n40 <- data40$MLE_gama2; MLEC_gama2_n40 <- data40$MLE_C_gama2
MLE_gama3_n40 <- data40$MLE_gama3; MLEC_gama3_n40 <- data40$MLE_C_gama3

MLE_beta1_n80 <- data80$MLE_beta1; MLEC_beta1_n80 <- data80$MLE_C_beta1
MLE_beta2_n80 <- data80$MLE_beta2; MLEC_beta2_n80 <- data80$MLE_C_beta2
MLE_beta3_n80 <- data80$MLE_beta3; MLEC_beta3_n80 <- data80$MLE_C_beta3
MLE_gama1_n80 <- data80$MLE_gama1; MLEC_gama1_n80 <- data80$MLE_C_gama1
MLE_gama2_n80 <- data80$MLE_gama2; MLEC_gama2_n80 <- data80$MLE_C_gama2
MLE_gama3_n80 <- data80$MLE_gama3; MLEC_gama3_n80 <- data80$MLE_C_gama3

MLE_beta1_n160 <- data160$MLE_beta1; MLEC_beta1_n160 <- data160$MLE_C_beta1
MLE_beta2_n160 <- data160$MLE_beta2; MLEC_beta2_n160 <- data160$MLE_C_beta2
MLE_beta3_n160 <- data160$MLE_beta3; MLEC_beta3_n160 <- data160$MLE_C_beta3
MLE_gama1_n160 <- data160$MLE_gama1; MLEC_gama1_n160 <- data160$MLE_C_gama1
MLE_gama2_n160 <- data160$MLE_gama2; MLEC_gama2_n160 <- data160$MLE_C_gama2
MLE_gama3_n160 <- data160$MLE_gama3; MLEC_gama3_n160 <- data160$MLE_C_gama3

MLE_beta1_n320 <- data320$MLE_beta1; MLEC_beta1_n320 <- data320$MLE_C_beta1
MLE_beta2_n320 <- data320$MLE_beta2; MLEC_beta2_n320 <- data320$MLE_C_beta2
MLE_beta3_n320 <- data320$MLE_beta3; MLEC_beta3_n320 <- data320$MLE_C_beta3
MLE_gama1_n320 <- data320$MLE_gama1; MLEC_gama1_n320 <- data320$MLE_C_gama1
MLE_gama2_n320 <- data320$MLE_gama2; MLEC_gama2_n320 <- data320$MLE_C_gama2
MLE_gama3_n320 <- data320$MLE_gama3; MLEC_gama3_n320 <- data320$MLE_C_gama3


#**************************SMLE estimates******************************#
SMLE_beta1_n40 <- data40$SMLE_beta1; SMLEC_beta1_n40 <- data40$SMLE_C_beta1
SMLE_beta2_n40 <- data40$SMLE_beta2; SMLEC_beta2_n40 <- data40$SMLE_C_beta2
SMLE_beta3_n40 <- data40$SMLE_beta3; SMLEC_beta3_n40 <- data40$SMLE_C_beta3
SMLE_gama1_n40 <- data40$SMLE_gama1; SMLEC_gama1_n40 <- data40$SMLE_C_gama1
SMLE_gama2_n40 <- data40$SMLE_gama2; SMLEC_gama2_n40 <- data40$SMLE_C_gama2
SMLE_gama3_n40 <- data40$SMLE_gama3; SMLEC_gama3_n40 <- data40$SMLE_C_gama3

SMLE_beta1_n80 <- data80$SMLE_beta1; SMLEC_beta1_n80 <- data80$SMLE_C_beta1
SMLE_beta2_n80 <- data80$SMLE_beta2; SMLEC_beta2_n80 <- data80$SMLE_C_beta2
SMLE_beta3_n80 <- data80$SMLE_beta3; SMLEC_beta3_n80 <- data80$SMLE_C_beta3
SMLE_gama1_n80 <- data80$SMLE_gama1; SMLEC_gama1_n80 <- data80$SMLE_C_gama1
SMLE_gama2_n80 <- data80$SMLE_gama2; SMLEC_gama2_n80 <- data80$SMLE_C_gama2
SMLE_gama3_n80 <- data80$SMLE_gama3; SMLEC_gama3_n80 <- data80$SMLE_C_gama3

SMLE_beta1_n160 <- data160$SMLE_beta1; SMLEC_beta1_n160 <- data160$SMLE_C_beta1
SMLE_beta2_n160 <- data160$SMLE_beta2; SMLEC_beta2_n160 <- data160$SMLE_C_beta2
SMLE_beta3_n160 <- data160$SMLE_beta3; SMLEC_beta3_n160 <- data160$SMLE_C_beta3
SMLE_gama1_n160 <- data160$SMLE_gama1; SMLEC_gama1_n160 <- data160$SMLE_C_gama1
SMLE_gama2_n160 <- data160$SMLE_gama2; SMLEC_gama2_n160 <- data160$SMLE_C_gama2
SMLE_gama3_n160 <- data160$SMLE_gama3; SMLEC_gama3_n160 <- data160$SMLE_C_gama3

SMLE_beta1_n320 <- data320$SMLE_beta1; SMLEC_beta1_n320 <- data320$SMLE_C_beta1
SMLE_beta2_n320 <- data320$SMLE_beta2; SMLEC_beta2_n320 <- data320$SMLE_C_beta2
SMLE_beta3_n320 <- data320$SMLE_beta3; SMLEC_beta3_n320 <- data320$SMLE_C_beta3
SMLE_gama1_n320 <- data320$SMLE_gama1; SMLEC_gama1_n320 <- data320$SMLE_C_gama1
SMLE_gama2_n320 <- data320$SMLE_gama2; SMLEC_gama2_n320 <- data320$SMLE_C_gama2
SMLE_gama3_n320 <- data320$SMLE_gama3; SMLEC_gama3_n320 <- data320$SMLE_C_gama3

#**************************MDPDE estimates******************************#

MDPDE_beta1_n40 <- data40$MDPDE_beta1; MDPDEC_beta1_n40 <- data40$MDPDE_C_beta1
MDPDE_beta2_n40 <- data40$MDPDE_beta2; MDPDEC_beta2_n40 <- data40$MDPDE_C_beta2
MDPDE_beta3_n40 <- data40$MDPDE_beta3; MDPDEC_beta3_n40 <- data40$MDPDE_C_beta3
MDPDE_gama1_n40 <- data40$MDPDE_gama1; MDPDEC_gama1_n40 <- data40$MDPDE_C_gama1
MDPDE_gama2_n40 <- data40$MDPDE_gama2; MDPDEC_gama2_n40 <- data40$MDPDE_C_gama2
MDPDE_gama3_n40 <- data40$MDPDE_gama3; MDPDEC_gama3_n40 <- data40$MDPDE_C_gama3

MDPDE_beta1_n80 <- data80$MDPDE_beta1; MDPDEC_beta1_n80 <- data80$MDPDE_C_beta1
MDPDE_beta2_n80 <- data80$MDPDE_beta2; MDPDEC_beta2_n80 <- data80$MDPDE_C_beta2
MDPDE_beta3_n80 <- data80$MDPDE_beta3; MDPDEC_beta3_n80 <- data80$MDPDE_C_beta3
MDPDE_gama1_n80 <- data80$MDPDE_gama1; MDPDEC_gama1_n80 <- data80$MDPDE_C_gama1
MDPDE_gama2_n80 <- data80$MDPDE_gama2; MDPDEC_gama2_n80 <- data80$MDPDE_C_gama2
MDPDE_gama3_n80 <- data80$MDPDE_gama3; MDPDEC_gama3_n80 <- data80$MDPDE_C_gama3

MDPDE_beta1_n160 <- data160$MDPDE_beta1; MDPDEC_beta1_n160 <- data160$MDPDE_C_beta1
MDPDE_beta2_n160 <- data160$MDPDE_beta2; MDPDEC_beta2_n160 <- data160$MDPDE_C_beta2
MDPDE_beta3_n160 <- data160$MDPDE_beta3; MDPDEC_beta3_n160 <- data160$MDPDE_C_beta3
MDPDE_gama1_n160 <- data160$MDPDE_gama1; MDPDEC_gama1_n160 <- data160$MDPDE_C_gama1
MDPDE_gama2_n160 <- data160$MDPDE_gama2; MDPDEC_gama2_n160 <- data160$MDPDE_C_gama2
MDPDE_gama3_n160 <- data160$MDPDE_gama3; MDPDEC_gama3_n160 <- data160$MDPDE_C_gama3

MDPDE_beta1_n320 <- data320$MDPDE_beta1; MDPDEC_beta1_n320 <- data320$MDPDE_C_beta1
MDPDE_beta2_n320 <- data320$MDPDE_beta2; MDPDEC_beta2_n320 <- data320$MDPDE_C_beta2
MDPDE_beta3_n320 <- data320$MDPDE_beta3; MDPDEC_beta3_n320 <- data320$MDPDE_C_beta3
MDPDE_gama1_n320 <- data320$MDPDE_gama1; MDPDEC_gama1_n320 <- data320$MDPDE_C_gama1
MDPDE_gama2_n320 <- data320$MDPDE_gama2; MDPDEC_gama2_n320 <- data320$MDPDE_C_gama2
MDPDE_gama3_n320 <- data320$MDPDE_gama3; MDPDEC_gama3_n320 <- data320$MDPDE_C_gama3

#********************MLE estimates of beta1****************#
boxplot(MLE_beta1_n40,MLE_beta1_n80,MLE_beta1_n160,MLE_beta1_n320,
MLEC_beta1_n40,MLEC_beta1_n80,MLEC_beta1_n160,MLEC_beta1_n320,
main=expression( paste("MLE for ", beta[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-0.5,1.4), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-0.35,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-0.35,"Presence of contamination",bty="n",cex=1.3)

#********************MLE estimates of beta2****************#
boxplot(MLE_beta2_n40,MLE_beta2_n80,MLE_beta2_n160,MLE_beta2_n320,
MLEC_beta2_n40,MLEC_beta2_n80,MLEC_beta2_n160,MLEC_beta2_n320,
main=expression( paste("MLE for ", beta[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.2,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-2.0,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-2.0,"Presence of contamination",bty="n",cex=1.3)

#********************MLE estimates of beta3****************#
boxplot(MLE_beta3_n40,MLE_beta3_n80,MLE_beta3_n160,MLE_beta3_n320,
MLEC_beta3_n40,MLEC_beta3_n80,MLEC_beta3_n160,MLEC_beta3_n320,
main=expression( paste("MLE for ", beta[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-1.95,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.8,"Presence of contamination",bty="n",cex=1.3)

#********************MLE estimates of gamma1****************#
boxplot(MLE_gama1_n40,MLE_gama1_n80,MLE_gama1_n160,MLE_gama1_n320,
MLEC_gama1_n40,MLEC_gama1_n80,MLEC_gama1_n160,MLEC_gama1_n320,
main=expression( paste("MLE for ", gamma[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.5,10), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(3.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.3,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.3,"Presence of contamination",bty="n",cex=1.3)

#********************MLE estimates of gamma2****************#
boxplot(MLE_gama2_n40,MLE_gama2_n80,MLE_gama2_n160,MLE_gama2_n320,
MLEC_gama2_n40,MLEC_gama2_n80,MLEC_gama2_n160,MLEC_gama2_n320,
main=expression( paste("MLE for ", gamma[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-6.85,6.6), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.7,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.7,"Presence of contamination",bty="n",cex=1.3)

#********************MLE estimates of gamma3****************#
boxplot(MLE_gama3_n40,MLE_gama3_n80,MLE_gama3_n160,MLE_gama3_n320,
MLEC_gama3_n40,MLEC_gama3_n80,MLEC_gama3_n160,MLEC_gama3_n320,
main=expression( paste("MLE for ", gamma[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-7.0,7.3), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.8,"Presence of contamination",bty="n",cex=1.3)

#**************************SMLE estimates******************************#

#********************SMLE estimates of beta1****************#
boxplot(SMLE_beta1_n40,SMLE_beta1_n80,SMLE_beta1_n160,SMLE_beta1_n320,
SMLEC_beta1_n40,SMLEC_beta1_n80,SMLEC_beta1_n160,SMLEC_beta1_n320,
main=expression( paste("SMLE for ", beta[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-0.5,1.4), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-0.35,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-0.35,"Presence of contamination",bty="n",cex=1.3)

#********************SMLE estimates of beta2****************#
boxplot(SMLE_beta2_n40,SMLE_beta2_n80,SMLE_beta2_n160,SMLE_beta2_n320,
SMLEC_beta2_n40,SMLEC_beta2_n80,SMLEC_beta2_n160,SMLEC_beta2_n320,
main=expression( paste("SMLE for ", beta[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.2,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-2.0,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-2.0,"Presence of contamination",bty="n",cex=1.3)

#********************SMLE estimates of beta3****************#
boxplot(SMLE_beta3_n40,SMLE_beta3_n80,SMLE_beta3_n160,SMLE_beta3_n320,
SMLEC_beta3_n40,SMLEC_beta3_n80,SMLEC_beta3_n160,SMLEC_beta3_n320,
main=expression( paste("SMLE for ", beta[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-1.95,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.8,"Presence of contamination",bty="n",cex=1.3)

#********************SMLE estimates of gamma1****************#
boxplot(SMLE_gama1_n40,SMLE_gama1_n80,SMLE_gama1_n160,SMLE_gama1_n320,
SMLEC_gama1_n40,SMLEC_gama1_n80,SMLEC_gama1_n160,SMLEC_gama1_n320,
main=expression( paste("SMLE for ", gamma[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.5,10), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(3.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.3,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.3,"Presence of contamination",bty="n",cex=1.3)

#********************SMLE estimates of gamma2****************#
boxplot(SMLE_gama2_n40,SMLE_gama2_n80,SMLE_gama2_n160,SMLE_gama2_n320,
SMLEC_gama2_n40,SMLEC_gama2_n80,SMLEC_gama2_n160,SMLEC_gama2_n320,
main=expression( paste("SMLE for ", gamma[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-6.85,6.6), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.7,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.7,"Presence of contamination",bty="n",cex=1.3)

#********************SMLE estimates of gamma3****************#
boxplot(SMLE_gama3_n40,SMLE_gama3_n80,SMLE_gama3_n160,SMLE_gama3_n320,
SMLEC_gama3_n40,SMLEC_gama3_n80,SMLEC_gama3_n160,SMLEC_gama3_n320,
main=expression( paste("SMLE for ", gamma[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-7.0,7.3), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.8,"Presence of contamination",bty="n",cex=1.3)

#**************************MDPDE estimates******************************#

#********************MDPDE estimates of beta1****************#
boxplot(MDPDE_beta1_n40,MDPDE_beta1_n80,MDPDE_beta1_n160,MDPDE_beta1_n320,
MDPDEC_beta1_n40,MDPDEC_beta1_n80,MDPDEC_beta1_n160,MDPDEC_beta1_n320,
main=expression( paste("MDPDE for ", beta[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-0.5,1.4), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-0.35,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-0.35,"Presence of contamination",bty="n",cex=1.3)

#********************MDPDE estimates of beta2****************#
boxplot(MDPDE_beta2_n40,MDPDE_beta2_n80,MDPDE_beta2_n160,MDPDE_beta2_n320,
MDPDEC_beta2_n40,MDPDEC_beta2_n80,MDPDEC_beta2_n160,MDPDEC_beta2_n320,
main=expression( paste("MDPDE for ", beta[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.2,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-2.0,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-2.0,"Presence of contamination",bty="n",cex=1.3)

#********************MDPDE estimates of beta3****************#
boxplot(MDPDE_beta3_n40,MDPDE_beta3_n80,MDPDE_beta3_n160,MDPDE_beta3_n320,
MDPDEC_beta3_n40,MDPDEC_beta3_n80,MDPDEC_beta3_n160,MDPDEC_beta3_n320,
main=expression( paste("MDPDE for ", beta[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-1.95,0.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(-1.2,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.8,"Presence of contamination",bty="n",cex=1.3)

#********************MDPDE estimates of gamma1****************#
boxplot(MDPDE_gama1_n40,MDPDE_gama1_n80,MDPDE_gama1_n160,MDPDE_gama1_n320,
MDPDEC_gama1_n40,MDPDEC_gama1_n80,MDPDEC_gama1_n160,MDPDEC_gama1_n320,
main=expression( paste("MDPDE for ", gamma[1])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-2.5,10), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(3.8,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-1.3,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-1.3,"Presence of contamination",bty="n",cex=1.3)

#********************MDPDE estimates of gamma2****************#
boxplot(MDPDE_gama2_n40,MDPDE_gama2_n80,MDPDE_gama2_n160,MDPDE_gama2_n320,
MDPDEC_gama2_n40,MDPDEC_gama2_n80,MDPDEC_gama2_n160,MDPDEC_gama2_n320,
main=expression( paste("MDPDE for ", gamma[2])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-6.85,6.6), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.7,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.7,"Presence of contamination",bty="n",cex=1.3)

#********************MDPDE estimates of gamma3****************#
boxplot(MDPDE_gama3_n40,MDPDE_gama3_n80,MDPDE_gama3_n160,MDPDE_gama3_n320,
MDPDEC_gama3_n40,MDPDEC_gama3_n80,MDPDEC_gama3_n160,MDPDEC_gama3_n320,
main=expression( paste("MDPDE for ", gamma[3])),
xlab="Sample size",ylab="",names=c(40,80,160,320,40,80,160,320),outline=T,
ylim=c(-7.0,7.3), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col=NULL)
abline(0.7,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,-5.8,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,-5.8,"Presence of contamination",bty="n",cex=1.3)










