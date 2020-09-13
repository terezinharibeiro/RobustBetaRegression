#***********************Optimal tuning values***************#
rm(list=ls())

data40 <- read.table("Qoptimal_Scenario1_n40.txt", h=T)
data80 <- read.table("Qoptimal_Scenario1_n80.txt", h=T)
data160 <- read.table("Qoptimal_Scenario1_n40.txt", h=T)
data320 <- read.table("Qoptimal_Scenario1_n40.txt", h=T)

SMLE_n40 <-  data40$qoptimal_SMLE;     SMLEC_n40 <-   data40$qoptimal_SMLE_C
MDPDE_n40 <- data40$qoptimal_MDPDE;    MDPDEC_n40 <-  data40$qoptimal_MDPDE_C
SMLE_n80 <-  data80$qoptimal_SMLE;     SMLEC_n80 <-   data80$qoptimal_SMLE_C
MDPDE_n80 <- data80$qoptimal_MDPDE;    MDPDEC_n80 <-  data80$qoptimal_MDPDE_C
SMLE_n160 <- data160$qoptimal_SMLE;    SMLEC_n160 <-  data160$qoptimal_SMLE_C
MDPDE_n160 <-data160$qoptimal_MDPDE;   MDPDEC_n160 <- data160$qoptimal_MDPDE_C
SMLE_n320 <- data320$qoptimal_SMLE;    SMLEC_n320 <-  data320$qoptimal_SMLE_C
MDPDE_n320 <-data320$qoptimal_MDPDE;   MDPDEC_n320 <- data320$qoptimal_MDPDE_C

#==================================================================
#==================================================================

par(mfrow=c(1,1))
boxplot(SMLE_n40,SMLE_n80,SMLE_n160,SMLE_n320,SMLEC_n40,SMLEC_n80,SMLEC_n160,SMLEC_n320
,main=expression( paste("Optimal ", q, " for SMLE")),
xlab="Sample size",ylab="",names=c(40,80, 160, 320,40,80, 160, 320),outline=T,
ylim=c(0.5,1.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col = NULL)
#abline(1.0,0,col=2,lty=2,lwd=2)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,0.54,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,0.54,"Presence of contamination",bty="n",cex=1.3)

boxplot(MDPDE_n40,MDPDE_n80,MDPDE_n160,MDPDE_n320,MDPDEC_n40,MDPDEC_n80,MDPDEC_n160,MDPDEC_n320,
main=expression( paste("Optimal ", q, " for MDPDE")),
xlab="Sample size",ylab="",names=c(40,80, 160, 320,40,80, 160, 320),outline=T,
ylim=c(0.5,1.0), cex=1.5,cex.lab=2.0,cex.axis=1.5,cex.main=2.0, pch=16, col = NULL)
abline(v=4.5,lty=3,lwd=2)
legend(0.0,0.54,"Absence of contamination",bty="n",cex=1.3)
legend(4.2,0.54,"Presence of contamination",bty="n",cex=1.3)







