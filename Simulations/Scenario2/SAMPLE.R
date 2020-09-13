#*************************************************************************Description*************************************************************************************************#
#Replicate X and Z of the smallest sample and generate the response variable y.
#***ARGUMENTS***#
# n - sample size. Must be 1, 2, 3, 4, 5, 6, 7, 8, 12, 16 or 32 times the smallest sample size.
# X -  regressor matrix for the mean submodel with number of rows equal to the smallest sample size.
# Z -  regressor matrix for the precision submodel with number of rows equal to the smallest sample size.
# size - smallest sample size.
# theta - vector of true parameters.
# linkmu - character specification of the link function in the mean submodel. Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Default is "logit".
# linkphi - character specification of the link function in the precision submodel. Currently, "identity", "log", "sqrt" are supported. Default is "log".
#**************************************************************************************************************************************************************************************#

FunSample <- function(n, mXini, mZini, size, theta, linkmu, linkphi){
              if(n == size)
	        {     
		mX <- mXini
		mZ <- mZini
		}
		else if (n == size*2)
		  {
		  mX <- rbind(mXini, mXini)
		  mZ <- rbind(mZini, mZini)
		  }
		  else if (n == size*3)
		  {
		  mX <- rbind(mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini)
		  }
		  else if (n == size*4)
		  {
                  mX <- rbind(mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini)
		  }	
		   else if (n == size*5)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini)
		  }	
		   else if (n == size*6)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini)
		  }
		   else if (n == size*7)
		  {
                  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini)
		  }
		   else if (n == size*8)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
		  }
		  else if (n == size*12)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
		  }
		  else if (n == size*16)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
		  }
		  else if (n == size*32)
		  {
		  mX <- rbind(mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini, mXini)
		  mZ <- rbind(mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini, mZini)
                   }	

 X <- mX; Z <- mZ;  kk1 <- ncol(X); kk2 <- ncol(Z)
 beta <- theta[1:kk1]
 gama <- theta[(kk1+1.0):(kk1+kk2)]                                                  
 eta <- as.vector(X%*%beta)	
 delta <- as.vector(Z%*%gama)
 if(linkmu == "logit") mu <- exp(eta)/(1.0+exp(eta)) #logit
 if(linkmu == "probit") mu <- pnorm(eta) #probit
 if(linkmu == "cloglog") mu <- 1.0 - exp(-exp(eta)) #cloglog
 if(linkmu == "log") mu <- exp(eta) #log
 if(linkmu == "loglog") mu <- exp(-exp(-eta)) #loglog
 if(linkmu == "cauchit") mu <- (pi^(-1))*atan(eta) + 0.5 #cauchit 
 if(linkphi == "log") phi <- exp(delta) #log
 if(linkphi == "identify") phi <- delta #identify
 if(linkphi == "sqrt") phi <- delta^2 #sqrt
 a <- mu*phi; b <- (1.0 - mu)*phi
 y <- rbeta(n, a, b)                                            
 results <- list(y=y, X=X, Z=Z, mu=mu, phi= phi)
 return(results)
}
