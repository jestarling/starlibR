#----------------------------------------------------------------
# Roxy package build comments:
#' ttestCorrSim(n,B,mu) A function to model effect of correlation on 2-sample t-test.
#'
#' This function calculates the power of an independent 2-sample t-test, to detect
#' various differences in pairs of means (mu), under the presence of a variety of 
#' correlation (rho) values.  Produces data and a plot.
#'
#' @param n The sample size of for each simulation.  Defaults to n=30.
#' @param B The number of simulations to perform for each pair of mean (mu) values. Defaults to 1000.
#' @param mu The matrix of pairs of mean values to compare.  Each pair of mu values is a row.
#'		Defaults to matrix(c(0,0,0,.25,0,.5,0,.75),byrow=T,nrow=4), which compares 0 to .25, .5, and .75.
#'
#' @keywords t-test correlation simulation power
#'
#' @return x A list, containing a data frame for each pair of means.
#' @return x[[m]] A data frame for the mth mean pair comparison. Contains means, rho, and power.
#' @return plot Displays a plot of the all simulations.
#'
#' @examples
#' ttestCorrSim(n=30,B=100,mu=matrix(c(0,0,0,.25,0,.5,0,.75),byrow=T,nrow=4))
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

ttestCorrSim = function(n=30,B=1000, mu=matrix(c(0,0,0,.25,0,.5,0,.75),byrow=T,nrow=4) ){
	
	library(mvtnorm)
	n=n 	#Sample size for each simulation (ex: 30 pairs).
	B=B		#The number of simulations to perform.
	rho <- seq(-0.5, 0.95, by = 0.05) #Vector of correlations to test.
	mu <- matrix(c(0,0,0,.25,0,.5,0,.75),byrow=T,nrow=4) #Vector of mu's to test

#Set up 2 vectors to hold p-values of 2-sample and paired t-tests for each sim.
	p_val_tt <- p_val_pt <- matrix(NA, nrow = length(rho), ncol = B)

	par(mfrow=c(nrow(mu),nrow(mu))) #Display mxm grid of the four plots.  (One plot per set of mu values.)
	
	x <- list() #Create list output to store values.

#Loop through four mu value pairs.
	for (m in 1:4){
	
	#Loop through each rho value.
	for(i in 1:length(rho)) {
  	
	# Define Sigma.
  		Sigma <- matrix(c(1, rho[i], rho[i], 1), nrow = 2)
  
		for(b in 1:B) {
    		# Simulate data.
    		y <- rmvnorm(n, mu[m,], Sigma) #Generate correlated data sample.

    		# Store p-value for Two-sample t-test.
    		p_val_tt[i, b] <- t.test(y[, 1], y[, 2], var.equal = TRUE)$p.value
    		# Store p-value for Paired t-test.
    		p_val_pt[i, b] <- t.test(y[, 1], y[, 2], paired = TRUE)$p.value
  		}
	}

	#Calculate proportion of H0 rejections for each row of rho values.
	#Proportion is calculated as (# Rejections / B Simulations)
	rej_tt <- ifelse(p_val_tt < 0.05,1,0)
	power_tt <- rowMeans(rej_tt)

	rej_pt <- ifelse(p_val_pt < 0.05,1,0)
	power_pt <- rowMeans(rej_pt)

	#Plot rho vs power for the specific mu value.
	plot(rho,power_tt,type='l',col='blue',xlim=c(-.5,1),ylim=c(0,1),xlab='rho',
		ylab='Proportion Times H0 Rejected',
		  main=paste('Proportion of Times H0: Equal Means Rejected, mu = ',mu[m,1]," & ",
		  mu[m,2]), cex.main=.8,cex.lab=.7)
	lines(rho,power_pt,type='l',col='red')	
	
	x[[m]] <- data.frame(mu1 = rep(mu[m,1],length(rho)),
						mu2 = rep(mu[m,2],length(rho)), 
						rho=rho, 
						ttestPower=power_tt)
	
	}

	return(x)
}




