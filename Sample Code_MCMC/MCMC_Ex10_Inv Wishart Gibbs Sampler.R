#SDS 383D - Exercise 4
#Function
#Jennifer Starling
#March 2017

#================================================================
# Housekeeping ==================================================
#================================================================
rm(list=ls())
library(mvtnorm)		#For sampling multivariate normal.
library(MCMCpack)		#For sampling inverse wishart.

#================================================================
# Generate Data =================================================
#================================================================
m = 5
n = 100

#Generate iid y_i's for i=1...100, with y_i ~ Nm(0,Im).
y = rmvnorm(n, rep(0,m), diag(m)) #Each y_i is a column.

#================================================================
# Inverse Wishart Gibbs Sampler Function ========================
#================================================================

gibbs = function(y,Omega,k,R,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Samler for math test hierarchical normal model.
	#			Unknown parameters: (tau.sq, sig.sq, mu)	
	#-------------------------------------------------------------
	#MODEL:		y_i ~ N(mu, Sigma) where y_i are (mx1) vectors.
	#			mu ~ N(0, Omega) where mu is a (mx1) vector.
	#			Sigma ~ InvWishart(k,R)
	#-------------------------------------------------------------
	#INPUTS: 	y = matrix of MVN realizations (each y_i is a column)
	#			Omega = hyperparameter value for variance of mu.
	#			k = df prior for Inv Wishart.
	#			R = scale matrix prior for Inv Wishart.
	#-------------------------------------------------------------
	#OUTPUTS:	theta.post = matrix of posterior theta_i samples (rows = samples)
	#			mu.post = vector of posterior mu samples
	#			sig.sq.post = vector of posterior sig2.sq samples
	#			tau.sq.post = vector of posterior tau2.sq samples
	#-------------------------------------------------------------
	
	n = ncol(y)			#Number of observations.
	m = nrow(y)			#Dimension of multivariate normal likelihood.
	ybar = rowMeans(y)	#The (mx1) vector of sample means.
	
	#Data structures to hold sampled values.
	mu = matrix(0,m,iter)	#Each column is a sampled mu.
	Sigma = array(0,dim=c(m,m,iter))

	#Initialize each element of chain.
	mu[,1] = rep(0,m)
	Sigma[,,1] = diag(m)
	
	#Precache inverse of Omega.
	Omega.inv = solve(Omega)
	
	#Iterate through sampler.
	for (i in 2:iter){
		
		#Update mu.
		Sigma.inv = solve(Sigma[,,i-1])	#Precache Sigma inverse.
		mu.var = solve(n * Sigma.inv + Omega.inv)		
		mu.mean = n * mu.var %*% Sigma.inv %*% ybar
		mu[,i] = rmvnorm(1,mu.mean,mu.var)
		
		#Update Sigma.
		Scale = R + (y - mu[,i]) %*% t(y - mu[,i])
		Sigma[,,i] = riwish(k+n,Scale)
	}
	
	#Burn beginning observations.
	if (burn > 0){
		mu = mu[,-burn]
		Sigma = Sigma[,,-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		mu = mu[,seq(1,ncol(mu),by=thin)]
		Sigma = Sigma[,,seq(1,dim(Sigma)[3],by=thin)]
	}
	
	#Return results.
	return(list(mu=mu,Sigma=Sigma))
}

#================================================================
# Set hyperparameters and run sampler. ==========================
#================================================================

#Hyperparameters
Omega = 100 * diag(m)
k = 6
R = diag(m)

#Run sampler.
output = gibbs(y,Omega,k,R,iter=11000,burn=1000,thin=2)

#Posterior mean of mu.
mu.post.mean = rowMeans(output$mu)
mu.post.mean

#Posterior mean of Sigma.
Sigma.post.mean = apply(output$Sigma,c(1,2),mean)
Sigma.post.mean


