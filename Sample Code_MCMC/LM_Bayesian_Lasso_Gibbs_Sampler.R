#SDS 387 - Homework 4
#Bayesian Lasso - Gibbs Sampler (Diabetes Data)
#Jennifer Starling
#March 2017

#================================================================
# Housekeeping ==================================================
#================================================================
rm(list=ls())

library(mvtnorm)	#For generating multivariate normal rvs.
library(SuppDists)	#For generating inverse gaussian rvs.
library(lars)   
library(plotrix)	#For credible intervals.   

#================================================================
# Load Data =====================================================
#================================================================

data(diabetes)  

y = diabetes$y
X = diabetes$x2  # 64 variables from all 10 main effects,
                  # two-way interactions and quadradics
X = scale(X,center=T)    # rescale so that all variables have mean 0 and sd 1
y = scale(y,center=T)

#Preview covariates.
X[1:10,1:5]

#================================================================
# Gibbs Sampler =================================================
#================================================================

bayesian.lasso.gibbs = function(X,y,r=1,delta=1.78,B=11000,thin=0,burn=1000){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for the Bayesian Lasso as
	#			described in Park & Casella (2008, JASA)	
	#-------------------------------------------------------------
	#INPUTS: 	X = design matrix of covariates (nxp)
	#			y = vectors of responses (nx1)
	#			r,d = the shape and rate hyperparameters
	#				for the LASSO penalty lambda.
	#			B = number of iterations for gibbs sampler.
	#			thin = number of obs to thin out to reduce cov.
	#			burn = number of first observations to burn.
	#-------------------------------------------------------------
	#OUTPUTS:	fx = vector of realizations from gaussian process.
	#-------------------------------------------------------------
	
	n = length(y)	#Number of observations.
	p = ncol(X)		#Number of covariates.

	#Precache X'X and X'y.
	XtX = crossprod(X)
	Xty = crossprod(X,y)

	#Initialize data structures to hold gibbs sampler realizations.
	beta = matrix(0,nrow=p,ncol=B)	#beta_1,...,beta_p in each col.
	sig2 = rep(0,B)
	tau2 = matrix(0,nrow=p,ncol=B) #tau1.sq, ...,taup.sq in each col.
	lambda = rep(0,B)	

	#Initialize first values to 1.
	beta[,1] = 0
	sig2[1] = 1
	tau2[,1] = 1
	lambda[1] = 1

	#Iterate gibbs sampler.
	for (b in 2:B){
	
		#Update beta full conditional.
		Dinv = diag(1/tau2[,b-1])
		Ainv = solve(XtX + Dinv)
		beta[,b] = rmvnorm(1,Ainv %*% Xty, sig2[b-1]*Ainv)
	
		#Update sig2 full conditional. (Uses new beta, old tau2.)
		sh = (n -1)/2 + p/2
		rt = .5 * (t(y - X %*% beta[,b]) %*% (y - X %*% beta[,b]) + t(beta[,b]) %*% Dinv %*% beta[,b])
		sig2[b] = 1/rgamma(1,shape=sh,rate=sc)
	
		#Update tau1...taup full conditionals. (Uses new beta, old lambda, new sig2.)
		nu = sqrt(sig2[b] * lambda[b-1]^2 / beta[,b]^2)
		tau2[,b] = 1 / rinvGauss(p,nu,lambda[b-1]^2)
		
		#Update lambda full conditional.
		lambda[b] = sqrt(rgamma(1,shape=p+r,rate=sum(tau2[,b])/2 + delta))
	}
	
	#Burn and thin.
	beta = beta[,-c(1:burn)]
	sig2 = sig2[-c(1:burn)]
	tau2 = tau2[,-c(1:burn)]

	if (thin>0){
		beta = beta[,seq(1,ncol(beta),by=thin)]
		sig2 = sig2[seq(1,length(sig2),by=thin)]
		tau2 = tau2[,seq(1,ncol(tau2),by=thin)]
	}
		
	#Return posterior samples from full conditionals.
	return(list(beta=beta,sig2=sig2,tau2=tau2,lambda=lambda))
	
} #End gibbs sampler code.

#================================================================
# Run Gibbs Sampler =============================================
#================================================================

#Set as the posterior mean from Park & Casella.
output = bayesian.lasso.gibbs(X,y,r=1,delta=1.78,B=11000,thin=0,burn=1000)

beta.post = output$beta
sig2.post = output$sig2
tau2.post = output$tau2
lambda.post = output$lambda

#Calculate posterior means.
beta.post.mean = rowMeans(beta.post)
sig2.post.mean = mean(sig2.post)
tau2.post.mean = rowMeans(tau2.post)
lambda.post.mean = mean(lambda.post)

#Posterior variances.
beta.post.var = apply(beta.post,1,var)
sig2.post.var = var(sig2.post)
tau2.post.var = apply(tau2.post,1,var)
lambda.post.var = var(lambda.post)

#Threshold all betas.
thresh = .05
selected.betas.idx = which(beta.post.mean > thresh) 
length(selected.betas.idx)

#Print names of selected betas.
cbind(colnames(X)[selected.betas.idx])

#================================================================
# Plot Results ==================================================
#================================================================

#Beta bayesian credible intervals.

pdf('/Users/jennstarling/UTAustin/2017S_Linear Models/Homework/HW 04/Figures/CIplots_beta.pdf')
t = qt(1-.05/2,n-p)

CI = matrix(cbind(
	beta.post.mean - t * sqrt(beta.post.var), 
	beta.post.mean, 
	beta.post.mean + t * sqrt(beta.post.var)),ncol=3,byrow=F)

plotCI(x=CI[1:10,2],li=CI[1:10,1],ui=CI[1:10,3],xlab='Variable Number',ylab='Standardized Coefficients',
	main=bquote("Bayesian Credible Intervals for First 10" ~ beta[j]))
abline(h=0,lty=2,col='grey',lwd=2)
dev.off()

