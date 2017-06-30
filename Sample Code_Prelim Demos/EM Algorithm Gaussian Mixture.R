# EM Algorithm for Gaussian Mixture
# May 2017
# Jennifer Starling

rm(list=ls())

#===========================================================================
#=== Genererate some simulated data from a 2-component Gaussian mixture. ===
#===========================================================================

true.weights = c(.3,.7)
n = 1000
true.mu = c(-1,5)
true.cov = matrix(c(.1,0,0,.1), nrow=2, byrow=T) #Two gaussians are independent.

y = c( rnorm(n*true.weights[1], true.mu[1], sqrt(true.cov[1,1])) , 
	   rnorm(n*true.weights[2], true.mu[2], sqrt(true.cov[2,2])) )

hist(y,breaks=30)

#===========================================================================
#=== EM Algorithm Functions ================================================
#===========================================================================

my.em = function(y,niter,mu0,cov0){
	#-----------------------------------------------------------------------
	# FUNCTION:		EM Algorithm function for 2-component gaussian mixture.
	#-----------------------------------------------------------------------
	# INPUTS:		y 		= data values
	#				niter 	= number of EM algorithm iterations	
	#				mu0 	= starting guess for vector of means
	#				cov0    = starting guess for covariance matrix  
	#-----------------------------------------------------------------------
	# OUTPUT:		w.j		=
	#				mu.j	=
	#				cov.j	=
	#				d.i		=
	#-----------------------------------------------------------------------
	
	#-----------------------------------------------------------------------
	# Setup
	#-----------------------------------------------------------------------
	
	# Initialize data structures.
	n 				= length(y)				#Number of observations.
	w.j 			= matrix(NA,2,niter)	#Vectors (cols) of component weights for each iter.
	rownames(w.j) 	= c('w1','w2')
	logl			= rep(NA,niter)			#Vector of logl for each iter.
	
	mu.j  = matrix(NA,2,niter)	#Vectors (cols) of component means for each iter.
	cov.j = list()				#List of cov matrices for each iter.
	
	# Initialize algorithm starting values.
	mu.j[,1] = rep(0,2)
	cov.j[[1]] = diag(1,2)
	
	#Initialize random group memberships.
	w.j[,1] = rep(.5,2)
	w.ij[,1] = sample(1:2,size=n,replace=T,prob=c(.5,.5))
	
	#-----------------------------------------------------------------------
	# EM algorithm updates.
	#-----------------------------------------------------------------------
	
	for (i in 2:niter){
		
		#----------------------------------------------------------------
		# E Step
		#----------------------------------------------------------------
		
		# Compute membership weights for each observation.
		probs = matrix(NA,n,2)
		
		for (j in 1:2){
			
			num 	= w.j[j,i-1] * dnorm(y,mu.j[j,i-1], sqrt(cov.j[[i-1]][j,j]))
			denom 	= w.j[1,i-1] * dnorm(y,mu.j[1,i-1], sqrt(cov.j[[i-1]][1,1])) + 
					  w.j[2,i-1] * dnorm(y,mu.j[2,i-1], sqrt(cov.j[[i-1]][2,2]))
		
			probs[,j] = num / denom
		}
		
		#Normalize probs.
		probs = probs / rowSums(probs) 
		
		#----------------------------------------------------------------
		# M Step
		#----------------------------------------------------------------
		
		# Compute new weights.
		Nk = colSums(probs) #Effective number of data points assigned to component k (sum of membership weights for kth component)
		w.j[,i] = Nk/n		#New weights.
		
		# Updated means.
		mu.1.new = (1/Nk[1]) * sum(w.j[,1] * y)
		mu.2.new = (1/Nk[2]) * sum(w.j[,2] * y)
		mu.j[,i] = c(mu.1.new, mu.2.new)
		
		# Updated covariance.
		cov.11.new = (1/Nk[1]) * sum(w.j[,1] * (y-mu.1.new) %*% t(y - mu.1.new))
		cov.22.new = (1/Nk[2]) * sum(w.j[,2] * (y-mu.2.new) %*% t(y - mu.2.new))
		
		cov.j[[i]] = matrix(c(cov.11.new, 0, 0, cov.22.new), nrow=2, byrow=T)
		
		#----------------------------------------------------------------
		# Update logl and check convergence.
		#----------------------------------------------------------------
		
		logl[i] = 0
		
		#for (k in 1:2){
			#logl[i] = logl[i] + w.j[k,i] * dnorm(y,mu.j[k,i],sqrt(cov.j[[i]][k,k]))
		#}
		
	}
	
	#-----------------------------------------------------------------------
	# Return output.
	#-----------------------------------------------------------------------
	
} # End my.em function



