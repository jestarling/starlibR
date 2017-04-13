# MCMC Exercise 12
# Gaussian Mixture Model with Infinite (M) Components
# April 2017
# Jennifer Starling

rm(list=ls())

#================================================================
# Gibbs Sampler =================================================
#================================================================

gibbs = function(y,alpha,Napp,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for infinite-component Gaussian mixture model with known weights.
	#-------------------------------------------------------------
	#MODEL:		
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of observations (data).
	#		a, b = hyperparams for sig.sq ~ IG(a,b).
	#		m, tau = hyperparams for mu.j ~ N(m,1/tau)
	#-------------------------------------------------------------
	#OUTPUTS:	
	#-------------------------------------------------------------
	require(Hmisc) #For generating random multinomials (vectorized).
	
	#Number of data points.
	n = length(y)
	
	#Hyperparameters for mu.j's.  mu.j ~ N(m,tau.sq)
	m = 0
	tau.sq = 100
	
	#Hyperparameters for lambda.  lambda ~ IG(a,b)
	a = .5
	b = .5
	
	#Calculate the deterministic (known) weight sequence.
	wi.seq = alpha^(1:Napp)
	
	#Data structures to hold values for each Gibbs iteration.
	#Using lists here, since dimension varies for each iteration.
	
	M = rep(0,iter)	#For tracking number of components in each sample.
	mu.j = list()		#mu components for each mu.j.
	w.j = list()		#Holds vector of weights for each component.
	d.i = list()		#Component assignments for each obs.
	u.i = list()		#Holds mu.i for each iteration for each obs.

	lambda = rep(.1,iter)	#Holds lambda values and initializes to .1.
	x.pred = rep(0,iter)	#Holds x from predictive for each iter.
	
	#Begin by selecting a random M>1, and randomizing some d.i values.
	M[1] = 5					#Indicates number of components in 1st Gibbs iter.
	d.i[[1]] = sample(1:M[1],n,replace=T)	#Indicates which component each obs belongs to.
	
	#Initialize weights based on initial M.
	w.j[[1]] = wi.seq[1:M[1]]
	
	#Initialize mu.j and u.i values based on initial M.
	mu.j[[1]] = rep(0,M[1])
	w.di = w.j[[1]][d.i[[1]]]
	u.i[[1]] = runif(n,0,w.di[[1]])
	
	#-------------------------------------------------------------
	#Iterate through sampler.
	for (k in 2:iter){
		
		### Progress report.
		if(k %% 100 == 0) {
			print(sprintf('Iteration %i out of %i...',k,iter))
		}
		
		#-------------------------------------------------------
		### Update lambda (1/sig.sq).
		
		#\sum_{i=1}^{n}(y.i - mu_{d.i})^2
		SS = sum((y - mu.j[[k-1]][d.i[[k-1]]])^2)
		
		lambda[[k]] = rgamma(1, n/2+a, b + SS/2)
		
		#-------------------------------------------------------
		### Update mu.j (normal-normal update within each group)
		nj = as.numeric(table(d.i[[k-1]]))
		ybar = aggregate(y,list(d.i[[k-1]]),mean)$x
		
		var = 1 / (1/tau.sq + nj*lambda[[k-1]])
		mean = var * ( (1/tau.sq)*m + nj*lambda[[k-1]]*ybar)
				
		mu.j[[k]] = rnorm(M[k-1],mean,sqrt(var))
		
		#-------------------------------------------------------
		### Update u.i.
		w.di = w.j[[k-1]][d.i[[k-1]]]	 #Weight for each obs i; w_{di}
		u.i[[k]] = runif(n,0,w.di)

		#-------------------------------------------------------
		### Update M, p(d.i=j) and d.i.
		
		M[k] = max(floor(log(u.i[[k]])/log(.5)))
		
			# NOTE: How to determine number of new components:
			# u.i < (1/2)^2 --> j < log(u.i)/log(.5)
			# Number of new components = max{ floor( log(u.i) / log(.5) ) }
		
		#Update weights based on new number of possible components.  
		#   (Known, just grabbing vector of right length.)
		w.j[[k]] = wi.seq[1:M[k]]
		
		#Adjust mu.j based on having > or < components in M+1 than M.
		if (M[k] > M[k-1]){
			#Fill in extra mu.j[[k]] components with samples from prior.
			new.mu = rnorm(M[k] - M[k-1], m, sqrt(tau.sq))
			mu.j[[k]] = c(mu.j[[k]], new.mu)
		}
		
		if (M[k] < M[k-1]){
			print('M[k] < M[k-1] IS TRUE')
			#Remove extra mu.j[[k]] components.
			mu.j[[k]] = mu.j[[k]][1:M[k]]
		}
		
		#Update vectors of P(d.i = j | ...) based on new weights. 
		probs = matrix(0,n,M[k])  #A (x by K) matrix to hold prob vecs for each obs.
		
		for (j in 1:M[k]){
			ind = u.i[[k]] < w.j[[k]][j]
			probs[,j] = ind * w.j[[k]][j] * dnorm(y,mu.j[[k]][j],sqrt(1/lambda[[k]]))
		}
		
		#Normalize probs for each obs i.
		probs = probs / rowSums(probs) 
		probs[is.na(probs)] = 0	#Weights with 0 in all cols will give NA.  Reset to zero.
	
		#Draw new P(d.i = j | ...) from multinomial.
		d.i[[k]] = as.numeric(rMultinom(probs,1))
			
		#-------------------------------------------------------------
		### GAPS UPDATE:
		#Calculate sample sizes in all categories.
		ncounts = as.numeric(table(factor(d.i[[k]],levels=1:M[[k]])))
		empty = which(ncounts==0)
		nonempty = which(ncounts>0)
		
		if (length(empty)>0){
			#Drop any zero categories from mu and w, and update M to include 
			#count of nonzero categories only.
			mu.j[[k]] = mu.j[[k]][-which(ncounts==0)]
			w.j[[k]] = w.j[[k]][-which(ncounts==0)]
		
			M[[k]] = sum(ncounts>0)
			
			#Shift d.i values to remove gaps.  
			#Ex: 2,4,5 as nonzeros become 1,2,3.
			shifted.di.list = sort(unique(d.i[[k]]))
			d.i[[k]] = match(d.i[[k]],shifted.di.list)	
		}
			
		#-------------------------------------------------------------
		### Sample a new x from predictive density.
		
		#First sample a component using the updated weights.
		new.component = sample(1:M[[k]], 1, prob=w.j[[k]])
		x.pred[[k]] = rnorm(1,mu.j[[k]][new.component],sqrt(1/lambda[k]))
			
	} #End Gibbs Sampler.
	
	#Burn beginning observations.
	if (burn > 0){
		x.pred = x.pred[-burn]
		lambda = lambda[-burn]
		M = M[-burn]
		mu.j = mu.j[-burn]
		d.i = d.i[-burn]
		w.j = w.j[-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		idx.keep = seq(1,length(x.pred),by=thin)
		x.pred = x.pred[idx.keep]
		lambda = lambda[idx.keep]
		M = M[idx.keep]
		mu.j = mu.j[idx.keep]
		d.i = d.i[idx.keep]
		w.j = w.j[idx.keep]
	}
	
	#Return results.
	return(list(x.pred=x.pred, lambda=lambda, M=M, mu.j=mu.j, d.i=d.i, w.j=w.j))
} #End function.


#================================================================
# Test Sampler on N(0,1) Single-Component Data ==================
#================================================================

# Generate data; n=100 obs from the standard normal density.
n=100
y = rnorm(n,0,1)

# Run Gibbs Sampler.
output = gibbs(y,alpha=.5,Napp=100,iter=6000,burn=1000,thin=2)

#Plot of predictive.
#pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-12/Figures/weights_known_stdNormal.pdf')
hist(y,breaks=30,freq=F,main='Predicted Density in Blue')
lines(density(output$x.pred),col='blue',lwd=2)
#dev.off()

#================================================================
# Run Sampler on 3-Class Well-Separated Data & Analyze Results ==
#================================================================

# Generate well-separated 2-class data.
y = c(rnorm(100,1,.5),rnorm(100,4,.5))
hist(y,breaks=50)

output = gibbs(y,alpha=.5,Napp=100,iter=6000,burn=1000,thin=2)

#Plot of predictive.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-12/Figures/inf_gaussian_posterior_pred_3class.pdf')
hist(output$x.pred,breaks=50,main='Plot of Posterior Predictive Density', xlab='Posterior Predictive',freq=F)
curve(.5*dnorm(x,1,.5) + .5*dnorm(x,4,.5),col='blue',add=T)
dev.off()