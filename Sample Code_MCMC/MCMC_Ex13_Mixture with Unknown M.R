# MCMC Exercise 13
# Mixture Model of M Exponential(mu)'s with Unknown M, P(M) prior is shifted poisson.
# April 2017
# Jennifer Starling

rm(list=ls())

#================================================================
# Function to Generate Data from g(x) ===========================
#================================================================

data.gen = function(n){
	#-------------------------------------------------------------
	# FUNCTION:	Generates samples from density g(x) = 24 / (2+x)^4 
	#			for x>0, using inverse cdf method.
	#-------------------------------------------------------------
	# INPUTS:	n = sample size.
	#-------------------------------------------------------------
	# OUTPUTS:  gx = sample of size n from density.
	#-------------------------------------------------------------
	u = runif(n,0,1)
	gx = (8 / (1-u))^(1/3) - 2
	return(gx)
}

#================================================================
# Gibbs Sampler with MH Step for (M,mu,w) =======================
#================================================================

gibbs = function(x, M.init = 5, iter=11000, burn=1000, thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for infinite-component Gaussian mixture model with unknown weights.
	#-------------------------------------------------------------
	#MODEL:		
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of observations (data).
	#			M.init = initial number of components.
	#			iter = Gibbs sampler iterations.
	#			burn = number of iterations to discard from beginning of chain.
	#			thin = thinning frequency to remove autocorrelation.
	#-------------------------------------------------------------
	#OUTPUTS:	x.pred = sample from predictive density.
	#-------------------------------------------------------------
	
	require(Hmisc) 		#For generating random multinomials (vectorized).
	require(MCMCpack)	#For Dirichlet.
	
	#Number of data points.
	n = length(x)
	
	#Track acceptance probs on MH step.
	accept.up = rep(NA,iter)
	accept.down = rep(NA,iter)
	
	#-------------------------------------------------------------
	# Set up data structures.
	#-------------------------------------------------------------	
	
	#Data structures to hold values for each Gibbs iteration.
	#Using lists here, since dimension varies for each iteration.
	
	M	 =rep(NA,iter)		#Tracks number of components in each iter.
	mu.j = list()		#mu components for each mu.j.
	w.j  = list()		#Holds vector of weights for each component.
	d.i  = list()		#Component assignments for each obs.
	mu.i = list()		#Holds mu.i for each iteration for each obs.

	lambda = rep(.1,iter)	#Holds lambda values and initializes to .1.
	x.pred = rep(0,iter)	#Holds x from predictive for each iter.
	
	#Begin by selecting a random M>1 (user input), and randomizing some d.i values.
	M[1] = M.init							#Indicates number of components in 1st Gibbs iter.
	d.i[[1]] = sample(1:M[1],n,replace=T)	#Indicates which component each obs belongs to.
	
	#Initialize weights based on initial M.
	w.j[[1]] = rdirichlet(1,rep(1,M[1]))
	
	#Initialize mu.j and u.i values based on initial M.
	mu.j[[1]] = rep(0,M[1])
	w.di 	  = w.j[[1]][d.i[[1]]]
	
	#-------------------------------------------------------------
	# Gibbs Sampler.
	#-------------------------------------------------------------
	
	#Iterate through sampler.
	for (k in 2:iter){
		
		### Progress report.
		if(k %% 100 == 0) {
			print(sprintf('Iteration %i out of %i...',k,iter))
		}
		
		#-------------------------------------------------------------
		# Update mu.j ('shuffle')
		#-------------------------------------------------------------
		nj 	  = table(d.i[[k-1]])	#Sample size for each component.
		alpha = nj + .5
		beta  = .5 + aggregate(x, list(d.i[[1]]), sum)$x
		
		mu.j[[k]] = rgamma(length(nj), alpha, beta)
		
		#-----------------------
		# Prepare for M+1 move.
		#-----------------------
		
		#Pre-draw extra mu.j in case we accept an M+1 move in MH step.
		#We will take the average of the current mu.j values.
		mu.new = mean(mu.j[[k]])			#New mu value.
		mu.j.Mp1 = c(mu.j[[k]], mu.new)		#Preache vector of updated mu values for M+1 move.
		
		#-----------------------
		# Prepare for M-1 move.
		#-----------------------
		
		#We will randomly select a weight and a mean to remove in case of an M-1 MH step.
		#Note: could also do this by choosing two weights and combining.
		remove 	 = sample(1:length(nj),1)
		mu.j.Mm1 =  mu.j[[k]][-remove]
		
		#-------------------------------------------------------------
		# Update weights w ('shuffle')
		#-------------------------------------------------------------
		nj 	     = table(d.i[[k-1]])	#Sample size for each component.
		w.j[[k]] = rdirichlet(1, nj+1) 
		
		#-----------------------
		# Prepare for M+1 move.
		#-----------------------
		
		#Pre-draw extra w.j in case we accept an M+1 move in MH step.
		#We will randomly pick a weight and split it in half.  
		wt.to.split  = sample(1:length(nj),1)
		
		#Precache vector of updated weights for M+1 move.
		w.new = w.j[[k]][wt.to.split]/2
		w.j.Mp1 = c(w.j[[k]], w.new)		
		w.j.Mp1[wt.to.split] = w.j.Mp1[wt.to.split]/2 #Split the value of the chosen weight also.
		
		#-----------------------
		# Prepare for M-1 move.
		#-----------------------
		
		#Remove the weight corresponding with the randomly chosen mean above, then renormalize.
		w.j.Mm1 = w.j[[k]][-remove]
		w.j.Mm1 = w.j.Mm1 / sum(w.j.Mm1)
	
		#-------------------------------------------------------------
		# MH Step to update (M, mu, w) jointly.
		#-------------------------------------------------------------
		
		# Propose a move to M+1 or M-1, with probs 1/2.
		# If M=1, propose a move to M+2 with prob 1.
		proposed.move = NA
		q1 = q2 = .5
		
		if(M[[k-1]]==1){
			proposed.move = 'up'
			q1 = 1
		} else{
			proposed.move = ifelse(rbinom(1,1,.5)==1,'up','down')
		}	
		
		#-----------------------
		# Proposed M+1 actions.
		#-----------------------
		
		#If M+1 proposed, calculate P(accept) and update M accordingly.
		if(proposed.move=='up'){
			
			#Calculate acceptance probability for M+1.
			p.accept = p.accept.up(y, M[[k-1]], mu = mu.j[[k]], w=w.j[[k]], mu.new=mu.j.Mp1, w.new=w.j.Mp1, q1=.5, q2=.5)
			
			#If p.accept > .5, select new M = M+1 and update mu and w, else set new M = M.
			if(p.accept > .5){
				
				#Update M.
				M[k] = M[[k-1]]+1
				
				#Save pre-cached mu and w as the new values.
				mu.j[[k]] = mu.j.Mp1
				w.j[[k]]  = w.j.Mp1
				
			} else{
				#Leave M at current value.  No change to mu and w.
				M[k] = M[[k-1]]
			}
			
			#Track move and acceptance.
			accept.up[k] = p.accept
		}
		
		#-----------------------
		# Proposed M-1 actions.
		#-----------------------
		
		#If M-1 proposed, calculate P(accept) and update M accordingly.
		if(proposed.move=='down'){
			
			#Calculate acceptance probability for M-1.
			p.accept = p.accept.down(y, M[[k-1]], mu = mu.j[[k]], w=w.j[[k]], mu.new=mu.j.Mm1, w.new=w.j.Mm1,remove, q1=.5, q2=.5)
			
			#If p.accept > .5, select new M = M-1 and update mu and w, else set new M = M.
			if(p.accept > .5){
				
				#Update M.
				M[k] = M[[k-1]]-1
				
				#Save pre-cached mu and w as the new values.
				mu.j[[k]] = mu.j.Mm1
				w.j[[k]]  = w.j.Mm1
				
			} else{
				#Leave M at current value.  No change to mu and w.
				M[k] = M[[k-1]]
			}
			
			#Track move and acceptance.
			accept.down[k] = p.accept
		}
		
		#-------------------------------------------------------------
		# Update P(di=j) and group memberships based on new M.
		#-------------------------------------------------------------
		
		#Update vectors of P(d.i = j | ...) based on new weights. 
		probs = matrix(0,n,M[k]) 	#A (x by K) matrix to hold prob vecs for each obs.
		
		for (j in 1:M[k]){
			probs[,j] = w.j[[k]][j] * dexp(x, mu.j[[k]][j])
		}
		
		#Normalize probs for each obs i.
		probs = probs / rowSums(probs) 
		probs[is.na(probs)] = 0	#Weights with 0 in all cols will give NA.  Reset to zero.
		
		#Draw new P(d.i = j | ...) from multinomial.  If M[[k]] = 1, all components are group 1.
		if(ncol(probs)>1){
			d.i[[k]] = as.numeric(rMultinom(probs,1))
		} else{
			d.i[[k]] = rep(1,n)
		}	
		
		#-------------------------------------------------------------
		# Gaps removal.
		#-------------------------------------------------------------
		
		#Calculate sample sizes in all categories.
		ncounts = as.numeric(table(factor(d.i[[k]],levels=1:M[k])))
		empty = which(ncounts==0)
		nonempty = which(ncounts>0)
		
		if (length(empty)>0){
			#Drop any zero categories from mu and w, and update M to include 
			#count of nonzero categories only.
			mu.j[[k]] = mu.j[[k]][-which(ncounts==0)]
			w.j[[k]] = w.j[[k]][-which(ncounts==0)]
		
			M[k] = sum(ncounts>0)
			
			#Shift d.i values to remove gaps.  
			#Ex: 2,4,5 as nonzeros become 1,2,3.
			shifted.di.list = sort(unique(d.i[[k]]))
			d.i[[k]] = match(d.i[[k]],shifted.di.list)	
		}
			
		#-------------------------------------------------------------
		### Sample a new x from predictive density.
		
		#First sample a component using the updated weights.
		new.component = sample(1:M[k], 1, prob=w.j[[k]])
		
		#Sample from exponential(mu.j) using selected component.
		x.pred[k] = rexp(1,mu.j[[k]][new.component])
			
	} #End Gibbs Sampler.
	
	#Thin and burn gibbs samples.
	x.pred 	= x.pred[seq(burn+1,iter,by=thin)]
	M 		= M[seq(burn+1,iter,by=thin)] 
	mu.j 	= mu.j[seq(burn+1,iter,by=thin)]
	d.i 	= d.i[seq(burn+1,iter,by=thin)]
	w.j 	= w.j[seq(burn+1,iter,by=thin)]
	
	accept.up 	= accept.up[seq(burn+1,iter,by=thin)]
	accept.down = accept.down[seq(burn+1,iter,by=thin)]

	#Remove NAs from acceptance probability trackers.
	accept.up 	= accept.up[which(!is.na(accept.up))]
	accept.down = accept.down[which(!is.na(accept.down))]
	
	#Return results.
	return(list(x.pred=x.pred, M=M, mu.j=mu.j, d.i=d.i, w.j=w.j))
	
} #End function.


#================================================================
# P(accept) for M+1 Proposal ====================================
#================================================================

p.accept.up = function(y, M, mu, w, mu.new, w.new, q1=.5, q2=.5){
	#-------------------------------------------------------------
	#FUNCTION: 	Calculate P(accept) for M+1 proposal.
	#-------------------------------------------------------------
	#INPUTS: 	y  		= vector of observations (data).
	#			M  		= current number of components.
	#			mu 		= vector of current mu.j values. 
	#			w  		= vector of current w.j values.  
	#			mu.new 	= vector of proposed M+1 current mu.j values, from precache. (Length M+1)
	#			w.new  	= vector of proposed M+1 current w.j values, from precahce.  (Length M+1)
	#			q1, q2 	= probabilities for MH step.
	#-------------------------------------------------------------
	#OUTPUTS:	p.accept = probability of accepting the proposed move.	
	#-------------------------------------------------------------
	
	#-------------------------------------------------------------
	# Set up prior density functions.
	#-------------------------------------------------------------
	
	#Prior for M, number of components.  Shifted Pois(5).
	p.M = function(M) {
		dpois(M-1,5)
	}
	
	#Prior for rates mu.
	p.mu = function(mu){
		dgamma(mu,.5,.5)
	}
	
	#Prior for weights.  Dirichlet with all components set to 1.
	p.w = function(w){
		ddirichlet(w,alpha=rep(1,length(w)))
	}
	
	#-------------------------------------------------------------
	# Calculate likelihoods.
	#-------------------------------------------------------------
	
	#Likelihoods for numerator and denominator.
	lhood.Mp1 = 0
	lhood.M = 0
	
	for (j in 1:M){
		lhood.M =+ w[j] * dexp(x, mu[j])
	}
	
	for (j in 1:(M+1)){
		lhood.Mp1 =+ w.new[j] * dexp(x,mu.new[j])
	}
	
	#Set lhoods to products of sums.
	lhood.Mp1 	= prod(lhood.Mp1)
	lhood.M 	= prod(lhood.M)
	
	#-------------------------------------------------------------
	# Calculate priors for M+1 and M.
	#-------------------------------------------------------------
	pM.Mp1  = p.M(M+1)
	pmu.Mp1 = prod(p.mu(mu.new)) 	#Since mu.j are indep.
	pw.Mp1  = p.w(w.new)
	
	pM.M	= p.M(M)
	pmu.M	= prod(p.mu(mu)) 		#Since mu.j are indep.
	pw.M	= p.w(w)
	
	#-------------------------------------------------------------
	# Calculate P(accept).
	#-------------------------------------------------------------
	
	#p(mu^(m+1) | mu^(m)) = draw from prior
	p.mu.p1 = p.mu(mu.new[M+1])
	
	#p(w^m | w^(m+1)) = 1/m * p(w) where p(w) can be uniform or beta.
	p.w.p1 = (1/M) * runif(1,0,1)
	
	numerator 	= pw.Mp1 * pM.Mp1 * pmu.Mp1 * lhood.Mp1 * p.w.p1
	denominator = pw.M * pM.M * pmu.M * lhood.M * p.mu.p1
	p.accept 	= min(1, (q1 * numerator) / (q2 * denominator))

	#Error handling for NAN.
	if(is.na(p.accept>0)){ p.accept=0}
	
	return(p.accept)
}


#================================================================
# P(accept) for M-1 Proposal ====================================
#================================================================

p.accept.down = function(y, M, mu, w, mu.new, w.new, remove, q1=.5, q2=.5){
	#-------------------------------------------------------------
	#FUNCTION: 	Calculate P(accept) for M+1 proposal.
	#-------------------------------------------------------------
	#INPUTS: 	y  		= vector of observations (data).
	#			M  		= current number of components.
	#			mu 		= vector of current mu.j values. 
	#			w  		= vector of current w.j values.  
	#			mu.new 	= vector of proposed M-1 current mu.j values, from precache. (Length M+1)
	#			w.new  	= vector of proposed M-1 current w.j values, from precahce.  (Length M+1)
	#			remove  = index of component which was selected for removal.
	#			q1, q2 	= probabilities for MH step.
	#-------------------------------------------------------------
	#OUTPUTS:	p.accept = probability of accepting the proposed move.	
	#-------------------------------------------------------------
	
	#-------------------------------------------------------------
	# Set up prior density functions.
	#-------------------------------------------------------------
	
	#Prior for M, number of components.  Shifted Pois(5).
	p.M = function(M) {
		dpois(M-1,5)
	}
	
	#Prior for rates mu.
	p.mu = function(mu){
		dgamma(mu,.5,.5)
	}
	
	#Prior for weights.  Dirichlet with all components set to 1.
	p.w = function(w){
		ddirichlet(w,alpha=rep(1,length(w)))
	}
	
	#-------------------------------------------------------------
	# Calculate likelihoods.
	#-------------------------------------------------------------
	
	#Likelihoods for numerator and denominator.
	lhood.Mm1 = 0
	lhood.M = 0
	
	for (j in 1:M){
		lhood.M =+ w[j] * dexp(x, mu[j])
	}
	
	for (j in 1:(M-1)){
		lhood.Mm1 =+ w.new[j] * dexp(x,mu.new[j])
	}
	
	#Set lhoods to products of sums.
	lhood.Mm1 	= prod(lhood.Mm1)
	lhood.M 	= prod(lhood.M)
	
	#-------------------------------------------------------------
	# Calculate priors for M+1 and M.
	#-------------------------------------------------------------
	pM.Mm1  = p.M(M+1)
	pmu.Mm1 = prod(p.mu(mu.new)) 	#Since mu.j are indep.
	pw.Mm1  = p.w(w.new)
	
	pM.M	= p.M(M)
	pmu.M	= prod(p.mu(mu)) 		#Since mu.j are indep.
	pw.M	= p.w(w)
	
	#-------------------------------------------------------------
	# Calculate P(accept).
	#-------------------------------------------------------------
	
	#p(mu^(m) | mu^(m-1)) = draw from prior using mu to be removed.
	p.mu.p1 = p.mu(mu[remove])
	
	#p(w^m | w^(m-1)) = 1/m * p(w) where p(w) can be uniform or beta.
	p.w.num   = (1/M) * runif(1,0,1)
	p.w.denom = 1/M
	
	numerator 	= pw.Mm1 * pM.Mm1 * pmu.Mm1 * lhood.Mm1 * p.w.num * p.mu.p1
	denominator = pw.M * pM.M * pmu.M * lhood.M * p.w.denom
	p.accept 	= min(1, (q1 * numerator) / (q2 * denominator))

	#Error handling for NAN.
	if(is.na(p.accept>0)){ p.accept=0}
	
	return(p.accept)
}

#================================================================
# Generate data and run gibbs sampler ===========================
#================================================================

n = 100
x = data.gen(n)

output = gibbs(x,iter=11000,burn=1000,thin=2)

par(mfrow=c(1,2))
hist(output$accept.up)
hist(output$accept.down)

#Plot of predictive.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-13/Figures/Exp_Mixture_Predictive.pdf')
hist(x, breaks=30, freq=F,main='Exponential Mixture, 1 Component')
lines(density(output$x.pred), col='blue',lwd=2)
dev.off()

#===========================================================================
# Try generating data from a mixture of three well-separated exponentials ==
#===========================================================================

x = c(rexp(300,.5), rexp(300,5), rexp(300,10))
hist(x,freq=F)

output = gibbs(x, M.init=5, iter=11000,burn=1000,thin=2)
output$M

#Plot predictive.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-13/Figures/Exp_Mixture_Predictive_3 Components.pdf')
hist(x, breaks=30, freq=F,main='Exponential Mixture, 3 Component')
lines(density(output$x.pred), col='blue',lwd=2)
dev.off()