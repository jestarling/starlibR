#MCMC Exercise 11 - Problem 1 & 2
#Mixture Model with Two Gaussian Components

rm(list=ls())

#================================================================
# Model: ========================================================
#================================================================
#	Model:  
#	g(x|theta) = w N(x|mu1,sig.sq) + (1-w) N(x|mu2,sig.sq)

#	Priors:
#	lambda = 1/sig.sq ~ Ga(1,1)
#	mu1 ~ N(0,100)
# 	mu2 ~ N(0,100)
#	w ~ Beta(1,1)	

#================================================================
# Gibbs Sampler =================================================
#================================================================

gibbs = function(y,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for two-component gaussian mixture model.
	#-------------------------------------------------------------
	#MODEL:		g(x|theta) = w N(x|mu1,sig.sq) + (1-w) N(x|mu2,sig.sq)
	#				lambda = 1/sig.sq ~ Ga(a,b) = Ga(1,1)
	#				mu1 ~ N(m1,v1) = N(0,100)
	# 				mu2 ~ N(m2,v2) = N(0,100)
	#				w ~ Beta(c,d) = Beta(1,1)
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of observed data.
	#			hyperparams = vector (m1,m2,v1,v2,a,b,c,d)
	#-------------------------------------------------------------
	#OUTPUTS:	mu1 = vector of posterior mu1 samples.
	#			mu2 = vector of posterior mu2 samples.
	#			lambda = vector of posterior lambda samples.
	#			w = vector of posterior w samples.
	#			d = matrix of d assignments.  
	#				d[i,k] = d value for obs i, gibbs sample k.
	#			y.pred = vector of draws from predictive distribution.
	#-------------------------------------------------------------
	### Data setup.
	
	#Hyperparameters: Prior means and variances for mu1 and mu2.
	mu0 = c(0,0)
	v0 = c(100,100)
	
	#Hyperparameters for lambda = 1/sig.sq ~ Ga(a,b).
	a=1
	b=1
	
	#Hyperparameters for w ~ Beta(c,d).
	c=1
	d=1
	
	n = length(y)	#Total observations.
	M = 2			#Number of components.
	
	#Set up structures to hold parameters.
	mu.j 	= array(NA, dim=c(M,iter))	#Each col is a mu.j.
	lambda 	= array(NA,dim=c(iter))		#Holds precisions.
	w.j 	= array(NA,dim=c(iter))		#Holds weights.  
	d.i 	= array(NA,dim=c(n,iter))	#Holds component memberships.
	y.pred	= array(NA,dim=c(iter))		#For posterior predictive draws.

	#-------------------------------------------------------------
	### Initialize first iteration values.
	mu.j[,1] = 0
	lambda[1] = .1
	w.j[1] = .5

	#Generating random 1's and 2's for d1...dn.
	d.i[,1] = ifelse(rbinom(n,1,w.j[1])==1,1,2) 	

	#Initialize first y value from predictive distribution.
	component = sample(1:M,1,prob=c(w.j[1],1-w.j[1]))
	y.pred[1] = rnorm(1,mu.j[component,1], 1 / sqrt(lambda[1]))
	
	#-------------------------------------------------------------
	#Gibbs sampler.
	for (k in 2:iter){ 
		
		### Progress report.
		if(k %% 100 == 0) {
			print(sprintf('Iteration %i out of %i...',k,iter))
		}
		
		#---------------------------------------------------------
		### Update lambda (1/sig.sq).
		
		#\sum_{i=1}^{n}(y.i - mu_{d.i})^2
		mu.di = mu.j[,k-1][d.i[,k-1]]
		SS = sum((y - mu.di)^2)
		
		lambda[[k]] = rgamma(1, n + a, b + SS/2)
		
		#---------------------------------------------------------
		### Update mu.j (normal-normal).
		nj = as.numeric(table(factor(d.i[,k-1],levels=1:M)))
		
		ybar = aggregate(y,list(d.i[,k-1]),mean,drop=F)$x

		#Fill in empty ybars with draws from prior.
		if(nj[1]==0){ 
			new.mu = rnorm(1,mu0[1],var[1])
			ybar = c(new.mu, ybar) 
		}
		if(nj[2]==0){ 
			new.mu = rnorm(1,mu0[2],var[2])
			ybar = c(ybar, new.mu) 
		}
		
		var = 1 / (1/v0 + nj * lambda[k])
		mean = var * ((1/v0)*mu0 + nj * lambda[k] * ybar)
		
		mu.j[,k] = rnorm(M,mean,sqrt(var))
		
		#---------------------------------------------------------
		### Update w.
		w.j[k] = rbeta(1,nj[1]+1,nj[2]+1)
		
		#---------------------------------------------------------
		### Update d.i and P(d.i=j).
		probs.d = w.j[k] * dnorm(y,mu.j[1,k], 1/sqrt(lambda[k])) / 
			(w.j[k] * dnorm(y,mu.j[1,k], 1/sqrt(lambda[k])) +
			(1-w.j[k]) * dnorm(y,mu.j[2,k], 1/sqrt(lambda[k])))
		
		
		#num = w.j[k] * dnorm(y,mu.j[1,k],1/sqrt(lambda[k]))
		#denom = w.j[k] * dnorm(y,mu.j[1,k],1/sqrt(lambda[k])) + (1-w.j[k]) * dnorm(y,mu.j[2,k],1/sqrt(lambda[k]))
		#probs.d = num / denom
		d.i[,k] = ifelse(rbinom(n,1,probs.d)==1,1,2) 	
		
		#Generate a y value from posterior predictive, using currently updated weight.
		component = sample(1:M,1,prob=c(w.j[k],1-w.j[k]))
		y.pred[k] = rnorm(1,mu.j[component,k], 1 / sqrt(lambda[k]))
	}
	
	#Burn beginning observations.
	if (burn > 0){
		mu.j = mu.j[,-burn]
		w.j = w.j[-burn]
		lambda = lambda[-burn]
		d.i = d.i[,-burn]
		y.pred = y.pred[-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		keep.idx = seq(1,length(lambda),by=thin)
		mu.j = mu.j[,keep.idx]
		w.j = w.j[keep.idx]
		lambda = lambda[keep.idx]
		d.i = d.i[,keep.idx]
		y.pred = y.pred[keep.idx]
	}
	
	#Return results.
	return(list(mu.j = mu.j, w.j = w.j, lambda = lambda, d.i = d.i, y.pred = y.pred))
}

#================================================================
# Generate Data & Run Sampler ===================================
#================================================================

#Generate data of size n=50, with 30 obs from N(-1,4) and 
#remaining 20 with mean (1,4).
y = c(rnorm(300,-1,2),rnorm(200,1,2))

output = gibbs(y,iter=1100,burn=100,thin=2)

#================================================================
# Investigate Results ===========================================
#================================================================

### Plot predictive.
par(mfrow=c(1,2))
hist(y,main='Posterior Predictive (blue) on Data Hist',breaks=50,freq=F)
lines(density(output$y.pred),lwd=2,col='blue')
hist(output$y.pred,breaks=20,freq=F,main='Histogram of Predictive')

### Check and see if the d values can be used to tell which obs came from each component.
true.d = c(rep(1,30),rep(2,20))
est.d = round(rowMeans(output$d.i),2)
cbind(true.d,est.d)

#Plot d values versus true d values.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-11/Figures/MCMC_Ex11_Prob1_d.pdf')
plot(true.d,col='blue',pch=19,main='True d1...dn (blue) vs Estimated')
points(est.d,col='red',pch=19)
dev.off()

#There are a few ds where we are really in the middle and don't know.
#Also, the ds with decisive 1 or 0 were flipped from the true.
#Components got switched.

#Histograms of posterior parameters (reasonableness checking).
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-11/Figures/MCMC_Ex11_Prob1_hist.pdf')
par(mfrow=c(2,2))
hist(output$mu.j[1,],xlab='True mu1: -1')				
hist(output$mu.j[2,],xlab='True mu2: 1')
hist(output$w,xlab='True w: .6')
hist(output$lambda,xlab='True lambda: .01')
dev.off()



#	Parameter estimates are generally poor. w was not bad.  Others have some issues:
#	1. mu1 and mu2 components have switched.  Like in lecture; cannot say which component is really which.
#	2. mu1 and mu2 estimates are poor, aside from flipping.  
#	3. lambda estimate is also poor.  
#			*More likely an issue due to Gamma prior, which weights away from zero values.
#			Choosing a Jeffreys Prior for sig.sq instead of IG may resolve this.

#================================================================
# What happens when change mu1 to -3, and mu2 to 3? =============
#================================================================

#Generate data of size n=50, with 30 obs from N(-3,4) and 
#remaining 20 with mean (3,4).
y2 = c(rnorm(30,-3,2),rnorm(20,3,2))
true.w = 30/50

output = gibbs(y,hyperparams=c(0,0,100,100,1,1,1,1),iter=11000,burn=1000,thin=2)


### Plot predictive.
hist(y,main='Posterior Predictive',breaks=20,freq=F)
lines(density(output$y.pred),lwd=2,col='blue')

# Plot parameters.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-11/Figures/MCMC_Ex11_Prob1_hist2.pdf')
par(mfrow=c(2,2))
hist(output2$mu1,xlab='True mu1: -3')				
hist(output2$mu2,xlab='True mu2: 3')
hist(output2$w,xlab='True w: .6')
hist(output2$lambda,xlab='True lambda: .01')
dev.off()

#Estimates of the mu1 and mu2 parameters improved.

#================================================================
# Problem 2:  ===================================================
#================================================================

#One way to do this has been added into the Gibbs Sampler.  Estimate using Gibbs output.
ghat.p.gibbs = mean(output$y.pred)	#Mean of predicted y values.
ghat.p.gibbs

#----------------------------------------------------------------
#Other way is as follows; ghat_p = 1/K * sum_{k=1}^{K} g(x | theta^{K}), 
#where theta^{K} are params at each Gibbs iteration.
ghat.p = 0
K = length(output$mu1)

for (k in 1:K){
	mu = c(output$mu1[k],output$mu2[k])
	sig.sq = sig.sq=c(1/output$lambda[k],1/output$lambda[k])
	w = output$w[k]
	d = output$d[,k]
	ghat.p = ghat.p + full.lhood(d,w,mu,sig.sq)
} 

ghat.p = ghat.p/K
ghat.p

