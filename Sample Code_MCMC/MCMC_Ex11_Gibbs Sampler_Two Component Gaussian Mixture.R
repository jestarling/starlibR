#MCMC Exercise 11 - Problem 1 & 2
#Mixture Model with Two Gaussian Components

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
# Likelihood Function ===========================================
#================================================================

full.lhood = function(d,w,mu,sig.sq){
	#-------------------------------------------------------------
	#FUNCTION: 	Generates full likelihood for two-component Gaussian mixture
	#			with latent indicators d=1, 0 to indicate component for each y.
	#-------------------------------------------------------------
	#INPUTS:	d = vector of indicators for components of y, length n.
	#			w = weight, between 0 and 1.
	#			mu = vector, (mu1,mu2) of means for the two components.
	#			sig.sq = vector, (sig.sq1,sig.sq2) of vars for the two components.
	#-------------------------------------------------------------
	#OUTPUT:	l = value of full likelihood function.
	#-------------------------------------------------------------
	n1 = sum(d==1)
	n2 = sum(d==0)
	
	comp.1 = w * rnorm(n1,mu[1],sqrt(sig.sq[1]))
	comp.2 = (1-w) * rnorm(n2,mu[2],sqrt(sig.sq[2]))
	
	l = prod(c(comp.1,comp.2))
	return(l)
}

#================================================================
# Gibbs Sampler =================================================
#================================================================

gibbs = function(y,hyperparams=c(0,0,100,100,1,1,1,1),iter=11000,burn=1000,thin=2){
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
	#			prob_d = array of posterior prob(d) samples.  
	#				d[i,j,k] access obs i, d level j, gibbs sample k.
	#			d = matrix of d assignments.  
	#				d[i,k] = d value for obs i, gibbs sample k.
	#-------------------------------------------------------------
	
	#Extract hyperparameters.
	m1 = hyperparams[1];	m2 = hyperparams[2];
	v1 = hyperparams[3];	v2 = hyperparams[4];
	a = hyperparams[5];		b = hyperparams[6];
	c = hyperparams[7];		d = hyperparams[8];
	
	n = length(y)	#Total observations.
	s = 2			#Number of components.
	
	#Set up structures to hold parameters.
	mu1 = rep(0,iter)
	mu2 = rep(0,iter)
	lambda = rep(0,iter)
	w = rep(0,iter)
	d = matrix(0,n,iter)
	y.pred = rep(0,iter)	#For posterior predictive, ie 
							#estimating \int g(y|theta)*f(theta|y1...yn) dtheta, with 
							#theta=(mu1,mu2,w,lambda)
	
	#Initialize first iteration values.
	mu1[1] = 0
	mu2[1] = 0
	lambda[1] = 1
	w[1] = .5					#Initial w = .5
	d[,1] = rbinom(n,1,w[1])	#Generating random 1's and 0's for d1...dn.
	
	#Initialize first y value from predictive distribution.
	y.pred[1] = full.lhood(d[,1],w[1],mu=c(mu1[1],mu2[1]),sig.sq=c(1/lambda[1],1/lambda[1]))
	
	#Iterate through sampler.
	for (i in 2:iter){
		
		#Set up observations in each group based on previous d1...dn.
		y1 = y[which(d[,i-1]==1)]
		y2 = y[which(d[,i-1]==0)]
		ybar1 = mean(y1)
		ybar2 = mean(y2)
		n1 = length(y1)
		n2 = length(y2)
		
		#Update mu1.
		var = 1 / (1/v1 + lambda[i-1]*n1)
		mean = var * lambda[i-1]*n1*ybar1
		mu1[i] = rnorm(1,mean,sqrt(var))
		
		#Update mu2.
		var = 1 / (1/v2 + lambda[i-1]*n2)
		mean = var * lambda[i-1]*n2*ybar2
		mu2[i] = rnorm(1,mean,sqrt(var))
		
		#Update w.
		w[i] = rbeta(1,n1+1,n2+1)
		
		#Update lambda.
		rt = .5 * ( sum((y1-mu1[i])^2) + sum((y2-mu2[i])^2)  )
		lambda[i] = rgamma(1,n,rt)
		
		#Update d1...dn.
		num = w[i] * dnorm(y,mu1[i],1/sqrt(lambda[i]))
		denom = w[i] * dnorm(y,mu1[i],1/sqrt(lambda[i])) + (1-w[i]) * dnorm(y,mu2[i],1/sqrt(lambda[i]))
		probs.d = num / denom
		d[,i] = rbinom(n,1,probs.d)
		
		#Generate a y value from predictive.
		y.pred[i] = full.lhood(d[,i],w[i],mu=c(mu1[i],mu2[i]),sig.sq=c(1/lambda[i],1/lambda[i]))

	}
	
	#Burn beginning observations.
	if (burn > 0){
		mu1 = mu1[-burn]
		mu2 = mu2[-burn]
		w = w[-burn]
		lambda = lambda[-burn]
		d = d[,-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		mu1 = mu1[seq(1,length(mu1),by=thin)]
		mu2 = mu2[seq(1,length(mu2),by=thin)]
		w = w[seq(1,length(w),by=thin)]
		lambda = lambda[seq(1,length(lambda),by=thin)]
		d = d[,seq(1,ncol(d),by=thin)]
	}
	
	#Return results.
	return(list(mu1=mu1,mu2=mu2,w=w,lambda=lambda,d=d,y.pred=y.pred))
}

#================================================================
# Generate Data =================================================
#================================================================

#Set hyperparameters.
hyperparameters = c(0,0,100,100,1,1,1,1)

#Generate data of size n=50, with 30 obs from N(-1,4) and 
#remaining 20 with mean (1,4).
y = c(rnorm(30,-1,2),rnorm(20,1,2))
true.w = 30/50

#================================================================
# Run Sampler ===================================================
#================================================================

K = 11000
output = gibbs(y,hyperparams=c(0,0,100,100,1,1,1,1),iter=K,burn=1000,thin=2)

#================================================================
# Investigate Results ===========================================
#================================================================

### Check and see if the d values can be used to tell which obs came from each component.
true.d = c(rep(1,30),rep(0,20))
est.d = round(rowMeans(output$d),2)
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
hist(output$mu1,xlab='True mu1: -1')				
hist(output$mu2,xlab='True mu2: 1')
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

output2 = gibbs(y,hyperparams=c(0,0,100,100,1,1,1,1),iter=11000,burn=1000,thin=2)

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

