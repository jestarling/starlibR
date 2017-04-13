#MCMC Homework 4
#Mixture Model with Two Gaussian Components
#Jennifer Starling
#March 2017

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
#	w ~ U(0,1)	

gx = function(di,mu,sig.sq){
	#-------------------------------------------------------------
	#FUNCTION: 	Generates a value of the function g(x|theta)
	#			where g(x|theta) = w*N(x|mu1,sig.sq) + (1-w)N(x|mu2,sig.sq)
	#-------------------------------------------------------------
	#INPUTS:	d = indicator for which component to use.
	#			mu = vector, (mu1,mu2) of means for the two components.
	#			sig.sq = vector, (sig.sq1,sig.sq2) of vars for the two components.
	#-------------------------------------------------------------
	#OUTPUT:	l = value of full likelihood function.
	#-------------------------------------------------------------
	gx = ifelse(di==1,
		rnorm(1,mu[1],sqrt(sig.sq[1])),
		rnorm(1,mu[2],sqrt(sig.sq[2]))
		)	
	
	return(gx)
}

#================================================================
# Gibbs Sampler =================================================
#================================================================

gibbs = function(y,m1,m2,v1,v2,a,b,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for two-component gaussian mixture model.
	#-------------------------------------------------------------
	#MODEL:		g(x|theta) = w N(x|mu1,sig.sq) + (1-w) N(x|mu2,sig.sq)
	#				lambda = 1/sig.sq ~ Ga(a,b) = Ga(1,1)
	#				mu1 ~ N(m1,v1) = N(0,100)
	# 				mu2 ~ N(m2,v2) = N(0,100)
	#				w ~ U(0,1)
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of observed data.
	#			hyperparams = {m1,m2,v1,v2,a,b}
	#-------------------------------------------------------------
	#OUTPUTS:	mu1 = vector of posterior mu1 samples.
	#			mu2 = vector of posterior mu2 samples.
	#			lambda = vector of posterior lambda samples.
	#			w = vector of posterior w samples.
	#			d = matrix of d assignments.  
	#				d[i,k] = d value for obs i, gibbs sample k.
	#			y.pred = vector of draws from predictive distribution.
	#-------------------------------------------------------------
	
	n = length(y)	#Total observations.
	
	#Set up structures to hold parameters.
	mu1 = rep(0,iter)
	mu2 = rep(0,iter)
	lambda = rep(0,iter)
	w = rep(0,iter)
	d = matrix(0,n,iter)
	d.pred = rep(0,iter)
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
	d.pred[1] = rbinom(1,1,w[1])
	y.pred[1] = gx(d.pred[1],mu=c(mu1[1],mu2[1]),sig.sq=c(1/lambda[1],1/lambda[1]))
	
	#Iterate through sampler.
	for (i in 2:iter){
		
		#Update d1...dn.
		prob.vec = w[i-1] * dnorm(y,mu1[i-1],1/sqrt(lambda[i-1])) / 
			(w[i-1] * dnorm(y,mu1[i-1],1/sqrt(lambda[i-1])) + (1-w[i-1]) * dnorm(y,mu2[i-1],1/sqrt(lambda[i-1])))	
		d[,i] = rbinom(n,rep(1,n),prob.vec)
		
		#Update observations in each group based on new d1...dn.
		y1 = y[which(d[,i]==1)]
		y2 = y[which(d[,i]==0)]
		ybar1 = mean(y1)
		ybar2 = mean(y2)
		n1 = length(y1)
		n2 = length(y2)
		
		#If there is no data in one of the groups, sample is from posterior only.
		#Since which() generates numeric() for empty set, handle as follows.
		#If no data from the posterior for a component, sample is prior only.
		if(length(y1)==0){
			n1=0
			ybar1=0
		}
		if(length(y2)==0){
			n2=0
			ybar2=0
		}
		
		#print(paste("n1 = ",n1))
		#print(paste("n2 = ",n2))
		
		#Update mu1.
		var = 1 / (1/v1 + n1*lambda[i-1])
		mean = var * ((1/v1)*m1 + n1*lambda[i-1]*ybar1)
		mu1[i] = rnorm(1,mean,sqrt(var))
		
		#Update mu2.
		var = 1 / (1/v2 + n2*lambda[i-1])
		mean = var * ((1/v2)*m2 + n2*lambda[i-1]*ybar2)
		mu2[i] = rnorm(1,mean,sqrt(var))
		
		#Update w.
		w[i] = rbeta(1,n1+1,n2+1)
		
		#Update lambda.
		RSS1 = sum((y1-mu1[i])^2)
		RSS2 = sum((y2-mu2[i])^2)
		lambda[i] = rgamma(1,n/2, RSS1/2 + RSS2/2 + 1) 
				
		#Generate a y value from posterior predictive, using currently updated weight.
		d.pred[i] = rbinom(1,1,w[i])
		y.pred[i] = gx(d.pred[i],mu=c(mu1[i],mu2[i]),sig.sq=c(1/lambda[i],1/lambda[i]))
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
# 1. Generate Data ==============================================
#================================================================

#Generate data of size n=100 from N(0,1)
y = rnorm(100,0,1)

#================================================================
# 1. Run Sampler ================================================
#================================================================

K = 11000
m1 = m2 = 0
v1 = v2 = 100
a = b = 1
c = d = 1

output = gibbs(y,m1,m2,v1,v2,a,b,iter=K,burn=1000,thin=2)

#================================================================
# 2. Plot predictive density & some of posteriors ===============
#================================================================

#Histogram of posterior predictive density.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 04/Figures/2_Posterior_Pred_hist.pdf')
hist(output$y.pred,breaks=50,freq=F,main='Posterior Predictive Density',xlim=c(-5,5))
dev.off()

#Histogram of posterior mu1, mu2, sigma and w.
pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 04/Figures/3_Posterior_hist.pdf')
par(mfrow=c(2,2))
hist(output$mu1,breaks=50,freq=F,main='Posterior of mu1',xlim=c(-5,5))
hist(output$mu2,breaks=50,freq=F,main='Posterior of mu1',xlim=c(-5,5))
hist(1/output$lambda,breaks=50,freq=F,main='Posterior of sig.sq')
hist(output$w,breaks=50,freq=F,main='Posterior of w')
dev.off()

#================================================================
# 3. Plot di values =============================================
#================================================================

post.di.means = colMeans(output$d)

pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 04/Figures/4_di.pdf')
hist(post.di.means,main='Posterior Mean of d1...dn values')
dev.off()

#================================================================
# 4. Evaluate the integral specified ============================
#================================================================

#Evaluate I = \int x * g_p(x) dx

Ihat = (1/K) * sum(output$y.pred)
Ihat
# 0.00901856