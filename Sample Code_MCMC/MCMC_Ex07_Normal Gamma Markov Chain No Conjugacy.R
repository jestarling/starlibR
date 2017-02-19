#MCMC Exercise 7
#Markov Chain for Normal-Gamma wo Nice Conjugacy
#Feb 19, 2017
#Jennifer Starling

### MODEL:
###	Likelihood: 	x|mu,lambda ~ Normal(mu,lambda) with lambda=precision.
###	Priors:			mu ~ N(0,k)	 with k=1/tau^2, also precision.
###					lambda ~ gamma(a,b)

#Length for chain.
B=11000

#Sample size n
n = 1000

#Generate some data using likelihood function
true.lambda = .5 
true.mu = 3
x = rnorm(n,true.mu,sqrt(1/true.lambda))

xbar = mean(x)
Sx = sum((x-xbar)^2)

#Set hyperparameters tau2,a,b.
#tau2 = .5
k = .05
a=5
b=5

#Predictors to hold mu,lambda. (Initialize chain.)
mu = 1
lambda = 1

#Run chain.
for (i in 2:B){
	
	#Update mu:
	 mu[i] = rnorm(1, (lambda[i-1]*n*xbar)/(lambda[i-1]*n + 2*k), 1/sqrt(lambda[i-1]*n+2*k) )
	
	#Update lambda
	lambda[i] = rgamma(1,(n+2*a)/2, .5 * (n*xbar^2 - 2*n*xbar*mu[i] + n*mu[i]^2 + Sx + 2*b))
}

#Burn in period of 1k.
mu = mu[1001:B]
lambda = lambda[1001:B]

#Thinning every 2nd obs to decrease autocorr.
mu = mu[seq(1,length(mu),by=2)]
lambda = lambda[seq(1,length(lambda),by=2)]

#Plot and print plot.
par(mfrow=c(2,1))
hist(mu)
hist(lambda)

#Output plot.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-07/Figures/Prob3_MarkovChain.pdf')
par(mfrow=c(2,1))
hist(mu)
hist(lambda)
dev.off()

#Posterior means.  Working - recovering original lambda and mu used to generate data.
mean(mu)
true.mu
mean(lambda)
true.lambda

