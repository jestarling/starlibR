#MCMC - Exercise 3 - Spring 2017

#Sample the gamma density using the ratio of uniform algorithm.

ratio_u_for_gammas = function(n,a,b){
	#Samples the gamma density via the ratio of uniforms.
	#INPUTS:	n = number of gamma realizations to generate.
	#			a = the alpha (shape) parameter
	#			b = the beta (rate) parameter
	#OUTPUTS:	fx = vector of realizations from the gamma density.
	
	#Accept x = u2/u1 as being from fx if u1 <= sqrt(f(u2/u1))
	fx = function(x,a,b) dgamma(x,a,b)
	
	#fx_max1 = First function to maximize; sqrt(f(x)).
	fx_max1 = function(x,a,b) -sqrt(dgamma(x,a,b))
	
	#fx_max2 = Second function to maximize; x*sqrt(f(x))
	fx_max2 = function(x,a,b)  -x * sqrt(dgamma(x,a,b))
	
		# Using negatives bc optim minimizes.  We want to maximize.
	
	#a is max_x for sqrt(f(x))
	a_bound = optim(par=.5,fn=fx_max1,a=a,b=b,lower=0,upper=Inf,method="L-BFGS-B")$par
	b_bound = optim(par=.5,fn=fx_max2,a=a,b=b,lower=0,upper=Inf,method="L-BFGS-B")$par
	
	#Generate n random uniform samples using calculated bounds.
	u1 = runif(n,0,a_bound)
	u2 = runif(n,0,b_bound)
	
	#Value of f(x) for all generated samples.
	fx_samps = u2/u1
	
	#Acceptance test.
	accept = u1 <= sqrt(fx(u2/u1,a,b))
	
	#Accepted samples from fx.
	fx_samps = fx_samps[accept]
	
	#Acceptance probability.
	prob_accept = sum(accept)/n
	
	return(list(fx_samps=fx_samps,prob_accept=prob_accept))
}


#Test function.
output = ratio_u_for_gammas(10000,a=2,b=1)
gamma_samps = output$fx_samps	#Generated gamma samples.
p_accept = output$prob_accept	#P(Acceptance)

#Histogram of sampled values.
hist(gamma_samps,freq=F)

#Overlay true gamma density.
lines(seq(0,10,by=.01),dgamma(seq(0,10,by=.01),a,b),type='l',col='blue') 
