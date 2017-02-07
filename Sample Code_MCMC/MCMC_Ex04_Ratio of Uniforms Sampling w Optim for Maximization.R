#MCMC Exercise 04
#Problem 2
#Jennifer Starling
#Feb 7, 2017

#Sample the following density using ratio of uniforms:
# f(x) /propto x^3 * exp(-2*x) * (1-exp(-2*x))^3

# To sample using ratio of uniforms, need two bounds:
#	a = sqrt(f(x)) at x value which maximizes sqrt of f(x).
#	b = x*sqrt(f(x)) at x value which maximizes x*sqrt(f(x)).

#Analytically intractible to find a, so will use optim to easily find max.

#Set up function fx and obtain bounds a and b.
f = function(x) (x^3 * exp(-2*x) * (1-exp(-2*x))^3)

#Call function and plot results.
output = ratio_of_unifs(f,10000)
f_samp = output$f_sample

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-04/Ratio_of_Unifs.pdf')
hist(f_samp,freq=F,main='Ratio of Uniforms Sample of f(x)',xlab='fx')
lines(seq(0,10,.001),3*f(seq(0,10,.001)),col='blue')
dev.off()

###############################################
###   Ratio of Uniforms Sampling Function   ###
###############################################

ratio_of_unifs = function(f,n){
	#INPUTS: 	f = the function to sample from.
	#			n = number of samples to generate
	#OUTPUTS:
	#			f_sample = sample of realizations from f.
	#			p_accept = acceptance probability.
	#			a = max_x sqrt(fx) 
	#			b = max_x x*sqrt(fx)
	#			n_accepted = number of accepted values.
	#			f = the original function input by the user.
	
	#Set up functions and bounds for sampling from uniform.
	fa = function(x) sqrt(f(x))
	fb = function(x) x*sqrt(f(x))

	a = optim(.1,fa,control=list(fnscale=-1),method='BFGS')$value 	#fnscale means to maximize instead of minimize.
	b = optim(.1,fb,control=list(fnscale=-1),method='BFGS')$value
	
	#Draw two independent uniform samples.
	u1 = runif(n,0,a)
	u2 = runif(n,0,b)
	
	#Calculate acceptance check vector. 
	#Accept if u1 <= sqrt(f(u2/u1))
	accept = u1 <= sqrt(f(u2/u1))
	
	#Keep only accepted values.
	f_sample = (u2/u1)[accept]
	n_accepted = length(f_sample)/n
	p_accept = n_accepted/n
	
	#Return output.
	return(list(f_sample=f_sample, n_accept=n_accepted, p_accept=p_accept,a=a,b=b,f=f))
}