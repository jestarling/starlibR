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