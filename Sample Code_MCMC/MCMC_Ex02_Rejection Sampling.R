#MCMC Exercise 2: Rejection Sampling

#Sample from f(x) proportional to: (-logx)^2 * x^3 * (1-x)^2 for 0<x<1.

#Let h(x) = 1/B(3,3) * x^2(1-x)^2, which has form of gamma(shape=3, rate=3).
#Let l(x) = B(3,3) * x(-logx)^2.

#Sample x~h(x) and y~U(0,1) independently.


#The upper bound for l(x) is M = 1/e^2
M = 1/e^2

#Set up lx function for convenience.
lx_fun = function(x){
	return((3/8)*(1-exp(-x))^2)
}

#------------------------------------------------------

#Repeat B times to generate a good sample of fx,
#and to approximate the acceptance probability.

B = 1000		#Number of iterations.
fx = NA			#Empty vector to hold fx realizations.
accepts = 0		#Hold number of accepts.

for (i in 1:1000){	
	M = 3/8				#Upper bound for function lx, found by taking limit as x->inf.
	y = runif(1,0,1)	#Generate a U(0,1) 
	hx = rgamma(1,shape=4,rate=2)	#Sample from gamma dist.
	lx = lx_fun(hx)		#Plug hx into lx.

	bound = lx / M		#Calculate bound for accept/reject.

	#If bound greater than the uniform realization, accept hx as a realization of fx:
	if(y <= bound) {
		fx = c(fx,hx)
		accepts = accepts + 1
	}
}

#Probability of acceptance.
p_accept = sum(accepts)/B
p_accept

#hx is the sample from fx.

#------------------------------------------------------
#PLOTTING:

fx_fun = function(x){
	return(x^3 * exp(-2*x) * (1-exp(-x))^2)
}

hist(fx,probability=T)
lines(x=seq(0,10,by=.01),y=4*fx_fun(seq(0,10,by=.01)),col='blue')
