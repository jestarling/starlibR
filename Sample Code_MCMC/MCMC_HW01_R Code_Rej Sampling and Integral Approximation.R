#MCMC Homework 01
#Jennifer Starling

#####################################
#   Problem 1: REJECTION SAMPLING   #
#####################################
#Sample from f(x) proportional to: x^3 * exp(-2x) * (1-exp(-x))^2.

#Let h(x) = x^3 * exp(-2x), which has form of gamma(shape=4, rate=2).
#Let l(x) = (3/8)*(1-exp(-x))^2.

#Sample x~h(x) and y~U(0,1) independently.

#M is the max for l(x) on the interval (0,1).
xmax = 1/exp(2)

#Set up lx function for convenience.
lx_fun = function(x){
	return(beta(3,3) * x * (-log(x))^2)
}

M = lx_fun(xmax)

#------------------------------------------------------
#Rejection sampling function.

rej_sample = function(N){
	
	y = runif(N,0,1)	#Generate N U(0,1) realizations.
	hx = rbeta(N,3,3)	#Sample N beta(3,3) realizations.
	lx = lx_fun(hx)		#Plug hx into lx.
	
	# If bound lx/M greater than the uniform realization, 
	# accept hx as a realization of fx:
	bound = lx/M
	accepts = sum(y<=bound)	#Number of acceptances.
	prob_accept = accepts/N	#Probability of acceptance.
	
	fx = hx[y <= bound]
	
	return(list(accepts=accepts,N=N,prob_accept = prob_accept,fx=fx))
}

#Run function.
output = rej_sample(N=10000)

#Probability of acceptance.
output$prob_accept

#hx is the sample from fx.

#------------------------------------------------------
#Repeat B times to generate a good sample of fx,
#and to approximate the acceptance probability.

B = 10000
p_accept_vec = rep(0,B)

for (b in 1:B){
	output = rej_sample(N=10000)
	p_accept_vec[b] = output$prob_accept
}

mean(p_accept_vec)

#------------------------------------------------------
#PLOTTING:

#Histogram of fx, plus sanity check that the function shape looks like f(x).
fx_fun = function(x){
	return( (-log(x))^2 * x^3 * (1-x)^2)
}

hist(output$fx,probability=T)
x = seq(0,1,by=.001)
lines(x,100*fx_fun(x),col='blue')

#########################################
#   Problem 2: INTEGRAL APPROXIMATION   #
#########################################


#Function for monte carlo approximation of integral.  
# INPUT: Number of samples to use in integral approximation.
# OUTPUT: Returns Ihat_N.
monte_carlo_integral = function(N){
	#N = number of monte carlo samples to approximate the integral.

	#Take N samples from the Beta(2,3/2) distribution, ie f(xi).
	fx = rbeta(N,2,7/2)

	#Calculate g(xi) for each sample as g(x) = beta(2,3/2) * x^2 * (-log(x))^2.
	gx = beta(2,7/2) * fx^2 * (-log(fx))^2

	#Ihat_N = 1/N * sum(gx)
	Ihat = (1/N) * sum(gx)
	return(Ihat)
}

#Call function and display results. This is for a single estimate of Ihat_N.
I_hat_N = monte_carlo_integral(N=1000)
I_hat_N	

#-----------------------------------------------------
#Graphical demonstration that method has worked.  (Try for vector of N values.)
N = seq(1,500000,by=100)	#Vector of N values to use in approximating each Ihat.
Ihat_vec = rep(0,length(N))		#Empty vector to hold ihat values.

for (b in 1:length(N)){
	Ihat_vec[b] = monte_carlo_integral(N[b])
}

plot(N,Ihat_vec,type='l',main='Monte Carlo Ihat for increasing N sample sizes',
	xlab='N',ylab='Ihat_N')
	
#-----------------------------------------------------
#For a fixed N, say N=1000, estimate variance of Ihat_N.
N = 1000
B = 10000
Ihat_vec = rep(0,B)

#Calculate B estimates of Ihat.
for (b in B){
	Ihat_vec[b] = monte_carlo_integral(N)
}	

#Estimate variance by calculating variance of Ihat estimates.
var(Ihat_vec)

###################################################################
#   Problem 3: Rejection Sampling, Std Normal w Cauchy Proposal   #
###################################################################

#Rejection sample from standard normal f(x) = (2pi)^(-.5) * exp(-.5x^2)
#Using h(x) = standard cauchy, so h(x) = 1 / pi*(1+x^2)
#Note: Sample x~h(x) and y~U(0,1) independently.

#M is the max for l(x) on the interval (0,1).
M = sqrt(2*pi/exp(1))

#Set up lx function for convenience.
lx_fun = function(x){
	return( (2*pi)^(-1/2) * pi * (1+x^2) * exp((-1/2) * x^2) )
}

#------------------------------------------------------
#Rejection sampling function.

rej_sample = function(N){
	
	y = runif(N,0,1)	#Generate N U(0,1) realizations.
	hx = rcauchy(N,0,1)	#Sample N Cauchy(0,1) realizations.
	lx = lx_fun(hx)		#Plug hx into lx.
	
	# If bound lx/M greater than the uniform realization, 
	# accept hx as a realization of fx:
	bound = lx/M
	accepts = sum(y<=bound)	#Number of acceptances.
	prob_accept = accepts/N	#Probability of acceptance.
	
	fx = hx[y <= bound]
	
	return(list(accepts=accepts,N=N,prob_accept = prob_accept,fx=fx))
}

#Run function.
output = rej_sample(N=10000)	#With p(accept) = .6ish, taking 10K samples
								#to make sure we have plenty of accepts to plot.

#Probability of acceptance.
output$prob_accept

#hx is the sample from fx.

#------------------------------------------------------
#Repeat B times to generate a good sample of fx,
#and to approximate the acceptance probability.

B = 10000
p_accept_vec = rep(0,B)

for (b in 1:B){
	output = rej_sample(N=10000)
	p_accept_vec[b] = output$prob_accept
}

mean(p_accept_vec)

#------------------------------------------------------
#PLOTTING:

#Randomly select 1000 of the accepted samples generated above.  
#Plot histogram.
fx = output$fx[1:1000]

hist(fx,main='Histogram of fx for 1000 accepted samples',freq=F)

#Sanity check that the function shape looks like f(x).
fx_fun = function(x){
	return( (2*pi)^(-1/2) * exp((-1/2) * x^2))
}

x = seq(-3,3,by=.001)
lines(x,fx_fun(x),col='blue')