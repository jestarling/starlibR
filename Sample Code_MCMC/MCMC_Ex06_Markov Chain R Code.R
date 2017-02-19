#MCMC Exercise 06
#Jennifer Starling
#Feb 16, 2017


#####################################################
###   Problem 1 Simulations: (Gamma stationary)   ###
#####################################################

#----------------------------------------------------------------
### Sample from stationary Gamma(a,1).  x>0,a>0.

fx = 1	#Vector to hold sample from f(x).
n=1000	#Sample size.
a=5		#Gamma parameter for Gamma(a,1).

for (i in 2:n){
	x0 = fx[i-1]
	
	# Not the V, W from iid U(0,1) way:
	u1 = runif(1,0,x0^(a-1)/gamma(a))
	u2 = runif(1,0,exp(-x0))
	
	x1 =  exp(log(u1*gamma(a))/(a-1)) 	#lower bound for x_{n+1}
	x2 = -log(u2)						#upper bound for x_{n+1}
	
	xnew = runif(1,x1,x2)				#sample new x_{n+1}
	fx = c(fx,xnew)
}

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_gamma.pdf')
hist(fx,freq=F,main='Gamma(a,1) Markov Chain',xlab='fx')
x = seq(.01,max(fx),by=.001)
lines(x,(1/gamma(a)) * x^(a-1) * exp(-x),col='firebrick3',lwd=2)
dev.off()

#Analytical & Sample AR-1 autocorr.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_gamma_cov.pdf')
acf(fx,main=paste('Corr for fx chain: AR1 = 1 - 1/2a, a=',a))
dev.off()
1-1/(2*a)

#Show autocorrelation.
acf(fx)
acf(fx,plot=F)

#Analytical AR-1 autocorr.
1-1/(2*a)

#----------------------------------------------------------------
### Sample from stationary density from Q1: Gamma(a,1).
###	Done the second way, with V1, W1, V2,W2 ~ U(0,1) and X_{n+1} = 

fx = 1	#Vector to hold sample from f(x).
n=1000	#Sample size.
a=2		#Gamma parameter for Gamma(a,1).

for (i in 2:n){
	x0 = fx[i-1]
	
	# Using V,W,Z ~ U(0,1)
	w = runif(1,0,1)
	y = runif(1,0,1)
	z = runif(1,0,1)
	
	xnew = x0 * w^(1/(a-1)) + y*(x0 - log(z) - x0 * w^(1/(a)))
	
	fx = c(fx,xnew)
}

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_gamma.pdf')
hist(fx,freq=F,main='Gamma(a,1) Markov Chain',xlab='fx')
x = seq(.01,max(fx),by=.001)
lines(x,(1/gamma(a)) * x^(a-1) * exp(-x),col='firebrick3',lwd=2)
dev.off()

#Analytical & Sample AR-1 autocorr.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_gamma_cov.pdf')
acf(fx,main=paste('Corr for fx chain: AR1 = 1 - 1/2a, a=',a))
dev.off()
1-1/(2*a)

###########################################################
###   Problem 2 Simulations: (Exponential stationary)   ###
###########################################################

#----------------------------------------------------------------
### Sample from stationary density from Q2: f(x) = exp(-x), x>0.
fx = 1
n=10000

for (i in 2:n){
	x0 = fx[i-1]
	u = runif(1,0,exp(-x0))
	x1 = runif(1,0,-log(u))
	fx = c(fx,x1)
}

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_exponential1.pdf')
hist(fx,freq=F,main='Exponential(1) Markov Chain',xlab='fx')
x = seq(0,max(fx),by=.001)
lines(x,exp(-x),col='firebrick3',lwd=2)
dev.off()

#Show autocorrelation.
acf(fx)
acf(fx,plot=F)
1/2	#Analytical expression for exponential(1) autocorr lag 1.



#----------------------------------------------------------------
### Sample from stationary density from Q2: f(x) = exp(-x), x>0.
###	Done the second way, with V, W ~ U(0,1) and X_{n+1} = W * (X_n - log(V)).
fx = 1
n = 10000

for (i in 2:n){
	x0 = fx[i-1]
	v = runif(1,0,1)
	w = runif(1,0,1)
	x1 = w * (x0 - log(v))
	fx = c(fx,x1)
}

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-06/Figures/mcmc_exponential1.pdf')
hist(fx,freq=F,main='Exponential(1) Markov Chain',xlab='fx')
x = seq(0,max(fx),by=.001)
lines(x,exp(-x),col='firebrick3',lwd=2)
dev.off()

#Show autocorrelation.
acf(fx)
acf(fx,plot=F)

