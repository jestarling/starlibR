# MCMC Exercise 14
# Kalman Filter Basic Example - Univariate Normal Example
# April 2017
# Jennifer Starling

rm(list=ls())

#Generate data.
data = data.gen(T=50, x0=1, sig.sq=1, tau.sq=1)
x = data$x
y = data$y

#NOTE: We only observe the y values.  The x values are latent; they are our target to estimate.

#Run Kalman Filter for t=1 to 50.
output.k = kalman.filter(y, sig.sq, tau.sq, m=0, v=1)

#Run Particle Filter for t=1 to 50.
output.pf = particle.filter(y, sig.sq, tau.sq, x0=1, N = 10000)

#Plot actual x values, with Kalman and Particle Filter estimates overlaid.
#pdf('/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-14/Kalman_vs_Particle.pdf',height=5,width=7)
	par(mar = c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

	plot(1:50, x, pch=19, col='black', main='Particle and Kalman Filter Comparison',xlab='time',ylab='x')
	lines(1:50, output.k$x.est, col='dodgerblue',lwd=2)
	lines(1:50, output.pf$m, col='green3',lwd=2)

	legend('topright', inset=c(-.3,0), legend = c('Kalman','Particle'), 
		   lty=c(1,1), lwd=c(2,2), col=c('dodgerblue','green3'), title='x 	Estimates')
#dev.off()

#================================================================
# Data-Generating Function ======================================
#================================================================

data.gen = function(T,x0,sig.sq,tau.sq){
	#------------------------------------------------------------
	# FUNCTION: Generates a time series of data according to the model.
	#------------------------------------------------------------
	# MODEL:	p(x_{t} | x_{t-1}) ~ N(x_t | x_{t-1} , tau.sq)
	#			p(y_{t} | x_{t}  ) ~ N(y_t | x_{t}   , sig.sq)
	#------------------------------------------------------------
	# INPUTS:	T = number of time points to generate.
	#			x0 = initial starting value for x.
	#			sig.sq = variance for system equation yt|xt.
	#			tau.sq = variance for state equation xt|xt-1.
	#-------------------------------------------------------------
	# OUTPUTS:	x = vector of generated x values.
	#			y = vector of generated y values.
	#-------------------------------------------------------------
	
	x = rep(NA,T)
	y = rep(NA,T)
	
	#Initialize first x and y value.
	x[1] = rnorm(1, x0, sqrt(tau.sq))
	y[1] = rnorm(1, x[1], sqrt(sig.sq))
	
	#Loop through remaining time points.
	for (t in 2:T){
		
		x[t] = rnorm(1, x[t-1], sqrt(tau.sq))
		y[t] = rnorm(1, x[t], sqrt(sig.sq))
	}
	
	return(list(x=x,y=y))
}

#================================================================
# Basic Univariate Gaussian Kalman Filter Function ==============
#================================================================

kalman.filter = function(y,sig.sq,tau.sq,m=0,v=1){
	#------------------------------------------------------------
	# FUNCTION: Basic univariate Kalman filter.
	#------------------------------------------------------------
	# MODEL:	Recursively estimates mu and psi.sq for:
	#			p(x_{t} | y_{1:t}) ~ N(mu_{t},psi.sq_{t}) where
	#			psi.sq_{t} = 1 / (1/sig.sq + 1/(tau.sq + psi.sq_{t-1}))
	#			mu_{t} = psi.sq_{t} * ( y_{t} / sig.sq + mu_{t-1} / (tau.sq + psi.sq_{t-1}))
	#-------------------------------------------------------------
	# INPUTS:	y 		= input vector of y values.
	#			sig.sq 	= variance for system equation yt|xt.
	#			tau.sq 	= variance for state equation xt|xt-1.
	#			m, v 	= prior mean and variance; starting values at t=0.
	#-------------------------------------------------------------
	# OUTPUTS:	x.est 	= predicted values of p(x_{T} | y_{1:T}) at each time t.
	#			mu	 	= recursive estimates of mean of p(x_{t} | y_{1:t}) at each time t.
	#			psi.sq 	= recursive estimates of var  of p(x_{t} | y_{1:t}) at each time t.
	#-------------------------------------------------------------
	
	#Total number of timepoints.
	T = length(y)
	
	#Set up data structures.
	mu = rep(NA,T)		
	psi.sq = rep(NA,T)
	x.est = rep(NA,T)
	
	#Initialize first iteration (using m and v for starting values).
	psi.sq[1] 	= 1 / ( (1/sig.sq) + 1/(tau.sq + v ))
	mu[1] 		= v * ( (1/sig.sq) * y[t] + 1/(tau.sq + v) * m)
	x.est[1] 	= rnorm(1, mu[1], sqrt(psi.sq[1]))
	
	#Loop through recursion.
	for (t in 2:T){
		
		# Update psi.sq.
		psi.sq[t] = 1 / ( (1/sig.sq) + 1/(tau.sq + psi.sq[t-1]) )
	
		# Update mu.
		mu[t] = psi.sq[t] * ( (1/sig.sq) * y[t] + 1/(tau.sq + psi.sq[t-1]) * mu[t-1])
		
		# Draw estimated x value.
		x.est[t] = rnorm(1, mu[t], sqrt(psi.sq[t]))
		
	} #End time loop.
	
	return(list(mu=mu, psi.sq=psi.sq, x.est=x.est))	
}

#================================================================
# Basic Univariate Particle Filter Function =====================
#================================================================

particle.filter = function(y,sig.sq,tau.sq,x0=NULL,N=100){
	#-------------------------------------------------------------
	# FUNCTION: 	Basic particle filter.
	#-------------------------------------------------------------
	# MODEL:		Proposal
	#-------------------------------------------------------------
	# INPUTS:	y 		= input vector of y values.
	#			sig.sq 	= variance for system equation yt|xt.
	#			tau.sq 	= variance for state equation xt|xt-1.
	#			m, v 	= prior mean and variance; starting values at t=0.
	#			N		= Number of iterations for filter.
	#-------------------------------------------------------------
	# OUTPUTS:	x.est 	= predicted values of p(x_{T} | y_{1:T}) at each time t.
	#			mu	 	= recursive estimates of mean of p(x_{t} | y_{1:t}) at each time t.
	#			psi.sq 	= recursive estimates of var  of p(x_{t} | y_{1:t}) at each time t.
	#-------------------------------------------------------------
	
	#Set up data structures.
	T = length(y)						#Total number of timepoints.
	x.pf = matrix(NA, nrow=T+1,ncol=N)		#Each col is a particle.  Rows are time points.
	rownames(x.pf) = paste('t',0:T,sep='')
	colnames(x.pf) = paste('N',1:N,sep='')
	
	m = rep(NA,T)						#Stores mean of the particle filter at each time.
	
	#-------------------------------------------------------------
	# 1. Initialization at t=0.
	#-------------------------------------------------------------
	
	#If x0 = NULL, then p(x0) ~ N(0,1).  If a specific hard-coded value provided, use it.
	#Note: Instead of drawing x0.i from p(x), here x0=1.
	x.pf[1,] = ifelse(is.null(x0), rnorm(N,0,1), x0 )
	
	#-------------------------------------------------------------
	# Iterate through times t=1 to T.
	#-------------------------------------------------------------
	
	for (t in 2:(T+1)){
		
		#-------------------------------------------------------------
		# 2. Importance Sampling Step.
		#-------------------------------------------------------------

		x.pf[t,] = rnorm(N, x.pf[t-1], sqrt(tau.sq))
		
		#Likelihood to calculate weights. p(y_t | x_t) = N(y_t | x_t, sig.sq)
		w.tilde = dnorm(y[t-1], mean = x.pf[t,], sd = sqrt(sig.sq))
		
		#Normalize weights.
		w = w.tilde / sum(w.tilde)
		
		#Compute mean of particle distribution for comparison with Kalman filter.  (Before resampling.)
		m[t-1] = sum(w * x.pf[t,])
		
		#-------------------------------------------------------------
		# 3. Resampling Step.
		#-------------------------------------------------------------

		#Draw indices to resample.
		samp.idx = sample(1:N, size=N, replace=T, prob=w)
		
		#Resample entire time path for each particle, not just x.pf[t,] single time point.
		x.pf = x.pf[,samp.idx]
	
	} #End time loop.
	
	return(list(x.pf=x.pf, m=m))
	
} #End function.
