#MCMC Homework 2
#Problem 2
#Feb 9, 2017
#Jennifer Starling

#Includes code for two ways of doing problem: 
#Uniform (1/2) stationary density, and bernoulli(q).

#------------------------------------------------------
#Uniform 1/2 stationary density.

#Desired chain length.
n = 1000

#Chain to hold realizations. Initialize first value P(X0=1) = .5
x = rbinom(1,1,.5)

#Transition probability matrix.
P = matrix(c(2/3,1/3,1/3,2/3),byrow=T,nrow=2)

#Markov chain.
for (i in 2:n){
	x0 = x[i-1]
	
	#X_{n+1} can be 0 or 1.  Draw a random uniform
	#based on probability of X_{n+1}
	tp = P[x0+1,2]
	xnew = rbinom(1,1,tp)
	
	#Select x1 with highest transition prob.
	#Append to new vector.
	#Note: possible x values corresond to P indices - 1.
	x = c(x,xnew)
}

#Plot resulting Markov Chain to show stationary.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 02/Figures/2_Markov_Chain.pdf')
par(mfrow=c(2,1))
plot(x,type='l',main='Markov Chain Trace')
plot(x,main='Monte Carlo Chain Sample')
dev.off()