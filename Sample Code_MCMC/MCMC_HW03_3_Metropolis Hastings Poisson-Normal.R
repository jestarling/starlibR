#MCMC Homework 3
#Problem 3
#Feb 23, 2017
#Jennifer Starling

#Want to sample posterior f(theta|data) \sim exp(theta*a) * exp(-n*exp(theta)) * exp(-.5*theta^2)
#This is the Poisson likelihood model with mean exp(theta), with a N(0,1) prior on theta.  (Centered on prev val.)
#Note: a>0, n = integer, -inf < theta < inf. 
#Find a Markov chain for sampling from f and implement it.

#------------------------------------------------------

target = function(theta,a,n){
	(exp(theta*a) * exp(-n*exp(theta)) * exp(-.5*theta^2)) 
}

q = function(theta){
	rnorm(1,theta,1)
}

mh = function(target,q,B,a,n,thin=NA,burn=NA){
	#FUNCTION: Performs Metropolis Hastings with the given
	#q with length of chain = B.
	#INPUTS:	B = length of chain.
	#			a = value greater than 0.
	#			n = integer value.
	
	theta = 0	#Initialize chain to 1.
	
	#Iterate chain.
	for (i in 2:B){
		#Save previous theta_n.
		theta.old = theta[i-1]
		
		#Sample from uniform and q.
		u = runif(1,0,1)
		theta.new = q(theta.old)
		
		#Calculate alpha.
		alpha = target(theta.new,a,n) / target(theta.old,a,n)
		
		#Acceptance check
		if(u < alpha){
			theta = c(theta,theta.new)
		} else{
			theta = c(theta,theta.old)
		}
	}
	
	#Thin and burn.
	if(sum(is.na(thin),is.na(burn))==0){
		theta = theta[-c(1:burn)]
		theta = theta[seq(1,length(theta),by=thin)]
	}
	
	return(theta)
}

#------------------------------------------------------
#Plot histogram of theta versus target function.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/3_Pois.pdf',width=12,height=6)
n=10
a=.5
par(mfrow=c(1,2))
theta= mh(target,q,B=20000,a=a,n=n)  	
hist(theta,freq=F,main=paste('Hist of theta, n= ',n,', a= ',a,sep=''),breaks=50)

x=seq(-10,10,by=.01)
f = target(x,a,n)
lines(x,42*f,type='l',col='red',lwd=2)

plot(theta,type='l',main='Trace of theta')
dev.off()

#------------------------------------------------------
### Try three combos of a and n to illustrate chain mixing issues.
theta= mh(target,q,B=5000,a=.05,n=1)  	#Best mixing.
theta1 = mh(target,q,B=5000,a=.5,n=5) 	#Better mixing but still not great.
theta2 = mh(target,q,B=5000,a=5,n=100)  	#Very chunky chain, gets stuck.

#Plot resulting Markov Chain to show stationary.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/3_Pois_Markov_Chain.pdf')
par(mfrow=c(3,1))
plot(theta,type='l',main='Chain Trace: theta, a=.05,n=1')
plot(theta1,type='l',main='Chain Trace: theta, a=.5,n=5')
plot(theta2,type='l',main='Chain Trace: theta, a=5,n=100')
dev.off()








