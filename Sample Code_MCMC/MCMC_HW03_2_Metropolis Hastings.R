#MCMC Homework 3
#Problem 2
#Feb 23, 2017
#Jennifer Starling

#Want to sample from pi(x) \propto exp(-.5x^2) / (1+x^2).
#using Metropolis Hastings with q(x'|x) = N(x'|x,sig2) with
#some sigma as proposal.

#------------------------------------------------------
### Simulate what happens if sigma too large/too small.

target = function(y){
	exp(-.5*y^2)/(1+y^2)
}

q = function(y,sigma){ #Proposal
	rnorm(1,y,sigma)
}

mh = function(target,q,n,sig){
	#------------------------------------------------------
	#FUNCTION: 	Performs Metropolis Hastings with the given
	#			q with length of chain = n.
	#------------------------------------------------------
	#INPUTS:	n = length of chain.
	#			sig = sd of normal proposal q.
	#			center = 'pos.x' to center q at x, else center at -x.
	#------------------------------------------------------
	#OUTPUTS:	x = sample from target density.
	#------------------------------------------------------
	x = 0			#Initialize chain to 1.
	p_accept = 0	#Initialize prob of acceptance.
	#Iterate chain.
	for (i in 2:n){
		#Save previous x_n.
		x.old = x[i-1]
		
		#Sample from uniform and q.
		u = runif(1,0,1)
		x.new = q(x.old,sig)
		
		#Calculate alpha.
		alpha = min(1,target(x.new)/target(x.old))
		
		#Acceptance check
		if(u < alpha){
			x = c(x,x.new)
			p_accept = p_accept + 1
		} else{
			x = c(x,x.old)
		}
	}
	print(paste('Probability of acceptance: ',p_accept / n,sep=''))
	return(x)
}

#------------------------------------------------------
### Try three possible sigmas, two which are extreme.
n = 10000
sig = c(.05,1,50)
fx = matrix(nrow=n,ncol=length(sig))	#To hold results for each sigma.
p_accepts = rep(0,length(sig))

for (s in 1:length(sig)){
	fx[,s] = mh(target,q,n,sig[s])
}

#------------------------------------------------------
### Plot resulting MH histograms and traces.

x = seq(-10,10,by=.01)
f = target(x)

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/2_Metropolis_Hastings.pdf',height=12,width=12)
par(mfrow=c(3,2))
for (i in 1:3){
	#Plot histogram for each sigma.
	hist(fx[,i],freq=F,main=paste('Sigma = ',sig[i],sep=''),breaks=30)
	lines(x,.6*f,lwd=2,col='red')
	
	#Plot trace for each sigma.
	plot(fx[,i],type='l',main=paste('Sigma = ',sig[i],sep=''))
	
}
dev.off()

#------------------------------------------------------
#A suitable sigma looks to be sig=1.  Use it to evaluate the integral (posterior mean).
fx = mh(target,q,n=10000,sig=1)
mean(fx)  #Posterior mean (evaluate the integral).

x = seq(-10,10,by=.01)
f = target(x)

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/2_Metropolis_Hastings_sig1.pdf')
par(mfrow=c(2,1))
hist(fx,freq=F,main=paste('Sigma = ',1,' Posterior Mean = ',round(mean(fx),4),sep=''),breaks=30)
lines(x,.6*f,lwd=2,col='red')
plot(fx,type='l',main=paste('Sigma = ',1,sep=''),xlab='Iterations',ylab='fx')
dev.off()

#------------------------------------------------------
#What happens if you use q(x'|x) = N(x'|-x,sig2) instead?

q = function(y,sigma){ #Proposal
	rnorm(1,-y,sigma)
}

### Try three possible sigmas, two which are extreme.
n = 10000
sig = c(.05,1,50)
fx = matrix(nrow=n,ncol=length(sig))	#To hold results for each sigma.
p_accepts = rep(0,length(sig))

for (s in 1:length(sig)){
	fx[,s] = mh(target,q,n,sig[s])
}

### Plot resulting MH histograms and traces.

x = seq(-10,10,by=.01)
f = target(x)

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/2_Metropolis_Hastings_NEGX.pdf',height=12,width=12)
par(mfrow=c(3,2))
for (i in 1:3){
	#Plot histogram for each sigma.
	hist(fx[,i],freq=F,main=paste('Sigma = ',sig[i],sep=''),breaks=30)
	lines(x,.6*f,lwd=2,col='red')
	
	#Plot trace for each sigma.
	plot(fx[,i],type='l',main=paste('Sigma = ',sig[i],sep=''))
	
}
dev.off()





