#MCMC Lecture 12, Example 2
#Feb 28, 2017
#Jennifer Starling

rm(list=ls())

#Want to sample the Poisson density.
#Target: f(x) \propto theta^x / x!, x \in {0,1,2...}

#Proposal q(x'|x) = 
#	\begin{cases}
#		When x >= 1:
#			x' = x + 1 with probablity 1/2.
#			x' = x - 1 with probability 1/2.
#		When x = 0:
#			q(1|0) = 1.
# 	\end{cases}

#This needs special case for alpha(1,0) and alpha(0,1).

target = function(x,theta){
	 theta^x / factorial(x)
}

q = function(x){
	if(x>=1){
		p = rbinom(1,1,.5)
		if (p ==1){
			return(x+1)
		} else{
			return(x-1)
		}
	} else{
		return(1)
	}
}

alpha = function(x.new,x.old,theta){
	
	tgt.new = theta^x.new / factorial(x.new)
	tgt.old = theta^x.old / factorial(x.old)
		
	if(x.new==1 && x.old==0){
		return(min(1,theta/2))
	} else if (x.new==0 && x.old==1){
		return(min(1,2/theta))
	} else{
		return(min(1,tgt.new/tgt.old))
	}
}

mh = function(target,alpha,q,n,theta){
	#------------------------------------------------------
	#FUNCTION: 	Performs Metropolis Hastings with the given
	#			q with length of chain = n.
	#------------------------------------------------------
	#INPUTS:	n = length of chain.
	#			target = target function.
	#			q = proposal function.
	#			theta = Poisson(theta) parameter.
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
		x.new = q(x.old)
		
		#Calculate alpha.
		alph = alpha(x.new,x.old,theta)
		
		#Acceptance check
		if(u < alph){
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
### Run function for a few values of theta.
n = 10000
theta = c(.5,1,10)
fx = matrix(nrow=n,ncol=length(theta))	#To hold results for each sigma.

for (s in 1:length(theta)){
	fx[,s] = mh(target,alpha,q,n,theta[s])
}

#------------------------------------------------------
### Plot resulting MH histograms and traces.

x = seq(0,50,by=1)

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-08a_Lecture12/Figures/Poisson.pdf',height=12,width=12)
par(mfrow=c(3,2))
for (i in 1:3){
	#Plot histogram for each sigma.
	plot(table(fx[,i])/length(fx[,i]))
	f = dpois(x,theta[i])
	lines(x,f,lwd=2,col='red')
	
	#Plot trace for each sigma.
	plot(fx[,i],type='l',main=paste('Theta = ',theta[i],sep=''))
	
}
dev.off()




