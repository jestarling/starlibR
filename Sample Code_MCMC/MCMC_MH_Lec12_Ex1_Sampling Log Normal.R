#MCMC Lecture 12, Example 2
#Feb 28, 2017
#Jennifer Starling

rm(list=ls())

#Want to sample the target density:
# f(x) \propto x^a * exp(-x), x>0, a>0

#Since support is only on positives for x,
#the normal proposal is problematic.
#Use log-normal instead.

target = function(x,a){
	 x^a * exp(-x)
}

q = function(x,sigma=1){
	#rlnorm(1,x,sigma)		#Built-in.
	exp(rnorm(1,x,sigma)) 	#Helps me remember log-normal.
}

alpha = function(x.new,x.old,a){
	min(1,(x.new/x.old)^(a+1) * exp(-(x.new-x.old)))
}

mh = function(target,alpha,q,n,a,sigma=1){
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
		x.new = q(x.old,sigma)
		
		#Calculate alpha.
		alph = alpha(x.new,x.old,a)
		
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
### Run function for a few values of a.
n = 10000
avals = c(.5,1,10)
fx = matrix(nrow=n,ncol=length(avals))	#To hold results for each sigma.

for (s in 1:length(avals)){
	fx[,s] = mh(target,alpha,q,n,avals[s],sigma=1)
}

#------------------------------------------------------
### Plot resulting MH histograms and traces.

x = seq(0,10,by=.01)

pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-08a_Lecture12/Figures/LogNorm.pdf',height=12,width=12)
par(mfrow=c(3,2))
for (i in 1:3){
	#Plot histogram for each sigma.
	hist(fx[,i],freq=F)
	f = target(x,avals[i])
	lines(x,f,lwd=2,col='red')
	
	#Plot trace for each sigma.
	plot(fx[,i],type='l',main=paste('a = ',avals[i],sep=''))
	
}
dev.off()




