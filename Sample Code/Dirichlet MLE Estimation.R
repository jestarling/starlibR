################################################################################
##   SDS 383C - Stats Modeling 1 - Fall 2016 - PROBLEM 3c:                    ##
################################################################################

setwd('/Users/jennstarling/UTAustin/starlib/data')

#i: Gives a scatter plot of the data.
data <- as.matrix(read.table('./dir1.txt'),ncol=3,byrow=T)

plot(data,main='Scatter Plot of Data') #Is a 2-simplex.

#------------------------------------------------------------------------------
#ii: Compute the MLE alpha_hat, and plot the log-likelihood as a function of 
#iteration. Briefly give a description of the algorithm you use.


#Dirichlet Log-likelihood function:
dir_logl <- function(data,alpha){	
	N <- nrow(data)						#Number of (j) observations.
	logxi <- colMeans(log(data))		#logxi = (1/N)*sum(log(x_ij)
	logl <- N* (lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha-1)*logxi))

	return(logl)	
}

#Inversion of digamma function using 
#"Estimating a Dirichlet Distribution" by Thomas P. Minka, 2000, Appendix 3 method.
inv_digamma <- function(y){
	
  #Initialize psi inverse iterations.	
  psiinv_y <- array(NA, dim=c(6,length(y)))

  #Set gamma to -Psi(1)
  gamma <- -digamma(1)
  
  #Initialize psiinv_y values.
  for (i in 1:length(y)){
    if (y[i] >= -2.22)
      psiinv_y[1,i] <- exp(y[i]) + 0.5
    else
      psiinv_y[1,i] <- -1/(y[i]+gamma)
  }
  
  #Iterate through Newton's algorithm to update psiinv_y values.
  #Only 6 iterations needed, per author.
  for (i in 2:6){
    psiinv_y[i,] <- ( psiinv_y[i-1,] - 
		(digamma(psiinv_y[i-1,]) - y)/(trigamma(psiinv_y[i-1,])) )
  }
  
  return (psiinv_y[6,])
}


#FIND MLE OF ALPHA VECTOR:
conv <- 1*10^(-5)	#Convergence tolerance
maxiter=1000		#Maximum iterations.
	
#Initialize matrix to hold gradients for each iteration.					
grad <- matrix(0,nrow=maxiter,ncol=ncol(data)) 	
	
#Initialize matrix to hold alphas for each iteration.
alphas <- matrix(0,nrow=maxiter+1,ncol=ncol(data)) 		
alphas[1,] <- rep(1,3)

N <- nrow(data)						#Number of (j) observations.
logxi <- colSums(log(data)) /N		#logxi = (1/N)*sum(log(x_ij)
	
#Initialize log-likelihood.
loglik <- rep(0,maxiter) 
loglik[1] <- dir_logl(data,alphas[1,])
	
#Maximize loglikelihood with iterative algorithm.	
for (i in 2:maxiter){
	
	#Update alpha.			
	alphas[i,] <- inv_digamma(logxi + digamma(sum(alphas[i-1,])))
	
	#Update loglikelihood.
	loglik[i] <- dir_logl(data,alphas[i,])
		
	#Convergence check.
	if (abs((loglik[i]-loglik[i-1]))/abs(loglik[i-1]+1E-3) < conv){
		print('Convergence reached.')
		break;
	}	
}

#Save alpha estimted values.		
alpha_hat <- alphas[i,] 

#Plot the log-likelihood function.
plot(1:i,loglik[1:i],type='l',col='blue',xlab='iteration',ylab='log-likelihood',
	main='Dirichlet Log-Likelihood Convergence')

#------------------------------------------------------------------------------
#3 iii: 

#The Dirichlet pdf function for a 2-simplex. (3 alphas)
f <- function(x1,x2,alpha){
	f <- ((gamma(sum(alpha)) / prod(gamma(alpha)) ) 
		* x1^(alpha[1]-1) * x2^(alpha[2]-1) * (1-x1-x2)^(alpha[3]-1))
	return(f)
}

#Create grid to evaluate Dirichlet distribution.
x1 <- x2 <- seq(.01,.99,by=.01)
z <- outer(x1,x2,f,alpha_hat)

#Contour plot of Dirichlet distribution with alpha_hat parameters, with data.
contour(z,col='red',main='Dirichlet with alpha_hat Parameters & Data',
	xlab='x1',ylab='x2')
points(data,col='black',pch=20)
