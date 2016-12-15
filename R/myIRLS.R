#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Iteratively Reweighted Least Squares Function
#(Binomial Distribution)
#	Inputs:
#		X = design matrix
#		y = vector of 1/0 response values
#		maxiter = maximum number of iterations allowed
#		conv = tolerance level for convergence check
#	Outputs:
#		loglik = loglikelihood function (should be maximized)
#		beta_hat = estimated coefficients
#		iter = iterations completed
#		converged = convergence indicator

irls = function(X,y,maxiter=500,conv=1E-14){
	
	#------------------------------------------------------------------	
	#Binomial Loglikelihood function. 
		#Inputs: 	X = design matrix 
		#			Y = vector of 1, 0 values. 
		#   		m = sample size vector m.
		#			w = probabilities vector.  (Diagonal of W matrix.)
		#Output: Binomial loglikelihood value.
	logl <- function(X,Y,beta,m,w){
		logl <- sum(Y*log(w+1E-4) + (m-Y)*log(1-w+1E-4)) #Calculate log-likelihood.
		return(logl)	
	}
	#------------------------------------------------------------------
	
	#Indicator for whether convergence met.
	converged <- 0			
	
	#1. Initialize matrix to hold beta vector for each iteration.
	betas <- matrix(0,nrow=maxiter+1,ncol=ncol(X)) 
	betas[1,] <- rep(0,ncol(X))	#Initialize beta vector to 0 to start.
	
	#2. Initialize values for log-likelihood.
	loglik <- rep(0,maxiter) 	#Initialize vector to hold loglikelihood fctn.
	loglik[1] <- logl(X,Y,betas[1,],m,w)

	#3. Perform Newton's method iterations.	
	for (i in 2:maxiter){
		
		#Calculate probabilities vector w_i.
		w = 1 / (1 + exp(-X %*% betas[i-1,]))
		
		#Set up W matrix, using w as diagonal elements.	
		W = diag(as.vector(w))
		
		#Update betas.
		betas[i,]  = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
		
		#Update loglikelihood.
		loglik[i] = logl(X,Y,betas[i,],m,w)

		#Check if convergence met: If yes, exit loop.
		if (i>2 && abs(loglik[i]-loglik[i-1])/abs(loglik[i-1]+1E-3) < conv ){
			converged=1;
			break;
		}
	} #end IRLS iterations.
	
	return(list(beta_hat=betas[i,],iter=i,logl = loglik[1:i]),converged=converged)
} #End IRLS function.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~