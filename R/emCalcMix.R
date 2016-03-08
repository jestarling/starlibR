#-------------------------------------------------------------
# Package: varstar
# Function: emCalcMix(x,densities,p,maxiter,conv)

#Purpose: Performs expectation-Maximization Algorithm to estimate mixing
#proportions for each feature in a data set.

#Performs EM algorithm for a mixture distribution 
#with K known densities with known parameters, and an unknown mixture
#proportion p.  Goal of function is to estimate p for specified feature.

#Inputs:
#    x = Vector of observed feature values for specified feature.
#    densities = List of k density function; one per class. Flexible length, which must 
#		equal k.
#    p = Vector of class starting proportion estimates.  Default is .2 for 5 classes.
#    maxiter = Maximum number of iterations of the EM algorithm allowed.
#    conv = Convergence criteria for change in likelihood function.  Change between
#           two iterations < conv indicates convergence.

#Output: emCalcMix list object with the following components.
#    emCalcMix$p = vector of 5 estimated mixing proportions for feature.
#    emCalcMix$iter = number of iterations before algorithm converged or was stopped.
#    emCalcMix$allp = (iter x 5) data frame containing all phat iterated estimates.
#    logl = Vector of log-likelihood values, with one log-lik value for each iteration.
	
#----------------------------------------------------------------
# Roxy package build comments:
#' emCalcMix(x,densities,p,maxiter,conv)
#'
#' This function performs the Expectation-Maximization Algorithm to estimate
#' a vector of class proportions for a single feature in a data set wih an 
#' unknown set of class proportions. Called by \code{\link{featureMixtureProportion}}.
#'
#' @param x A vector of observed feature values for specified feature. Defaults to NULL.
#' @param densities List of k density function; one per class. Flexible length, which must 
#'		equal k.
#' @param p Vector of class starting proportion estimates.  Default is .2 for 5 classes.
#' @param maxiter Maximum number of iterations of the EM algorithm allowed. Default is 2000.
#' @param conv Convergence criteria for change in likelihood function.  Change between
#       two iterations < conv indicates convergence. Default is 0.00001.
#'
#' @keywords EM expectation-maximization
#'
#' @return A list containing the following components:
#' @return p = vector of 5 estimated mixing proportions for feature.
#' @return iter = number of iterations before algorithm converged or was stopped.
#' @return allp = (iter x 5) data frame containing all phat iterated estimates.
#' @return logl = Vector of log-likelihood values, with one log-lik value for each iteration.
#'
#' @examples
#'
#' ## Define ctrl object.
#' em <- emCalcMix(f_sample,densities,p=rep(.2,5),maxiter=1000,conv=.00001)
#'
#' ##Sample density function:
#' f_class1 = df[df$class=="eb",f_name]
#'    	kd1  = density(f_class1)
#'    	dens1 = function(x)
#'		{
#'			min_den <- quantile(kd1$y,.05)
#'			de <- approx(kd1$x,kd1$y,x)$y
#'			de[is.na(de)] <- min_den
#'			return(de)
#'		}
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

emCalcMix <- function(x,densities,p=rep(.2,5),maxiter=2000,conv=.00001){
	
	#Argument Definitions:
	# x = vector of values
	# densities = list of density functions, containing k functions
	# p = vector of starting probabilities.  length(p) must equal k.
	# maxiter = maximum number of EM algorithm iterations allowed.
	# conv = change in the log-likelihood from one iteration to next
	#			to consider convergence to be met.
	
	K = length(densities) #Defines number of densities provided.
		
	## --------------------------------------
	##Error check sum(probs = 1):
	if (sum(p)!=1) stop("Sum of probabilities must equal 1.")
	
	##Error check: Number of starting probs p = number of densities provided in list
	if (length(densities)  != length(p)){  
		stop("Provide same number of densities and starting probabilities")}
	
	## Error check density inputs:
	for (i in 1:K){
		try(densities[[i]](x),silent=T)
		if(substr(geterrmessage(),1,26)=='Error in densities[[i]](x)'){
			stop("All densities must be in valid format. \n
				Example: dens1 = function(x) {return(dnorm(x,mean=0,sd=1))}")}
	}
	## --------------------------------------
	
	#Set up placeholders for p estimates, log-likelihood estimates, and iter count.
	len <- length(p)
	
	#Vector to store est props for each iter.
	phat <- matrix(rep(-10,len*(maxiter+1)),nrow=maxiter+1,byrow=T)
		
	phat[1,] <- p #Set first phat estimate row to starting value.
	l <- rep(0,maxiter+1)
	iterations <- 1
	epsilon <- 0.000001
	
	#Loop until convergence is met, or max iterations is met.
	for (i in 1:maxiter){
		
		#EXPECTATION STEP
		#Calculate estimated responsibility for each obs.
		denom <- 0 #Reset denominator before beginning step.
		ghat <- matrix(0,nrow=K,ncol=length(x)) #Store ghat values.
		
		#Calculate denominator for each ghat_i value.
		for (j in 1:K) { 
			denom <- denom + phat[i,j]*densities[[j]](x,j) 
		}

		#Calculate estimated responsibilities
		for (j in 1:K){
			ghat[j,] <- phat[i,j]*densities[[j]](x,j) / denom
		}
		
		#MAXIMIZATION STEP
		phat[i+1,] <- rowMeans(ghat)
		
		#CALCULATE LOG-LIKELIHOOD FUNCTION
		#(Is actually minimizing the neg-log-likelihood fctn.)
		for (j in 1:K){
			l[i+1] <- l[i+1] + sum(log(phat[i+1,j]*densities[[j]](x,j)+epsilon))
		}
		
		iterations <- i+1 #Increment iteration count.
		
		#CHECK IF CONVERGENCE MET
		if(abs(l[i+1]-l[i])<=conv) {break}
	}
	
	#Data cleanup.
	phat <- phat[1:iterations,] #Get rid of any unused phat placeholders
	l <- l[1:iterations] #Get rid of any unused l placeholders
	
	#Function output
	return(list(p=phat[i+1,],iter=iterations,allp=phat,logl=l))	
}