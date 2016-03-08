#-------------------------------------------------------------
# Package: varstar
# Function: emAllFeatures(x,densities,p,maxiter,conv)

#Purpose: Performs expectation-Maximization Algorithm to estimate mixing
#proportions for the entire data set.  Uses all features at once.

#Performs EM algorithm for a data set with K unknown classes, to obtain
#an estimate of the unknown class proportions.  Uses a data set with known
#class proportions to estimate densities for the EM algorithm.

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
#' emAllFeatures(known,unknown,p=rep(.2,5),maxiter=2000,conv=.00001)
#'
#' This function performs the Expectation-Maximization Algorithm to estimate
#' a vector of class proportions for a single feature in a data set wih an 
#' unknown set of class proportions. Called by \code{\link{featureMixtureProportion}}.
#'
#' @param known A data set with class as the first column.  This data set has known 
#'		class proportions.
#' @param unknown A data set whose proportions will be estimated.  This data set may have 
#'	class as the first column, to return an MSE to compare the estimates to the actual proportions.
#' @param p Vector of class starting proportion estimates.  Default is .2 for 5 classes.
#' @param maxiter Maximum number of iterations of the EM algorithm allowed. Default is 2000.
#' @param conv Convergence criteria for change in likelihood function.  Change between
#       two iterations < conv indicates convergence. Default is 0.00001.
#'
#' @keywords EM expectation-maximization
#'
#' @return A list containing the following components:
#' @return p = vector of estimated mixing proportions for the unknown data set.
#' @return iter = number of iterations before EM algorithm converged or was stopped.
#' @return allp = (iter x 5) data frame containing all phat iterated estimates.
#' @return logl = Vector of log-likelihood values, with one log-lik value for each iteration.
#' @return mse = MSE of actual unknown class proportions (if in data set) versus estimated.
#'
#' @examples
#'
#' ## Define ctrl object.
#' em <- emAllFeatures(df1,df2,p=rep(.2,5),maxiter=1000,conv=.00001)
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------
emAllFeatures <- function(known,unknown,p=rep(.2,5),maxiter=2000,conv=.00001){
	
	### Argument Definitions:
	# known = A data frame containing class as the first column.
	# unknown = A data frame where classes are unknown (may or may not contain class col).
	# p = vector of starting probabilities.  length(p) must equal # classes in known set.
	# maxiter = maximum number of EM algorithm iterations allowed.
	# conv = change in the log-likelihood from one iteration to next
	#			to consider convergence to be met.
	
	df = unknown			#Initializes unknown data frame.
	known = known			#Initializes known data frame.
	K = length(p) 			#Defines number of class proportions.
	J = ncol(known) - 1 	#Defines number of features (in known set).
	classes = sort(unique(known$class)) #Defines classes in known set.
	
	#Known data frame is used to create estimated densities for each class/feature combo.
	#Unknown data frame is the data whose class proportions are being estimated.
		
	#--------------------------------------
	### Error Checking inputs:
	# Error check sum(probs = 1):
	if (sum(p)!=1) stop("Sum of probabilities must equal 1.")
		
	#Error check: First column of known must be 'class'.
	if (colnames(known)[1] != 'class'){  
		stop("First column of known data frame must be 'class'.")}
		
	#Error check: Length of p must match number of classes in known data set.
	if (K != length(classes)){  
		stop("Starting proportion vector length must = number of classes in known data set.")}		
		
	#--------------------------------------
	### Create a JxK matrix of estimated density functions, using known class data.
	
	#Setup matricies to hold f_class, kd, and density objects.
	#Note: To use a matrix of lists, must reference each element as densities[[i,j]]
	j <- sort(rep(1:J,times=K)) #List vector of j indices to loop over.
	k <- rep(1:K,times=J)		#List vector of k indices to loop over.
	indices <- data.frame(j=j, k=k) #Set up data frame of all possible index pairs.
	
	#Create f_class and kd estimators for each class/feature combo.
	f_class <- mapply(function(j,k) known[known$class==classes[k],j+1], j,k)
	kd <- matrix( lapply(f_class,density), nrow=J, ncol=K, byrow=T)
	
	#Iterate over features and classes to create matrix of density functions.
	densities <- matrix(list(),nrow=J,ncol=K)
	for (i in 1:nrow(indices)){
		j = indices[i,1]; k = indices[i,2]
		
		#Create density function for class/feature.
		densities[[j,k]] = function(x){
			min_den <- quantile(kd[[j,k]]$y, 0.05)
			de <- approx(kd[[j,k]]$x,kd[[j,k]]$y,x)$y
			de[is.na(de)] <- min_den
			return(de)
			}
	}
	  
	#--------------------------------------
	### Set up f_jk(x_i) vectors, ie each function evaluated at 
	#Create a JxK matrix, each cell contains evaluated density for j/k combo, ie f_jk(x_i)	
	densities_eval <- matrix(0,ncol=nrow(indices),nrow=nrow(df))
	
	for (i in 1:nrow(indices)){
		j = indices[i,1]; k = indices[i,2]
		densities_eval[,i] <- densities[[j,k]](df[ ,j+1])
	}
	
	#--------------------------------------
	### EM Algorithm variables setup:
	
	#Vector to store est props for each iter.
	phat <- matrix(rep(-10,K*(maxiter+1)),nrow=maxiter+1,byrow=T)
	phat[1,] <- p #Set first phat estimate row to starting value.
	
	l <- rep(0,maxiter+2) 	#Store log-likelihoods for each iter.
	iter <- 1				#Initialize iteration tracker to 1.
	epsilon <- 0.000001 	#Initialize epsilon value (error prevention).

	#--------------------------------------
	### EM Algorithm Execution - Loop until convergence met, or max iterations met:
	
	for (i in 1:maxiter){
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		### EXPECTATION STEP:
		
		#Calculate estimated responsibility for each obs:
		# g_k(x_i) = {p_k * PRODUCT_j=1toJ(f_jk(x_i))}  / 
		#	{SUM_k=1toK[ p_k * PRODUCT_j=1toJ(f_jk(x_i)) ]}
		
		#1. Numerator calculation
		
		#Store k ghat(x_i) numerator values for ea class.
		ghat_num <- matrix(0,nrow=nrow(df),ncol=K) 
			
		#Calculate product term by multiplying all J densities together for the class.
		for (k in 1:K){ #Loop through the K classes.				
			cols <- which(indices[,2]==k) #Select the columns of densities_eval for each k val.
			ghat_num[,k] <- apply(densities_eval[,cols],1,prod) * phat[i,k]
		} #End class loop.
		
		#2. Denominator calculation
		ghat_denom <- rowSums(ghat_num)
		
		#3. Put numerator/denominator together to get ghat_k(x_i) values.
		ghat <- ghat_num / ghat_denom
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
		### MAXIMIZATION STEP:
		#p_k = (1/n) * SUM_i=1ton[g_k(x_i)]
		phat[i+1,] <- colMeans(ghat,na.rm=T)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		### CALCULATE LOG-LIKELIHOOD FUNCTION:
		
		#The log-likelihood function is:
		# L = log ( prod_1toK(phat_k * {prod_1toJ(density_kj(y_j))} ))
		# Take the log, and simplify to get:
		# logL = sum_1toK(log(phat_k)) + sum_1toK( sum_1toJ(log(density_kj(y_j))) )
		
		#This is using the evaluated density values from earlier in the function.
		l[i+1] <- sum(log(phat[i+1,])) + sum(log(densities_eval+epsilon))
			
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		### CHECK IF CONVERGENCE MET:
		iterations <- i+1 #Increment iteration count.
		
		if(is.na(abs(l[i+1]-l[i]))){diff <- 0} else {diff <- abs(l[i+1]-l[i])}
		if(diff <= conv) {break}
		
	} #End maxiter loop.
	
	#Data cleanup.
	phat <- phat[1:iterations,] #Get rid of any unused phat placeholders
	l <- l[1:iterations] #Get rid of any unused l placeholders
	
	#Calculate MSE of phat estimate, compared to unknown class proportions.
		#(Execute only if class is first col in unknown data set.)
	if (colnames(unknown)[1] == 'class'){
		#Calculate class proportions.
		unknownProps <- rep(0,K)	
		for (k in 1:K){unknownProps[k] <- sum(unknown$class==classes[k]) / length(unknown$class)}
		mse <- (1/K) * sum((unknownProps-phat[i+1,])^2)
	}
	if (colnames(unknown)[1] != 'class'){ mse <- 'unknown class not provided'}
	
	#Function output
	return(list(p=phat[i+1,],iter=iterations,allp=phat,logl=l,mse=mse))	
	
} #End function.


