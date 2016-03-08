#-------------------------------------------------------------
# Package: varstar
# Function: threeClassSim(nA,nB,nC)

#Purpose: Creates a simulated data set with three classes (A,B,C).

#Inputs:
#   nA:  Number of observations to generate of class A.
#	nB:  Number of observations to generate of class B.
#	nC:  Number of observations to generate of class C.
#   wellSep: T/F value to indicate whether distributions should be well-separated for the classes.

#Output: 
#	x, a list containing two objects:
#	x$data: The simulated data set.
#	x$key: The distribution key for the simulated data set.

#----------------------------------------------------------------
# Roxy package build comments:
#' threeClassSim(nA,nB,nC)
#'
#' Creates a simulated data set with three classes (A,B,C).
#'
#' @param nA Number of observations to generate from class A. Defaults to NULL.
#' @param nB Number of observations to generate from class B. Defaults to NULL.
#' @param nC Number of observations to generate from class C. Defaults to NULL.
#' @param wellSep Value to indicate whether distributions of three classes should be well separated.  Defaults to TRUE.
#'
#' @return A list containing the following components:
#' @return x$data The simulated data set.
#' @return x$key The distribution key for the simulated data set.
#'
#' @examples
#'
#' ## Define ctrl object.
#' datasim <- threeClassSim(100,200,300)
#'
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

threeClassSim <- function(nA,nB,nC,wellSep=TRUE){
	
	nA <- nA; nB <- nB; nC <- nC; #initialize variables
	
	x <- list() #Empty list for function to return.

	#Create vectors of random values to be used as means and
	#standard deviations for the variables.  THere is some overlap
	#if WellSep = FALSE.
	
	if (wellSep==TRUE){
		#WELL_SEPARATED:
		A.means <- runif(50, 4.5, 5.0)
		A.sd <- runif(50, 0.5, 1.0)
		B.means <- runif(50, 0.5, 1.5)
		B.sd <- runif(50, 0.5, 1.0)
		C.means <- runif(50, 3, 3.5)
		C.sd <- runif(50, 0.5, 1.0)
	}

	if (wellSep==FALSE){
		A.means <- runif(50, 0.5, 1.0)
		A.sd <- runif(50, 0.5, 1.0)
		B.means <- runif(50, 0.5, 1.5)
		B.sd <- runif(50, 0.5, 1.0)
		C.means <- runif(50, 1.0, 1.5)
		C.sd <- runif(50, 0.5, 1.0)
	}

	#Create empty matrices for each class 
	A <- matrix(NA, nrow=nA, ncol=50)
	B <- matrix(NA, nrow=nB, ncol=50)
	C <- matrix(NA, nrow=nC, ncol=50)

	#Sample from normal distributions with the above
	#randomly generated parameters
	for(i in 1:50){
		A[,i] <- rnorm(nA, A.means[i], A.sd[i])
		B[,i] <- rnorm(nB, B.means[i], B.sd[i])
		C[,i] <- rnorm(nC, C.means[i], C.sd[i])
	}

	#Combine the above matrices and label with classes.
	data <- rbind(A,B,C)
	data <- data.frame(data)
	data <- cbind(c(rep("A", nA), rep("B", nB), rep("C", nC)), data)
	colnames(data)[1] <- "class"

	#Create a key, with the underlying distributions
	A.dists <- NULL; B.dists <- NULL; C.dists <- NULL
	for(i in 1:50){
		A.dists[i] <- paste("N(", round(A.means[i], 2), 
                    ", ", round(A.sd[i], 2), ")", sep="")
		B.dists[i] <- paste("N(", round(B.means[i], 2),
                    ", ", round(B.sd[i], 2), ")", sep="")
		C.dists[i] <- paste("N(", round(C.means[i], 2),
                    ", ", round(C.sd[i], 2), ")", sep="")
	}
	dists.table <- cbind(A.dists, B.dists, C.dists)

	x$data <- data
	x$key <- dists.table
	
	return(x)

}