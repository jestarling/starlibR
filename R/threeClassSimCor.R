#-------------------------------------------------------------
# Package: varstar
# Function: threeClassSimCor(nA,nB,nC)

#Purpose: Creates a simulated data set with three classes (A,B,C), with some corrupted features.

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
#' threeClassSimCor(nA,nB,nC)
#'
#' Creates a simulated data set with three classes (A,B,C).
#' Some features are corrupted, to simulate real-world data.
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
#' datasim <- threeClassSimCor(100,200,300)
#'
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

threeClassSimCor <- function(nA,nB,nC,corFeatures,wellSep=TRUE){
	
	
	nA <- nA; nB <- nB; nC <- nC; #initialize variables
	wellSep <- wellSep 
	
	x <- list() #Empty list for function to return.
	
	nCor <- corFeatures
	
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

	#Create a corrupted data set

	#Randomly select the columns to be corrupted
	#cor.index <- sample(1:50, 25)
	cor.index <- sample(1:50,nCor)

	#Generate random means for the corrupted variables.  
	#(Not all will be used)
	cor.means.A <- runif(50, 0, 1.25)
	cor.means.B <- runif(50, 0, 1.25)
	cor.means.C <- runif(50, 0, 1.25)


	#Assign the newly-generated, corrupted means to the variables that were 
	#selectd for corruption.
	A.means.c <- A.means
	B.means.c <- B.means
	C.means.c <- C.means
	for(i in cor.index){
  		A.means.c[i] <- cor.means.A[i]
  		B.means.c[i] <- cor.means.B[i]
  		C.means.c[i] <- cor.means.C[i]
	}

	#Create a new data set, as above, including the corrupted feature.
	A.c <- matrix(NA, nrow=nA, ncol=50)
	B.c <- matrix(NA, nrow=nB, ncol=50)
	C.c <- matrix(NA, nrow=nC, ncol=50)

	for(i in 1:50){
  		A.c[,i] <- rnorm(nA, A.means.c[i], A.sd[i])
  		B.c[,i] <- rnorm(nB, B.means.c[i], B.sd[i])
  		C.c[,i] <- rnorm(nC, C.means.c[i], C.sd[i])
	}
	cor.data <- rbind(A.c, B.c, C.c)
	cor.data <- data.frame(cor.data)
	cor.data <- cbind(c(rep("A", nA), rep("B", nB), rep("C", nC)), cor.data)
	colnames(cor.data)[1] <- "class"

	##Distibutions of the Corrupted data (Key)
	cor.A.dists <- NULL; cor.B.dists <- NULL; cor.C.dists <- NULL
	for(i in 1:50){
  		cor.A.dists[i] <- paste("N(", round(A.means.c[i], 2), 
                      ", ", round(A.sd[i], 2), ")", sep="")
  		cor.B.dists[i] <- paste("N(", round(B.means.c[i], 2),
                      ", ", round(B.sd[i], 2), ")", sep="")
  		cor.C.dists[i] <- paste("N(", round(C.means.c[i], 2),
                      ", ", round(C.sd[i], 2), ")", sep="")
	}

	cor.dists.table <- cbind(cor.A.dists, cor.B.dists, cor.C.dists)
	colnames(cor.dists.table) <- c("Corrupted A", "Corrupted B", "Corrupted C")

	x$data <- cor.data	
	x$key <- cor.dists.table

	return(x)

}