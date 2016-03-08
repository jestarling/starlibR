#-------------------------------------------------------------
# Package: varstar
# Function: calcProps(data)

#Purpose: Calculates the class proportions of the 'class' column
# for the data set.  Data set must contain a column named class.

#Inputs:
#    data = A data frame with a 'class' column, containing various 
#		class values for each row of data.

#Output: 
#	classProps, a vector of class proportions for the data$class column.

	
#----------------------------------------------------------------
# Roxy package build comments:
#' emCalcMix(x,densities,p,maxiter,conv)
#'
#' This function calculates teh class proportions of the 'class' column
#' for the data set.  Data set must contain a column named class.
#'
#' @param data A data frame with a 'class' column, containing class values
#'	for each row. Defaults to NULL.
#'
#' @keywords 
#'
#' @return classProps A vector of class proportions for the data$class column.
#'
#' @examples
#' ##data is a data frame containing a column named 'class'.
#' p <- classProps(data)
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------
calcProps <- function(data){
	
	#Extract all individual classes, in alphabetical order.
	classes <- sort(unique(data$class))
	
	#Number of classes.
	k <- length(classes)
	
	#Total sample size (length of class_vector).
	n <- length(data$class)
	
	#Calculate class proportions.
	classProps <- rep(0,k)
	
	for (i in 1:k){
		classProps[i] <- sum(data$class==classes[i]) / n
	}
	
	return(classProps)
}