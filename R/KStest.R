#----------------------------------------------------------------------------
#FUNCTION 5: KS Function to compare two densities. (For feature selection)
#  KS-Test function for EM Algorithm
#  This function takes the data input into our code to produce the mixture proportions 
#  and compares it to the cdf computed from the EM algorithm.
#  Inputs:
#  data = the data input into the EM algorthim
#  simData = simulation of the data based on the proportions of the EM algorithm
#
#  Outputs:
#  p-value for the KS-Test, in a list
#  Test statistic for the KS-Test, in a list
#  plot of the cdfs of the data and estimates, in png format

#----------------------------------------------------------------
# Roxy package build comments:
#' KStest: A function to compare two estimated densities.
#'
#' This function takes two vectors of values and compares them, to determine
#' if they come from the same distribution.
#'
#' @param x1 A vector of values. Defaults to NULL.
#' @param x2 A vector of values.  Defaults to NULL.
#'
#' @keywords KS test
#'
#' @return A list object, containing the following list elements.
#' @return ks$pvalue A pvalue indicating the confidence that the data are the same distribution.
#' @return ks$statistics The KS test statistic for comparing the two data sets.
#'
#' @examples
#' ## Define ctrl object.
#' x <- rnorm(1000,0,1)
#' myProbHist(x)
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------

KStest = function(x1,x2){
	
	library(ks) 			#For K-S feature selection.
	
	ks = NULL
	test = ks.test(x1,x2)
	ks$pvalue = as.numeric(test$p.value)
	ks$statistic = as.numeric(test$statistic)
	
	return(ks)
}