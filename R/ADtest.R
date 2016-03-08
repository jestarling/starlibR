#----------------------------------------------------------------------------
#FUNCTION 6: AD Function to compare two densities. (For feature selection)
#  AD-Test function for EM Algorithm
#  This function takes the data input into our code to produce the mixture proportions 
#  and compares it to the cdf computed from the EM algorithm.
#  Inputs:
#  data = the data input into the EM algorthim
#  simData = simulation of the data based on the proportions of the EM algorithm
#
#  Outputs:
#  p-value for the AD-Test, in a list
#  Test statistic for the AD-Test, in a list
#  plot of the cdfs of the data and estimates, in png format

#----------------------------------------------------------------
# Roxy package build comments:
#' ADtest: A function to compare two estimated densities.
#'
#' This function takes two vectors of values and compares them, to determine
#' if they come from the same distribution.
#'
#' @param x1 A vector of values. Defaults to NULL.
#' @param x2 A vector of values.  Defaults to NULL.
#'
#' @keywords AD test Anderson-Darling
#'
#' @return A list object, containing the following list elements.
#' @return ad$pvalue A pvalue indicating the confidence that the data are the same distribution.
#' @return ad$statistics The AD test statistic for comparing the two data sets.
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

ADtest = function(x1,x2){
	
	library(kSamples)		#For A-D feature selection.
	
	ad = NULL
	test = ad.test(x1,x2,method='asymptotic')
	ad$pvalue = as.numeric(test$ad[1,3]) #Version 1 p-val. Use [2,3] for v2 pval.
	ad$statistic = as.numeric(test$ad[1,2]) #Version 2 TS.  use [2,2] for v2 TS.
	
	return(ad)
}
