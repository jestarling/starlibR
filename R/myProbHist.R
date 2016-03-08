#This function adjusts the y-labels of an R histogram, to show a y-axis range of 
#probability values from 0 to 1.

#Source: http://stackoverflow.com/questions/13615666/histogram-does-not-show-densities

#----------------------------------------------------------------
# Roxy package build comments:
#' myProbHist: A probability histogram function.
#'
#' This function plots the probability histogram for the input vector of data (x).
#' The y-axis is scaled from 0 to 1, to reflect the probability of each bin.
#'
#' @param x A vector of values. Defaults to NULL.
#'
#' @keywords probability histogram
#'
#' @return Plots a probability histogram of the x data.
#'
#' @examples
#' ## Define ctrl object.
#' x <- rnorm(1000,0,1)
#' myProbHist(x)
#'
#' @author Jennifer Starling
#'
#' @export

myProbHist <- function(x){
	tmp <- hist(x, yaxt='n',ylab='Proportion') #yaxt='n' suppresses y axis labels in plot
	tmp2 <- pretty( tmp$counts/sum(tmp$counts) ) #Calculates new y-axis labels.
	
	#Extend new y-labels from 0 to 1, to show whole y-axis.
	tmp2_breaks <- tmp2[2] - tmp2[1] #Calculates the width of each y-label increment.
	tmp3 <- round(seq(0,1,by=tmp2_breaks),2) #Fills in y-values from 0 to 1.
	
	axis(2, at=tmp3*sum(tmp$counts), labels=tmp3) #Adds new y-axis labels to hist plot.
}