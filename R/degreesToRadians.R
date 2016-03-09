#----------------------------------------------------------------
# Roxy package build comments:
#' degreesToRadians(degrees, minutes)
#'
#' This function converts degrees to radians.
#'
#' @param degrees A degrees value.  Defaults to 45.
#' @param minutes A minutes value.  Defaults to NULL. 
#'
#' @keywords trigonometry
#'
#' @return radians The radians equivalent of the input degrees value.
#'
#' @examples
#' radiansToDegrees(45,0)	##Example entering degrees and minutes.
#' radiansToDegrees()		##Defaults to 45 degrees, 0 minutes.
#' radiansToDegrees(45)  	##No need to input a minutes value.
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

degreesToRadians<-function(degrees=45,minutes=0)
{
	if(!is.numeric(minutes)) stop("Please enter a numeric value for minutes!\n")
	if(!is.numeric(degrees)) stop("Please enter a numeric value for degrees!\n")
	decimal<-minutes/60
	c.num<-degrees+decimal
	radians<-c.num*pi/180
	radians
}