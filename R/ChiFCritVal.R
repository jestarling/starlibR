#Chisq & F-statistic critical value calculator:

#----------------------------------------------------------------
# Roxy package build comments:
#' ChiFCritVal(conf,n,p) A critical value function for Z and T distributions.
#'
#' This function returns two tables of critical values (one-sided and two-sided)
#' for the chi-squared and F (central) distributions.
#'
#' @param conf A vector of confidence levels (1-alpha) to calculate critical values.  
#' 	Defaults to c(.99, .985, .975, .95, .9)
#' @param n A sample size for calculating degrees of freedom (df1 for F).  Defaults to 1000.
#' @param p The number of parameters being estimated for degrees of freedom (df1 for F).  Defaults to 1.
#' @param n2 A sample size for calculating degrees of freedom (df2 for F).  Defaults to 1000.
#' @param p2 The number of parameters being estimated for degrees of freedom (df2 for F).  Defaults to 1.
#'
#' @keywords critical values hypothesis testing
#'
#' @return Returns a list object x containing the following values.
#' @return x$n The sample size input by the user.
#' @return x$p The number of estimated parameters input by the user.
#' @return x$n2 The second sample size input by the user.
#' @return x$p2 The second number of estimated parameters input by the user.
#' @return x$twoSided The table of 2-sided critical values for the chi-squred and F dists.
#' @return x$oneSided The table of 1-sided critical values (right tail) for the chi-squred and F dists.
#'
#' @examples
#' ## Input specific parameters.
#' ChiFCritVal(conf=c(.99,.95), n=1000, p=1)
#' ## Run with default parameters.
#' ChiFCritVal()
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------

ChiFCritVal = function(conf=c(.99,.98,.975,.95,.9),n=1000,p=1,n2=1000,p2=1){
	
	conf=conf #Takes a user input for a vector of 1+ conf (1-alpha) values to use.
		#Default value is c(.99, .98, .975, .95, .9)
	n=n #Sets up n value for df calculation for chi-sq, and df1 for F.
	p=p #Sets up p parameter for df calculation for chi-sq, and df1 for F.
	n2=n2 #Sets up n2 value for df calculation for df2 for F.
	p2=p2 #Sets up p2 parameter for df calculation for df2 for F.
	
	a<-1-conf #Set up vector of alpha values.
	
	#Calculate 2-sided critical values. (Both tail probabilities)
	c <- qchisq(1-a/2,df=n-p-1) 			#Chisq dist.
	f <- qf(1-a/2,df1=n-p-1,df2=n2-p2-1)	#F dist.
	
	#Calculate 1-sided critical values. (Right tail probabilities)
	cU <- qchisq(1-a,df=n-p-1) 				#Chisq dist.
	fU <- qf(1-a,df1=n-p-1,df2=n2-p2-1)		#F dist.
	
	x = list()
	x$n <- n
	x$p <- p
	x$n2 <- n2
	x$p2 <- p2
	x$twoSided <- cbind(conf,c,f)
	x$oneSided <- cbind(conf,cU,fU)
	return(x)
}


