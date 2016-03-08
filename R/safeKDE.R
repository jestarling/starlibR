#----------------------------------------------------------------
# Roxy package build comments:
#' safeKDE(x,classes,kdeClass)
#'
#' Function purpose is to create kernel density estimates for a given vector of data.
#' Function handles errors seen when the data has a variance of zero.  A vector of 
#' categorical 'group' levels can be provided (grps input), along with a level
#'
#' @param x A vector of values for which a KDE will be created (either for the entire x vector,
#'		or for just the specified kdeClass, depending on inputs.)
#' @param classes A vector of categorical class values for vector x.  Length must equal length(x).
#' @param kdeClass The class within classes for which a kde should be created.  If classes and
#' 	kdeclass are not populated, a kde is created for the entire x vector. 
#'
#' @keywords kde kernel density estimation
#'
#' @return ks A list object with three components.
#' @return ks$eval.points The x values for the kde estimate.
#' @return ks$estimate The y values for the kde estimate.
#' @return ks$h The bandwidth for the kde estimate.
#'
#' @examples
#' ##data is a data frame containing a column named 'class'.
#' p <- classProps(data)
#'
#' @author Jennifer Starling
#'
#' @export


#----------------------------------------------------------------

safeKDE = function(x,classes=NULL,kdeClass=NULL){

	library(ks) #Loads ks library for kde() function.
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#ERROR CHECKING:
	
	#If calculating KDE for just a particular class of x,
	#must provide both a class vector and the kdeClass, ie the specific class for
	#which to calculate kde.	
	
	if (is.null(classes) & !is.null(kdeClass)) stop("If populating classes, must provide both classes and kdeClass")
	
	if (!is.null(classes) & is.null(kdeClass)) stop("If populating classes, must provide both classes and kdeClass")
	
	if (!is.null(classes) & (length(classes) != length(x))) stop("x and classes must be same length")
	
	if (!is.null(classes) & !is.null(kdeClass) & (length(kdeClass) != 1)) stop("If populating classes, kdeClass must be length 1.")
	
	if (!is.null(kdeClass) & !(kdeClass %in% classes)) stop("kdeClass must be a value in classes.")

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	
	if(is.null(classes)) {x_subset <- x}
	else {x_subset <- x[which(classes==kdeClass)]}
	
    #Bandwidth 1 (bw1) as std dev for the class.
    bw1 <- sd(x_subset,na.rm=T) 
    
    #In order to create bw2, must check that kde() for all of x doesn't error.
    bw2 <- NA #Initialize bw2.
    if (class( try(kde(x),silent=T) )!='try-error')
	{ 
    	     bw2 <- kde(x)$h #Bw for all of x, not just x in the specified group.
	}
	
	#1. If default kde() does not result in error, use default kde().
	if (class( try(kde(x_subset),silent=T) )!='try-error')
	{ 
		result <- kde(x_subset)
		print("default kde")
	}
	
	#2. If kde() with entire feature bandwidth errors, try bw = sdev for class.
	else if (class( try(kde(x_subset,h=bw1),silent=T) )!='try-error')
	{ 
		result <- kde(x_subset,h=bw1)
		print("bw=sdev kde")
	}
	
	#3. If default kde() errors, first try kde() using bandwidth for entire feature.
	#Only proceed with this method if bw2 was successfully created above.
	else if (class( try(kde(x_subset,h=bw2),silent=T) )!='try-error' & is.numeric(bw2))
	{ 
		result <- kde(x_subset,h=bw2)
		print("bw-for-all-of-x kde")
	}	
	
	#4. If all previous kde() attempts error, use density() fctn and cast as kde().
	else{
		#Placeholder kde() to put density() values into.
		result <- kde(rnorm(length(x_subset)))
		
		#Density function.
		temp <- density(x_subset)
		
		#Set kde object values to density object values.
		result$eval.points <- temp$x
		result$estimate <- temp$y
		result$h <- temp$bw	
		
		print("density-cast kde")
	}
	
	return (result)
}	

