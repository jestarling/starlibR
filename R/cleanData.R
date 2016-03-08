#-------------------------------------------------------------
# Package: varstar
# FUNCTION: cleanData(train,test)

# Function purpose is to eliminate features which cannot be used for classification
# due to problematic or unclean data.  

#Inputs:  The two data sets, test and train (example: hip and linear).

#Outputs: A list containing two cleaned data sets.  
#    cleanData[[1]] = training set. 
#    cleanData[[2]] = test set.

#----------------------------------------------------------------
# Roxy package build comments:
#' cleanData(train,test)
#'
#' This function cleans a test and a training data set, to eliminate features
#' which cannot be used for classification due to problematic or unclean data.
#' Eliminates columns (features) where stdev=NA or 0.  Eliminates rows with
#' infinite values. 
#'
#' @param train A data set whose class values are known. Defaults to NULL.
#' @param test A data set whose class values are unknown.  Defaults to NULL.
#'
#' @keywords
#'
#' @return A list containing the following components:
#' @return cleanData[[1]] = Cleaned training data set.
#' @return cleanData[[2]] = Cleaned test training set.
#'
#' @examples
#' clean_data <- cleanData(train=dataset1, test=dataset2)
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------

cleanData = function(train,test){
	# Get rid of variables with stdev == NA
	train_pred = train[,2:ncol(train)]
	drop_sd_na=colnames(train_pred[,is.na(colSds(as.matrix(train_pred)))])
	train=train[,!(names(train) %in% drop_sd_na)]
	test=test[,!(names(test) %in% drop_sd_na)]

	# Get rid of variables with stdev == 0
	train_pred = train[,2:ncol(train)]
	drop_sd_zero=colnames(train_pred[,0==(colSds(as.matrix(train_pred)))])
	train=train[,!(names(train) %in% drop_sd_zero)]
	test=test[,!(names(test) %in% drop_sd_zero)]

	# Get rid of rows that are not finite numbers
	train_row_ok = apply(train[,2:ncol(train)], 1, function(x){all(is.finite(x))})
	train=subset(train, train_row_ok)

	test_row_ok = apply(test[,2:ncol(test)], 1, function(x){all(is.finite(x))})
	test=subset(test, test_row_ok)
	
	return(list(train,test))
}