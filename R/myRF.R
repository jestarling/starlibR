#-------------------------------------------------------------
# Package: varstar
# Function: myRF(known, unknown, ctrl, grid, keeps)

#Purpose: 
# 1. Builds a random forest model using the caret package, on the known data set.
# 2. Predicts class values for unknown set, using the trained rf model.
# 4. Returns error rates and confusion matrices for the known and unknown data sets.

#Input data set must have a variable called class.

#Inputs: 
#  known = data set with known labels.
#  unknown = data set with unknown labels.
#  ctrl = the trainControl statement (caret).
#  grid = the grid of .mtry values to use.  Optional, defaults to NULL.
#  keeps = the vector of feature names (including class) to consider in model.

#Output: A list object x, with components as follows.
#  x$model = Random forest model object (created using caret pckg).
#  x$classPred = Predicted class values for unknown data set.
#  x$conf_matrix_known = Confusion matrix for cross-validated model (on training set).
#  x$result = Accuracy for training model.
#  x$unknown.error = Error rate for applying model to unknown data.
#  x$conf_matrix_unknown = Confusion matrix for applying model to unknown data.

#----------------------------------------------------------------
# Roxy package build comments:
#' myRF: A Random Forest function
#'
#' This function builds a random forest classifier, predicts class values for the
#' unknown data set, and returns error rates and confusion matrices for the 
#' known and unknown data sets.
#'
#' @param known A data set with known classes, used to classify the unknown data set.
#' Defaults to NULL.
#' @param unknown A data set whose classes are considered unknown.  Classes will be 
#' predicted for this data set.  Defaults to NULL.
#' @param ctrl A trainControl statement from the caret package.  Defaults to NULL.
#' @param grid A grid for the tuneGrid parameter in the train (caret) function. Defaults to NULL.
#' @param keeps A vector of feature names to consider in the model (must include 'class'). 
#' Defaults to NULL.
#' @param samps A vector of sample sizes by class for the sampsize random forest argument.
#'
#' @keywords random forest variable star classification
#'
#' @return A list containing the following components:
#' @return x$model = Random forest model object (created using caret pckg).
#' @return x$classPred = Predicted class values for unknown data set.
#' @return x$conf_matrix_known = Confusion matrix for cross-validated model (on training set).
#' @return x$result = Accuracy for training model.
#' @return x$unknown.error = Error rate for applying model to unknown data.
#' @return x$conf_matrix_unknown = Confusion matrix for applying model to unknown data.
#'
#' @examples
#' ## Define ctrl object.
#' c <- trainControl(method='cv',number=5,classProbs=F)
#' ## Define list of features to keep, including 'class' as the first feature.
#' features <- c('class','feature1','feature2','feature3')
#' ## Known and Unknown data sets must contain a 'class' column.
#' model <- myRF(known=labeled_data_set, unknown=unlabeled_data_set, ctrl=c,keeps=features)
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

myRF = function(known, unknown, ctrl, grid, keeps, sampsize=NULL)
{	
	x = list() #Set up empty list to hold function output.
	grid=grid #Initialize grid input to avoid errors.
	
	#If no samps value provided, set equal to rows in known set.
	if (is.null(sampsize)) sampsize=nrow(known)


	#Train model using known data.
	x$model <- train(class ~ ., 
	                 data=known[,keeps],
			             method='rf',
			             metric='Accuracy',
			             trControl= ctrl,
			             tuneGrid=grid,
			             sampsize=samps)

	#Predict class values for unknown data.
	x$classPred <- predict(x$model, newdata=unknown)						

	#Calculate error rates and confusion matrices.

	#1. Test model confusion matrix (for cross-validated model):
	x$conf_matrix_known <- confusionMatrix(predict(x$model), known$class) 
	x$result <- x$model$result
	
	#2. Error rate & confusion matrix for model applied to unknown set.
	x$unknown.error <- mean(1-(x$classPred==unknown$class))
	x$conf_matrix_unknown <- confusionMatrix(x$classPred, unknown$class)

	#Return all of the pieces calculated previously.
	return(x)
}