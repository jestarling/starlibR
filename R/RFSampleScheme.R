
# RFSampleScheme
#
# Compute sample levels for each class in a training data set.
# Returns a vector for randomforest's sampsize parameter.
# Initial sample sizes are driven by estimated class proportions from the unlabeled data set
# and are then bounded to be no lower than by a desired given class minimum 
# and no larger than the class levels in the labeled data set.
#
# totalObs       : # observations desired
# minObsPerClass : minimum # observations desired per class
# nKnown         : # observations in the data set with known class labels
# propClassKnown : class proportions in the data set with known class labels
# nUnknown       : # observations in the data set with unknown class labels
# propClassUnknown : estimated class proportions in the data set with unknown class labels
#

#----------------------------------------------------------------
# Roxy package build comments:
#' RFSampleScheme: A Random Forest function.
#'
#' This function builds a vector of sample sizes to be used as the sampsize parameter
#' in the random forest function.  The vector contains sample sizes for each class level.
#'
#' @param totalObs The number of observatiosn desired for the entire data set, all classes together.
#' Defaults to NULL.
#' @param minObsPerClass The minimum number of observations allowed per class.  Defaults to NULL.
#' @param nKnown Number of observataions in the data set with known class labels. Defaults to NULL.
#' @param nUnknown Number of observataions in the data set with unknown class labels. Defaults to NULL.
#' @param propClassUnknown Estimated class proportiosnin the data set with unknown class labels.
#' Defaults to NULL.
#'
#' @keywords random forest variable star classification
#'
#' @return A vector containing sampsize class proportions to use in the RandomForest model.
#'
#' @author Robert Saw
#'
#' @export

#----------------------------------------------------------------

RFSampleScheme <- function(.totalObs, .minObsPerClass, .nKnown, .propClassKnown, .nUnknown, .propClassUnknown) {

  .nClassKnown <- round(.nKnown * .propClassKnown)

  .sample <- round(.totalObs * .propClassUnknown)
  
  for (.i in 1:length(.sample)) {
    if ( .sample[.i] < .minObsPerClass) {
      .sample[.i] = .minObsPerClass
    }
    if ( .sample[.i] > .nClassKnown[.i]) {
      .sample[.i] = .nClassKnown[.i]
    }
  }
  .sample
}

