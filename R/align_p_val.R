#----------------------------------------------------------------
# Roxy package build comments:
#' align_p_val
#'
#' This function computes a p-value for the alignment of two sequences.
#' Length(sequence1) must be <= length(sequence2) in the arguments.
#'
#' @param seq_1 A sequence of DNA or amino acids.  Ex: "ASEDLTI" or "GACT"
#' @param seq_2 A second sequence of DNA or amino acids.  Ex: "AEDFGI" or "GACCT"
#' @param scoringMat A scoring matrix for the alignment. See ?pairwiseAlignment in library Biostrings for 	
#'		details.
#' @param gapOpen A numeric value for the penalty for opening a gap. See ?pairwiseAlignment in 
#' 		library Biostrings for details.
#' @param gapExtend A numeric value for continuing a gap.  See ?pairwiseAlignment in library 
#'		Biostrings for details.
#' @param B The number of bootstrap samples used in the p-value bootstrap calculation.
#' @param type The text 'global' or 'local'.  See ?pairwiseAlignment in library 
#'		Biostrings for details.
#'
#' @keywords sequence alignment genomics
#' @return tt_0 The t-statistic for the original alignment.
#' @return tt_B A vector of t-statistics for B randomly generated alignments.
#' @return p_val The p-value, ie the proportion of randomization-based alignment scores that are
#'	great than or equal to the original observed alignment score.  mean(tt_B >= tt_0).
#'
#' @examples
#' ##EXAMPLE 1: AMINO ACID SEQUENCE ALIGNMENT
#' ##Set up two sequences:
#' seq_1 <- "ASEDLTI"
#' seq_2 <- "AEEDFGI"
#' ##Set up gap parameters:
#' gapOpen <- 0
#' gapExtend <- -2
#' ##Set up scoring matrix: (Biostrings contains protein scoring matrices.  
#' ##Can also manually set up your own scoring matrix; be sure to name rows and columns.)
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("Biostrings")
#' library(Biostrings) 
#' data(PAM30) ##load PAM30 scoring matrix
#' myScoringMat <- "PAM30"
#' ##Perform alignment and obtain p-value:
#' myAlignment <- pairwiseAlignment(seq_1, seq_2, substitutionMatrix = myScoringMat, 
#'  	gapOpening = gapOpen, gapExtension = gapExtend, type = "global", scoreOnly = FALSE)
#' pval_output <- align_p_val(seq_1, seq_2, myScoringMat, gapOpen, gapExtend, B = 1000, type='global')
#'
#'
#' ##EXAMPLE 2: DNA SEQUENCE ALIGNMENT:
#' seq_1 <- "ACT"
#' seq_2 <- "GCAT"
#' myScoringMat <- matrix(c(3,-2,-2,-1,-2,2,0,-1,-2,0,2,-2,-1,-1,-2,1),nrow=4,byrow=T,
#'	dimnames=list(c('A','C','T','G'),c('A','C','T','G')))
#' gapOpen <- 0
#' gapExtend <- -1
#' ##Perform global alignment:
#' myAlignment <- pairwiseAlignment(seq_1, seq_2, substitutionMatrix = myScoringMat, 
#'  	gapOpening = gapOpen, gapExtension = gapExtend, type = "global", scoreOnly = FALSE)
#' pval_output <- align_p_val(seq_1, seq_2, myScoringMat, gapOpen, gapExtend, B = 1000, type='global')
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

# FUNCTION DEFINITION:  Function for computing a p-value for the alignment of 
# two sequences. Length(seq_1) must be <= Length(seq_2) in the arguments.
# Type must be 'local' or 'global', in quotes.

align_p_val <- function(seq_1, seq_2, scoringMat, gapOpen, gapExtend, B, type) {
  ## Compute alignment score for original pair of sequences.
  tt_0 <- pairwiseAlignment(seq_1, seq_2, substitutionMatrix = scoringMat, 
    gapOpening = gapOpen, gapExtension = gapExtend, type = type, scoreOnly = TRUE)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  		# FUNCTION DEFINITION to generate random permutations of a sequence.
		gen_random_seqs <- function(sq, B) {
 		## Break input sequence into vector of characters.
 		sq <- strsplit(sq, "")[[1]]
  		n <- length(sq)

  		## Compute sample proportions of each observed letter.
  		sq_tbl <- table(sq)
  		sq_letters <- names(sq_tbl)
  		n_letters <- length(sq_letters)
  
  		pp <- numeric()
  		for(i in 1:length(sq_letters)) {
  			pp[i] <- sq_tbl[i] / n
  		}
    	
  		## Generate B random sequences by sampling with replacement using a multinomial model.
  		sqs <- numeric(B)
  		for(i in 1:B){
  			sqs[i] <- paste(sample(sq_letters, n, rep = TRUE, prob = pp), collapse = "")
  		} 
  
  		return(sqs)
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Generate B random versions of first sequence.
  random_seqs <- gen_random_seqs(seq_1, B)
  
  ## For each random sequence, re-run the alignment and store the resulting score.
  tt_B <- numeric(B)
  for(i in 1:B){
  	tt_B[i] <- pairwiseAlignment(random_seqs[i], seq_2, substitutionMatrix = scoringMat, 
    	gapOpening = gapOpen, gapExtension = gapExtend, type = type, scoreOnly = TRUE)
  } 


  ## Compute a p-value as the proportion of randomization-based alignment scores that are 
  ## equal to or greater than our original observed score.
  p_val <- mean(tt_B >= tt_0)
  
  ## Compute an estimated density function of randomized-alignment scores.
  dens <- density(tt_B)
  
  return(list("tt_0" = tt_0, "tt_B" = tt_B, "p_val" = p_val,'density' = dens))
}