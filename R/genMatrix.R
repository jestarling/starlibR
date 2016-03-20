#----------------------------------------------------------------
# Roxy package build comments:
#' genMatrix(type,n=5)
#'
#' This function generates examples of various square matrix types.
#'
#' @param type The type of the example matrix to generate.  Required field; default is 'identity'.
#' @param type='zero' to generate a zero matrix.
#' @param type='identity' to generate an identity matrix.
#' @param type='ut' to generate an upper triangular matrix.
#' @param type='uut' to generate a unit upper triangular matrix.
#' @param type='ust' to generate a strictly upper triangular matrix.
#' @param type='lt' to generate a lower triangular matrix.
#' @param type='lut' to generate a unit lower triangular matrix.
#' @param type='lst' to generate a strictly lower triangular matrix.
#' @param n The number of rows and colums of the square matrix (nxn).  Default is n=5.
#'
#' @keywords matrix
#'
#' @return mat A matrix of the type and dimension specified.
#'
#' @examples
#' genMatrix('identity',5)
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------
genMatrix <- function(type='identity',n=5){
	
	#types: zero, identity, diagonal,uTriangular,lTriangular,
	#	uuTriangular,luTriangular,usTriangular,lsTriangular
	
	mat = NULL #Generate NULL value to hold matrix result.
	
	#Zero matrix.
	if(type == 'zero') { mat = matrix(rep(0,n*n),nrow=n,byrow=T) }
	
	#Identity matrix.
	if(type == 'identity') { mat = diag(n) }
	
	#Diagonal matrix (random normal diagonal).
	if(type == 'diagonal') { mat = diag(rnorm(n)) }
	
	#Upper triangular
	if(type == 'ut') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[lower.tri(mat)] <- 0
	}
	
	#Unit Upper triangular
	if(type == 'uut') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[lower.tri(mat)] <- 0
		diag(mat) <- 1
	}
	
	#Strictly Upper triangular
	if(type == 'ust') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[lower.tri(mat,diag=T)] <- 0
	}
	
	#Lower triangular
	if(type == 'lt') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[upper.tri(mat)] <- 0
	}
	
	#Unit Lower triangular
	if(type == 'lut') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[upper.tri(mat)] <- 0
		diag(mat) <- 1
	}
	
	#Strictly Lower triangular
	if(type == 'lst') { 
		mat = matrix(rnorm(n*n),nrow=n)
		mat[upper.tri(mat,diag=T)] <- 0
	}
	
	return(mat)
}