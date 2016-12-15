#------------------------------------------------------------------
#QR Solver Function:
qr_solver <- function(A,b){
	#Solves linear system Ax=b.
	
	#Obtain QR decomposition of matrix A.  Extract components.
	QR <- qr(A)
	Q <- qr.Q(QR)
	R <- qr.R(QR)
	
	#Backsolve for x.
	x <- qr.solve(A,b)
	return(x)
}