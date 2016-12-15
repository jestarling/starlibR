#include <RcppEigen.h>
#include <algorithm>    // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;

typedef MapMatd::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

//#################################################################
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

VectorXd gauss_seidel_solver(MatrixXd A, VectorXd b, VectorXd x0, int iters=100){
	
	int n = A.rows();	//Number of observations.
	VectorXd x = x0;	//Vector to hold x responses.
	
	//Gauss-Seidel iterations.
	for (int k=0; k<iters; k++){	//Loop through number of iterations.
		for(int i=0; i<n; i++){	//Loop through data set once for each iteration.
			//Update each x_i element.
			x(i) = (1/A(i,i)) * (b(i) - A.row(i)*x + A(i,i)*x(i));
		} //end 1:n inner for loop.
	} //end outer iterations (k) for loop.
	return x;
}
