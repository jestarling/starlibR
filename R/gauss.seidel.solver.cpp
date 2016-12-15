#include <RcppArmadillo.h>
//#include <cmath.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
vec gauss_seidel_solver(mat A, vec b, vec x0, int iters=100){
	
	int n = A.n_rows;
	vec x = x0;
	
	for (int k=0; k<iters; k++){
		for(int i=0; i<n; i++){
			x(i) = (1/A(i,i)) * (b(i) - A.row(i)*x + A(i,i)*x(i));
		}
	}	
	return x;
};

/*


function [x] = gauss_seidel(A, b, x0, iters)
    n = length(A);
    x = x0;
    for k = 1:iters
        for i = 1:n
            x(i) = (1/A(i, i))*(b(i) - A(i, 1:n)*x + A(i, i)*x(i));
        end
    end
end

*/

