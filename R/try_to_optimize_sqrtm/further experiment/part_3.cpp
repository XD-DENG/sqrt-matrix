#include <Rcpp.h>

using namespace Rcpp;



// [[Rcpp::export]]
ComplexMatrix part_3(int n, int k, List R_index, NumericMatrix S, List Sch_x_EValues, ComplexMatrix X){

	if(n-k > 1){
		for（int j=1; j<n-k; j++）{
			List ij = R_index[j];  // ij may be of length 1 or 2. So here use "List" to assign its type firstly
			for (int i=j-2; i>=0; i--){
				List ii = R_index[i];  // the same reason as "ij" to use "List"
				double sumU = 0;
				if(ij.size() == 1 && ii.size() == 1){
					if (j - i > 1){
						for ( int l = i+1; l<j-1, l++){  // !!!!!!!!in doubt
							List il = R_index[l];
							if(il.size() == 2){
								sumU = sumU + X(ii, il[0]) * x(il[0], ij) + X(ii, il[1]) * x(il[1], ij) 
							}else{
								sumU = sumU + X(ii, il) * X(il, ij)
							}
						}
						X[ii, ij] = 
					}
				}
			}
		}

	}

}

  