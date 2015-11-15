#include <Rcpp.h>
#include <complex>   // to compute the root square of negative number

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List part_1_plus_2(int n, int k, List R_index, NumericMatrix S, List Sch_x_EValues){ // here I use "List" type for Sch_x_EValues since NumericVector type will cause warning :"Warning message: imaginary parts discarded in coercion"
  

  int l = 0;
  int i = 0;

  while (i < n-1) {
    if(S(i+1, i) == 0){   // in Rcpp, we use ( instead of [ for index of matrix
      R_index[l] = i+1;
    }
    else{
      NumericVector temp(2);
      temp[0] = i+1;
      temp[1] = i+2;
      R_index[l] = temp;
      i = i + 1;
    }
    i = i + 1;
    l = l + 1;        
  }


  List NULL_check =R_index[n-k-1];
  if(NULL_check.size() == 0){
    R_index[n-k-1] = n;
  } 

  ComplexMatrix X(n, n);  

  for (int j = 0; j<(n-k); j++){
    List ij = R_index[j];
    int length_ij = ij.size();

    if(length_ij == 1){
      double  ij = R_index[j];
      float dot_s = S(ij-1, ij-1);
      if(dot_s < 0){
        complex<double> temp_1 (dot_s, 0.0);
        complex<double> temp_result = sqrt(temp_1);
        X(ij-1, ij-1).r = real(temp_result);
        X(ij-1, ij-1).i = imag(temp_result);
      }else{
        X(ij-1, ij-1).r = real(sqrt(dot_s));
        X(ij-1, ij-1).i = 0;
      }
    } else {   // the length of ij =2
      double temp=ij[0];
      complex<double> ev1 = Sch_x_EValues[temp-1];
      double r1 = real(sqrt(ev1));  
      double ij_1 = ij[0];
      double ij_2 = ij[1];
      X(ij_1-1, ij_1-1).r = r1 + 1/(2 * r1) * (S(ij_1-1, ij_1-1) - real(ev1));
      X(ij_1-1, ij_2-1).r = 1/(2 * r1) * (S(ij_1-1, ij_2-1));
      X(ij_2-1, ij_1-1).r = 1/(2 * r1) * (S(ij_2-1, ij_1-1));
      X(ij_2-1, ij_2-1).r = r1 + 1/(2 * r1) * (S(ij_2-1, ij_2-1) - real(ev1));

    }
  }

  List result(2);
  result[0] = R_index;
  result[1] = X;

  return result;


}