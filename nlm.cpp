/** @brief Function for applying non-local means (NLM) to 1D signal in RCPP
* for papaer, see: https://hal.science/hal-00271141/document
*
* @param [in] x input noisy signal
* @param [in] h filter parameter
* @param [in] alpha neighborhood parameter, bandwidth of Gaussian kernel
* @param [in] pad indicator whether signal should be symmertically mirrored along edges
* @return yhat denoised signal after applying non local means
*/


#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector nlf(NumericVector x, double h, double alpha, bool pad) {
	double h2 = h*h;
	double a2 = alpha*alpha;
	int N = x.size();
	int b = alpha*6;
	NumericVector yhat(N), galpha(2*b+1);
	
	// gaussian kernel values for neighborhood similarity
	for (int l=0; l<b+1; l++){
		galpha[l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
		galpha[-l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
	}
	
	int ki, kj; // indexes for padding
	for(int i=0; i<N; i++){
		
		double tmp_sum_y = 0;
		double sum_wij = 0;
		for (int j=0; j<N; j++){
			
			double tmp_sum_j = 0;
			for (int l=-b; l<b+1; l++){
				if (i+l<0 or j+l<0 or i+l>=N or j+l>=N){
					if (pad) {		
						
						// padded indexes
						if (i+l<N) {
							ki = abs(i+l);
						} else {
							ki = (N-1)-abs((N-1)-(i+l));
						}
						if (j+l<N) {
							kj = abs(j+l);
						} else {
							kj = (N-1)-abs((N-1)-(j+l));
						}
						
						tmp_sum_j += (x[ki]-x[kj])*(x[ki]-x[kj])*galpha[l+b]; //sum of neighborhood similarity
					} else {
						continue;
					}
				} else {
					tmp_sum_j += (x[i+l]-x[j+l])*(x[i+l]-x[j+l])*galpha[l+b]; //sum of neighborhood similarity
				}
			}
			double wij = exp(-tmp_sum_j/h2);
			tmp_sum_y += wij*x[j];
			sum_wij += wij;  // sum for normalizing factor
		}
		yhat[i] = tmp_sum_y/sum_wij;
	}
	
	return yhat;
}



