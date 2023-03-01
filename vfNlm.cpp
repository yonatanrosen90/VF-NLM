/** @brief Function for applying variance factored non-local means (VF-NLM) to 1D signal in RCPP
* VF-NLM generalizes NLM for spatially varying noise
* for papaer, see: --- 
* Note that sig is a vector of the same size as x
* Applying VF-NLM with constnat noise reduces to the standard NLM
* For applying standard NLM set sig to a constant vector
*
* @param [in] x input noisy signal
* @param [in] h filter parameter
* @param [in] alpha neighborhood parameter, bandwidth of Gaussian kernel
* @param [in] sig vector of the standard deviation of the noise at each point
* @param [in] pad indicator whether signal should be symmertically mirrored along edges
* @return yhat denoised signal after applying non local means
*/



#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector vf_nlf(NumericVector x, double h, double alpha, NumericVector sig, bool pad) {
	double h2 = h*h;
	double a2 = alpha*alpha;
	int N = x.size();
	int b = alpha*6;
	NumericVector yhat(N), galpha(2*b+1);
	NumericVector s2 = sig*sig; // vector of variances of each point
	
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
					// padded indexes
					if (pad) {
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
						tmp_sum_j += (x[ki]-x[kj])*(x[ki]-x[kj])*galpha[l+b]/(s2[ki] + s2[kj]);
					} else {
						continue;
					}
				} else {
					tmp_sum_j += (x[i+l]-x[j+l])*(x[i+l]-x[j+l])*galpha[l+b]/(s2[i+l] + s2[j+l]);
				}
			}
			double wij = exp(-tmp_sum_j/h2)/s2[j];
			tmp_sum_y += wij*x[j];
			sum_wij += wij; // sum for normalizing factor
		}
		yhat[i] = tmp_sum_y/sum_wij;
	}
	
	return yhat;
}



