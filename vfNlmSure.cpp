/** @brief Function for applying variance factored non-local means (VF-NLM) with Stein's Unbiased Risk Estimate (SURE) for the MSE for 1D signal in RCPP
* VF-NLM generalizes NLM for spatially varying noise
* for paper, see: --- 
* for paper on SURE, see: https://www.jstor.org/stable/2240405
* Note that sig is a vector of the same size as x
* Applying VF-NLM with constnat noise reduces to the standard NLM
* For applying standard NLM with SURE set sig to a constant vector
*
* @param [in] x input noisy signal
* @param [in] h filter parameter
* @param [in] alpha neighborhood parameter, bandwidth of Gaussian kernel
* @param [in] sig vector of the standard deviation of the noise at each point
* @param [in] pad indicator whether signal should be symmertically mirrored along edges
* @return res List containing: yhat denoised signal after applying non local means, sure - SURE for the MSE of yhat
*/

#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]


List vf_nlf_sure(NumericVector x, double h, double alpha, NumericVector sig, bool pad) {
	List res;
	double h2 = h*h;
	double a2 = alpha*alpha;
	int N = x.size();
	int b = alpha*6;
	NumericVector yhat(N);
	NumericVector s2 = sig*sig; // vector of variances of each point
	
	// gaussian kernel values for neighborhood similarity
	NumericVector galpha(N);	
	for (int l=0; l<b+1; l++){
		galpha[l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
		galpha[-l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
	}
	
	int ki,kj; // indexes for padding
	double sure = 0;
	for(int i=0; i<N; i++){
		
		double tmp_sum_i = 0;
		double sum_wij = 0;
		
		// initialize sums for sure
		double div_tx = 0;
		double div_ty = 0;
		double term_x = 0;
		double term_y = 0;
		for (int j=0; j<N; j++){
			
			double tmp_sum_j = 0;
			for (int l=-b; l<b+1; l++){
				if (i+l<0 or j+l<0 or i+l>=N or j+l>=N){
					// if pad - find indexes, else - skip
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
			double w_ij = exp(-tmp_sum_j/h2)/s2[j];
			tmp_sum_i += w_ij*x[j];
			sum_wij += w_ij; // sum for normalizing factor
			
			// calc derivative dy1_dx1
			term_x += w_ij * ( (x[j] - x[i]) / (s2[i] + s2[j]) ) * x[j] / sqrt(2 * M_PI * a2);
			term_y += w_ij * ( (x[j] - x[i]) / (s2[i] + s2[j]) ) / sqrt(2 * M_PI * a2);
			
			int tmp_t = i-j;
			if (tmp_t>=-b and tmp_t<b+1){
				// if pad - find indexes, else and out of bounds - skip
				if ((!pad) and (i+tmp_t >= N or i+tmp_t<0) ){
					continue;
				}
				if (i+tmp_t<N) {
					ki = abs(i+tmp_t);
				} else {
					ki = (N-1)-abs((N-1)-(i+tmp_t));
				}
				
				// sums for sure
				double divals = w_ij * ( (x[ki]-x[i]) / ( s2[ki]+s2[i] ) ) * galpha[tmp_t+b];
				div_tx += divals * x[j];
				div_ty += divals;
			}
			
		}
		yhat[i] = tmp_sum_i/sum_wij;
		
		// sure calculation
		double dy_dx = 1/sum_wij * ( 1/s2[i] + (2/h2)*term_x - (2/h2)*term_y*yhat[i] + (2/h2)*div_tx - (2/h2)*div_ty*yhat[i] );
		sure += (yhat[i]-x[i])*(yhat[i]-x[i]) - s2[i] + 2*s2[i]*dy_dx;
	}
	sure = sure/N;
	
	res["yhat"] = yhat;
	res["sure"] = sure;
	
	return res;
}


