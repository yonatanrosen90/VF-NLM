#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]


List nlf_sd_pad(NumericVector x, double h, double alpha, bool pad) {
	List res;
	double h2 = h*h;
	double a2 = alpha*alpha;
	int N = x.size();
	int b = alpha*6;
	NumericVector yhat(N), galpha(2*b+1), wi_vec(N), sd(N);
		
	for (int l=0; l<b+1; l++){
		galpha[l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
		galpha[-l+b] = exp(-l*l/(2*a2)) / sqrt(2 * M_PI * a2);
	}
	
	int ki,kj;
	for(int i=0; i<N; i++){
		
		double tmp_sum_y = 0;
		double sum_wij = 0;
		for (int j=0; j<N; j++){
			
			double tmp_sum_j = 0;
			for (int l=-b; l<b+1; l++){
				if (i+l<0 or j+l<0 or i+l>=N or j+l>=N){
					if (pad){
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
						tmp_sum_j += (x[ki]-x[kj])*(x[ki]-x[kj])*galpha[l+b];
						
					} else {
						continue;
					}
				} else {
					tmp_sum_j += (x[i+l]-x[j+l])*(x[i+l]-x[j+l])*galpha[l+b];
				}
			}
			double wij = exp(-tmp_sum_j/h2);
			tmp_sum_y += wij*x[j];
			sum_wij += wij;
			wi_vec[j] = wij;
		}
		yhat[i] = tmp_sum_y/sum_wij;
		
		sd[i] = 0;
		
		for (int j=0; j<N; j++){
			sd[i] += (x[j]-yhat[i])*(x[j]-yhat[i])*wi_vec[j] / sum_wij;
		}
	}
	
	res["yhat"] = yhat;
	res["sd"] = sd;
	
	return res;
}



