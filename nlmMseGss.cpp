#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
//#include "/sci/home/rosen90/thesis/chrom21/nlmMse.cpp"
#include "C:/Users/yonatan/Documents/university/thesis/rcode/chrom21/nlmMse.cpp"

using namespace Rcpp;

// [[Rcpp::export]]

//golden section search algorithm for nlm

List nlf_mse_gss(NumericVector x, NumericVector true_y,
					double min_init_h, double max_init_h,
					double min_init_a, double max_init_a, 
					double tol_params, double tol_mse, 
					int MaxIter_mse, int MaxIter_params) {
	
	int N = x.size();
	List res;

	int iter=0;
	
	double mse = 0;
	double mse_prev = 0;
	
	double prev_h;
	double prev_a;
	
	NumericVector yhat;
	
	double rho = ( 3 - sqrt(5) ) / 2;
	
	double best_a = max_init_a - (max_init_a - min_init_a) * rho;
	double best_h = 0;
	
	List best_res;
	
	do {
		iter += 1;
		
		mse_prev = mse;
		mse=0;
		
		prev_a = best_a;
		prev_h = best_h;
		
		// optimize h
		
		int h_iter = 0;
		
		double min_h_rg = min_init_h;
		double max_h_rg = max_init_h;
		//double mid_h1 = best_h;
		double mid_h1 = min_init_h + (max_init_h - min_init_h) * rho;
		double mid_h2 = max_init_h - (max_init_h - min_init_h) * rho;
		double mid_h1_val = nlf_mse(x, mid_h1, best_a, true_y)["mse"];
		double mid_h2_val = nlf_mse(x, mid_h2, best_a, true_y)["mse"];
		
		while ( abs(max_h_rg-min_h_rg)>tol_params and h_iter<MaxIter_params ) {
			h_iter += 1;
			if (mid_h2_val > mid_h1_val){
				max_h_rg = mid_h2;
				mid_h2 = mid_h1;
				mid_h2_val = mid_h1_val;
				mid_h1 = min_h_rg + (max_h_rg - min_h_rg)*rho;
				mid_h1_val = nlf_mse(x, mid_h1, best_a, true_y)["mse"];
			} else {
				min_h_rg = mid_h1;
				mid_h1 = mid_h2;
				mid_h1_val = mid_h2_val;
				mid_h2 = max_h_rg - (max_h_rg - min_h_rg)*rho;
				mid_h2_val = nlf_mse(x, mid_h2, best_a, true_y)["mse"];
			}
			
		}
		best_h = min_h_rg * 0.5 + max_h_rg * 0.5;
		
		
		// optimize a
		
		int a_iter = 0;
		
		double min_a_rg = min_init_a;
		double max_a_rg = max_init_a;
		double mid_a1 = min_init_a + (max_init_a - min_init_a) * rho;
		//double mid_a2 = best_a;
		double mid_a2 = max_init_a - (max_init_a - min_init_a) * rho;
		double mid_a1_val = nlf_mse(x, best_h, mid_a1, true_y)["mse"];
		double mid_a2_val = nlf_mse(x, best_h, mid_a2, true_y)["mse"];
		
		while ( abs(max_a_rg-min_a_rg)>tol_params and a_iter<MaxIter_params) {
			a_iter += 1;
			if (mid_a2_val > mid_a1_val){
				max_a_rg = mid_a2;
				mid_a2 = mid_a1;
				mid_a2_val = mid_a1_val;
				mid_a1 = min_a_rg + (max_a_rg - min_a_rg)*rho;
				mid_a1_val = nlf_mse(x, best_h, mid_a1, true_y)["mse"];
			} else {
				min_a_rg = mid_a1;
				mid_a1 = mid_a2;
				mid_a1_val = mid_a2_val;
				mid_a2 = max_a_rg - (max_a_rg - min_a_rg)*rho;
				mid_a2_val = nlf_mse(x, best_h, mid_a2, true_y)["mse"];
			}
			
		}
		best_a = min_a_rg * 0.5 + max_a_rg * 0.5;
		
		best_res = nlf_mse(x, best_h, best_a, true_y);
		mse = best_res["mse"];
	}
	while ( abs(mse_prev-mse)>tol_mse and iter<MaxIter_mse );
	
	res["yhat"] = best_res["yhat"];
	res["niter"] = iter;
	res["h"] = best_h;
	res["a"] = best_a;
	res["prev_h"] = prev_h;
	res["prev_a"] = prev_a;
	res["mse"] = mse;
	res["mse_prev"] = mse_prev;
	
	return res;
}



