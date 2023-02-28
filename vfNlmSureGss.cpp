#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include <math.h>
//#include "/sci/home/rosen90/thesis/chrom21/vfNlmSure.cpp"
#include "C:/Users/yonatan/Documents/university/thesis/rcode/chrom21/vfNlmSure.cpp"

using namespace Rcpp;

// [[Rcpp::export]]

//golden section search algorithm for nlm

List vf_nlf_sure_gss(NumericVector x, NumericVector sig, 
					double min_init_h, double max_init_h,
					double min_init_a, double max_init_a, 
					double tol_params, double tol_sure, 
					int MaxIter_sure, int MaxIter_params) {
	
	int N = x.size();
	List res;

	int iter=0;
	int h_iter = 0;
	int a_iter=0;
	
	double sure = 0;
	double sure_prev = 0;
	
	double prev_h = 0;
	double prev_a = 0;
	NumericVector yhat;
	
	double rho = ( 3 - sqrt(5) ) / 2;
	
	double best_a = max_init_a - (max_init_a - min_init_a) * rho;
	double best_h = 0;
	
	List best_res;
	
	do {
		iter += 1;
		
		sure_prev = sure;
		sure=0;
		
		prev_a = best_a;
		prev_h = best_h;
		
		// optimize h
		
		h_iter = 0;
		
		double min_h_rg = min_init_h;
		double max_h_rg = max_init_h;
		//double mid_h1 = best_h;
		double mid_h1 = min_init_h + (max_init_h - min_init_h) * rho;
		double mid_h2 = max_init_h - (max_init_h - min_init_h) * rho;
		double mid_h1_val = vf_nlf_sure(x, mid_h1, best_a, sig)["sure"];
		double mid_h2_val = vf_nlf_sure(x, mid_h2, best_a, sig)["sure"];
		
		while ( abs(max_h_rg-min_h_rg)>tol_params and h_iter<MaxIter_params ) {
			h_iter += 1;
			
			if (mid_h2_val > mid_h1_val){
				max_h_rg = mid_h2;
				mid_h2 = mid_h1;
				mid_h2_val = mid_h1_val;
				mid_h1 = min_h_rg + (max_h_rg - min_h_rg)*rho;
				mid_h1_val = vf_nlf_sure(x, mid_h1, best_a, sig)["sure"];
			} else {
				min_h_rg = mid_h1;
				mid_h1 = mid_h2;
				mid_h1_val = mid_h2_val;
				mid_h2 = max_h_rg - (max_h_rg - min_h_rg)*rho;
				mid_h2_val = vf_nlf_sure(x, mid_h2, best_a, sig)["sure"];
			}
			
		}
		best_h = min_h_rg * 0.5 + max_h_rg * 0.5;
		
		// optimize a
		
		a_iter = 0;
		
		double min_a_rg = min_init_a;
		double max_a_rg = max_init_a;
		double mid_a1 = min_init_a + (max_init_a - min_init_a) * rho;
		//double mid_a2 = best_a;
		double mid_a2 = max_init_a - (max_init_a - min_init_a) * rho;
		double mid_a1_val = vf_nlf_sure(x, best_h, mid_a1, sig)["sure"];
		double mid_a2_val = vf_nlf_sure(x, best_h, mid_a2, sig)["sure"];
		
		while ( abs(max_a_rg-min_a_rg)>tol_params and a_iter<MaxIter_params) {
			a_iter += 1;
			if (mid_a2_val > mid_a1_val){
				max_a_rg = mid_a2;
				mid_a2 = mid_a1;
				mid_a2_val = mid_a1_val;
				mid_a1 = min_a_rg + (max_a_rg - min_a_rg)*rho;
				mid_a1_val = vf_nlf_sure(x, best_h, mid_a1, sig)["sure"];
			} else {
				min_a_rg = mid_a1;
				mid_a1 = mid_a2;
				mid_a1_val = mid_a2_val;
				mid_a2 = max_a_rg - (max_a_rg - min_a_rg)*rho;
				mid_a2_val = vf_nlf_sure(x, best_h, mid_a2, sig)["sure"];
			}
			
		}
		best_a = min_a_rg * 0.5 + max_a_rg * 0.5;
		
		best_res = vf_nlf_sure(x, best_h, best_a, sig);
		sure = best_res["sure"];
		
	}
	while ( abs(sure_prev-sure)>tol_sure and iter<MaxIter_sure );
	
	res["yhat"] = best_res["yhat"];
	res["niter"] = iter;
	res["h_niter"] = h_iter;
	res["a_niter"] = a_iter;
	res["h"] = best_h;
	res["a"] = best_a;
	res["prev_h"] = prev_h;
	res["prev_a"] = prev_a;
	res["sure"] = sure;
	res["sure_prev"] = sure_prev;
	
	return res;
}



