/* Stata function for 
Ryo Okui and Takahide Yanagi. Kernel Estimation for Panel Data       
   with Heterogeneous Dynamics. 2019.
Ryo Okui and Takahide Yanagi. Panel Data Analysis with Heterogeneous 
   Dynamics. 2019. */
   
/*
Note : 
    1. Data should be xtset.
    2. Data should be strongly balanced.

Contents :
    1. Empirical CDF Estimaton    : phecdf
	2. Moment Estimation          : phmoment
	3. Kernel Density Estimation  : phkd

*/

/// 1. Empirical CDF Estimation

capture prog drop phecdf
program define phecdf, eclass
    version 14.0
	syntax varlist(numeric) [if] [in] [, method(string) acov_order(integer 0) acor_order(integer 1) ]
	
	capture drop mean_ecdf acov_ecdf acor_ecdf 
	capture drop mean_grid acov_grid acor_grid
	capture graph drop meanecdf acovecdf acorecdf
	
	quietly xtset
	marksample touse
	
	if ("`r(balanced)'" != "strongly balanced"){
	display as error "Error: The given data are not strongly balanced."
	exit
	}
	
	mata: data = st_data(., "`varlist'")
	mata: T = (`r(tmax)' - `r(tmin)')/`r(tdelta)' + 1
	mata: data = colshape(data, T)
	mata: acov_order = `acov_order'
	mata: acor_order = `acor_order'
	
	if ("`method'" == "hpj"){
	    mata: m_hpjecdf(data, acov_order, acor_order)
	}
	else if ("`method'" == "toj"){
	    mata: m_tojecdf(data, acov_order, acor_order)
	}
	else {
	    mata: m_neecdf(data, acov_order, acor_order)
	}
	
    graph twoway line mean_ecdf mean_grid, ytitle("Cumulative Distribution") ///
	                  xtitle("Mean") title("Empirical CDF Estimation for Mean")
	graph rename meanecdf
	
	graph twoway line acov_ecdf acov_grid, ytitle("Cumulative Distribution") ///
	                  xtitle("Autocovariance") title("Empirical CDF Estimation for Autocovariance of Order `acov_order'")
	graph rename acovecdf
	
	graph twoway line acor_ecdf acor_grid, ytitle("Cumulative Distribution") ///
	                  xtitle("Autocorrelation") title("Empirical CDF Estimation for Autocorrelation of Order `acor_order'")
	graph rename acorecdf
	
end
