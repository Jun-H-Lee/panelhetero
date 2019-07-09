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

/// 3. Kernel Density Estimation

capture prog drop phkd
program define phkd, eclass
    version 14.0
	syntax varlist(numeric) [if] [in] [, method(string) acov_order(integer 0) acor_order(integer 1) ]
	
	capture drop mean_dest acov_dest acor_dest 
	capture drop mean_grid acov_grid acor_grid
	capture graph drop meandest acovdest acordest
	
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
	    mata: m_hpjkd(data, acov_order, acor_order)
	}
	else if ("`method'" == "toj"){
	    mata: m_tojkd(data, acov_order, acor_order)
	}
	else {
	    mata: m_nekd(data, acov_order, acor_order)
	}
	
    graph twoway line mean_dest mean_grid, ytitle("Density") ///
	                  xtitle("Mean") title("Kernel Density Estimation for Mean")
	graph rename meandest
	
	graph twoway line acov_dest acov_grid, ytitle("Density") ///
	                  xtitle("Autocovariance") title("Kernel Density Estimation for Autocovariance of Order `acov_order'")
	graph rename acovdest
	
	graph twoway line acor_dest acor_grid, ytitle("Density") ///
	                  xtitle("Autocorrelation") title("Kernel Density Estimation for Autocorrelation of Order `acor_order'")
	graph rename acordest
	
end
