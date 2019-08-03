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
	
	capture drop mean_dest
	capture drop acov_dest
	capture drop acor_dest
	capture drop mean_dest
	capture drop acov_dest
	capture drop acor_dest
	capture graph drop meandest
	capture graph drop acovdest
	capture graph drop acordest
	
	local obs = _N
	
	quietly xtset
	marksample touse
	
	if ("`r(balanced)'" != "strongly balanced"){
	    display as error "error: The given data is not xtset or strongly balanced."
	    exit
	}
	
	capture mata: mm_kern("e",0)
	if (_rc != 0) {
	    display as error "error: Cannot find the package MOREMATA." _newline `"The package can be installed with the command "ssc install moremata"."'
	    exit
	}
	
	capture mata: kdens_bw_dpi(runiform(10,1), level = 2)
	if (_rc != 0) {
		display as error "error: Cannot find the package KDENS." _newline `"The package can be installed with the command "ssc install kdens"."'
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
	else if ("`method'" == "sim"){
	    mata: m_nekd(data, acov_order, acor_order)
	}
	else{
	    display as error "error: The name of method is not correctly specified."
		exit
	}
	
    graph twoway line mean_dest mean_grid, ytitle("Density") ///
	                  xtitle("Mean") title("Kernel Density Estimation for Mean")
	graph rename meandest
	
	if (`acov_order' == 0){
	    graph twoway line acov_dest acov_grid, ytitle("Density") ///
	                  xtitle("Variance") title("Kernel Density Estimation for Variance")
	    graph rename acovdest
	}
	else{
	    graph twoway line acov_dest acov_grid, ytitle("Density") ///
	                  xtitle("Autocovariance") title("Kernel Density Estimation for Autocovariance of Order `acov_order'")
	    graph rename acovdest
	}
	
	graph twoway line acor_dest acor_grid, ytitle("Density") ///
	                  xtitle("Autocorrelation") title("Kernel Density Estimation for Autocorrelation of Order `acor_order'")
	graph rename acordest
	
	drop mean_dest acov_dest acor_dest mean_grid acov_grid acor_grid
	quietly keep in 1/`obs'
	
end
