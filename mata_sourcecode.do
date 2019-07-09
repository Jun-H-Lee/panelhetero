/* Mata function library for 
Ryo Okui and Takahide Yanagi. Kernel Estimation for Panel Data       
   with Heterogeneous Dynamics. 2019.
Ryo Okui and Takahide Yanagi. Panel Data Analysis with Heterogeneous 
   Dynamics. 2019. */

/* Contents
1. Basic
  1.1. Autocovariance
  1.2. Autocorrelation
2. Empirical CDF Estimation
  2.1. Naive Estimation
  2.2. Half-Panel-Jackknife Estimation
  2.3. Third-Order-Jackknife Estimation
3. Moment Estimation
  3.1. Naive Estimation
  3.2. Half-Panel-Jackknife Estimation
  3.3. Third-Order-Jackknife Estimation
4. Kernel Density Estimation
  4.1. Naive Estimation
  4.2. Half-Panel-Jackknife Estimation
  4.3. Third-Order-Jackknife Estimation
*/

version 14.0
mata:
mata clear

///  1. Basic

// 1.1. Autocovariance

function mataacov (y, acov_order){
	k = acov_order
	s = length(y)
	mean_est = sum(y) / s
	y1 = y[1::s-k]
	y2 = y[k+1::s]
	acov_est = sum((y1 :- mean_est) :* (y2 :- mean_est)) / (s - k)
	return(acov_est)
}

// 1.2. Autocorrelation

function mataacor (y, acov_order){
	acor_est = mataacov(y, acov_order) / mataacov(y, 0)
	return(acor_est)
}

///  2. Empirical CDF estimation

function m_freq(x, X){
    est = mean(X:<=x)
    return(est)
}

// 2.1. Naive Estimation 
function m_neecdf(data, acov_order, acor_order) {
    N = rows(data)
	S = cols(data)
	grid = 400
	
	mean_est = J(N, 1, 0)
	acov_est = J(N, 1, 0)
	acor_est = J(N, 1, 0)
	
	for (i = 1; i <= N; i++) {
	    mean_est[i] = mean(data[i,.]')
		acov_est[i] = mataacov(data[i,.], acov_order)
		acor_est[i] = mataacor(data[i,.], acor_order)
	}
	
	mean_lim = (min(mean_est), max(mean_est))
    acov_lim = (min(acov_est), max(acov_est))
    acor_lim = (min(acor_est), max(acor_est))
    
	mean_grid = rangen(mean_lim[1], mean_lim[2], grid)
	acov_grid = rangen(acov_lim[1], acov_lim[2], grid)
	acor_grid = rangen(acor_lim[1], acor_lim[2], grid)
	
	mean_ecdf = J(grid, 1, 0)
	acov_ecdf = J(grid, 1, 0)
	acor_ecdf = J(grid, 1, 0)

	for (i=1; i<=grid; i++) {
	    mean_ecdf[i] = m_freq(mean_grid[i], mean_est)
	    acov_ecdf[i] = m_freq(acov_grid[i], acov_est)
		acor_ecdf[i] = m_freq(acor_grid[i], acor_est)
	}
	
	temp=st_addvar("double", "mean_ecdf")
    temp=st_addvar("double", "acov_ecdf")
    temp=st_addvar("double", "acor_ecdf")
    temp=st_addvar("double", "mean_grid")
	temp=st_addvar("double", "acov_grid")
	temp=st_addvar("double", "acor_grid")
	
    st_addobs(max((0,grid  - st_nobs())))
    st_store(.,"mean_ecdf", mean_ecdf\J(st_nobs()-rows(mean_ecdf),1,.))
    st_store(.,"acov_ecdf", acov_ecdf\J(st_nobs()-rows(acov_ecdf),1,.))
    st_store(.,"acor_ecdf", acor_ecdf\J(st_nobs()-rows(acor_ecdf),1,.))
	
    st_store(.,"mean_grid", mean_grid\J(st_nobs()-rows(mean_grid),1,.))
	st_store(.,"acov_grid", acov_grid\J(st_nobs()-rows(acov_grid),1,.))
    st_store(.,"acor_grid", acor_grid\J(st_nobs()-rows(acor_grid),1,.))
}


// 2.2. Half-Panel-Jackknife Estimation

function hpjecdfest1(x, X, X1, X2) {
    est = m_freq(x, X)
	est1 = m_freq(x, X1)
    est2 = m_freq(x, X2)
	
	hpjest = 2 * est - (est1 + est2) / 2
	
	if (hpjest < 0) {
	    hpjest = 0
	}
	if (hpjest > 1){
	    hpjest = 1
	}
    return(hpjest)
}

function hpjecdfest2(x, X, X1, X2, X3, X4) {
    est = m_freq(x, X)
	est1 = m_freq(x, X1)
    est2 = m_freq(x, X2)
	est3 = m_freq(x, X3)
    est4 = m_freq(x, X4)
	
	
	hpjest = 2 * est - (est1 + est2 + est3 + est4) / 4
	
	if (hpjest < 0) {
	    hpjest = 0
	}
	if (hpjest > 1){
	    hpjest = 1
	}
    return(hpjest)
}

function m_hpjecdf(data, acov_order, acor_order) {
    N = rows(data)
	S = cols(data)
	grid = 400

	if (mod(S,2)==0) {
	    mean_est = J(N, 1, 0)
	    acov_est = J(N, 1, 0)
	    acor_est = J(N, 1, 0)
	
	    data1 = data[,1::(S/2)]
	    data2 = data[,(S/2+1)::S]
		
	    mean_est1 = J(N, 1, 0)
	    mean_est2 = J(N, 1, 0)
		
	    acov_est1 = J(N, 1, 0)
	    acor_est1 = J(N, 1, 0)
	    
		acov_est2 = J(N, 1, 0)
	    acor_est2 = J(N, 1, 0)
		
		for(i=1; i<=N ;i++){
		    mean_est[i] = mean(data[i,.]')
			acov_est[i] = mataacov(data[i,.], acov_order)
			acor_est[i] = mataacor(data[i,.], acor_order)
			
		    mean_est1[i] = mean(data1[i,.]')
			acov_est1[i] = mataacov(data1[i,.], acov_order)
			acor_est1[i] = mataacor(data1[i,.], acor_order)
			
			mean_est2[i] = mean(data2[i,.]')
			acov_est2[i] = mataacov(data2[i,.], acov_order)
			acor_est2[i] = mataacor(data2[i,.], acor_order)
		}
		
		mean_lim = (min(mean_est), max(mean_est))
        acov_lim = (min(acov_est), max(acov_est))
        acor_lim = (min(acor_est), max(acor_est))
    
	    mean_grid = rangen(mean_lim[1], mean_lim[2], grid)
	    acov_grid = rangen(acov_lim[1], acov_lim[2], grid)
	    acor_grid = rangen(acor_lim[1], acor_lim[2], grid)
	
	    mean_ecdf = J(grid, 1, 0)
	    acov_ecdf = J(grid, 1, 0)
	    acor_ecdf = J(grid, 1, 0)
		
		for (i=1; i<= grid; i++){
		    mean_ecdf[i] = hpjecdfest1(mean_grid[i], mean_est, mean_est1, mean_est2)
			acov_ecdf[i] = hpjecdfest1(acov_grid[i], acov_est, acov_est1, acov_est2)
			acor_ecdf[i] = hpjecdfest1(acor_grid[i], acor_est, acor_est1, acor_est2)
		}
		mean_ecdf = sort(mean_ecdf, 1)
		acov_ecdf = sort(acov_ecdf, 1)
		acor_ecdf = sort(acor_ecdf, 1)
	}	
	else{
	    mean_est = J(N, 1, 0)
	    acov_est = J(N, 1, 0)
	    acor_est = J(N, 1, 0)
		
		data1 = data[., 1::floor(S/2)]
		data2 = data[., (floor(S/2)+1)::S]
		data3 = data[.,1::(floor(S/2)+1)]
		data4 = data[.,(floor(S/2)+2)::S]
		
		mean_est1 = J(N, 1, 0)
		mean_est2 = J(N, 1, 0)
		mean_est3 = J(N, 1, 0)
		mean_est4 = J(N, 1, 0)
		
		acov_est1 = J(N, 1, 0)
		acov_est2 = J(N, 1, 0)
		acov_est3 = J(N, 1, 0)
		acov_est4 = J(N, 1, 0)
		
		acor_est1 = J(N, 1, 0)
		acor_est2 = J(N, 1, 0)
		acor_est3 = J(N, 1, 0)
		acor_est4 = J(N, 1, 0)
		
		for(i=1; i<=N; i++){
		    mean_est[i] = mean(data[i,.]')
			mean_est1[i] = mean(data1[i,.]')
			mean_est2[i] = mean(data2[i,.]')
			mean_est3[i] = mean(data3[i,.]')
			mean_est4[i] = mean(data4[i,.]')
			
			acov_est[i] = mataacov(data[i,.], acov_order)
			acov_est1[i] = mataacov(data1[i,.], acov_order)
			acov_est2[i] = mataacov(data2[i,.], acov_order)
			acov_est3[i] = mataacov(data3[i,.], acov_order)
			acov_est4[i] = mataacov(data4[i,.], acov_order)
	
			acor_est[i] = mataacor(data[i,.], acor_order)
			acor_est1[i] = mataacor(data1[i,.], acor_order)
			acor_est2[i] = mataacor(data2[i,.], acor_order)
			acor_est3[i] = mataacor(data3[i,.], acor_order)
			acor_est4[i] = mataacor(data4[i,.], acor_order)
		} 
		
		mean_lim = (min(mean_est), max(mean_est))
        acov_lim = (min(acov_est), max(acov_est))
        acor_lim = (min(acor_est), max(acor_est))
    
	    mean_grid = rangen(mean_lim[1], mean_lim[2], grid)
	    acov_grid = rangen(acov_lim[1], acov_lim[2], grid)
	    acor_grid = rangen(acor_lim[1], acor_lim[2], grid)
	
	    mean_ecdf = J(grid, 1, 0)
	    acov_ecdf = J(grid, 1, 0)
	    acor_ecdf = J(grid, 1, 0)
		
		for (i=1; i<= grid; i++){
		    mean_ecdf[i] = hpjecdfest2(mean_grid[i], mean_est, mean_est1, mean_est2, mean_est3, mean_est4)
			acov_ecdf[i] = hpjecdfest2(acov_grid[i], acov_est, acov_est1, acov_est2, acov_est3, acov_est4)
			acor_ecdf[i] = hpjecdfest2(acor_grid[i], acor_est, acor_est1, acor_est2, acor_est3, acov_est4)
		}
		
		mean_ecdf = sort(mean_ecdf, 1)
		acov_ecdf = sort(acov_ecdf, 1)
		acor_ecdf = sort(acor_ecdf, 1)
	}
	
	temp=st_addvar("double", "mean_ecdf")
    temp=st_addvar("double", "acov_ecdf")
    temp=st_addvar("double", "acor_ecdf")
    temp=st_addvar("double", "mean_grid")
	temp=st_addvar("double", "acov_grid")
	temp=st_addvar("double", "acor_grid")
	
    st_addobs(max((0,grid  - st_nobs())))
    st_store(.,"mean_ecdf", mean_ecdf\J(st_nobs()-rows(mean_ecdf),1,.))
    st_store(.,"acov_ecdf", acov_ecdf\J(st_nobs()-rows(acov_ecdf),1,.))
    st_store(.,"acor_ecdf", acor_ecdf\J(st_nobs()-rows(acor_ecdf),1,.))
	
    st_store(.,"mean_grid", mean_grid\J(st_nobs()-rows(mean_grid),1,.))
	st_store(.,"acov_grid", acov_grid\J(st_nobs()-rows(acov_grid),1,.))
    st_store(.,"acor_grid", acor_grid\J(st_nobs()-rows(acor_grid),1,.))
}

// 2.3. Third-Order-Jackknife Estimation

// To be done later.


///  3. Moment Estimation

// 3.1. Naive Estimation

function momentest (quantity, indices) {
	mean_est = quantity[indices, 1]
	acov_est = quantity[indices, 2]
	acor_est = quantity[indices, 3]
	
	mean_mean = sum(mean_est) / length(mean_est)
	acov_mean = sum(acov_est) / length(acov_est)
	acor_mean = sum(acor_est) / length(acor_est)
	
	mean_var = variance(mean_est)
	acov_var = variance(acov_est)
	acor_var = variance(acor_est)
	
	mean_acov_cor = correlation((mean_est, acov_est))[2,1]
	mean_acor_cor = correlation((mean_est, acor_est))[2,1]
	acov_acor_cor = correlation((acov_est, acor_est))[2,1]
	
	estimate = (mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
	
	return(estimate)
}

function m_nemoment (data, acov_order, acor_order, B) {
	N = rows(data)
	S = cols(data)
	level = 0.95
	mean_est = J(N, 1, 0)
	acov_est = J(N, 1, 0)
	acor_est = J(N, 1, 0)
	
	for (i = 1; i <= N; i++) {
	    mean_est[i] = mean(data[i,.]')
		acov_est[i] = mataacov(data[i,.], acov_order)
		acor_est[i] = mataacor(data[i,.], acor_order)
	}
	equantity = (mean_est, acov_est, acor_est)
	
	estimate_value = momentest(equantity, 1::N)
	number_par = length(estimate_value)
	estimate_boot = J(B, number_par, 0)
	for (b = 1; b <= B; b++) {
		index_boot = rdiscrete(N, 1, J(N, 1, 1/N))
		estimate_boot[b,.] = momentest(equantity, index_boot)
	}
	
	se = sqrt(diagonal(variance(estimate_boot)))
	quantile_boot_1 = mm_quantile(estimate_boot, 1, 0.975)
	quantile_boot_2 = mm_quantile(estimate_boot, 1, 0.025)
	ci_1 = 2 * estimate_value - quantile_boot_1
	ci_2 = 2 * estimate_value - quantile_boot_2
	ci = (ci_1 \ ci_2)
	result =  (estimate_value', se, ci')
	
	st_matrix("est", estimate_value')
	st_matrix("se", se)
	st_matrix("ci", ci')
	
	printf("\n")
    printf("Estimates for Moments.\n")
    printf("Parameters                                              Estimate          \n")
    printf("________________________________________________________________ \n")
    printf("Mean of Mean                                            %f\n",estimate_value[1])
	printf("Mean of Autocovariance                                  %f\n",estimate_value[2])
	printf("Mean of Autocorrelation                                 %f\n",estimate_value[3])
	printf("Variance of Mean                                        %f\n",estimate_value[4])
	printf("Variance of Autocovariance                              %f\n",estimate_value[5])
	printf("Variance of Autocorrelation                             %f\n",estimate_value[6])
	printf("Correlation between Mean and Autocovariance             %f\n",estimate_value[7])
	printf("Correlation between Mean and Autocorrelation            %f\n",estimate_value[8])
	printf("Correlation between Autocovariance and Autocorelation   %f\n",estimate_value[9])
	
	printf("\n")
    printf("%f %% Confidence Intervals for Moments.\n", level*100)
    printf("Parameters                                              Low                   High\n")
    printf("__________________________________________________________________________________\n")
    printf("Mean of Mean                                            %f        %f\n",ci_1[1],ci_2[1])
	printf("Mean of Autocovariance                                  %f        %f\n",ci_1[2],ci_2[2])
	printf("Mean of Autocorrelation                                 %f        %f\n",ci_1[3],ci_2[3])
	printf("Variance of Mean                                        %f        %f\n",ci_1[4],ci_2[4])
	printf("Variance of Autocovariance                              %f        %f\n",ci_1[5],ci_2[5])
	printf("Variance of Autocorrelation                             %f        %f\n",ci_1[6],ci_2[6])
	printf("Correlation between Mean and Autocovariance             %f        %f\n",ci_1[7],ci_2[7])
	printf("Correlation between Mean and Autocorrelation            %f        %f\n",ci_1[8],ci_2[8])
	printf("Correlation between Autocovariance and Autocorelation   %f        %f\n",ci_1[9],ci_2[9])
	
	printf("\n")
    printf("Standard Errors for Moments.\n")
    printf("Parameters                                              Stanadard Errors          \n")
    printf("________________________________________________________________________\n")
    printf("Mean of Mean                                            %f\n",se[1])
	printf("Mean of Autocovariance                                  %f\n",se[2])
	printf("Mean of Autocorrelation                                 %f\n",se[3])
	printf("Variance of Mean                                        %f\n",se[4])
	printf("Variance of Autocovariance                              %f\n",se[5])
	printf("Variance of Autocorrelation                             %f\n",se[6])
	printf("Correlation between Mean and Autocovariance             %f\n",se[7])
	printf("Correlation between Mean and Autocorrelation            %f\n",se[8])
	printf("Correlation between Autocovariance and Autocorelation   %f\n",se[9])
	
}

// 3.2. Half-Panel-Jackknife Moment Estimation

function hpjmomentest1(quantity, indices){
	mean_est = quantity[indices, 1]
	mean_est1 = quantity[indices, 2]
	mean_est2 = quantity[indices, 3]
	
	acov_est = quantity[indices, 4]
	acov_est1 = quantity[indices, 5]
	acov_est2 = quantity[indices, 6]
	
	acor_est = quantity[indices, 7]
	acor_est1 = quantity[indices, 8]
	acor_est2 = quantity[indices, 9]
	
	mean_mean = sum(mean_est) / length(mean_est)
	mean_mean1 = sum(mean_est1) / length(mean_est1)
	mean_mean2 = sum(mean_est2) / length(mean_est2)

	acov_mean = sum(acov_est) / length(acov_est)
	acov_mean1 = sum(acov_est1) / length(acov_est1)
	acov_mean2 = sum(acov_est2) / length(acov_est2)

	acor_mean = sum(acor_est) / length(acor_est)
	acor_mean1 = sum(acor_est1) / length(acor_est1)
	acor_mean2 = sum(acor_est2) / length(acor_est2)
	
	mean_var = variance(mean_est)
	mean_var1 = variance(mean_est1)
	mean_var2 = variance(mean_est2)
	
	acov_var = variance(acov_est)
	acov_var1 = variance(acov_est1)
	acov_var2 = variance(acov_est2)
	
	acor_var = variance(acor_est)
	acor_var1 = variance(acor_est1)
	acor_var2 = variance(acor_est2)
	
	mean_acov_cor = correlation((mean_est, acov_est))[2,1]
	mean_acov_cor1 = correlation((mean_est1, acov_est1))[2,1]
	mean_acov_cor2 = correlation((mean_est2, acov_est2))[2,1]
	
	mean_acor_cor = correlation((mean_est, acor_est))[2,1]
	mean_acor_cor1 = correlation((mean_est1, acor_est1))[2,1]
	mean_acor_cor2 = correlation((mean_est2, acor_est2))[2,1]
	
	acov_acor_cor = correlation((acov_est, acor_est))[2,1]
	acov_acor_cor1 = correlation((acov_est1, acor_est1))[2,1]
	acov_acor_cor2 = correlation((acov_est2, acor_est2))[2,1]
	
	estimate_1 = (mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
	estimate_2 = (mean_mean1, acov_mean1, acor_mean1, mean_var1, acov_var1, acor_var1, mean_acov_cor1, mean_acor_cor1, acov_acor_cor1)
	estimate_3 = (mean_mean2, acov_mean2, acor_mean2, mean_var2, acov_var2, acor_var2, mean_acov_cor2, mean_acor_cor2, acov_acor_cor2)

	hpjestimate = 2 * estimate_1 - (estimate_2 + estimate_3) / 2
	return(hpjestimate)
}

function hpjmomentest2(quantity, indices){
	mean_est = quantity[indices, 1]
	mean_est1 = quantity[indices, 2]
	mean_est2 = quantity[indices, 3]
	mean_est3 = quantity[indices, 4]
	mean_est4 = quantity[indices, 5]
	
	acov_est = quantity[indices, 6]
	acov_est1 = quantity[indices, 7]
	acov_est2 = quantity[indices, 8]
	acov_est3 = quantity[indices, 9]
	acov_est4 = quantity[indices, 10]
	
	acor_est = quantity[indices, 11]
	acor_est1 = quantity[indices, 12]
	acor_est2 = quantity[indices, 13]
	acor_est3 = quantity[indices, 14]
	acor_est4 = quantity[indices, 15]
	
	mean_mean = sum(mean_est)/length(mean_est)
	mean_mean1 = sum(mean_est1)/length(mean_est1)
	mean_mean2 = sum(mean_est2)/length(mean_est2)
	mean_mean3 = sum(mean_est3)/length(mean_est3)
	mean_mean4 = sum(mean_est4)/length(mean_est4)

	acov_mean = sum(acov_est)/length(acov_est)
	acov_mean1 = sum(acov_est1)/length(acov_est1)
	acov_mean2 = sum(acov_est2)/length(acov_est2)
	acov_mean3 = sum(acov_est3)/length(acov_est3)
	acov_mean4 = sum(acov_est4)/length(acov_est4)

	acor_mean = sum(acor_est)/length(acor_est)
	acor_mean1 = sum(acor_est1)/length(acor_est1)
	acor_mean2 = sum(acor_est2)/length(acor_est2)
	acor_mean3 = sum(acor_est3)/length(acor_est3)
	acor_mean4 = sum(acor_est4)/length(acor_est4)
	
	mean_var = variance(mean_est)
	mean_var1 = variance(mean_est1)
	mean_var2 = variance(mean_est2)
	mean_var3 = variance(mean_est3)
	mean_var4 = variance(mean_est4)
	
	acov_var = variance(acov_est)
	acov_var1 = variance(acov_est1)
	acov_var2 = variance(acov_est2)
	acov_var3 = variance(acov_est3)
	acov_var4 = variance(acov_est4)
	
	acor_var = variance(acor_est)
	acor_var1 = variance(acor_est1)
	acor_var2 = variance(acor_est2)
	acor_var3 = variance(acor_est3)
	acor_var4 = variance(acor_est4)
	
	mean_acov_cor = correlation((mean_est, acov_est))[2,1]
	mean_acov_cor1 = correlation((mean_est1, acov_est1))[2,1]
	mean_acov_cor2 = correlation((mean_est2, acov_est2))[2,1]
	mean_acov_cor3 = correlation((mean_est3, acov_est3))[2,1]
	mean_acov_cor4 = correlation((mean_est4, acov_est4))[2,1]
	
	mean_acor_cor = correlation((mean_est, acor_est))[2,1]
	mean_acor_cor1 = correlation((mean_est1, acor_est1))[2,1]
	mean_acor_cor2 = correlation((mean_est2, acor_est2))[2,1]
	mean_acor_cor3 = correlation((mean_est3, acor_est3))[2,1]
	mean_acor_cor4 = correlation((mean_est4, acor_est4))[2,1]
	
	acov_acor_cor = correlation((acov_est, acor_est))[2,1]
	acov_acor_cor1 = correlation((acov_est1, acor_est1))[2,1]
	acov_acor_cor2 = correlation((acov_est2, acor_est2))[2,1]
	acov_acor_cor3 = correlation((acov_est3, acor_est3))[2,1]
	acov_acor_cor4 = correlation((acov_est4, acor_est4))[2,1]
	
	estimate_1 = (mean_mean, acov_mean, acor_mean, mean_var, acov_var, acor_var, mean_acov_cor, mean_acor_cor, acov_acor_cor)
	estimate_2 = (mean_mean1, acov_mean1, acor_mean1, mean_var1, acov_var1, acor_var1, mean_acov_cor1, mean_acor_cor1, acov_acor_cor1)
	estimate_3 = (mean_mean2, acov_mean2, acor_mean2, mean_var2, acov_var2, acor_var2, mean_acov_cor2, mean_acor_cor2, acov_acor_cor2)
	estimate_4 = (mean_mean3, acov_mean3, acor_mean3, mean_var3, acov_var3, acor_var3, mean_acov_cor3, mean_acor_cor3, acov_acor_cor3)
	estimate_5 = (mean_mean4, acov_mean4, acor_mean4, mean_var4, acov_var4, acor_var4, mean_acov_cor4, mean_acor_cor4, acov_acor_cor4)

	hpjestimate = 2 * estimate_1 - (estimate_2 + estimate_3 + estimate_4 + estimate_5) / 4
	return(hpjestimate)
}

function m_hpjmoment(data, acov_order, acor_order, B){
	N = rows(data)
	S = cols(data)
	level = 0.95
	
    if (mod(S,2)==0) {
	    mean_est = J(N, 1, 0)
	    acov_est = J(N, 1, 0)
	    acor_est = J(N, 1, 0)
	
	    data1 = data[,1::(S/2)]
	    data2 = data[,(S/2+1)::S]
		
	    mean_est1 = J(N, 1, 0)
	    mean_est2 = J(N, 1, 0)
		
	    acov_est1 = J(N, 1, 0)
	    acor_est1 = J(N, 1, 0)
	    
		acov_est2 = J(N, 1, 0)
	    acor_est2 = J(N, 1, 0)
		
		for(i=1; i<=N ;i++){
		    mean_est[i] = mean(data[i,.]')
			acov_est[i] = mataacov(data[i,.], acov_order)
			acor_est[i] = mataacor(data[i,.], acor_order)
			
		    mean_est1[i] = mean(data1[i,.]')
			acov_est1[i] = mataacov(data1[i,.], acov_order)
			acor_est1[i] = mataacor(data1[i,.], acor_order)
			
			mean_est2[i] = mean(data2[i,.]')
			acov_est2[i] = mataacov(data2[i,.], acov_order)
			acor_est2[i] = mataacor(data2[i,.], acor_order)
		}
		equantity2 = (mean_est, mean_est1, mean_est2, acov_est, acov_est1, acov_est2, acor_est, acor_est1, acor_est2)
		estimate_value = hpjmomentest1(equantity2, 1::N)
	
		number_par = length(estimate_value)
		estimate_boot = J(B, number_par, 0)
		
		for (b = 1; b <= B; b++) {
			index_boot = rdiscrete(N, 1, J(N, 1, 1/N))
			estimate_boot[b,] = hpjmomentest1(equantity2, index_boot)
		}
		
		se = sqrt(diagonal(variance(estimate_boot)))
		quantile_boot_1 = mm_quantile(estimate_boot, 1, 0.975)
		quantile_boot_2 = mm_quantile(estimate_boot, 1, 0.025)	
		ci_1 = 2 * estimate_value - quantile_boot_1
		ci_2 = 2 * estimate_value - quantile_boot_2
		ci = (ci_1 \ ci_2)
		result =  (estimate_value', se, ci')
	}
		
	else{
	    mean_est = J(N, 1, 0)
	    acov_est = J(N, 1, 0)
	    acor_est = J(N, 1, 0)
		
		data1 = data[., 1::floor(S/2)]
		data2 = data[., (floor(S/2)+1)::S]
		data3 = data[.,1::(floor(S/2)+1)]
		data4 = data[.,(floor(S/2)+2)::S]
		
		mean_est1 = J(N, 1, 0)
		mean_est2 = J(N, 1, 0)
		mean_est3 = J(N, 1, 0)
		mean_est4 = J(N, 1, 0)
		
		acov_est1 = J(N, 1, 0)
		acov_est2 = J(N, 1, 0)
		acov_est3 = J(N, 1, 0)
		acov_est4 = J(N, 1, 0)
		
		acor_est1 = J(N, 1, 0)
		acor_est2 = J(N, 1, 0)
		acor_est3 = J(N, 1, 0)
		acor_est4 = J(N, 1, 0)
		
		for(i=1; i<=N; i++){
		    mean_est[i] = mean(data[i,.]')
			mean_est1[i] = mean(data1[i,.]')
			mean_est2[i] = mean(data2[i,.]')
			mean_est3[i] = mean(data3[i,.]')
			mean_est4[i] = mean(data4[i,.]')
			
			acov_est[i] = mataacov(data[i,.], acov_order)
			acov_est1[i] = mataacov(data1[i,.], acov_order)
			acov_est2[i] = mataacov(data2[i,.], acov_order)
			acov_est3[i] = mataacov(data3[i,.], acov_order)
			acov_est4[i] = mataacov(data4[i,.], acov_order)
	
			acor_est[i] = mataacor(data[i,.], acor_order)
			acor_est1[i] = mataacor(data1[i,.], acor_order)
			acor_est2[i] = mataacor(data2[i,.], acor_order)
			acor_est3[i] = mataacor(data3[i,.], acor_order)
			acor_est4[i] = mataacor(data4[i,.], acor_order)
		} 
		
		equantity3 = (mean_est,mean_est1,mean_est2,mean_est3,mean_est4, ///
		              acov_est,acov_est1,acov_est2,acov_est3,acov_est4, ///
				      acor_est,acor_est1,acor_est2,acor_est3,acor_est4)
	
		estimate_value = hpjmomentest2(equantity3, 1::N)
	
		number_par = length(estimate_value)
		estimate_boot = J(B, number_par, 0)
		for (b = 1; b <= B; b++) {
			index_boot = rdiscrete(N, 1, J(N, 1, 1/N))
			estimate_boot[b,.] = hpjmomentest2(equantity3, index_boot)
		}
		se = sqrt(diagonal(variance(estimate_boot)))
		quantile_boot_1 = mm_quantile(estimate_boot, 1, 0.975)
		quantile_boot_2 = mm_quantile(estimate_boot, 1, 0.025)
		ci_1 = 2 * estimate_value - quantile_boot_1
		ci_2 = 2 * estimate_value - quantile_boot_2
		ci = (ci_1 \ ci_2)
	    result =  (estimate_value', se, ci')
	}
	
	st_matrix("est", estimate_value')
	st_matrix("se", se)
	st_matrix("ci", ci')
	
	printf("\n")
    printf("Estimates for Moments.\n")
    printf("Parameters                                              Estimate          \n")
    printf("________________________________________________________________ \n")
    printf("Mean of Mean                                            %f\n",estimate_value[1])
	printf("Mean of Autocovariance                                  %f\n",estimate_value[2])
	printf("Mean of Autocorrelation                                 %f\n",estimate_value[3])
	printf("Variance of Mean                                        %f\n",estimate_value[4])
	printf("Variance of Autocovariance                              %f\n",estimate_value[5])
	printf("Variance of Autocorrelation                             %f\n",estimate_value[6])
	printf("Correlation between Mean and Autocovariance             %f\n",estimate_value[7])
	printf("Correlation between Mean and Autocorrelation            %f\n",estimate_value[8])
	printf("Correlation between Autocovariance and Autocorelation   %f\n",estimate_value[9])
	
	printf("\n")
    printf("%f %% Confidence Intervals for Moments.\n", level*100)
    printf("Parameters                                              Low                   High\n")
    printf("__________________________________________________________________________________\n")
    printf("Mean of Mean                                            %f        %f\n",ci_1[1],ci_2[1])
	printf("Mean of Autocovariance                                  %f        %f\n",ci_1[2],ci_2[2])
	printf("Mean of Autocorrelation                                 %f        %f\n",ci_1[3],ci_2[3])
	printf("Variance of Mean                                        %f        %f\n",ci_1[4],ci_2[4])
	printf("Variance of Autocovariance                              %f        %f\n",ci_1[5],ci_2[5])
	printf("Variance of Autocorrelation                             %f        %f\n",ci_1[6],ci_2[6])
	printf("Correlation between Mean and Autocovariance             %f        %f\n",ci_1[7],ci_2[7])
	printf("Correlation between Mean and Autocorrelation            %f        %f\n",ci_1[8],ci_2[8])
	printf("Correlation between Autocovariance and Autocorelation   %f        %f\n",ci_1[9],ci_2[9])
	
	printf("\n")
    printf("Standard Errors for Moments.\n")
    printf("Parameters                                              Stanadard Errors          \n")
    printf("________________________________________________________________________\n")
    printf("Mean of Mean                                            %f\n",se[1])
	printf("Mean of Autocovariance                                  %f\n",se[2])
	printf("Mean of Autocorrelation                                 %f\n",se[3])
	printf("Variance of Mean                                        %f\n",se[4])
	printf("Variance of Autocovariance                              %f\n",se[5])
	printf("Variance of Autocorrelation                             %f\n",se[6])
	printf("Correlation between Mean and Autocovariance             %f\n",se[7])
	printf("Correlation between Mean and Autocorrelation            %f\n",se[8])
	printf("Correlation between Autocovariance and Autocorelation   %f\n",se[9])
}

// 3.2. Third-Order-Jackknife Moment Estimation

// To be done later.


///  4. Kernel Density Estimation

// 4.1. Naive Estimation of Kernel Density
function kdest (x, X, h) {
    est = sum(normalden((x :- X)/h))/(length(X)*h)
    return(est)
}

function m_nekd (data, acov_order, acor_order){
    N = rows(data)
	S = cols(data)
	grid = 400
	
	mean_est = J(N,1,0)
	acov_est = J(N,1,0)
	acor_est = J(N,1,0)
	for (i=1 ; i<=N ; i++){
	    mean_est[i] = mean(data[i,.]')
		acov_est[i] = mataacov(data[i,.], acov_order) 
		acor_est[i] = mataacor(data[i,.], acor_order)
	}
	
	mean_bw = kdens_bw_dpi(mean_est, level = 2)
	acov_bw = kdens_bw_dpi(acov_est, level = 2)
	acor_bw = kdens_bw_dpi(acor_est, level = 2)
	
	mean_lim = (min(mean_est), max(mean_est))
    acov_lim = (min(acov_est), max(acov_est))
    acor_lim = (min(acor_est), max(acor_est))
    
	mean_grid = rangen(mean_lim[1], mean_lim[2], grid)
	acov_grid = rangen(acov_lim[1], acov_lim[2], grid)
	acor_grid = rangen(acor_lim[1], acor_lim[2], grid)
	
	mean_dest = J(grid, 1, 0)
	acov_dest = J(grid, 1, 0)
	acor_dest = J(grid, 1, 0)
	
	for (i = 1; i <= grid; i++) {
	    mean_dest[i] = kdest(mean_grid[i], mean_est, mean_bw)
		acov_dest[i] = kdest(acov_grid[i], acov_est, acov_bw)
		acor_dest[i] = kdest(acor_grid[i], acor_est, acor_bw)
	}
	
	temp=st_addvar("double", "mean_dest")
    temp=st_addvar("double", "acov_dest")
    temp=st_addvar("double", "acor_dest")
    temp=st_addvar("double", "mean_grid")
	temp=st_addvar("double", "acov_grid")
	temp=st_addvar("double", "acor_grid")
	
    st_addobs(max((0,grid  - st_nobs())))
    st_store(.,"mean_dest", mean_dest\J(st_nobs()-rows(mean_dest),1,.))
    st_store(.,"acov_dest", acov_dest\J(st_nobs()-rows(acov_dest),1,.))
    st_store(.,"acor_dest", acor_dest\J(st_nobs()-rows(acor_dest),1,.))
	
    st_store(.,"mean_grid", mean_grid\J(st_nobs()-rows(mean_grid),1,.))
	st_store(.,"acov_grid", acov_grid\J(st_nobs()-rows(acov_grid),1,.))
    st_store(.,"acor_grid", acor_grid\J(st_nobs()-rows(acor_grid),1,.))
}

// 4.2. Half-Panel-Jackknife 

function hpjkdest1 (x, X, X1, X2, h) {
    N = length(X)
	
	est = sum(normalden((x :- X)/h))/(N * h)
	est1 = sum(normalden((x :- X1)/h))/(N * h)
	est2= sum(normalden((x :- X2)/h))/(N * h)
    
	hpjest = 2 * est - (est1 + est2) / 2
	
	if (hpjest < 0) {
	    hpjest = 0
	}
	return(hpjest)
}

function hpjkdest2 (x, X, X1, X2, X3, X4, h) {
    N = length(X)
	
	est = sum(normalden((x :- X)/h))/(N * h)
	est1 = sum(normalden((x :- X1)/h))/(N * h)
	est2= sum(normalden((x :- X2)/h))/(N * h)
    est3 = sum(normalden((x :- X3)/h))/(N * h)
	est4= sum(normalden((x :- X4)/h))/(N * h)
 
    hpjest = 2 * est - (est1 + est2 + est3 + est4) / 4
	
	if (hpjest < 0) {
	    hpjest = 0
	}
	return(hpjest)
}


function m_hpjkd (data, acov_order, acor_order) {
    N = rows(data)
	S = cols(data)
	grid = 400
	
	mean_est = J(N,1,0)
	acov_est = J(N,1,0)
	acor_est = J(N,1,0)
	for (i=1 ; i<=N ; i++){
	    mean_est[i] = mean(data[i,.]')
		acov_est[i] = mataacov(data[i,.], acov_order) 
		acor_est[i] = mataacor(data[i,.], acor_order)
	}

	mean_bw = kdens_bw_dpi(mean_est, level = 2)
	acov_bw = kdens_bw_dpi(acov_est, level = 2)
	acor_bw = kdens_bw_dpi(acor_est, level = 2)
	
	mean_lim = (min(mean_est), max(mean_est))
    acov_lim = (min(acov_est), max(acov_est))
    acor_lim = (min(acor_est), max(acor_est))
    
	mean_grid = rangen(mean_lim[1], mean_lim[2], grid)
	acov_grid = rangen(acov_lim[1], acov_lim[2], grid)
	acor_grid = rangen(acor_lim[1], acor_lim[2], grid)
	
	mean_dest = J(grid, 1, 0)
	acov_dest = J(grid, 1, 0)
	acor_dest = J(grid, 1, 0)
	
	if (mod(S,2) == 0) {
	    data1 = data[., 1::(S / 2)]
        data2 = data[., (S / 2 + 1)::S]
	    
		mean_est1 = J(N,1,0)
	    acov_est1 = J(N,1,0)
	    acor_est1 = J(N,1,0)
		
		mean_est2 = J(N,1,0)
	    acov_est2 = J(N,1,0)
	    acor_est2 = J(N,1,0)
		
	    for (i=1 ; i<=N ; i++){
	        mean_est1[i] = mean(data1[i,.]')
		    acov_est1[i] = mataacov(data1[i,.], acov_order) 
		    acor_est1[i] = mataacor(data1[i,.], acor_order)
			
			mean_est2[i] = mean(data2[i,.]')
		    acov_est2[i] = mataacov(data2[i,.], acov_order) 
		    acor_est2[i] = mataacor(data2[i,.], acor_order)
	    }
	    for (i = 1; i <= grid; i++) {
	        mean_dest[i] = hpjkdest1(mean_grid[i], mean_est, mean_est1, mean_est2, mean_bw)
		    acov_dest[i] = hpjkdest1(acov_grid[i], acov_est, acov_est1, acov_est2, acov_bw)
		    acor_dest[i] = hpjkdest1(acor_grid[i], acor_est, acor_est1, acor_est2, acor_bw)
	    }    
	}
	else {
	    data1 = data[., 1::floor(S / 2)]
        data2 = data[., (floor(S / 2) + 1)::S]
        data3 = data[., 1::ceil(S / 2)]
        data4 = data[., (ceil(S / 2) + 1)::S]
	    
		mean_est1 = J(N,1,0)
	    acov_est1 = J(N,1,0)
	    acor_est1 = J(N,1,0)
		
		mean_est2 = J(N,1,0)
	    acov_est2 = J(N,1,0)
	    acor_est2 = J(N,1,0)
		
		mean_est3 = J(N,1,0)
	    acov_est3 = J(N,1,0)
	    acor_est3 = J(N,1,0)
		
		mean_est4 = J(N,1,0)
	    acov_est4 = J(N,1,0)
	    acor_est4 = J(N,1,0)
		
		for (i=1 ; i<=N ; i++){
	        mean_est1[i] = mean(data1[i,.]') 
		    acov_est1[i] = mataacov(data1[i,.], acov_order) 
		    acor_est1[i] = mataacor(data1[i,.], acor_order)
			
			mean_est2[i] = mean(data2[i,.]') 
		    acov_est2[i] = mataacov(data2[i,.], acov_order) 
		    acor_est2[i] = mataacor(data2[i,.], acor_order)
			
			mean_est3[i] = mean(data3[i,.]') 
		    acov_est3[i] = mataacov(data3[i,.], acov_order) 
		    acor_est3[i] = mataacor(data3[i,.], acor_order)
			
			mean_est4[i] = mean(data4[i,.]')
		    acov_est4[i] = mataacov(data4[i,.], acov_order) 
		    acor_est4[i] = mataacor(data4[i,.], acor_order)
	    }
		for (i = 1; i <= grid; i++) {
	        mean_dest[i] = hpjkdest2(mean_grid[i], mean_est, mean_est1, mean_est2, mean_est3, mean_est4, mean_bw)
		    acov_dest[i] = hpjkdest2(acov_grid[i], acov_est, acov_est1, acov_est2, acov_est3, acov_est4, acov_bw)
		    acor_dest[i] = hpjkdest2(acor_grid[i], acor_est, acor_est1, acor_est2, acor_est3, acor_est4, acor_bw)
	    }  
	}
    temp=st_addvar("double", "mean_dest")
    temp=st_addvar("double", "acov_dest")
    temp=st_addvar("double", "acor_dest")
    temp=st_addvar("double", "mean_grid")
	temp=st_addvar("double", "acov_grid")
	temp=st_addvar("double", "acor_grid")
	
    st_addobs(max((0,grid  - st_nobs())))
    st_store(.,"mean_dest", mean_dest\J(st_nobs()-rows(mean_dest),1,.))
    st_store(.,"acov_dest", acov_dest\J(st_nobs()-rows(acov_dest),1,.))
    st_store(.,"acor_dest", acor_dest\J(st_nobs()-rows(acor_dest),1,.))
	
    st_store(.,"mean_grid", mean_grid\J(st_nobs()-rows(mean_grid),1,.))
	st_store(.,"acov_grid", acov_grid\J(st_nobs()-rows(acov_grid),1,.))
    st_store(.,"acor_grid", acor_grid\J(st_nobs()-rows(acor_grid),1,.))
}

// 4.3. Third-Order-Jackknife

// To be done later

mata mlib create lpanelhetero, dir(PERSONAL) replace
mata mlib add lpanelhetero *()
mata mlib index
end
