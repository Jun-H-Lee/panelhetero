/// Exercise Code for panelhetro.ado


/// 0. Initialization

clear 

webuse pig
xtset id week

/// 1. Empiricl CDF Estimaton

phecdf weight, method("sim") acov_order(0) acor_order(1)

/// 2. Moment Estimation

phmoment weight, method("hpj") boot(200) acov_order(0) acor_order(1)

ereturn list
matrix list e(ci)
matrix list e(se)
matrix list e(est)

/// 3. Kernel Density Estimation

phkd weight, method("hpj") acov_order(0) acor_order(1)
