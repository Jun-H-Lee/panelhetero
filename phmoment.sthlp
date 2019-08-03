{smcl}

{title:Title}

{p2colset 9 23 20 2}{...}
{p2col :{opt phmoment} {hline 2}}Moments Estimation for Heterogeneous Panel Data{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{opt phmoment} {it:panelvar} {ifin}[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Order}
{synopt :{opt acov_order(#)}}set order of the autocovariance; default is 0.{p_end}
{synopt :{opt acor_order(#)}}set order of the autocorrelation; default is 1.{p_end}
{synopt :{opt boot(#)}}set number of bootstrap replication; default is 200.{p_end}

{syntab:Method}
{synopt :{opth method(string)}}{it:string} must be one of three estimation method {it:"naive", "hpj", "toj"}.{p_end}
{synoptline}

{p 4 6 2}{it:panelvar} must be {help xtset} and strongly balanced.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:phmoment} performs estimation of 9 moments({it:mean of mean, mean of autocovariance, mean of autocorrelation,}
{it:variance of mean, variance of autocovariance, variance of autocorrelation,}
{it:correlation between mean and autocovariance, correlation between mean and autocorrelation and correlation between autocovariance and autocorrelation})
when the panel data exhibits heterogeneity across its cross-sectional units.


{marker options}{...}
{title:Options}

{dlgtab:Order}

{phang}
{opt acov_order} non-negative integer {it:k} for the order of autocovariance. The default is 0. 

{phang}
{opt acor_order} positive integer {it:k} for the order of autocorrelation. The default is 1. 

{phang}
{opt boot} positive interger {it:k} for the number of bootstrap replication. The default is 200.

{dlgtab:Method}

{phang}
{opth method:(strings:string)} specifies how the densities of moments are estimated. 
{it:"naive"} stands for naive estimation without bias-correction, {it:"hpj"} for half panel jackknife and {it:"toj"} for third order jackknife.

{marker results}
{title:Stored results} 

{pstd}
{cmd:phmoment} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(ci)}} 95% confidence intervals for the moments based on cross-sectional bootstrap.{p_end}
{synopt:{cmd:e(se)}} standard errors for the estimators based on cross-sectional bootstrap. {p_end}
{synopt:{cmd:e(est)}} estimates for the moments.{p_end}

{pstd}
All these are ordered by {it:mean of mean, mean of autocovariance, mean of autocorrelation, variance of mean, variance of autocovariance, variance of autocorrelation,}
{it:correlation between mean and autocovariance, correlation between mean and autocorrelation and correlation between autocovariance and autocorrelation.}{p_end}


{marker example}{...}
{title:Examples:  moments estimation}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse pig}{p_end}
{phang2}{cmd:. xtset id week}{p_end}

{pstd}Estimate the moments of the variable {it:weight} about mean, autocovariance of order 2 and autocorrelation of order 3 using Half Panel Jackknife with 300 bootstrap replications.{p_end}
{phang2}{cmd:. phmoment weight, method("hpj") boot(300) acov_order(2) acor_order(3)}{p_end}


{marker reference}{...}
{title:Reference}

{marker DM1993}{...}
{phang}
Ryo Okui. and Takahide Yanagi. 2019.
{browse "https://www.sciencedirect.com/science/article/pii/S0304407619301022?via%3Dihub":{it:Panel Data Analysis with Heterogeneous Dynamics}.}
{it:Journal of Econometrics}.
{p_end}
