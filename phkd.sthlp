{smcl}

{title:Title}

{p2colset 9 23 20 2}{...}
{p2col :{opt phkd} {hline 2}}Kernel Density Estimation for Heterogeneous Panel Data{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{opt phkd} {it:panelvar} {ifin}[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Order}
{synopt :{opt acov_order(#)}}set order of the autocovariance; default is 0.{p_end}
{synopt :{opt acor_order(#)}}set order of the autocorrelation; default is 1.{p_end}

{syntab:Method}
{synopt :{opth method(string)}}{it:string} must be one of three estimation method {it:"naive", "hpj", "toj"}.{p_end}
{synoptline}

{p 4 6 2}{it:panelvar} must be {help xtset} and strongly balanced.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:phkd} performs kernel density estimation when the panel data exhibits heterogeneity across its cross-sectional units.


{marker options}{...}
{title:Options}

{dlgtab:Order}

{phang}
{opt acov_order} non-negative integer {it:k} for the order of autocovariance. The default is 0. 

{phang}
{opt acor_order} positive integer {it:k} for the order of autocorrelation. The default is 1. 


{dlgtab:Method}

{phang}
{opth method:(strings:string)} specifies how the densities of moments are estimated. 
{it:"naive"} stands for naive estimation without bias-correction, {it:"hpj"} for half panel jackknife and {it:"toj"} for third order jackknife.

{marker results}
{title:Results} 

{pstd}{cmd:phkd} gives three plots of densities for mean, autocovariance and autocorrelation.

{marker example}{...}
{title:Examples:  kernel density estimation}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse pig}{p_end}
{phang2}{cmd:. xtset id week}{p_end}

{pstd}Kernel Estimate the variable {it:weight} about mean, autocovariance of order 2 and autocorrelation of order 3 using Half Panel Jackknife{p_end}
{phang2}{cmd:. phkd weight, method("hpj") acov_order(2) acor_order(3)}{p_end}


{marker references}{...}
{title:References}

{marker OY2019}{...}
{phang}
Ryo Okui. and Takahide Yanagi. 2019.
{browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3128885":{it:Kernel Estimation for Panel Data with Heterogeneous Dynamics}.}
{it:Working Paper}.

{marker DM1993}{...}
{phang}
Ryo Okui. and Takahide Yanagi. 2019.
{browse "https://www.sciencedirect.com/science/article/pii/S0304407619301022?via%3Dihub":{it:Panel Data Analysis with Heterogeneous Dynamics}.}
{it:Journal of Econometrics}.
{p_end}
