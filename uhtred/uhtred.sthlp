{smcl}
{* *! version 1.0.0. 27aug2024}{...}
{vieweralsosee "uhtred model description options" "help uhtred_models"}{...}
{vieweralsosee "uhtred estimation options" "help uhtred_estimation"}{...}
{vieweralsosee "uhtred reporting options" "help uhtred_reporting"}{...}
{vieweralsosee "uhtred postestimation" "help uhtred_postestimation"}{...}
{viewerjumpto "Syntax" "uhtred##syntax"}{...}
{viewerjumpto "Description" "uhtred##description"}{...}
{viewerjumpto "Options" "uhtred##options"}{...}
{viewerjumpto "Examples" "uhtred##examples"}{...}
{viewerjumpto "Stored results" "uhtred##results"}{...}
{title:Title}

{p2colset 5 15 19 2}{...}
{p2col:{bf:uhtred} {hline 2}}Flexible multivariate mixed-effects survival analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:uhtred} {help uhtred_models:{it:models}} {ifin}
[{cmd:,} {it:options}]

{pstd}
where
{it:models} are the model specifications; see {helpb uhtred_models:uhtred models}.{p_end}

{synoptset 30}{...}
{synopthdr:options}
{synoptline}
{synopt :{help uhtred_model_options:{it:model_description_options}}}fully
define, along with {it:models}, the model to be fit{p_end}

{synopt :{help uhtred_estimation:{it:estimation_options}}}method
used to obtain estimation results, including specifying initial values{p_end}

{synopt :{help uhtred_reporting:{it:reporting_options}}}reporting
of estimation results{p_end}
{synoptline}
{p 4 6 2}
Also see {helpb uhtred_postestimation:uhtred postestimation} for features
available after estimation.
{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:uhtred} fits univariate and multivariate survival models, each of which 
could be repeatedly measured (recurrent events), with any number of levels, and 
with any number of random effects at each level. Interval censoring is suported 
for univariate models, and left truncation is supported for single level models. 

{pstd}
{cmd:uhtred} provides a flexible predictor syntax, allowing the user to define 
variables, random effects, splines, and interaction between each of them. 
Non-linear and time-dependent effects are seamlessly incorporated into the 
predictor.


{marker options}{...}
{title:Options}

{phang}
{it:model_description_options}
describe the model to be fit.  The model to be fit is fully specified by
{it:models} -- which appear immediately after {cmd:uhtred} -- and the option 
{opt covariance()}.  See {helpb uhtred_model_options:uhtred model description options} and 
{helpb uhtred_models:uhtred model notation}.

{phang}
{it:estimation_options}
control how the estimation results are obtained.  These options control how
the standard errors (VCE) are obtained and control technical issues
such as choice of estimation method.  See 
{helpb uhtred_estimation:uhtred estimation options}.

{phang}
{it:reporting_options}
control how the results of estimation are displayed.  See 
{helpb uhtred_reporting:uhtred reporting options}.


{marker examples}{...}
{title:Examples}

{phang}
These examples are intended for quick reference.  For detailed examples, see the 
{bf:{browse "https://reddooranalytics.se/software/uhtred":uhtred homepage}}.

{phang}
{ul:{bf:Example 1: Royston-Parmar model}}
{p_end}

{phang2}Setup{p_end}
{phang3}{cmd:. webuse brcancer}{p_end}

{phang2}Flexible parametric model{p_end}
{phang3}{cmd:. uhtred (rectime hormon, family(rp, df(3) failure(censrec)))}{p_end}

{phang}
{ul:{bf:Example 2: Royston-Paramr model with a random intercept}}
{p_end}

{phang2}Setup{p_end}
{phang3}{cmd:. use http://fmwww.bc.edu/repec/bocode/s/stmixed_example1}{p_end}

{phang2}Flexible parametric model with random intercept{p_end}
{phang3}{cmd:. uhtred (rectime stime M1[centre]@1, family(rp, df(3) failure(event)))}{p_end}


{title:Author}

{p 5 12 2}
{bf:Michael J. Crowther}{p_end}
{p 5 12 2}
Red Door Analytics AB{p_end}
{p 5 12 2}
Stockholm, Sweden{p_end}
{p 5 12 2}
michael@reddooranalytics.se{p_end}


{title:References}

{phang}
{bf:Crowther MJ}. Extended multivariate generalised linear and non-linear mixed effects models. 
{browse "https://arxiv.org/abs/1710.02223":https://arxiv.org/abs/1710.02223}
{p_end}

{phang}
{bf:Crowther MJ}. merlin - a unified framework for data analysis and methods development in Stata. {browse "https://journals.sagepub.com/doi/pdf/10.1177/1536867X20976311":{it:Stata Journal} 2020;20(4):763-784}.
{p_end}

{phang}
{bf:Crowther MJ}. Multilevel mixed effects parametric survival analysis: Estimation, simulation and application. {browse "https://journals.sagepub.com/doi/abs/10.1177/1536867X19893639?journalCode=stja":{it:Stata Journal} 2019;19(4):931-949}.
{p_end}
