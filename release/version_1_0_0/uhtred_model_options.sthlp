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

{p2colset 5 17 19 2}{...}
{p2col:{helpb uhtred} {hline 2}}Flexible multivariate mixed-effects survival analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:uhtred} {it:{help uhtred_models:models}} ...{cmd:,} ...
    {it:model_description_options}

{synoptset 28 tabbed}{...}
{synopthdr:model_description_options}
{synoptline}
{synopt:{opt family()}, ...}see {helpb uhtred_models##family:uhtred families}{p_end}
{synopt:{opt cov:ariance()}}notation for treatment of covariances{p_end}
{synopt :{opt const:raints()}}specify constraints{p_end}
{synopt :{opt from()}}specify starting values{p_end}
{synopt :{opt nogen}}do not generate the component and element variables; see details{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{it:{help uhtred_models:models}} and the options above describe the model to be fit by {cmd:uhtred}.


{marker options}{...}
{title:Options}

{phang}
{cmd:family()} specifies the distribution, such as {cmd:family(poisson)}. See
      {helpb uhtred_models##family:uhtred family}.
	  
{phang}
{opt covariance(struct_list)} specifies the covariance structure for the random effects at each level. The structures of the 
random effects should be defined in order from highest level to lowest. Available structures include

{phang2}{opt diag:onal} unique variances for each random effect, all covariances are 0; the default

{phang2}{opt iden:tity} all variances equal, all covariances are 0

{phang2}{opt un:structured} all variances and covariances uniquely estimated

{phang2}{opt ex:changeable} all variances equal, all covariances equal

{phang}
{opt constraints()} specifies parameter constraints you wish to impose on your
model.

{phang}
{opt from()} specifies the starting values to be used in the optimization
process.

{phang}
{opt nogen} stops {cmd:uhtred} from generating Stata variables, representing those created in the component and element 
specifications. By default, they are created and indexed by model, component and number.

{marker examples}{...}
{title:Examples}

{phang}Setup{p_end}
{phang2}{cmd:. webuse brcancer}{p_end}

{phang}Flexible parametric model{p_end}
{phang2}{cmd:. uhtred (rectime hormon, family(rp, df(3) failure(censrec)))}{p_end}


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
{bf:Crowther MJ}. uhtred - a unified framework for data analysis and methods development in Stata. {browse "https://journals.sagepub.com/doi/pdf/10.1177/1536867X20976311":{it:Stata Journal} 2020;20(4):763-784}.
{p_end}

{phang}
{bf:Crowther MJ}. Multilevel mixed effects parametric survival analysis: Estimation, simulation and application. {browse "https://journals.sagepub.com/doi/abs/10.1177/1536867X19893639?journalCode=stja":{it:Stata Journal} 2019;19(4):931-949}.
{p_end}
