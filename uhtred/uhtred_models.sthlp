{smcl}
{* *! version 1.0.0. 27aug2024}{...}
{vieweralsosee "uhtred model description options" "help uhtred_models"}{...}
{vieweralsosee "uhtred estimation options" "help uhtred_estimation"}{...}
{vieweralsosee "uhtred reporting options" "help uhtred_reporting"}{...}
{vieweralsosee "uhtred postestimation" "help uhtred_postestimation"}{...}
{title:Title}

{p2colset 5 15 17 2}{...}
{p2col:{helpb uhtred} {hline 2}}Command syntax for model specification{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:uhtred} ({it:model1}) [({it:model2})] [...] [, {bf:covariance()}]

{pstd}
where the syntax of a {it:model} is

{p 8 12 2}
[{depvar}] [{it:component1}] [{it:component2}] [...] {ifin} [, {it:{help uhtred_models##model_options:model_options}}]

{pmore}
where the syntax of a {it:component} is

{pmore2}
{it:element1}[{bf:#}[{bf:#}]{it:element2}][{bf:#}[{bf:#}]{it:element3}][...][@{it:real}]

{pmore2}
and each {it:elementn} can take one of the forms described in {it:{help uhtred_models##elements:element}}, 
with {it:@{it:real}} described in {it:{help uhtred_models##elements_details:elements details}}.


{synoptset 35}{...}
{marker model_options}{...}
{synopthdr:model_options}
{synoptline}
{synopt :{cmd: {ul:f}amily({it:{help uhtred_models##family:family}})}}distributional family{p_end}
{synopt :{opt nocons:tant}}omit the constant term{p_end}
{synoptline}

{synoptset 35}{...}
{marker family}{...}
{synopthdr:family}
{synoptline}
{synopt :{cmd: rp [, {it:{help uhtred_models##survival:survival}} {it:{help uhtred_models##rpopts:rpopts}}]}}Royston-Parmar model on the log cumulative hazard scale{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 35}{...}
{marker survival}{...}
{synopthdr:survival}
{synoptline}
{synopt :{opth fail:ure(varname)}}indicator for failure event{p_end}
{synopt :{opth li:nterval(varname)}}lower interval time for interval censored observations{p_end}
{synopt :{opth lt:runcated(varname)}}entry time for left-truncated/delayed-entry model{p_end}
{synopt :{opth bh:azard(varname)}}expected mortality rate at event times, invokes a relative survival model{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 35}{...}
{marker rpopts}{...}
{synopthdr:rpopts - family(rp)}
{synoptline}
{synopt :{opth df(#)}}degrees of freedom for the baseline log cumulative hazard function (does not include the intercept){p_end}
{synopt :{opt knots(knots_list)}}knot locations for the baseline log cumulative hazard function - includes 
boundary knots, should be in increasing order.{p_end}
{synopt :{opt noorth:og}}turns off the default orthogonalisation of the spline terms{p_end}
{synoptline}
{p2colreset}{...}


{synoptset 35}{...}
{marker elements}{...}
{synopthdr:element}
{synoptline}
{synopt :{opt {varname}}}a variable in the dataset{p_end}
{synopt :{bf:M}#{cmd:[}{it:levelvar...}{cmd:]}}a random effect at the specified level; see {it:{help uhtred_models##elements_details:{it:details}}}{p_end}
{synopt :{cmd: rcs(}{it:varname}, {help uhtred_models##rcs_opts:{it:rcs_opts}}{cmd:)}}a restricted cubic spline function of {it:varname}{p_end}
{synoptline}

{synoptset 35}{...}
{marker rcs_opts}{...}
{synopthdr:rcs_opts}
{synoptline}
{synopt :{opth df(#)}}degrees of freedom for the spline function; see details{p_end}
{synopt :{opt knots(knots_list)}}knot locations for the spline function, including boundary knots{p_end}
{synopt :{opt orthog}}use Gram-Schmidt orthogonalisation on the spline variables{p_end}
{synopt :{opt event}}when using {cmd:df()}, calculate internal knot locations based on centiles of the observations that had an event (i.e. for survival models){p_end}
{synopt :{opt log}}calculate splines of the log of {it:varname}, rather than {it:varname}{p_end}
{synopt :{opth off:set(varname)}}to add before the spline function is calculated{p_end}
{synopt :{opth moff:set(varname)}}adds the negative of {it:varname} before the spline function is calculated{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
The command syntax of {cmd:uhtred} is fully specified by the {it:models} and the {bf:covariance()} option. 

{pstd}
For full details and many tutorials, take a look at the accompanying website: 
{browse "https://reddooranalytics.se/software/uhtred":reddooranalytics.se/software/uhtred}


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}{cmd: family({it:{help uhtred_models##family:family}})} distributional family{p_end}

{phang}{opt noconstant} suppresses the constant (intercept) term in the extended linear predictor.


{marker elements_details}{...}
{dlgtab:Elements}

{phang}Elements are the fundamental flexibility of {helpb uhtred}. Each {it:component} can have any number of {it:elements}. 
Elements within the same {it:component} are split using a {cmd:#} or {cmd:##}, as usual. Each {it:element} can be one of the following:

{phang2}{opt {varname}} a variable in the dataset. Any of Stata's {cmd:.} notation is not currently supported.

{phang2}{bf:M}#{cmd:[}{it:levelvar...}{cmd:]} a random effect at the specified level.

{p 12 12 2}{bf:M}# defines a random effect, which must begin with {cmd:M}, and be followed by a positive integer

{p 12 12 2}{it:levelvar} defines the cluster variables, for example, {cmd:M1[level1]} defines a random effect called {cmd:M1} 
at {cmd:level1}, and {cmd:M2[level1>level2]} defines a random effect called {cmd:M2} at {cmd:level2} which is nested within {cmd:level1}.
{* {phang2}{cmd: {f(&t)}} a user-defined function of time {bf:&t}. The function must be written in Mata code using colon {cmd::} notation }
{* (for element by element operations), and enclosed within braces {cmd:{}}. {cmd:timevar()} is also required.}

{phang2}{cmd:rcs(}{it:varname}, {help uhtred_models##rcsopts:{it:rcs_opts}}{cmd:)} a restricted cubic spline function of {it:varname}, where {it:rcs_opts} are:

{p 12 16 2}{cmd:df(#)} specifies the degrees of freedom for the spline function. Boundary knots are assumed to be the min and max of 
the {it:varname}. Internal knots are placed at equally spaced centiles.

{p 12 16 2}{cmd:knots(numlist)} specifies the knot locations, which includes the boundary knots. Must be in ascending order.

{p 12 16 2}{cmd:log} use splines of log {it:varname} rather than {it:varname}, the default.

{p 12 16 2}{cmd:orthog} orthogonalise the spline terms.

{p 12 16 2}{cmd:event} can be used in combination with {cmd:df()} and specifies that the internal knots are calculated 
based on centiles of event times. Only valid when the {cmd:family()} is a survival model.

{phang3}{opt offset(varname)} defines an offset to be added to the {cmd:rcs()} variable prior to the spline function being derived.

{phang3}{opt moffset(varname)} defines a negative offset to be taken away from the {cmd:rcs()} variable prior to the spline function being derived.

{phang2}{bf:@}{it:real} specifies a constraint to be applied to all elements of the component. This is generally used in combination with random 
effects, i.e to specify a random intercept, one would type {cmd:M1[id]@1}, which would constrain the coefficient of the random effect 
{cmd:M1} to be 1. Without the {cmd:@1}, {cmd:uhtred} would attempt to estimate a coefficient.


{dlgtab:Familys}

{phang}{opt family(family)} specifies the distributional family. The inbuilt distributions listed above are self explanatory. 


{dlgtab:Survival}

{phang}{opth failure(varname)} specifies the censoring/failure event. Should be coded 0 for right-censored observations, 
1 for exactly observed events, or 2 for interval-censored observations.

{phang}{opth linterval(varname)} specifies the lower interval for interval censored observations. The upper interval 
should be specified in the response variable for the model. Observations that are right censored or events should 
be coded as missing.

{phang}{opth ltruncated(varname)} specifies the time at which observations become at risk of the event. Allows 
fitting a delayed-entry survival model. If there are random effects in the associated survival model, then 
the likelihood is calculated by dividing through by the marginal survival function at the entry times, which results 
in a second set of numerical integration.

{phang}{opth bhazard(varname)} invokes a relative survival (excess hazard) model, by specifying the expected mortality (event) rate in the reference population at the observed event times.


{dlgtab:Royston-Parmar options}

{phang}{opt df(#)} degrees of freedom for the baseline log cumulative hazard function, i.e. number of restricted cubic 
spline terms. Internal knots are placed at centiles of the event times. Boundary knots are placed at the minimum and maximum event times.

{phang}{opt knots(knots_list)} defines the knot locations for the spline functions used to model the baseline 
log cumulative hazard function. Must include boundary knots. Knots should be specified in increasing order.


{title:Examples}

{phang}
For detailed examples, see {bf:{browse "https://reddooranalytics.se/software/uhtred":reddooranalytics.se/software/uhtred}}.


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
