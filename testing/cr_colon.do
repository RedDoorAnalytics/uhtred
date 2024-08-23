//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"
clear all

tr:do ./build/buildmlib.do
mata mata clear


use ./data/colon, clear
 
gen female = (sex==2)

// Expand data
// This creates 2 copies of each observation
set seed 978
// keep if runiform()<0.1
expand 2

// Recode and set up data for competing risk analysis
bysort id: gen cause=_n		// cause =1 for cause 1, cause =2 for cause 2

gen cancer=(cause==1)		// indicator for observation for cancer
gen other=(cause==2)		// indicator for observation for other

// Create dummy variables for age
tab agegrp, gen(ag)

// Categorize age and create interactions with cause
forvalues i = 0/3 {
	gen age`i'can=(agegrp==`i' & cancer==1) 
	gen age`i'oth=(agegrp==`i' & other==1) 
}

// Allow different effect of sex for cancer and other */
gen fem_can = female*cancer
gen fem_other = female*other

// Event indicator
gen event=(cause==status)  // status=1 death due to cancer, =2 death due to other

stset surv_mm, failure(event) scale(12) exit(time 120.5)

timer clear
timer on 1
merlin 	(_t 	female ag2 ag3 ag4				///
		female#rcs(_t , df(3) log event orthog)		///
		ag2#rcs(_t , df(3) log event orthog)		///
		ag3#rcs(_t , df(3) log event orthog)		///
		ag4#rcs(_t , df(3) log event orthog)		///
		if cause==1					///
		, family(rp, failure(event) df(4)) timevar(_t))	///
	(_t 	female ag1 ag2 ag3				///
		if cause==2					///
		, family(rp, failure(event) df(4)))		///
		, evaltype(gf1)
timer off 1
timer on 2
uhtred 	(_t 	female ag2 ag3 ag4				///
		c.female#rcs(_t , df(3) log event orthog)	///
		c.ag2#rcs(_t , df(3) log event orthog)		///
		c.ag3#rcs(_t , df(3) log event orthog)		///
		c.ag4#rcs(_t , df(3) log event orthog)		///
		if cause==1					///
		, family(rp, failure(event) df(4)))		///
		
		
	(_t 	female ag1 ag2 ag3				///
		if cause==2					///
		, family(rp, failure(event) df(4)))		///
		, evaltype(gf1)
timer off 2
timer list
