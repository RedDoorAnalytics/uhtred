//source paths
local drive /Users/michael/Library/CloudStorage/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"
clear all

tr:do ./build/buildmlib.do
mata mata clear


set seed 98798
clear
set obs 100
gen id 	= _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.1))
gen u1 	= rnormal(0,sd1)
gen age = rnormal()
expand 100
sort id 
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.1 u1 1) maxt(5) 
stset stime, f(died)

timer clear
timer on 1
mestreg trt age || id:, dist(weib) //evaltype(gf0) intmethod(gh)
predict s2, surv marginal
// predict s5, surv cond(ebmeans)
// predict d1, density marginal
// predict r1, reffects
// merlin (stime trt age M1[id]@1, family(weib, failure(died)))	///
// 	, 
// merlin (stime trt age  M1[id]@1, family(rp, df(1) failure(died)))
// timer off 1
// predict s01, surv fitted
// predict h01, haz fitted

// predict s012, surv marg
// predict h012, haz marg
timer on 2
uhtred (stime trt age M1[id]@1, family(rp, df(1) failure(died))) , ///
	restartvalues(M1 0.01)
timer off 2

timer on 3
// stuhtred trt age || id: , dist(rp) df(1)
timer off 3

// timer list
// predict s1, surv fitted
// predict h1, haz fitted
//
// predict s12, surv marg
// predict h12, haz marg
//
// // scatter s0 s1
