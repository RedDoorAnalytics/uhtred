//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"
clear all

tr:do ./build/buildmlib.do
mata mata clear
clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)
gen t0 = runiform() * 2
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) //ltruncated(t0)
stset stime, f(died) //enter(t0)
merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, knots(-4 1.5 2.3) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
		, 
		
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, knots(-4 1.5 2.3) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
		, 
