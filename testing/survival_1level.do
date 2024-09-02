//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"
clear all

tr:do ./build/buildmlib.do
mata mata clear

pr drop _all

clear 
set seed 725
set obs 50000
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
gen bh = 0.02

timer clear
timer on 1
// stpm2 trt bmi x1 x2 x3 , scale(hazard) df(3)
timer off 1
timer on 2
merlin (_t 	trt bmi x1 x2 x3 		///
		trt#rcs(_t, log orthog event df(2)) ///
		, family(rp, df(3) failure(died))) 	///
                , evaltype(gf2)
timer off 2

predict h0, hazard ci
predict s0, survival ci


timer on 3
uhtred (_t 	trt bmi x1 x2 x3 		///
		c.trt#rcs(_t, log orthog event df(2)) ///
		, family(rp, df(3) failure(died) /*ltruncated(t0)*/)) 	///
                , evaltype(gf2) 
timer off 3
timer list

predict h1, hazard ci
predict s1, survival ci

