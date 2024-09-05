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
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)

gen t0 = 0
replace t0 = runiform() * 2 //if runiform()<0.9

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		/*tde(trt 0.01) tdefunc(log({t}))*/	///
		maxt(10) ltruncated(t0)

stset stime, f(died) enter(t0)
gen bh = 0.02

timer clear
timer on 1
stpm2 trt /*bmi x1 x2 x3*/ , scale(hazard) df(2)
timer off 1
timer on 2
// merlin (_t 	trt bmi x1 x2 x3 		///
// 		, family(rp, df(2) failure(died) ltruncated(t0))) 	///
//                 , evaltype(gf0)
timer off 2

predict h0, hazard 
predict s0, survival timevar(_t)

mat b = -0.5,-0.7,0.6,0

timer on 3
uhtred (_t 	trt /*bmi x1 x2 x3 */		///
		, family(rp, df(2) failure(died) ltruncated(t0))) 	///
                , evaltype(gf0)  //from(b)
timer off 3
timer list

predict h1, hazard 
predict s1, survival timevar(_t)



// streg trt bmi x1 x2 x3 , dist(weib) nohr
// replace _t0 = 0
// predict s2, surv
//
//
// scatter s1 s2
