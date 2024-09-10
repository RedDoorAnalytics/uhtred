
clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) 

stset stime, f(died) 


uhtred (_t trt bmi x1 x2 x3, family(rp, df(3) failure(died)))


//testing syntax errors
cap drop h
rcof "predict h, hj " == 198
rcof "predict h, h jhf" == 198
rcof "prodict h, hazard" == 199
rcof "predict h, hazard at(trtr 1)" == 111
gen h=.
rcof "predict h, hazard at(trtr 1)" == 110
rcof "predict h, hazard timevar(jsdk)" == 111