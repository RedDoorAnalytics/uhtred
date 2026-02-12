//source paths
local drive /Users/michael.crowther/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
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

merlin (stime trt trt#rcs(stime, df(1) orthog) age M1[id]@1, ///
	family(rp, df(3) failure(died)) timevar(stime))

cap drop t1
range t1 0 5 100
predict s1, surv fitted timevar(t1)

uhtred (stime trt c.trt#rcs(stime, df(1) orthog) age M1[id]@1, ///
	family(rp, df(3) failure(died)) timevar(stime))
predict s2, surv fitted timevar(t1)

assert reldif(s1,s2)<1E-05 if _n<=100
