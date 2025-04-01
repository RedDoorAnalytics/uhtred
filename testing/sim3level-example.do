//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

//build mlib
clear all
tr:do ./build/buildmlib.do
mata mata clear
set seed 7254
clear
 
//simulate hospitals
set obs 50
gen hosp = _n
gen age = runiform()
gen sd1 = exp(log(0.1))
gen u1 = rnormal(0,sd1)
//simulate surgeons in each hospital
expand 30
bys hosp: gen surg = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(0.1))
gen u2 = rnormal(0,sd2)
//simulate patients within each surgeon & hospital
expand 10
gen id = _n
sort hosp surg id
 
//simulate survival times from a weibull model with random effects at
//the hospital and surgeon level
// log hazard ratios are specified on trt and age
survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(0.8)    ///
     cov(trt -0.5 age 0.02 u1 1 u2 1)                 ///
     maxt(10)
    
stset stime1, f(dead1)
 
//fit the true model with mestreg
mestreg trt age || hosp: || surg: , dist(weib) nohr
predict refs*, reffects
su refs*
predict s1, surv cond(ebmeans)
 
//compare to uhtred (merlin2)
uhtred     (stime1 trt age M2[hosp>surg]@1 M1[hosp]@1 , ///
     family(rp, df(1) noorthog failure(dead1))) , intpoints(15)
predict urefs*, reffects
su urefs*

predict s2, surv fitted
