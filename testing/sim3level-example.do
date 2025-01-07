//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software

cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

clear all
tr:do ./build/buildmlib.do
mata mata clear
/*
save the uhtred files to a folder of your choosing, 
and upodate the below local in drive
*/

// local drive /Users/michael/uhtred
// adopath ++ "`drive'"
// pr drop _all
// mata mata mlib index
//
// log using sim3level-comparison, replace text
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
expand 100
gen id = _n
sort hosp surg id

//simulate survival times from a weibull model with random effects at 
//the hospital and surgeon level
// log hazard ratios are specified on trt and age
survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) 	///
	cov(trt -0.5 age 0.02 u1 1 u2 1) 			///
	maxt(10)
	
stset stime1, f(dead1)

//fit the true model with mestreg
timer clear		
timer on 1
// mestreg trt age || hosp: || surg: , dist(weib) nohr
timer off 1

//compare to merlin with royston-parmar and df=1 (equivalent to a weibull)
timer on 2
// merlin 	(stime1 trt age M2[hosp>surg]@1 M1[hosp]@1 , ///
// 	family(rp, df(1) failure(dead1))) , adaptopts(log) trace
// equivalent to:	
// stmixed trt age || id1: || id2: , dist(rp) df(1)

/*
-- Iteration 0:   Adapted log likelihood = -329277.43
-- Iteration 1:   Adapted log likelihood =  -329429.1
-- Iteration 2:   Adapted log likelihood = -329381.89
-- Iteration 3:   Adapted log likelihood = -329336.85
-- Iteration 4:   Adapted log likelihood = -329300.86
-- Iteration 5:   Adapted log likelihood =  -329273.1
-- Iteration 6:   Adapted log likelihood = -329253.92
-- Iteration 7:   Adapted log likelihood = -329246.52
-- Iteration 8:   Adapted log likelihood = -329245.47

*/

timer off 2

//compare to uhtred (merlin2)
timer on 3
uhtred 	(stime1 trt age M2[hosp>surg]@1 M1[hosp]@1 , ///
	family(rp, df(1) failure(dead1))) , ///
	adaptopts(log)
	//trace restartvalues(M1 0.01 M2 0.01)
timer off 3

//compare timings
timer list
log close
