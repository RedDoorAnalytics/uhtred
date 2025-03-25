/*
save the uhtred files to a folder of your choosing, 
and update the below local in drive
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
mestreg trt age || hosp: || surg: , dist(weib) nohr
timer off 1

//compare to uhtred (merlin2)
timer on 2
uhtred 	(stime1 trt age M2[hosp>surg]@1 M1[hosp]@1 , ///
	family(rp, df(1) failure(dead1))) ,
timer off 2

//compare timings
timer list

//more complex baseline
uhtred 	(stime1 trt age M2[hosp>surg]@1 M1[hosp]@1 , ///
	family(rp, df(5) failure(dead1))) ,

predict refs*, reffects //posterior means
predict serefs*, reses //standard errors of posterior means
predict s1, surv fitted //survival conditional on posterior means
