
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

//============================================================================//
//error checks

//new var doesn't exist
rcof "uhtred (stime newvar trt bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 111
rcof "uhtred (stime  trt i.bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog) , family(rp, df(3) failure(died)))" == 452

//tde misspecified
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#ecs(stime, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 198
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(time, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 111
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(stime, df(1) lrg orthog) 	, family(rp, knots(3) failure(died)))" == 3598
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(stime, df(1) log prthog) 	, family(rp, knots(3) failure(died)))" == 3598
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(stime, df(-1) log orthog) 	, family(rp, knots(3) failure(died)))" ==  1986
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(stime, df(1 5) log orthog) 	, family(rp, knots(3) failure(died)))" == 1986
rcof "uhtred (stime trt bmi x1 x2 x3 trt#rcs(stime, df(1) knots(2 4) log orthog) 	, family(rp, knots(3) failure(died)))" == 1986
rcof "uhtred (stime trt bmi x1 x2 x3 trt#rcs(stime,  knots(2 2) log orthog) 	, family(rp, knots(3) failure(died)))" == 3598

//family options misspecifed
rcof "uhtred (stime trt bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog) 	, family(fam, knots(3) failure(died)))" == 198
rcof "uhtred (stime  trt bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog) 	, family(rp, knots(3) failure(newvar)))" == 111
rcof "uhtred (stime  trt bmi x1 x2 x3 	, family(rp, knots(3)))" == 198
*rcof "uhtred (stime trt bmi x1 x2 x3 , family(rp,    failure(died))" == 198

//other incorrect syntax
rcof "uhtred (stime trt bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog 	, family(rp, knots(3) failure(died)))" == 198
rcof "uhtred (stime trt bmi x1 x2 x3 trt#rcs(stime, df(1) log orthog)	 family(rp, knots(3) failure(died)))" == 198
rcof "uhtred stime trt bmi x1 x2 x3, family(rp, df(3) failure(died)))" == 198
rcof "uhtred stime trt bmi x1 x2 x3, family(rp, df(3 failure(died)))" == 198


//left truncated syntax errors
rcof "uhtred (stime trt bmi age x1, family(rp, df(2) failure(died) ltruncated(t02)))" == 111
rcof "uhtred (stime trt bmi age x1, family(rp, df(2) failure(died) letruncated(t0)))" ==3598
rcof "uhtred (stime trt bmi age x1, family(rp, df(2) failure(died)) ltruncated(t0))" == 198
