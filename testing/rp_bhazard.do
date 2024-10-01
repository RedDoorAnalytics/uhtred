/*============================================================================
program: rp_bhazard
description: testing for relative survival RP models 
created: by Hannah Bower

tests: 
	 
=============================================================================*/


// simulate data
set seed 72549

clear
set obs 5000

//exp mortailty
local truebhaz = 0.1
gen bhaz = exp(rnormal(log(`truebhaz'),0.5))
gen lbhaz = log(bhaz)

//excess mortality
local lambda = 0.1
local gamma  = 1.2
local beta = -0.5
gen trt = rbinomial(1,0.5)


survsim stime dead , hazard(`truebhaz' :+ `lambda':*`gamma':*{t}:^(`gamma':-1):*exp(`beta':*trt)) maxt(5)
stset stime, f(dead)


//now fit stpm2 model
stpm2 trt, scale(h) df(1) bhaz(bhaz) 

//merlin
stmerlin trt, dist(rp) df(1) bhazard(bhaz)

//uhtred
uhtred (stime 	trt 						///
        , family(rp, df(1) failure(dead) bhazard(bhaz)))

