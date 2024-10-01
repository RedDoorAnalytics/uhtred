/*============================================================================
program: rp_2level_int_slope
description: testing mixed models with random intercept
created: by Hannah Bower

tests: 	compare random intercept models to mestreg and merlin

		check interactions, index notation, df >1, user-defined knots, 
		covariance and intmethod() options
		
		Add error checks for invalid syntax
	
=============================================================================*/

//simulate data using survsim 
set seed 725488
clear
set obs 1000
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,0.25\0.25,1)
drawnorm u1 u2, means(0 0) sds(1 0.5) corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt
gen age = rnormal() + u2

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(age 0.01 trtui 1 u1 1) maxt(5) 
stset stime, f(dead)


//============================================================================
//random intercept and random effect 

//lets try and replicate mestreg with merlin & uhtred
//first random intercepts
mestreg trt age || id1:trt, distribution(weibull) cov(unstr) nohr
est store mestreg


merlin (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 
est store merlin	

//getting some mata error from uhtred, let's just ignore for now and write out the 
// tests
uhtred (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 
est store uhtred

est table mestreg merlin uhtred
est drop mestreg merlin uhtred


//can check the fixed effect covariates and the variances of the random effects/intercept
// corr=cov/sd1sd2
/*test*/
mestreg trt age || id1:trt, distribution(weibull) cov(unstr) nohr

merlin (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 

uhtred (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 

	
uhtred (stime trt rcs(age, df(1))#M1[id1]@1 M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	intmethod(gh) 	///
	cov(unstr) 
	
	
//note the random intercept is constrained at 1 (@1 notation)- ask MC why we do this?
//can add constraints to the mestreg too using constraint(1) and like the following, but not sure what to constrain
* constraint 1 [consump]wagepriv = [consump]wagegovt




//lets try and replicate mestreg with merlin & uhtred
//now random effect of trt
mestreg trt age || id1:trt, distribution(weibull) cov(unstr)
est store mestreg


merlin (stime trt age M1[id1] trt#M2[id1], ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 
	
est store merlin	

merlin (stime trt age#M1[id1]@1 M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 


uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 
est store uhtred

est table mestreg merlin uhtred
