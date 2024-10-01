/*============================================================================
program: rp_2level_pred
description: testing predictions from 2 level models 
created: by Hannah Bower

tests: 	
	
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
gen age = rnormal(55,5)

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(age 0.01 trtui 1 u1 1) maxt(5) 
stset stime, f(dead)

xtile agecat = age, nq(4)
tab agecat, gen(agecat)


//============================================================================
//test code- mestreg
mestreg trt agecat? || id1:, distribution(weibull) cov(unstr) nohr

//different predictions we can use
//loop over marginal predictions: marginalise over random intercepts
//predictions for fixed part of the model not working- not sure what's going on?
//fixedonly seems to be an option?
foreach p in surv haz dens eta {
	predict `p'_mef, `p' fixed
	predict `p'_mem, `p' marginal 
}

//reffects calculates estimates of the random effects using empirical Bayes predictions
//reses for standard errors
predict reffect_mem, reffects reses(reses_mem)

//testcode for merlin 
merlin (stime trt agecat1 agecat2 agecat3 M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
	
foreach p in surv haz dens eta {
	predict `p'_mf, `p' fixedonly ci
	predict `p'_mm, `p' marginal ci
	predict `p'_mat, `p' at(trt 1) ci
}
predict reffect_m, reffects 
predict reses_m, reses

foreach v in surv haz dens /*eta*/ reffect reses {
	di in red "`v'"
	assert abs(`v'_mem -`v'_mm)< 1E-5
	if "`v'"!= "reffect" | "`v'"!= "reses" {
	di in red "`v' 2 ""	
		assert abs(`v'_mef -`v'_mf)< 1E-5
	}
}




//note eta is different, why? 
//now check with uhtred - assuming same/similar syntax to merlin 
uhtred (stime trt agecat1 agecat2 agecat3 M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
	
//error message for 
foreach p in surv haz dens eta {
	di in red "`p'"
	predict `p'_uf, `p' fixedonly
	predict `p'_um, `p' marginal

}
predict reffect_u, reffects 
predict reses_u, reses

foreach v in surv haz dens eta reffect reses {
	di in red "`v'"
	assert abs(`v'_mf -`v'_uf)< 1E-5
	assert abs(`v'_mm -`v'_um)< 1E-5

}






//============================================================================
//error messages- syntax breaks



//============================================================================//
// end of file
//============================================================================//
