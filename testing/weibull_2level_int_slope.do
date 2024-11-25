/*============================================================================
program: weibull_2level_int_slope
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


mestreg trt age || id1:trt, distribution(weibull) cov(unstr) nohr intpoints(10)
est store mestreg


merlin (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) intpoints(10)
est store merlin	

/*
uhtred (stime trt age M1[id1]@1 c.trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) intpoints(10)
est store uhtred

est table mestreg merlin uhtred

*/


//save estimates in locals for testing
est restore mestreg
foreach v in trt age {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
foreach v in trt _cons {
	local mestreg_sd`v' = `=sqrt(_b[/:var(`v'[id1])])'
}
local mestreg_corr = `= _b[/:cov(trt[id1],_cons[id1])] /(sqrt(_b[/:var(trt[id1])] * _b[/:var(_cons[id1])] ))'

estimate restore merlin
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local j 1
foreach v in _cons trt {
	local mer_sd`v' = `=exp(_b[lns1_`j':_cons])'
	local j `=`j'+1'
}
local merlin_corr = `= tanh(_b[art1_1_2:_cons])'

/*
est restore uhtred
*/

//check estimates 
set tracedepth 3
set trace on 
foreach v in trt age {
	/*
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-3
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-3
	*/
	
	assert abs(`mestreg_b_`v''- `mer_b_`v'')< 1E-3
	assert abs(`mestreg_se_`v''- `mer_se_`v'')< 1E-3		
}
	
foreach v in _cons trt {
	/*
	assert abs(`mestreg_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
	assert abs(`merlin_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
	*/
	
	assert abs(`mestreg_sd`v''- `mer_sd`v'')< 1E-2
}
/*
assert abs(`mestreg_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
assert abs(`merlin_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
*/
	
assert abs(`mestreg_corr'- `merlin_corr')< 1E-2
	

*est drop mestreg merlin uhtred






//============================================================================
//check covariance structures

//loop over all covariances- note diagonal not available for mestreg  
//note2 exchange doesn't seem to be working for merlin or uhtred right now
foreach cov in diag iden un /*ex*/ {
	if "`cov'" != "diag" {
		mestreg trt age || id1:trt, distribution(weibull) cov(`cov') nohr intpoints(10)
		foreach v in trt age {
			local mestreg_b_`v' =_b[_t:`v']
			local mestreg_se_`v' =_se[_t:`v']
		}
		foreach v in trt _cons {
			local mestreg_sd`v' = `=sqrt(_b[/:var(`v'[id1])])'
		}
		local mestreg_corr = `= _b[/:cov(trt[id1],_cons[id1])] /(sqrt(_b[/:var(trt[id1])] * _b[/:var(_cons[id1])] ))'
	}
	
	//merlin 
	merlin (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') intpoints(10)
	local j 1
	foreach v in trt age {
		local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
		local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
		local j `=`j'+1'
	}
	local j 1
	foreach v in _cons trt {
		local mer_sd`v' = `=exp(_b[lns1_`j':_cons])'
		local j `=`j'+1'
	}
	local merlin_corr = `= tanh(_b[art1_1_2:_cons])'

	//uhtred
	uhtred (stime trt age M1[id1]@1 c.trt#M2[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') intpoints(10)
		
	//check estimates 
	foreach v in trt age {
		if "`cov'" != "diag" {
			assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-3
			assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-3
		}		
		assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-3
		assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-3
	}
	foreach v in _cons trt {
		if "`cov'" != "diag" {
			assert abs(`mestreg_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
		}
	assert abs(`merlin_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
	}
	
	assert abs(`mestreg_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
	assert abs(`merlin_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
	
	
}

	
	
//============================================================================
//check inmethods- ghermite


mestreg trt age || id1:trt, distribution(weibull) cov(unstr) nohr ///
	 intmethod(ghermite)
est store mestreg

merlin (stime trt age M1[id1]@1 trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr)  intmethod(ghermite)
est store merlin	


uhtred (stime trt age M1[id1]@1 c.trt#M2[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr)  intmethod(ghermite)
est store uhtred

est table mestreg merlin uhtred




//save estimates in locals for testing
est restore mestreg
foreach v in trt age {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
foreach v in trt _cons {
	local mestreg_sd`v' = `=sqrt(_b[/:var(`v'[id1])])'
}
local mestreg_corr = `= _b[/:cov(trt[id1],_cons[id1])] /(sqrt(_b[/:var(trt[id1])] * _b[/:var(_cons[id1])] ))'

estimate restore merlin
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local j 1
foreach v in _cons trt {
	local mer_sd`v' = `=exp(_b[lns1_`j':_cons])'
	local j `=`j'+1'
}
local merlin_corr = `= tanh(_b[art1_1_2:_cons])'


est restore uhtred


//check estimates 
set tracedepth 3
set trace on 
foreach v in trt age {
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-3
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-3
}
	
foreach v in _cons trt {
	assert abs(`mestreg_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
	assert abs(`merlin_sd`v''- `=exp(_b[lns1_1:`v'])')< 1E-3
}
assert abs(`mestreg_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
assert abs(`merlin_corr'- `= tanh(_b[art1_1_2:_cons])')< 1E-3
	

est drop mestreg merlin uhtred

	
	




//============================================================================//
// end of file
//============================================================================//

