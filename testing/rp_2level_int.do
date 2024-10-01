/*============================================================================
program: rp_2level_int
description: testing mixed models with random intercept
created: by Hannah Bower

tests: 	compare random intercept models to mestreg and merlin

		check interactions, factor variables, df >1, user-defined knots, 
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
//random intercepts only

//mestreg model
mestreg trt age || id1:, distribution(weibull) cov(unstr)

est store mestreg
foreach v in trt age {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
local mestreg_b_cons = _b[/:var(_cons[id1])]

//merlin 
merlin (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
	
est store merlin	
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_cons = `=exp(_b[lns1_1:_cons])'


//uhtred
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
	
est store uhtred

//check by eye
est table mestreg merlin uhtred, se
est drop mestreg merlin uhtred

//check coefficient estimates 

//check estimates 
foreach v in trt age {
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-5
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}

//note variance of random intercept is a bit diff from mestreg
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1
assert abs(`mer_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-5
	
//make assert for cscript	
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 

mkassert eclass


//============================================================================
//check covariance structures

//check different covariance structures - independent and exchangeable not working 
//for uhtred or merlin?
foreach cov in /*independent exchangeable*/  identity unstr {

	mestreg trt age || id1:, distribution(weibull) cov(`cov') nohr

	uhtred (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') noorthog
	
	merlin (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') noorthog
		
}

//loop over all covariances- note diagonal not available for mestreg & exchange 
//doesn't seem to be working for merlin or uhtred right now
foreach cov in diag iden un /*ex*/ {
	if "`cov'" != "diag" {
		mestreg trt age || id1:, distribution(weibull) cov(`cov')
		foreach v in trt age {
			local mestreg_b_`v' =_b[_t:`v']
			local mestreg_se_`v' =_se[_t:`v']
		}
		local mestreg_b_cons = _b[/:var(_cons[id1])]
	}
	//merlin 
	merlin (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') noorthog
	local j 1
	foreach v in trt age {
		local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
		local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
		local j `=`j'+1'
	}
	local mer_b_cons = `=exp(_b[lns1_1:_cons])'

	//uhtred
	uhtred (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(`cov') noorthog
		
	//check estimates 
	foreach v in trt age {
		if "`cov'" != "diag" {
			assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-5
			assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-5
		}		
		assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
		assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
	}
	if "`cov'" != "diag" {
		assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1
	}
	assert abs(`mer_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-5
}


//make asserts 

uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(diag) 

mkassert eclass	

uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(iden) 
	
mkassert eclass	

	
/*
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(ex) 
	
mkassert eclass	
*/

//============================================================================
//check intmethods- integration of latent variables

//intmethods
foreach int in mvaghermite ghermite {
	mestreg trt age || id1:, distribution(weibull) ///
		cov(unstr) nohr ///
		intmethod(`int') 
	
	uhtred (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(unstr) ///
		intmethod(`int') noorthog
	
	merlin (stime trt age M1[id1]@1, ///
		family(rp, df(1) failure(dead))), 	///
		cov(unstr) ///
		intmethod(`int') noorthog
		
}

//can just check ghermite since the default is mvaghermite
mestreg trt age || id1:, distribution(weibull) intmethod(ghermite)
foreach v in trt age {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
local mestreg_b_cons = _b[/:var(_cons[id1])]

//merlin 
merlin (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	intmethod(ghermite) noorthog
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_cons = `=exp(_b[lns1_1:_cons])'

//uhtred
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	intmethod(ghermite) noorthog
	
//check estimates 
foreach v in trt age {
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-5
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1
assert abs(`mer_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-5

//mkassert
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	intmethod(ghermite)

mkassert eclass


//============================================================================
//let's test some other stuff in the random intercept models- factor variables

//categorical variables
xtile agecat = age, nq(4)

mestreg trt i.agecat || id1:, distribution(weibull) cov(unstr) nohr
est store mestreg
foreach v in trt {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
forvalues i=2/4 {
	local mestreg_b_age`i' =_b[_t:`i'.agecat]
	local mestreg_se_age`i' =_se[_t:`i'.agecat]
}
local mestreg_b_cons = _b[/:var(_cons[id1])]


uhtred (stime trt i.agecat M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
est store uhtred

est table mestreg uhtred
est drop mestreg uhtred
	
//check estimates are similar 
assert abs(`mestreg_b_trt'- `=_b[xb1:trt]')< 1E-5
assert abs(`mestreg_se_trt'- `=_se[xb1:trt]')< 1E-5
forvalues i=2/4 {
	assert abs(`mestreg_b_age`i''- `=_b[xb1:`i'.agecat]')< 1E-5
	assert abs(`mestreg_se_age`i''- `=_se[xb1:`i'.agecat]')< 1E-5
}

//note variance of random intercept is a bit diff from mestreg
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1

	
	
	
//============================================================================
//let's test some other stuff in the random intercept models- interactions

//cts-cts interaction
mestreg c.trt##c.age || id1:, distribution(weibull) cov(unstr) nohr
foreach v in trt age {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
local mestreg_b_int =_b[_t:c.trt#c.age]
local mestreg_se_int =_se[_t:c.trt#c.age]
local mestreg_b_cons = _b[/:var(_cons[id1])]


uhtred (stime c.trt##c.age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog

foreach v in trt age {
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-5	
}
assert abs(`mestreg_b_int'- `=_b[xb1:c.trt#c.age]')< 1E-5
assert abs(`mestreg_se_int'- `=_se[xb1:c.trt#c.age]')< 1E-5	
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1
	
	
	

//cat interaction
mestreg c.trt##i.agecat || id1:, distribution(weibull) cov(unstr) nohr
foreach v in trt {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
forvalues i=2/4 {
	local mestreg_b_age`i' =_b[_t:`i'.agecat]
	local mestreg_se_age`i' =_se[_t:`i'.agecat]
	local mestreg_b_int`i' =_b[_t:`i'.agecat#c.trt]
	local mestreg_se_int`i' =_se[_t:`i'.agecat#c.trt]
}
local mestreg_b_cons = _b[/:var(_cons[id1])]


uhtred (stime c.trt##i.agecat M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) noorthog
assert abs(`mestreg_b_trt'- `=_b[xb1:trt]')< 1E-5
assert abs(`mestreg_se_trt'- `=_se[xb1:trt]')< 1E-5
forvalues i=2/4 {
	assert abs(`mestreg_b_age`i''- `=_b[xb1:`i'.agecat]')< 1E-5
	assert abs(`mestreg_se_age`i''- `=_se[xb1:`i'.agecat]')< 1E-5
	assert abs(`mestreg_b_int`i''- `=_b[xb1:`i'.agecat#c.trt]')< 1E-5
	assert abs(`mestreg_se_int`i''- `=_se[xb1:`i'.agecat#c.trt]')< 1E-5
}
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1




//============================================================================
//let's test some other stuff in the random intercept models- splines for cts covariate

//splines for age covariate
cap drop agercs*
rcsgen age, gen(agercs) noorthog df(2)
mestreg trt agercs* || id1:, distribution(weibull) cov(unstr) nohr
foreach v in trt agercs1 agercs2 {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
local mestreg_b_cons = _b[/:var(_cons[id1])]


uhtred (stime trt rcs(age, df(2) noorthog) M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr)  noorthog
	
assert abs(`mestreg_b_trt'- `=_b[xb1:trt]')< 1E-5
assert abs(`mestreg_se_trt'- `=_se[xb1:trt]')< 1E-5	
assert abs(`mestreg_b_agercs1'- `=_b[xb1:_rcs1_2_1_1]')< 1E-5
assert abs(`mestreg_se_agercs1'- `=_se[xb1:_rcs1_2_1_1]')< 1E-5	
assert abs(`mestreg_b_agercs2'- `=_b[xb1:_rcs1_2_1_2]')< 1E-5
assert abs(`mestreg_se_agercs2'- `=_se[xb1:_rcs1_2_1_2]')< 1E-5	
assert abs(`mestreg_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-1
	

	
//============================================================================
//let's test some other stuff in the random intercept models->1 df and user-defined knots

//check higher df- now checks are against merlin	
merlin (stime trt age M1[id1]@1, ///
	family(rp, df(3) failure(dead))), 	///
	cov(unstr) 	 noorthog
matrix A=e(rcsmat_1)
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_cons = `=exp(_b[lns1_1:_cons])'


uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(3) failure(dead))), 	///
	cov(unstr)  noorthog
matrix B=e(rcsmat_1)
assert mreldif( A , B ) < 1E-0
matrix drop A B	
foreach v in trt age {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
assert abs(`mer_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-5
	
	
//user-defined knots	
merlin (stime trt age M1[id1]@1, ///
	family(rp, knots(-6.056272 .2149434 1.022322 1.608266) failure(dead))), 	///
	cov(unstr)  noorthog
matrix A=e(rcsmat_1)
local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_cons = `=exp(_b[lns1_1:_cons])'
	
	
uhtred (stime trt age M1[id1]@1, ///
	family(rp, knots(-6.056272 .2149434 1.022322 1.608266) failure(dead))), 	///
	cov(unstr) noorthog
matrix B=e(rcsmat_1)
assert mreldif( A , B ) < 1E-0
matrix drop A B
foreach v in trt age {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
assert abs(`mer_b_cons'- `=exp(_b[lns1_1:_cons])')< 1E-5	



//============================================================================
//error messages- syntax breaks
//2 level random intercept
rcof "uhtred (stime trt age M1[id1]@1, family(rp, df(1) failure(dead))), cov(bla) " == 198
rcof "uhtred (stime trt age M1[id1]@1, family(rp, df(1) failure(dead))), intmethod(bla) " == 198
rcof "uhtred (stime trt age M[id1]@1, family(rp, df(1) failure(dead))) " == 1986
rcof "uhtred (stime trt age M1[id1]@trt, family(rp, df(1) failure(dead))) " == 1986
rcof "uhtred (stime trt age M1[bla]@1, family(rp, df(1) failure(dead))) " == 3598


	
//============================================================================//
// end of file
//============================================================================//

	