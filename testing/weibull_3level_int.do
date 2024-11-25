/*============================================================================
program: weib_3level_int
description: testing mixed models with random intercept (3level)
created: by Hannah Bower

tests: 	
	
=============================================================================*/

clear 
set seed 7254

set obs 100
gen id1 = _n
gen age = runiform()
gen sd1 = exp(log(0.1))
gen u1 = rnormal(0,sd1)
expand 100
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(0.1))
gen u2 = rnormal(0,sd2)
expand 10
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.02 u1 1 u2 1) maxt(10)
stset stime1, f(dead1)

replace id2 = _n

//fit models 
mestreg trt age || id1: || id2: , dist(weib) 
est store mestreg

*gsem (stime1 <- trt age M2[id1>id2]@1 M1[id1]@1 , family(weib, failure(dead1))) , //intmethod(gh)
*est store gsem
	
merlin 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 ///
	, family(rp, df(1) failure(dead1))) ///
			, 
est store merlin

uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 ///
	, family(rp, df(1) failure(dead1))) ///
			, 
est store uhtred



//check estimates against each other
est restore merlin	

local j 1
foreach v in trt age {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
forvalues v=1/2 {
	local mer_b_lns`v' =_b[lns`v'_1:_cons]
	local mer_se_lns`v' =_se[lns`v'_1:_cons]
}

est restore uhtred

	
foreach v in trt age {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
forvalues v=1/2{
	assert abs(`mer_b_lns`v''- `=_b[lns`v'_1:_cons]')< 1E-5
	assert abs(`mer_se_lns`v''- `=_se[lns`v'_1:_cons]')< 1E-5
}



//mkassert
uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 ///
	, family(rp, df(1) failure(dead1))) ///
			, 

mkassert eclass


/*============================================================================
testing with jobhistory dataset
=============================================================================*/


clear
webuse jobhistory


stset tend, origin(tstart) fail(failure)

sort birthyear id
order birthyear id, first


//models 
//note here stmixed is a wrapper for merlin so they are the same
mestreg education njobs || birthyear: || id:, 	///
	distribution(weib) nohr intpoints(12) 
est store mestreg
	
//test against gsem
*gsem (_t <- education njobs M2[birthyear>id]@1 M1[birthyear]@1 , family(weib, failure(failure))), intpoints(15) 	
	
/*	
stmixed education njobs ///
	|| birthyear: ///
	|| id: 	///
	, dist(weib) intpoints(15) 
est store stmixed
*/

merlin 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(weib, failure(_d))), intpoints(12)  intmethod(gh)	
est store merlin

uhtred 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(rp, df(1) failure(_d))), 	intpoints(12) intmethod(gh)
est store uhtred

est table mestreg  merlin uhtred




//check b and se
est restore mestreg
foreach v in education njobs {
	local mestreg_b_`v' =_b[_t:`v']
	local mestreg_se_`v' =_se[_t:`v']
}
local mestreg_b_cons1 = `=log(sqrt(_b[/:var(_cons[birthyear])]))'
local mestreg_b_cons2 =`=log(sqrt(_b[/:var(_cons[birthyear>id])]))'


est restore merlin	

local j 1
foreach v in education njobs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
forvalues i=1/2 {
	local mer_b_cons`i' = _b[lns`i'_1:_cons]
}
est restore uhtred

foreach v in education njobs {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-2
	
	assert abs(`mestreg_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`mestreg_se_`v''- `=_se[xb1:`v']')< 1E-2
	
}
forvalues v=1/2{
	assert abs(`mer_b_cons`v''- `=_b[lns`v'_1:_cons]')< 1E-2
	assert abs(`mestreg_b_cons`v''- `=_b[lns`v'_1:_cons]')< 1E-2
}

est drop mestreg merlin uhtred




/*============================================================================
end of file
=============================================================================*/
