/*============================================================================
program: rp_1level_lt
description: testing for Royston-Parmer 1 level models left truncation
created: by Hannah Bower

tests: PH and nonPH models
	 
 
=============================================================================*/

local drive C:\Users\Hannah Bower\Documents\GitHub
cd "`drive'\uhtred"
adopath ++ "`drive'\uhtred"
adopath ++ "`drive'\uhtred/uhtred"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear


set seed 72549

clear
set obs 500
gen id 	= _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen bmi = rnormal(30,3)


survsim stime died, dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1) maxt(10)

gen t0 = 0
replace t0 = runiform() * 2 //if _n>200

drop if stime<t0


stset stime, enter(t0) f(died)

//============================================================================//
//PH model

//start by checking weibull since this is what is simulated- change obs to higher number
uhtred (stime trt age bmi x1, family(rp, df(1) failure(died) ltruncated(t0)))



stpm2 trt bmi age x1, df(3) scale(h)
est store stpm2
matrix A = e(R_bh)

foreach v in trt bmi age x1 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}



merlin (stime trt bmi age x1, family(rp, df(3) ltruncated(t0) failure(died)) timevar(stime)), level(70)
est store merlin
		
local j 1
foreach v in trt bmi age x1 {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
matrix B = e(rcsrmat_1)


uhtred (stime trt bmi age x1, family(rp, df(3) failure(died) ltruncated(t0)))
est store uhtred

est table stpm2 merlin uhtred
est drop stpm2 merlin uhtred


//check orthog matrices and estimates of b and se
assert mreldif(e(rcsrmat_1) , A) < 1E-1
assert mreldif(e(rcsrmat_1) , B) < 1E-1
matrix drop A B

foreach v in trt bmi age x1 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-3

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-3
}


//make assert
uhtred (stime trt bmi age x1, family(rp, df(3) failure(died) ltruncated(t0)))
mkassert

qui {
	mat A = J(4,4,0)
	mat A[1,1] = .53387312
	mat A[1,2] = -2.6951611
	mat A[1,3] = -1.1767272
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = .67730034
	mat A[2,3] = .32481165
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .03234029
	mat A[3,4] =  0

	mat A[4,1] =  2.0129868
	mat A[4,2] = -8.5758075
	mat A[4,3] =  -3.6244006
	mat A[4,4] =  1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-5
matrix drop A


//============================================================================//
//non- PH test
//sim nonPH data
set seed 72549

clear
set obs 500
gen id 	= _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen bmi = rnormal(30,3)


survsim stime died, dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1) ///
	tde(trt 0.03) tdefunction(log({t}))maxt(10)

gen t0 = 0
replace t0 = runiform() * 2 //if _n>200

drop if stime<t0


stset stime, enter(t0) f(died)

//start by checking weibull since this is what is simulated- change obs to larger number
uhtred (stime trt age bmi x1 ///
		c.trt#rcs(stime, df(1) log), family(rp, df(1) failure(died) ltruncated(t0)))

//now check against other models for >> df
stpm2 trt bmi age x1, df(2) scale(h) tvc(trt) dftvc(2)
est store stpm2

matrix A=e(R_bh)
matrix Ar=e(R_trt)

foreach v in trt bmi age x1 _rcs_trt1 _rcs_trt2 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (stime trt bmi age x1 trt#rcs(stime, df(2) log orthog event), family(rp, df(2) ltruncated(t0) failure(died)) ), nolog
est store merlin

matrix B=e(rcsrmat_1)
matrix Br=e(rmat_1_5_2)

local j 1
foreach v in trt bmi age x1 trtrcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_trtrcs2 =_b[_cmp_1_5_2:_cons]
local mer_se_trtrcs2 =_se[_cmp_1_5_2:_cons]



uhtred (stime trt bmi age x1 c.trt#rcs(stime, df(2) log orthog event), family(rp, df(2) failure(died) ltruncated(t0)))
est store uhtred

est table stpm2 merlin uhtred
est drop stpm2 merlin uhtred

*assert mreldif(e(rcsrmat_1) , A) < 1E-2
assert mreldif(e(rcsrmat_1) , B) < 1E-2

*assert mreldif(e(rmat_1_5_2) , Ar) < 1E-2
assert mreldif(e(rmat_1_5_2) , Br) < 1E-2
matrix drop A Ar B Br


//check estimates of b and se
foreach v in trt bmi age x1 {
*	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-5
*	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-4
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-4
}
/*	
assert abs(`stpm2_b__rcs_trt1'- `=_b[tb1:c.trt#c._rcs1_5_2_1]')< 1E-3
assert abs(`stpm2_se__rcs_trt1'- `=_se[tb1:c.trt#c._rcs1_5_2_1]')< 1E-5
assert abs(`stpm2_b__rcs_trt2'- `=_b[tb1:c.trt#c._rcs1_5_2_2]')< 1E-3
assert abs(`stpm2_se__rcs_trt2'- `=_se[tb1:c.trt#c._rcs1_5_2_2]')< 1E-5
*/	
assert abs(`mer_b_trtrcs'- `=_b[tb1:c.trt#c._rcs1_5_2_1]')< 1E-5
assert abs(`mer_se_trtrcs'- `=_se[tb1:c.trt#c._rcs1_5_2_1]')< 1E-5
assert abs(`mer_b_trtrcs2'- `=_b[tb1:c.trt#c._rcs1_5_2_2]')< 1E-5
assert abs(`mer_se_trtrcs2'- `=_se[tb1:c.trt#c._rcs1_5_2_2]')< 1E-5

//mkassert
uhtred (stime trt bmi age x1 c.trt#rcs(stime, df(2) log orthog event), family(rp, df(2) failure(died) ltruncated(t0)))
mkassert

qui {
	mat A = J(3,3,0)
	mat A[1,1] = .53358255
	mat A[1,2] = -1.8181829
	mat A[1,3] = 0

	mat A[2,1] = 0
	mat A[2,2] = .48823343
	mat A[2,3] = 0

	mat A[3,1] =  2.0097647
	mat A[3,2] = -5.6294308 
	mat A[3,3] = 1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-5
matrix drop A

qui {
	mat B = J(3,3,0)
	mat B[1,1] = .53358255
	mat B[1,2] = -1.8181829
	mat B[1,3] = 0

	mat B[2,1] = 0
	mat B[2,2] = .48823343 
	mat B[2,3] = 0

	mat B[3,1] = 2.0097647
	mat B[3,2] = -5.6294308  
	mat B[3,3] =  1
}
assert mreldif( e(rmat_1_5_2) , B ) < 1E-5
matrix drop B




//============================================================================//
//error checks syntax

rcof "uhtred (_t trt bmi age x1, family(rp, df(2) failure(died) ltruncated(t02)))" == 111
rcof "uhtred (_t trt bmi age x1, family(rp, df(2) failure(died) letruncated(t0)))" ==3598
rcof "uhtred (_t trt bmi age x1, family(rp, df(2) failure(died)) ltruncated(t0))" == 198

/*============================================================================//
 end of file
//============================================================================*/




