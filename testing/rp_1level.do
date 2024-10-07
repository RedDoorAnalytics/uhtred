/*============================================================================
program: rp_1level
description: testing for Royston-Parmer 1 level models
created: by Hannah Bower

tests: 
	 RP PH model
	 RP nonPH model
	 Factor variables
	 user defined knots (PH and nonPH)
	 interactions
 
=============================================================================*/

//============================================================================//
//sim data- PH

clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		/*tde(trt 0.01) tdefunc(log({t}))*/	///
		maxt(10) 

stset stime, f(died) 

//============================================================================//
//RP PH model- simple Weib

//testing
stpm2 trt bmi x1 x2 x3, scale(h) df(1)
est store stpm2

foreach v in trt bmi x1 x2 x3 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(1) failure(died))) 	///
		,
est store merlin
		
local j 1
foreach v in trt bmi x1 x2 x3 {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}

uhtred (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(1) failure(died))) 	///
		,
est store uhtred

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
			
//mkassert
uhtred (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(1) failure(died))) 	///
		,

mkassert eclass




//============================================================================//
//RP PH model- simple 3df

//testing
stpm2 trt bmi x1 x2 x3, scale(h) df(3)
est store stpm2
matrix A = e(R_bh)
matrix list A


foreach v in trt bmi x1 x2 x3 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,
est store merlin
		
local j 1
foreach v in trt bmi x1 x2 x3 {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}

uhtred (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,
est store uhtred
matrix B = e(rcsrmat_1)
matrix list B


est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
			
//mkassert
uhtred (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,

mkassert eclass

assert colsof(A)==colsof(B)
assert rowsof(A)==rowsof(B)
assert mreldif( A , B ) < 1E-0

matrix drop A B

// make cscript by hand
uhtred (_t 	trt bmi x1 x2 x3   			///
		, family(rp, df(3) failure(died))) 	///
		,

mkassert eclass
		
qui {
	mat A = J(4,4,0)
	mat A[1,1] = .79717296
	mat A[1,2] = -12.329397
	mat A[1,3] = -5.0923271
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = 2.6519554
	mat A[2,3] = 1.1812266
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .09176479 
	mat A[3,4] =  0

	mat A[4,1] =  1.8714958
	mat A[4,2] = -47.20906
	mat A[4,3] =  -19.07834
	mat A[4,4] =  1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-0
matrix drop A



//============================================================================//
//RP PH model factor variables- note doesn't work in merlin

stpm2 trt bmi x1 x2 x3 i.agecat, scale(h) df(3)
est store stpm2

foreach v in trt bmi x1 x2 x3 agecat {
	if `v'== agecat {
		forvalues j=1/4 {
			local stpm2_b_`v'`j' =_b[xb:`j'.`v']
			local stpm2_se_`v'`j' =_se[xb:`j'.`v']
		}
	}
	else{
		local stpm2_b_`v' =_b[xb:`v']
		local stpm2_se_`v' =_se[xb:`v']
	}
}

//merlin no factor variables
cap drop agecat?
tab agecat, gen(agecat)
merlin (_t	trt bmi x1 x2 x3 agecat2 agecat3 agecat4 agecat5			///
		, family(rp, df(1) failure(died))) 	///
		,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
forvalues i=1/4 {
	local j `=`i'+5'
	local mer_b_agecat`i' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_agecat`i' =_se[_cmp_1_`j'_1:_cons]
}


//test uhtred without factor variables too
uhtred (_t	trt bmi x1 x2 x3 agecat2 agecat3 agecat4 agecat5	///
		, family(rp, df(1) failure(died))) 	///
		,
est store uhtred_nof

//now with 
uhtred (_t	trt bmi x1 x2 x3 i.agecat			///
		, family(rp, df(1) failure(died))) 	///
		,
est store uhtred


est table stpm2 merlin uhtred*
estimates drop stpm2 uhtred*


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-4
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}

forvalues i=1/4 {
	assert abs(`stpm2_b_agecat`i''- `=_b[xb1:`i'.agecat]')< 1E-2
	assert abs(`stpm2_se_agecat`i''- `=_se[xb1:`i'.agecat]')< 1E-5
	
	assert abs(`mer_b_agecat`i''- `=_b[xb1:`i'.agecat]')< 1E-5
	assert abs(`mer_se_agecat`i''- `=_se[xb1:`i'.agecat]')< 1E-5
}

			
//mkassert
uhtred (_t	trt bmi x1 x2 x3 i.agecat			///
		, family(rp, df(1) failure(died))) 	///
		,

mkassert eclass


//============================================================================//
//user-defined knots PH model

stpm2 trt bmi x1 x2 x3, knots(-3.5 1 1.5 2) knscale(log) scale(h) 
est store stpm2

foreach v in trt bmi x1 x2 x3  {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3  {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}


uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
est store uhtred

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-2
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}


//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
mkassert eclass


//============================================================================//
//user-defined BOUNDARY knots PH model

stpm2 trt bmi x1 x2 x3, df(3) bknots(-3.5  2) knscale(log) scale(h) 
est store stpm2

foreach v in trt bmi x1 x2 x3  {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3  {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}


uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
est store uhtred

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred

//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}


//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-3.5 1 1.5 2) failure(died))) 	///
,
mkassert eclass


//============================================================================//
//interaction term in PH - cat#cat

stpm2 i.trt bmi x1 x2 x3 i.trt#i.agecat, scale(h) df(3) 
est store stpm2

local stpm2_b_trt = _b[xb:1.trt]
local stpm2_se_trt = _se[xb:1.trt]

foreach v in bmi x1 x2 x3 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}
foreach v in 0b 1 {
	forvalues i=1/4 {
		local stpm2_b_`v'`i'age =_b[xb:`v'.trt#`i'.agecat]
		local stpm2_se_`v'`i'age =_se[xb:`v'.trt#`i'.agecat]
	}
}


uhtred (_t 	i.trt bmi x1 x2 x3 			///
		i.trt#i.agecat 	///
		, family(rp, df(3) failure(died))) 	///
		,
est store uhtred
est table stpm2 uhtred
est drop stpm2 uhtred

//check b and se
assert abs(`stpm2_b_trt'- `=_b[xb1:1.trt]')< 1E-3
assert abs(`stpm2_se_trt'- `=_se[xb1:1.trt]')< 1E-5
foreach v in bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5
}
foreach v in 0b 1 {
	forvalues i=1/4 {
		assert abs(`stpm2_b_`v'`i'age'- `=_b[xb1:`v'.trt#`i'.agecat]')< 1E-3
		assert abs(`stpm2_se_`v'`i'age'- `=_se[xb1:`v'.trt#`i'.agecat]')< 1E-5
	}
}

uhtred (_t 	i.trt bmi x1 x2 x3 			///
		i.trt#i.agecat 	///
		, family(rp, df(3) failure(died))) 	///
		,

mkassert eclass


//============================================================================//
//interaction term in PH - cts cts

stpm2 trt bmi c.x1#c.x2 x3 , scale(h) df(3) 
est store stpm2

local stpm2_b_int = _b[xb:c.x1#c.x2]
local stpm2_se_int = _se[xb:c.x1#c.x2]

foreach v in bmi trt bmi x3 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


uhtred (_t 	trt bmi c.x1#c.x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,
est store uhtred
est table stpm2 uhtred
est drop stpm2 uhtred

//check b and se
assert abs(`stpm2_b_int'- `=_b[xb1:c.x1#c.x2]')< 1E-3
assert abs(`stpm2_se_int'- `=_se[xb1:c.x1#c.x2]')< 1E-5
foreach v in trt bmi x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5
}


uhtred (_t 	trt bmi c.x1#c.x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,
mkassert eclass

//============================================================================//
//check orthog matrices


//higher df
stpm2 trt bmi x1 x2 x3, scale(h) df(5) 
matrix A=e(R_bh)


uhtred (_t 	trt bmi x1 x2 x3   			///
		, family(rp, df(5) failure(died))) 	///
		,
matrix B= e(rcsrmat_1)

assert colsof(A) == colsof(B)
assert rowsof(A) == rowsof(B)
assert mreldif( A , B ) < 1E-0
matrix drop A B

// make cscript by hand

uhtred (_t 	trt bmi x1 x2 x3   			///
		, family(rp, df(5) failure(died))) 	///
		,

qui {
	mat A = J(6,6,0)
	mat A[1,1] = .79717296
	mat A[1,2] = -15.764215
	mat A[1,3] = -10.66939
	mat A[1,4] = -6.3882997
	mat A[1,5] = -2.8095331
	mat A[1,6] = 0

	mat A[2,1] = 0
	mat A[2,2] = 3.1802201
	mat A[2,3] = 2.3450493
	mat A[2,4] = 1.4630702
	mat A[2,5] = .65478289
	mat A[2,6] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .1654318 
	mat A[3,4] = .16638174
	mat A[3,5] = .08792099 
	mat A[3,6] = 0

	mat A[4,1] = 0
	mat A[4,2] = 0
	mat A[4,3] = 0
	mat A[4,4] = .02533339
	mat A[4,5] = .02089929
	mat A[4,6] = 0

	mat A[5,1] = 0
	mat A[5,2] = 0
	mat A[5,3] = 0
	mat A[5,4] = 0
	mat A[5,5] = .00455911
	mat A[5,6] = 0

	mat A[6,1] =  1.8714958
	mat A[6,2] = -61.529876
	mat A[6,3] = -40.570331 
	mat A[6,4] =  -23.991209
	mat A[6,5] =  -10.496074
	mat A[6,6] =  1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-0
matrix drop A







//============================================================================//
//non-PH model 
//simulate nonPH data and then check agains Weib first

clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) 

stset stime, f(died) 

//testing
stpm2 trt bmi x1 x2 x3, df(1) dftvc(1) tvc(trt) scale(h) knscale(log) orthog
est store stpm2

foreach v in trt bmi x1 x2 x3 _rcs_trt1 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(1) failure(died))) 	///
                , 
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 trtrcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}


uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(1) failure(died))) 	///
                , 
est store uhtred
est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
	
assert abs(`stpm2_b__rcs_trt1'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`stpm2_se__rcs_trt1'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
	
assert abs(`mer_b_trtrcs'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`mer_se_trtrcs'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5

	
//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(1) failure(died))) 	///
		,
		
mkassert eclass

//============================================================================//
//non-PH model  3 df RP


//testing
stpm2 trt bmi x1 x2 x3, df(3) dftvc(1) tvc(trt) scale(h) knscale(log) orthog
est store stpm2
matrix A1 = e(R_bh)
matrix A2 = e(R_trt)
matrix list A2


foreach v in trt bmi x1 x2 x3 _rcs_trt1 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
                , 
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 trtrcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}


uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
                , 
est store uhtred
matrix B1 = e(rcsrmat_1)
matrix B2 = e(rmat_1_6_2)
matrix list B2

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred

assert mreldif(A1 , B1) < 1E-0
assert mreldif(A2 , B2) < 1E-0
matrix drop A1 A2 B1 B2



//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
	
assert abs(`stpm2_b__rcs_trt1'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`stpm2_se__rcs_trt1'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
	
assert abs(`mer_b_trtrcs'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`mer_se_trtrcs'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5

	
//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
		,
		
mkassert eclass


qui {
	mat A = J(4,4,0)
	mat A[1,1] = .79717296
	mat A[1,2] = -12.329397
	mat A[1,3] = -5.0923271
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = 2.6519554
	mat A[2,3] = 1.1812266
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .09176479 
	mat A[3,4] = 0

	mat A[4,1] = 1.8714958
	mat A[4,2] = -47.20906
	mat A[4,3] = -19.07834
	mat A[4,4] = 1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-0

qui {
	mat B = J(2,2,0)
	mat B[1,1] = .79717296
	mat B[1,2] = 0

	mat B[2,1] = 1.8714958 
	mat B[2,2] = 1
}
assert mreldif( e(rmat_1_6_2) , B ) < 1E-0
matrix drop A B






//============================================================================//
//non-PH model >1 df in tde
		
//testing
stpm2 trt bmi x1 x2 x3, df(3) dftvc(2) tvc(trt) scale(h) knscale(log) orthog
est store stpm2


foreach v in trt bmi x1 x2 x3 _rcs_trt1 _rcs_trt2 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, df(2) log orthog event) 	///
		, family(rp, df(3) failure(died))) 	///
                , 
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 trtrcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_trtrcs2 =_b[_cmp_1_6_2:_cons]
local mer_se_trtrcs2 =_se[_cmp_1_6_2:_cons]


uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(2) log orthog event) 	///
		, family(rp, df(3) failure(died))) 	///
                , 
est store uhtred
est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


//check estimates of b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
	
assert abs(`stpm2_b__rcs_trt1'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-3
assert abs(`stpm2_se__rcs_trt1'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`stpm2_b__rcs_trt2'- `=_b[tb1:c.trt#c._rcs1_6_2_2]')< 1E-3
assert abs(`stpm2_se__rcs_trt2'- `=_se[tb1:c.trt#c._rcs1_6_2_2]')< 1E-5
	
assert abs(`mer_b_trtrcs'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`mer_se_trtrcs'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`mer_b_trtrcs2'- `=_b[tb1:c.trt#c._rcs1_6_2_2]')< 1E-5
assert abs(`mer_se_trtrcs2'- `=_se[tb1:c.trt#c._rcs1_6_2_2]')< 1E-5

	
//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(2) log orthog event) 	///
		, family(rp, df(3) failure(died))) 	///
		,
		
mkassert eclass








//============================================================================//
//user-defined knots TDE

stpm2 trt bmi x1 x2 x3, knots(1 1.5) bknots(-4 2) knotstvc(trt 1.5) bknotstvc(trt -4 2.3) tvc(trt) ///
	scale(h) knscale(log) orthog
est store stpm2

foreach v in trt bmi x1 x2 x3 _rcs_trt1 _rcs_trt2 {
	local stpm2_b_`v' =_b[xb:`v']
	local stpm2_se_`v' =_se[xb:`v']
}


merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, knots(-4 1.5 2.3) log orthog) 	///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
		,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 rcs1 {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}
local mer_b_rcs2 =_b[_cmp_1_6_2:_cons]
local mer_se_rcs2 =_se[_cmp_1_6_2:_cons]


uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, knots(-4 1.5 2.3) log orthog ) 	///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
		,
est store uhtred
est table stpm2 merlin uhtred
est drop stpm2 merlin uhtred

//check b and se
foreach v in trt bmi x1 x2 x3 {
	assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-5

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}

forvalues v=1/2 {
	assert abs(`stpm2_b__rcs_trt`v''- `=_b[tb1:c.trt#_rcs1_6_2_`v']')< 1E-3
	assert abs(`stpm2_se__rcs_trt`v''- `=_se[tb1:c.trt#_rcs1_6_2_`v']')< 1E-5

	assert abs(`mer_b_rcs`v''- `=_b[tb1:c.trt#_rcs1_6_2_`v']')< 1E-5
	assert abs(`mer_se_rcs`v''- `=_se[tb1:c.trt#_rcs1_6_2_`v']')< 1E-5
}

uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, knots(-4 1.5 2.3) log orthog ) 	///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
		,

mkassert eclass




//============================================================================//
//throw error messages


//new var doesn't exist
rcof "uhtred (_t newvar trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 111
rcof "uhtred (_t  trt i.bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog) , family(rp, df(3) failure(died)))" == 452



//tde misspecified
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#ecs(_t, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 198
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(time, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 111
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(_t, df(1) lrg orthog) 	, family(rp, knots(3) failure(died)))" == 3598
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(_t, df(1) log prthog) 	, family(rp, knots(3) failure(died)))" == 3598
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(_t, df(-1) log orthog) 	, family(rp, knots(3) failure(died)))" ==  1986
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(_t, df(1 5) log orthog) 	, family(rp, knots(3) failure(died)))" == 1986
rcof "uhtred (_t trt bmi x1 x2 x3 trt#rcs(_t, df(1) knots(2 4) log orthog) 	, family(rp, knots(3) failure(died)))" == 1986
rcof "uhtred (_t trt bmi x1 x2 x3 trt#rcs(_t,  knots(2 2) log orthog) 	, family(rp, knots(3) failure(died)))" == 3598

//family options misspecifed
rcof "uhtred (_t trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog) 	, family(fam, knots(3) failure(died)))" == 198
rcof "uhtred (_t  trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog) 	, family(rp, knots(3) failure(newvar)))" == 111
rcof "uhtred (_t  trt bmi x1 x2 x3 	, family(rp, knots(3)))" == 198
*rcof "uhtred (_t trt bmi x1 x2 x3 , family(rp,    failure(died))" == 198

//other incorrect syntax
rcof "uhtred (_t trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog 	, family(rp, knots(3) failure(died)))" == 198
rcof "uhtred (_t trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog)	 family(rp, knots(3) failure(died)))" == 198
rcof "uhtred _t trt bmi x1 x2 x3, family(rp, df(3) failure(died)))" == 198
rcof "uhtred _t trt bmi x1 x2 x3, family(rp, df(3 failure(died)))" == 198

				
//============================================================================//
// end of file
//============================================================================//

