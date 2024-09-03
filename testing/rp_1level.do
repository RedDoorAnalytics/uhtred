
//============================================================================//
//sim data

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
		cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) //ltruncated(t0)

stset stime, f(died) //enter(t0)


//============================================================================//
//RP PH model- simple


//testing
stpm2 trt bmi x1 x2 x3, scale(h) df(3)
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
foreach v in trt bmi x1 x2 x3 {
	*assert abs(`stpm2_b_`v''- `=_b[xb1:`v']')< 1E-8
	*assert abs(`stpm2_se_`v''- `=_se[xb1:`v']')< 1E-8

	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-8
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-8

}
			
//mkassert
uhtred (_t	trt bmi x1 x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
		,

mkassert eclass

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



uhtred (_t	trt bmi x1 x2 x3 i.agecat			///
		, family(rp, df(1) failure(died))) 	///
		,
est store uhtred


est table stpm2 merlin uhtred*
estimates drop stpm2 uhtred*


foreach v in trt bmi x1 x2 x3 {
	
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}

forvalues i=1/4 {
	assert abs(`mer_b_agecat`i''- `=_b[xb1:`i'.agecat]')< 1E-5
	assert abs(`mer_se_agecat`i''- `=_se[xb1:`i'.agecat]')< 1E-5

}
	
			
//mkassert
uhtred (_t	trt bmi x1 x2 x3 i.agecat			///
		, family(rp, df(1) failure(died))) 	///
		,


mkassert eclass


//============================================================================//
//non-PH model 

//testing
stpm2 trt bmi x1 x2 x3, df(3) dftvc(1) tvc(trt) scale(h) knscale(log) orthog
est store stpm2

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
est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred



foreach v in trt bmi x1 x2 x3 {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
	
assert abs(`mer_b_trtrcs'- `=_b[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5
assert abs(`mer_se_trtrcs'- `=_se[tb1:c.trt#c._rcs1_6_2_1]')< 1E-5

	
	
//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
		,
		
mkassert eclass



//============================================================================//
//non-PH model >1 df in tde

				
//testing
stpm2 trt bmi x1 x2 x3, df(3) dftvc(2) tvc(trt) scale(h) knscale(log) orthog
est store stpm2

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



foreach v in trt bmi x1 x2 x3 {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
	
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
//user-defined knots PH model
stpm2 trt bmi x1 x2 x3, knots(-4 1 1.5 2) knscale(log) scale(h) 
est store stpm2


merlin (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,
est store merlin
local j 1
foreach v in trt bmi x1 x2 x3  {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}



uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,
est store uhtred

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


foreach v in trt bmi x1 x2 x3 {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}


//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,
mkassert eclass

//============================================================================//
//user-defined BOUNDARY knots PH model
stpm2 trt bmi x1 x2 x3, df(3) bknots(-4  2) knscale(log) scale(h) 
est store stpm2


merlin (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2)failure(died))) 	///
,
est store merlin
local j 1
foreach v in trt bmi x1 x2 x3  {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}



uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,
est store uhtred

est table stpm2 merlin uhtred
estimates drop stpm2 merlin uhtred


foreach v in trt bmi x1 x2 x3 {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-5
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}


//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,
mkassert eclass

//============================================================================//
//user-defined knots TDE

/*
//testing
stpm2 trt bmi x1 x2 x3, df(3) knotstvc(trt 1.5) bknotstvc(trt -4 2.3) tvc(trt) scale(h) knscale(log) orthog
est store stpm2

merlin (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, knots(-4 1.5 2.3) log orthog) 	///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
		,
est store merlin


uhtred (_t 	trt bmi x1 x2 x3 			///
		trt#rcs(_t, knots(-4 1.5 2.3) log orthog ) 	///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
		,
est store uhtred
est table merlin uhtred

*/





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

