/*============================================================================
program: rp_cr
description: testing for Royston-Parmer with competing risks
created: by Hannah Bower

tests: 
	 RP PH model
	 RP nonPH model
	 Factor variables
	 user defined knots (PH and nonPH)
	 interactions
 
=============================================================================*/

//data setup
use ./data/colon, clear
 
gen female = (sex==2)
set seed 978
expand 2

bysort id: gen cause=_n		

gen cancer=(cause==1)		
gen other=(cause==2)		

tab agegrp, gen(ag)
forvalues i = 0/3 {
	gen age`i'can=(agegrp==`i' & cancer==1) 
	gen age`i'oth=(agegrp==`i' & other==1) 
}

gen fem_can = female*cancer
gen fem_other = female*other

gen event=(cause==status)  

stset surv_mm, failure(event) scale(12) exit(time 120.5)


//checks for diff in stpm2 and merlin/uhtred- just look at one cause to check 
stpm2 female age1can age2can age3can if cause==1, df(3) scale(h) 
merlin (_t female age1can age2can age3can if cause==1, family(rp, df(3) failure(event)))
uhtred (_t female age1can age2can age3can if cause==1, family(rp, df(3) failure(event)))


//stpm2 models for competing risks
stpm2 female ag2 ag3 ag4 if cause==1, df(3) scale(h) 
est store stpm2_c1
stpm2 female ag2 ag3 ag4 if cause==2, df(3) scale(h) 
est store stpm2_c2
/*getting different vals for stpm2- what have I done?*/


//ignore differenes with stpm2 for now 
//ph models - merlin
merlin 	(_t 	female age1can age2can age3can	///
		if cause==1								///
		, family(rp, failure(event) df(3)))		///
		(_t 	female age1oth age2oth age3oth	///
		if cause==2								///
		, family(rp, failure(event) df(3)))		///
		,
		
est store merlin 
matrix A1=e(rcsrmat_1)
matrix A2=e(rcsrmat_2)

local list1 female age1can age2can age3can
local list2 female age1oth age2oth age3oth
forvalues cr=1/2 {	
	local j 1
	foreach v in `list`cr'' {
		local mer_b_`v'`cr' =_b[_cmp_`cr'_`j'_1:_cons]
		local mer_se_`v'`cr' =_se[_cmp_`cr'_`j'_1:_cons]
		local j `=`j'+1'
	}
}

//ph model- uhtred
uhtred 	(_t 	female age1can age2can age3can	///
		if cause==1								///
		, family(rp, failure(event) df(3)))		///
		(_t 	female age1oth age2oth age3oth	///
		if cause==2								///
		, family(rp, failure(event) df(3)))		///
		,
		
est store uhtred
matrix B2=e(rcsrmat_2)
matrix B1=e(rcsrmat_1)


est table stpm2_* merlin uhtred
est drop stpm2_* merlin uhtred

//check b and se
forvalues cr=1/2 {
	foreach v in `list`cr''{
		assert abs(`mer_b_`v'`cr''- `=_b[xb`cr':`v']')< 1E-3
		assert abs(`mer_se_`v'`cr''- `=_se[xb`cr':`v']')< 1E-3
	}
}

//check orthog matrices
assert mreldif(A1 , B1 ) < 1E-5
assert mreldif(A2 , B2 ) < 1E-5
matrix drop A1 A2 B1 B2

//=============================================================================
//nonPH models 
merlin 	(_t 	female ag2 ag3 ag4						///
			female#rcs(_t , df(2) log event orthog)		///
			if cause==1									///
			, family(rp, failure(event) df(3)) timevar(_t))	///
		(_t 	female ag2 ag3 ag4						///
			if cause==2									///
			, family(rp, failure(event) df(3)))			///
			, 
est store merlin
matrix A1=e(rcsrmat_1)
matrix A2=e(rcsrmat_2)		


forvalues cr=1/2 {	
	local j 1
	foreach v in female ag2 ag3 ag4 {
		local mer_b_`v'`cr' =_b[_cmp_`cr'_`j'_1:_cons]
		local mer_se_`v'`cr' =_se[_cmp_`cr'_`j'_1:_cons]
		local j `=`j'+1'
	}
}
local mer_b_frcs21 =_b[_cmp_1_5_2:_cons]
local mer_se_frcs21 =_se[_cmp_1_5_2:_cons]

			
uhtred 	(_t 	female ag2 ag3 ag4						///
			c.female#rcs(_t , df(2) log event orthog)	///
			if cause==1									///
			, family(rp, failure(event) df(3)) )		///
		(_t 	female ag2 ag3 ag4						///
			if cause==2									///
			, family(rp, failure(event) df(3)))			///
			,  
est store uhtred
matrix B1=e(rcsrmat_1)
matrix B2=e(rcsrmat_2)			


est table merlin uhtred
est drop  merlin uhtred



forvalues cr=1/2 {
	foreach v in female ag2 ag3 ag4 {
		assert abs(`mer_b_`v'`cr''- `=_b[xb`cr':`v']')< 1E-5
		assert abs(`mer_se_`v'`cr''- `=_se[xb`cr':`v']')< 1E-5
	}
}
assert abs(`mer_b_frcs11'- `=_b[tb1:c.female#c._rcs1_5_2_1]')< 1E-5
assert abs(`mer_b_frcs21'- `=_b[tb1:c.female#c._rcs1_5_2_2]')< 1E-5
assert abs(`mer_se_frcs11'- `=_se[tb1:c.female#c._rcs1_5_2_1]')< 1E-5
assert abs(`mer_se_frcs21'- `=_se[tb1:c.female#c._rcs1_5_2_2]')< 1E-5






//=============================================================================



