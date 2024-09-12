global drive C:\Users\Hannah Bower\Documents\GitHub

//============================================================================//

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
gen id1	= _n
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
stpm2 trt bmi age x1, df(3) scale(h)
est store stpm2


uhtred (stime trt bmi age x1, family(rp, df(3) failure(died) ltruncated(t0)))
est store uhtred

//predictions from stpm2
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density

local list3 h s xb

//ci option not available for cumh and dens in stpm2_pred
est restore stpm2
foreach v in `list1' {
	predict `v', `v' 
	predict `v'_z, `v' zeros 
	predict `v'_at, `v' at(trt 1 bmi 20) 
}

est restore stpm2
foreach v in `list3' {
	cap drop `v' `v'_z `v'_at
	predict `v', `v' ci
	predict `v'_z, `v' zeros ci
	predict `v'_at, `v' at(trt 1 bmi 20) ci

}

est restore uhtred
local j 1
foreach v in `list2' {
	tokenize `list1'
	predict u_``j'', `v' ci
	predict u_``j''_z, `v' zeros ci
	predict u_``j''_at, `v' at(trt 1 bmi 20) ci


	local j `=`j'+1'
}


//check predictions are the same - timevar doesn't seem to work

foreach v in `list3' {
	assert abs(`v' == u_`v') < 1E-5
	assert abs(`v'_lci == u_`v'_lci) < 1E-5
	assert abs(`v'_uci == u_`v'_uci) < 1E-5
	foreach s in z at {
		assert abs(`v'_`s' == u_`v'_`s') < 1E-5
		assert abs(`v'_`s'_lci == u_`v'_`s'_lci) < 1E-5
		assert abs(`v'_`s'_uci == u_`v'_`s'_uci) < 1E-5
	}
}

foreach v in cumh dens {
		assert abs(`v' == u_`v') < 1E-5
		assert abs(`v'_z == u_`v'_z) < 1E-5
		assert abs(`v'_at == u_`v'_at) < 1E-5
}





//============================================================================//
//make cscript stuff (dataset with predicted values)


set seed 72549

clear
set obs 500
gen id1	= _n
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

uhtred (stime trt bmi age x1, family(rp, df(3) failure(died) ltruncated(t0)))

foreach v in hazard survival chazard logchazard density {
	predict `v', `v' ci
	predict `v'_at, `v' at(trt 1 bmi 20) ci
	predict `v'_z, `v' zeros ci
	*predict `v'_tv, `v' timevar(timevar) ci	
}


local keeplist id1
foreach i in hazard survival chazard logchazard density {
	local keeplist `keeplist' `i' `i'_lci `i'_uci
	local keeplist `keeplist' `i'_z `i'_z_lci `i'_z_uci
	local keeplist `keeplist' `i'_at `i'_at_lci `i'_at_uci
	*local keeplist `keeplist' `i'_tv `i'_tv_lci `i'_tv_uci

}

preserve
keep `keeplist'
save "${drive}\uhtred\cert\cert-data\rp_1level_lt_pred.dta", replace
restore

//============================================================================//
// end of file
//============================================================================//

