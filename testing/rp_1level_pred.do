global drive C:\Users\Hannah Bower\Documents\GitHub

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
		maxt(10)

stset stime, f(died) 


//============================================================================//
//going to loop over all models so we can save them with different suffixes

//PH RP
local ms1 stpm2 trt bmi x1 x2 x3, scale(h) df(3)
local mu1 uhtred (_t trt bmi x1 x2 x3, family(rp, df(3) failure(died)))

//RP PH factor
local ms2 stpm2 trt bmi x1 x2 x3 i.agecat, scale(h) df(3)
local mu2 uhtred (_t	trt bmi x1 x2 x3 i.agecat, family(rp, df(3) failure(died))) 	
		
//nPH models
local ms3 stpm2 trt bmi x1 x2 x3, df(3) dftvc(1) tvc(trt) scale(h) knscale(log) orthog
local mu3 uhtred (_t 	trt bmi x1 x2 x3 c.trt#rcs(_t, df(1) log orthog) , family(rp, df(3) failure(died))) 

//nPH models dftvc>1
local ms4 stpm2 trt bmi x1 x2 x3, df(3) dftvc(2) tvc(trt) scale(h) knscale(log) orthog
local mu4 uhtred (_t 	trt bmi x1 x2 x3 c.trt#rcs(_t, df(2) log orthog event), family(rp, df(3) failure(died))) 	

//user-defined knots nonPH
local ms5 stpm2 trt bmi x1 x2 x3, df(3) knotstvc(trt 1.5) bknotstvc(trt -4 2.3) tvc(trt) scale(h) knscale(log) orthog
local mu5 uhtred (_t 	trt bmi x1 x2 x3 c.trt#rcs(_t, knots(-4 1.5 2.3) log orthog ) , family(rp, knots(-4 1 1.5 2) failure(died))) 	

forvalues i=1/5 {
	qui `ms`i''
	est store s`i'
	qui `mu`i''
	est store u`i'
	
}


//predictions from stpm2
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density

forvalues i=1/5 {
	est restore s`i'
	foreach v in `list1' {
		predict `v'`i', `v' //timevar(timevar)
	}
	est restore u`i'
	local j 1
	foreach v in `list2' {
		tokenize `list1'
		predict u_``j''`i', `v' //timevar(timevar)
		local j `=`j'+1'
	}
}


//check predictions are the same - timevar doesn't seem to work
forvalues i=1/5 {
	foreach v in `list1' {
		assert abs(`v'`i' == u_`v'`i') < 1E-5
	}
}

//============================================================================//
//test ci 


//ci- can't compare for chazard and density to stpm2 since they don't work for stpm2
// look like they work for uhtred, but we compare to stpm2 for the other options
local list1 h s xb 
local list2 hazard survival logchazard

est restore s1
foreach v in `list1' {
	cap drop `v' 
	cap drop `v'_lci 
	cap drop `v'_uci
	predict `v', `v' ci
}

//predictions from uhtred
est restore u1
local i 1
foreach v in `list2' {
	cap drop u_`i' u_`i'_lci u_`i'_uci
	tokenize `list1'
	predict u_``i'', `v' ci
	local i `=`i'+1'
}


//check predictions are the same - timevar doesn't seem to work
foreach i in lci uci {
	foreach v in `list1' {
		assert abs(`v'_`i' == u_`v'_`i') < 1E-5
	}
}


//test ci for model with tde
est restore s3
foreach v in `list1' {
	cap drop `v'* 
	predict `v', `v' ci
}

est restore u3
local i 1
foreach v in `list2' {
	cap drop u_``i''* 
	tokenize `list1'
	predict u_``i'', `v' ci
	local i `=`i'+1'
}

foreach i in lci uci {
	foreach v in `list1' {
		assert abs(`v'_`i' == u_`v'_`i') < 1E-5
	}
}


//============================================================================//
//test at

//predictions from stpm2
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density


est restore s1
foreach v in `list1' {
	predict `v'_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1)
}

//predictions from uhtred
est restore u1
local i 1
foreach v in `list2' {
	tokenize `list1'
	predict u_``i''_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1)
	local i `=`i'+1'
}


//check predictions are the same 
foreach v in `list1' {
	assert abs(`v'_at == u_`v'_at) < 1E-5
}

//test in nPH model
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density


est restore s3
foreach v in `list1' {
	cap drop `v'_at
	predict `v'_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1)
}

est restore u3
local i 1
foreach v in `list2' {
	tokenize `list1'
	cap drop u_``i''_at
	predict u_``i''_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1)
	local i `=`i'+1'
}


foreach v in `list1' {
	assert abs(`v'_at == u_`v'_at) < 1E-5
}


//check at with ci 
local list1 h s xb 
local list2 hazard survival  logchazard 


est restore s3
foreach v in `list1' {
	cap drop `v'_at*
	predict `v'_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1) ci
}

est restore u3
local i 1
foreach v in `list2' {
	tokenize `list1'
	cap drop u_``i''_at*
	predict u_``i''_at, `v' at(trt 1 bmi 20 x1 0.5 x2 1 x3 1) ci
	local i `=`i'+1'
}

foreach i in lci uci {
	foreach v in `list1' {
		assert abs(`v'_at_`i' == u_`v'_at_`i') < 1E-5
	}
}

//============================================================================//
//test zeros

//predictions from stpm2
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density


est restore s1
foreach v in `list1' {
	predict `v'_z, `v' zeros
}

//predictions from uhtred
est restore u1
local i 1
foreach v in `list2' {
	tokenize `list1'
	predict u_``i''_z, `v' zeros
	local i `=`i'+1'
}


//check predictions are the same 
foreach v in `list1' {
	assert abs(`v'_z== u_`v'_z) < 1E-5
}



//test in nonPH models- CHECK
/*
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density


est restore s2
foreach v in `list1' {
	cap drop `v'_z
	predict `v'_z, `v' zeros
}



est restore u2
local i 1
foreach v in `list2' {
	tokenize `list1'
	cap drop u_``i''_z
	predict u_``i''_z, `v' zeros
	local i `=`i'+1'
}


foreach v in `list1' {
	assert abs(`v'_z== u_`v'_z) < 1E-5
}
*/



//test with cis 
local list1 h s  xb 
local list2 hazard survival  logchazard 


est restore s1
foreach v in `list1' {
	cap drop `v'_z*
	predict `v'_z, `v' zeros ci
}

//predictions from uhtred
est restore u1
local i 1
foreach v in `list2' {
	tokenize `list1'
	cap drop u_``i''_z*
	predict u_``i''_z, `v' zeros ci
	local i `=`i'+1'
}


//check predictions are the same 
foreach i in lci uci {
	foreach v in `list1' {
		assert abs(`v'_z_`i'== u_`v'_z_`i') < 1E-5
	}
}



//============================================================================//
//test timevar-not working right now
/*
cap drop timevar
range timevar 0 10 200

//predictions from stpm2
local list1 h s cumh xb dens
local list2 hazard survival chazard logchazard density


est restore s1
foreach v in `list1' {
	predict `v'_tv, `v' timevar(timevar)
}

//predictions from uhtred
est restore u1
local i 1
foreach v in `list2' {
	tokenize `list1'
	predict u_``i''_tv, `v' timevar(timevar)
	local i `=`i'+1'
}


//check predictions are the same - timevar doesn't seem to work
foreach v in `list1' {
	assert abs(`v'_tv== u_`v'_tv) < 1E-5
}

//TO TEST cis with timevar and in nonPH
*/

//============================================================================//
//check to see that all other predictions actually run 
//not much thought here, just syntax checks for now
uhtred (_t trt bmi x1 x2 x3, family(rp, df(3) failure(died)))

/*
//none of the following seem to work
foreach v in mu eta totalsurv cif rmst timelost totaltimelost {
	predict `v', `v'	
}
*/

//check the difference and ratio predictions for hazard and survival (which work by themselves)
cap drop hdiff*
cap drop myhdiff
predict hdiff0, h at(trt 0) 
predict hdiff1, h  at(trt 1)

predict hdiff, hdiff at1(trt 0) at2(trt 1)
gen myhdiff=hdiff1-hdiff0
assert abs(abs(hdiff)-abs(myhdiff)) < 1E-5

cap drop sdiff* 
cap drop mysdiff
predict sdiff0, surv at(trt 0) 
predict sdiff1, surv  at(trt 1)

predict sdiff, sdiff at1(trt 0) at2(trt 1)
gen mysdiff=sdiff1-sdiff0
assert abs(abs(sdiff)-abs(mysdiff)) < 1E-5

//============================================================================//
//make cscript stuff (dataset with predicted values)


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
		maxt(10)

stset stime, f(died) 

uhtred (_t trt bmi x1 x2 x3, family(rp, df(3) failure(died)))

foreach v in hazard survival chazard logchazard density {
	predict `v', `v' ci
	predict `v'_at, `v' at(trt 1 bmi 20) ci
	predict `v'_z, `v' zeros ci
	*predict `v'_tv, `v' timevar(timevar) ci	
}


uhtred (_t 	trt bmi x1 x2 x3 c.trt#rcs(_t, df(1) log orthog) , family(rp, df(3) failure(died))) 
foreach v in hazard survival chazard logchazard density {
	predict `v'nph, `v' ci
	predict `v'nph_at, `v' at(trt 1 bmi 20) ci
	*predict `v'nph_tv, `v' timevar(timevar) ci	
}

local keeplist id1
foreach i in hazard survival chazard logchazard density {
	local keeplist `keeplist' `i' `i'_lci `i'_uci
	local keeplist `keeplist' `i'_z `i'_z_lci `i'_z_uci
	local keeplist `keeplist' `i'nph `i'nph_lci `i'nph_uci	

	foreach j in at /*tv z*/ {
		local keeplist `keeplist' `i'_`j' `i'_`j'_lci `i'_`j'_uci
		local keeplist `keeplist' `i'nph_`j' `i'nph_`j'_lci `i'nph_`j'_uci
	}
}

preserve
keep `keeplist'
save "${drive}\uhtred\cert\cert-data\rp_1level_pred.dta", replace
restore

//============================================================================//
//throw error messages
uhtred

//testing syntax errors
cap drop h
rcof "predict h, hj " == 198
rcof "predict h, h jhf" == 198
rcof "prodict h, hazard" == 199
rcof "predict h, hazard at(trtr 1)" == 111
gen h=.
rcof "predict h, hazard at(trtr 1)" == 110
rcof "predict h, hazard timevar(jsdk)" == 111


				
//============================================================================//
// end of file
//============================================================================//

