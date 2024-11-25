/*============================================================================
program: rp_1level_mts
description: testing for Royston-Parmer 1 level models with multiple timescales
created: by Hannah Bower

tests: moffset and offset options
	 
 
=============================================================================*/

//============================================================================//
//sim data

clear
set seed 249857
set obs 10000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen agec = age - 50
gen year = 1990 + floor(20*runiform())
gen yearc = year - 2000

gen x1 = rnormal()
gen bmi = rnormal(30,3)


survsim stime died, maxtime(5) cov(trt -0.5) 	///
hazard(	0.1:*1.2:*{t}:^0.2 :*			///
        exp(					///
                        0.1 :* (agec :+ {t}) 	///
                        :- 0.1 :* (yearc :+ {t}) ///
                )				///
        )



stset stime, f(died)

//============================================================================//
//check offset & moffset (against merlin)

merlin (stime 	trt bmi x1   			///
	trt#rcs(stime df(1) log orthog event offset(age)) ///
	, family(rp, df(3) failure(died))) 	///
	,
est store merlin


uhtred (stime 	trt 						///
        rcs(stime, df(1) offset(agec)) 				///
        rcs(stime, df(1) offset(yearc)) 			///
        , family(rp, df(3) failure(died)))



//============================================================================//
//old code
merlin (_t 	trt bmi x1 x2 x3   			///
	trt#rcs(_t, df(1) log orthog event offset(age)) ///
	, family(rp, df(3) failure(died))) 	///
	,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 rcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}

uhtred (_t 	trt bmi x1 x2 x3   			///
	c.trt#rcs(_t, df(1) log orthog event offset(age)) ///
	, family(rp, df(3) failure(died))) 	///
	,
estimate store uhtred

estimates table merlin uhtred
estimates drop merlin uhtred
foreach v in trt bmi x1 x2 x3 {
	assert abs(`mer_b_`v''- `=_b[xb1:`v']')< 1E-3
	assert abs(`mer_se_`v''- `=_se[xb1:`v']')< 1E-5
}
assert abs(`mer_b_rcs'- `=_b[tb1:c.trt#c.rcs1_6_2_1]')< 1E-3
assert abs(`mer_se_rcs'- `=_se[tb1:c.trt#c.rcs1_6_2_1]')< 1E-5

	
	
	
//moffset	
cap drop negage
gen negage=0-age

merlin (_t 	trt bmi x1 x2 x3   			///
	trt#rcs(_t, df(1) log orthog event moffset(negage)) ///
	, family(rp, df(3) failure(died))) 	///
	,
est store merlin

local j 1
foreach v in trt bmi x1 x2 x3 rcs {
	local mer_b_`v' =_b[_cmp_1_`j'_1:_cons]
	local mer_se_`v' =_se[_cmp_1_`j'_1:_cons]
	local j `=`j'+1'
}

uhtred (_t 	trt bmi x1 x2 x3   			///
	c.trt#rcs(_t, df(1) log orthog event moffset(negage)) ///
	, family(rp, df(3) failure(died))) 	///
	,
estimate store uhtred

estimates table merlin uhtred



uhtred (_t 	trt bmi x1 x2 x3   			///
		c.trt#rcs(_t, df(1) log orthog event moffset(negage)) ///
		, family(rp, df(3) failure(died))) 	///
		,
		
*these two are the same models, but which command can I use to check this against?

merlin (_t 	trt bmi x1 x2 x3   			///
		trt#rcs(_t, df(1) log orthog event offset(age)) ///
		, family(rp, df(3) failure(died))) 	///
		,
merlin (_t 	trt bmi x1 x2 x3   			///
		trt#rcs(_t, df(1) log orthog event moffset(negage)) ///
		, family(rp, df(3) failure(died))) 	///
		,

		
//not fitting the same model?		
stmt trt bmi x1 x2 x3,  time1(df(3)) time2(tvc(trt) dftvc(1) start(age) df(1)) nohr



//============================================================================//
//throw error messages


*rcof "uhtred (_t newvar trt bmi x1 x2 x3 trt#rcs(_t, df(1) log orthog) 	, family(rp, knots(3) failure(died)))" == 111


//============================================================================//
// end of file
//============================================================================//

