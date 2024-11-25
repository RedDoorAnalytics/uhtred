
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

merlin (stime trt bmi x1 x2 x3, family(rp, df(3) failure(died)))

foreach v in hazard survival chazard logchazard density {
	predict `v'n, `v' ci
	predict `v'n_at, `v' at(trt 1 bmi 20) ci
	predict `v'n_z, `v' zeros ci
}


uhtred (stime 	trt bmi x1 x2 x3 c.trt#rcs(stime, df(1) log orthog) , family(rp, df(3) failure(died))) 
foreach v in hazard survival chazard logchazard density {
	predict `v'nnph, `v' ci
	predict `v'nnph_at, `v' at(trt 1 bmi 20) ci
}

merge 1:1 id1 using ./cert/cert-data/rp_1level_pred, keep(match) nogen

foreach i in hazard survival chazard logchazard density {
	assert abs(`i'- `i'n)< 1E-5
	assert abs(`i'_lci- `i'n_lci)< 1E-5
	assert abs(`i'_uci- `i'n_uci)< 1E-5

	assert abs(`i'nph- `i'nnph)< 1E-5
	assert abs(`i'nph_lci- `i'nnph_lci)< 1E-5
	assert abs(`i'nph_uci- `i'nnph_uci)< 1E-5
	
	assert abs(`i'_z- `i'n_z)< 1E-5
	assert abs(`i'_z_lci- `i'n_z_lci)< 1E-5
	assert abs(`i'_z_uci- `i'n_z_uci)< 1E-5

	
	foreach j in at /*tv z*/ {
		assert abs(`i'_`j'- `i'n_`j')< 1E-5
		assert abs(`i'_`j'_lci- `i'n_`j'_lci)< 1E-5
		assert abs(`i'_`j'_uci- `i'n_`j'_uci)< 1E-5

		assert abs(`i'nph_`j'- `i'nnph_`j')< 1E-5
		assert abs(`i'nph_`j'_lci- `i'nnph_`j'_lci)< 1E-5	
		assert abs(`i'nph_`j'_uci- `i'nnph_`j'_uci)< 1E-5	
	}
}

foreach i in hazard survival chazard logchazard density {
	cap drop `i'*
}

//left truncated data predictions
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
replace t0 = runiform() * 2
drop if stime<t0

uhtred (stime trt bmi age x1, family(rp, df(3) failure(died) ltruncated(t0)))

foreach v in hazard survival chazard logchazard density {
	predict `v'n, `v' ci
	predict `v'n_at, `v' at(trt 1 bmi 20) ci
	predict `v'n_z, `v' zeros ci
}

merge 1:1 id1 using ./cert/cert-data/rp_1level_lt_pred, keep(match) nogen


foreach i in hazard survival chazard logchazard density {
	assert abs(`i'- `i'n)< 1E-5
	assert abs(`i'_lci- `i'n_lci)< 1E-5
	assert abs(`i'_uci- `i'n_uci)< 1E-5

	foreach j in at z /*tv*/ {
		assert abs(`i'_`j'- `i'n_`j')< 1E-5
		assert abs(`i'_`j'_lci- `i'n_`j'_lci)< 1E-5
		assert abs(`i'_`j'_uci- `i'n_`j'_uci)< 1E-5

	}
}

foreach i in hazard survival chazard logchazard density {
	cap drop `i'*
}


