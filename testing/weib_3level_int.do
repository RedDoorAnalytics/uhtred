//source paths
local drive C:\Users\Hannah Bower\Documents\GitHub
cd "`drive'\uhtred"
adopath ++ "`drive'\uhtred"
adopath ++ "`drive'\uhtred\uhtred"

clear all
tr:do .\build\buildmlib.do
mata mata clear

set seed 7254
clear 

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

timer clear
// mestreg trt age || id1: || id2: , dist(weib) //intmethod(gh)
// gsem (stime1 <- age M2[id1>id2]@1 M1[id1]@1 , family(weib, failure(dead1))) , //intmethod(gh)
			
timer on 1
merlin 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , family(rp, df(1) failure(dead1))) ///
			, evaltype(gf0) //intmethod(gh) intpoints(15) 
timer off 1
est store merlin
timer on 2
uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , family(rp, df(1) failure(dead1))) ///
			, evaltype(gf0) //intmethod(gh) //intpoints(15) 
timer off 2
timer list
est store uhtred
// uhtred 	(stime1 age M2[id1>id2]@1 M1[id1]@1 , family(cox, failure(dead1))) ///
// 			, intmethod(gh) intpoints(7) devcode5(294820)
			

//check b and se
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

uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , family(rp, df(1) failure(dead1))) ///
			, evaltype(gf0) //intmethod(gh) //intpoints(15) 

mkassert eclass












