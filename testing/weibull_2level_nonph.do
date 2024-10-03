/*============================================================================
program: weibull_2level_nonph
description: testing mixed models with nonPH
created: by Hannah Bower

tests: 	Test 2 level mixed models with random intercept and slope + TDE
		Test 3 level mixed models with random intercept and slope + TDE
		

	
=============================================================================*/


	

//simulate random intercept data for two level model - with TDE
set seed 98798
clear
set obs 10000
gen id 	= _n

local coeftrt -0.5
local coefage 0.1
local coeftde 0.03
local sdref1 0.1

gen trt = runiform()>0.5	
gen sd1 = exp(log(`=`sdref1''))
gen u1 	= rnormal(0,sd1)
gen age = rnormal()

expand 100
sort id 
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt `=`coeftrt'' age `=`coefage'' u1 1) tde(trt `=`coeftde'') tdefunction(log({t}))  maxt(5) 
	
stset stime, f(died)	
	
uhtred (stime trt age c.trt#rcs(stime, df(1) log) ///
		M1[id]@1,  ///
	family(rp, df(1) failure(died))) 			///
	,

//ensure estimates are close to the simulated data
set tracedepth 3
set trace on
assert abs(`=_b[xb1:trt]' - `=`coeftrt'' )< 1E-1
assert abs(`=_b[xb1:age]' - `=`coefage'' )< 1E-1
assert abs(`=_b[tb1:c.trt#c._rcs1_3_2_1]' - `=`coeftde'' )< 1E-1
assert abs(`=exp(_b[lns1_1:_cons])' - `=`sdref1'' )< 1E-1
set trace off
	

		
//simulate two level random intercept and effect data with TDE
//run after bug fixes for random slopes
/*
clear
set obs 20000
gen id1 = _n
expand 5

local coeftrt -0.5
local coefage 0.1
local coeftde 0.03
local sdref1 1
local sdref2 0.5
local corr 0.25

bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,`corr' \ `corr',1)
drawnorm u1 u2, means(0 0) sds(`=`sdref1'' `=`sdref1'') corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (`=`coeftrt''+u2) * trt
gen age = rnormal() + u2


survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(age `=`coefage'' trtui 1 u1 1) ///
	tde(trt `=`coeftde'') tdefunction(log({t})) maxt(5) 


stset stime, f(dead)


//testing 
uhtred (stime trt age c.trt#rcs(stime, df(1) log) ///
		M1[id1]@1 c.trt#M2[id1]@1 ,  ///
	family(rp, df(1) failure(dead))) 			///
	,
	
//ensure estimates are close to the simulated data
assert abs(`=_b[xb1:trt]' - `=`coeftrt'' )< 1E-1
assert abs(`=_b[xb1:age]' - `=`coefage'' )< 1E-1
assert abs(`=_b[tb1:c.trt#c._rcs1_3_2_1]' - `=`coeftde'' )< 1E-1

assert abs(`=exp(_b[lns1_1:_cons]' - `=`sdref1'' )< 1E-1
assert abs(`=exp(_b[lns1_2:_cons]' - `=`sdref2'' )< 1E-1
assert abs(`=tanh(_b[art1_1_2:_cons]' - `=`corr'' )< 1E-1
*/


	
	
	

//simulate three level random intercept models with TDE
clear 
set obs 100
gen id1 = _n

local coeftrt -0.5
local coefage 0.1
local coeftde 0.03
local sdref1 0.1
local sdref2 0.1

gen age = runiform()
gen sd1 = exp(log(0.1))
gen u1 = rnormal(0,`=`sdref1'')
expand 100
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(0.1))
gen u2 = rnormal(0,`=`sdref1'')
expand 10
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt `=`coeftrt'' age `=`coefage'' u1 1 u2 1) ///
	tde(trt `=`coeftde'') tdefunction(log({t})) maxt(10) 
	
stset stime1, f(dead1)

replace id2 = _n

//testing 
uhtred (stime1 trt age c.trt#rcs(stime1, df(1) log)  ///
		M1[id1]@1 M2[id1>id2]@1 ,  ///
	family(rp, df(1) failure(dead))) 			///
	,

assert abs(`=_b[xb1:trt]' - `=`coeftrt'' )< 1E-1
assert abs(`=_b[xb1:age]' - `=`coefage'' )< 1E-1
assert abs(`=_b[tb1:c.trt#c._rcs1_3_2_1]' - `=`coeftde'' )< 1E-1
assert abs(`=exp(_b[lns1_1:_cons])' - `=`sdref1'' )< 1E-1
assert abs(`=exp(_b[lns2_1:_cons])' - `=`sdref2'' )< 1E-1





//simulate three level random intercept & effect models with TDE
//run after bug fixes for random slopes
/*
clear
set obs 20000
gen id1 = _n
expand 5

local coeftrt -0.5
local coefage 0.1
local coeftde 0.03
local sdref1 1
local sdref2 0.5
local corr 0.25

bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,`corr' \ `corr',1)
drawnorm u1 u2, means(0 0) sds(`=`sdref1'' `=`sdref1'') corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (`=`coeftrt''+u2) * trt
gen age = rnormal() + u2
expand 10
gen id3 = _n
sort id1 id2 id3


survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt `=`coeftrt'' age `=`coefage'' u1 1 u2 1) ///
	tde(trt `=`coeftde'') tdefunction(log({t})) maxt(10) 
	
stset stime1, f(dead1)

//testing 
uhtred (stime1 trt age c.trt#rcs(stime1, df(1) log)  ///
		M1[id1]@1 c.trt#M2[id1]@1 M3[id1>id2]@1 ,  ///
	family(rp, df(1) failure(dead))) 			///
	,

*/
	
//============================================================================//
// end of file
//============================================================================//

	