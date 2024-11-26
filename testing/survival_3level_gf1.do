//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

clear all
tr:do ./build/buildmlib.do
mata mata clear

set seed 725
clear 

set obs 100
gen id1 = _n
gen age = runiform()
gen u1 = rnormal(0,0.5)
expand 30
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen u2 = rnormal(0,0.5)
expand 10
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.02 u1 1 u2 1) maxt(10)
stset stime1, f(dead1)

timer clear
// mestreg trt age || id1: || id2: , dist(weib) intmethod(gh)
// gsem (stime1 <- age M2[id1>id2]@1 M1[id1]@1 , family(weib, failure(dead1))) , //intmethod(gh)

timer on 99
// merlin 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , ///
// 		family(rp, df(1) failure(dead1))) ///
// 		, evaltype(gf0) intmethod(mvagh) //
timer off 99			
timer on 3
uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , ///
		family(rp, df(1) failure(dead1))) ///
		, evaltype(gf1) intmethod(mvagh) //
timer off 3
timer on 2
// uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , ///
// 		family(rp, df(1) failure(dead1))) ///
// 		, evaltype(gf1) intmethod(mvagh) //
timer off 2

// timer on 3			
// stuhtred trt age || id1: || id2:, dist(rp) df(1)
// timer off 3
timer list
