//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software

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
gen u1 = rnormal(0,1)
expand 30
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen u2 = rnormal(0,0.5)
expand 10
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(`=exp(-2.4)') gamma(`=exp(-0.188)') cov(u1 1 u2 1) maxt(10)
stset stime1, f(dead1)

gsem (_t <- M1[id1]@1 M2[id1>id2]@1, family(weib,fail(_d))), 
// mestreg trt age || id1: || id2: , dist(weib) nohr //intmethod(mcagh)
predict refs*, latent ebmeans 
// predict crefs*, reffects ebmodes

gen res1 = .
bys id1 : gen flag = _n==1

uhtred 	(stime1 M2[id1>id2]@1 M1[id1]@1 , ///
		family(rp, df(1) noorthog failure(dead1))) ///
		, restartvalues(M1 0.5 M2 0.25) 

predict urefs*, reffects
//
// su refs* urefs*
