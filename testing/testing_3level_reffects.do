//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

//build mlib
clear all
tr:do ./build/buildmlib.do
mata mata clear

webuse jobhistory
stset tend, origin(tstart) fail(failure)
sort birthyear id

mestreg education njobs 		///
	|| birthyear: || id:		///
	, distribution(weibull) nohr //intmethod(gh) intpoints(25)

predict m1, median cond(ebmeans)

// stmixed education njobs 		///
// 	|| birthyear: || id:		///
// 	, distribution(weibull)	showmerlin
//	
	
uhtred 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(rp, df(1) failure(_d))) 	///
	, intmethod(mvagh) intpoints(7) 

// predict rfs*, reffects

predict m2, median fitted
