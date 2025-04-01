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
//
// mestreg education njobs 		///
// 	|| birthyear: || id:		///
// 	, distribution(weibull) nohr intpoints(25)
//
// predict refs*, reffects

// stmixed education njobs 		///
// 	|| birthyear: || id:		///
// 	, distribution(weibull)	showmerlin
//	
// sort birthyear id
	
uhtred 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(rp, noorthog df(1) failure(_d))) 	///
	, 

predict rfs*, reffects 
su rfs*
