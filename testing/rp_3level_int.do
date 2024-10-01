/*============================================================================
program: rp_3level_int
description: testing mixed models with random intercept (3level)
created: by Hannah Bower

tests: 	
	
=============================================================================*/
clear
webuse jobhistory
stset tend, origin(tstart) fail(failure)

sort birthyear id
order birthyear id, first



mestreg education njobs || birthyear: || id:, 	///
	distribution(weib) adaptopts(log) nohr
est store mestreg
	
stmixed education njobs ///
	|| birthyear: ///
	|| id: 	///
	, dist(weib) adaptopts(log)
est store stmixed


merlin 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(weib, failure(_d))) 	

est store merlin

uhtred 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(rp, df(1) failure(_d))) 	


/*============================================================================
end of file
=============================================================================*/
