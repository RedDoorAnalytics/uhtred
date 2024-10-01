
clear
set seed 249857
set obs 10000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen agec = age - 50
gen year = 1990 + floor(20*runiform())
gen yearc = year - 2000

survsim stime died, maxtime(5) cov(trt -0.5) 	///
hazard(	0.1:*1.2:*{t}:^0.2 :*			///
        exp(					///
                        0.1 :* (agec :+ {t}) 	///
                        :- 0.1 :* (yearc :+ {t}) ///
                )				///
        )



stset stime, f(died)

//============================================================================//
//RP PH model- simple weibull

