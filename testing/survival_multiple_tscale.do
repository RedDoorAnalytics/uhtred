//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"
clear all

tr:do ./build/buildmlib.do
mata mata clear

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

       
timer clear
timer on 1
merlin (stime 	trt 						///
        rcs(stime, df(1) offset(agec)) 				///
        rcs(stime, df(1) offset(yearc)) 			///
        , family(rp, df(3) failure(died)) timevar(stime)), 
timer off 1

timer on 2
uhtred (stime 	trt 						///
        rcs(stime, df(1) offset(agec)) 				///
        rcs(stime, df(1) offset(yearc)) 			///
        , family(rp, df(3) failure(died)))

timer off 2


       
stset stime, f(died)

timer on 4
stmt trt, 	time1(df(1)) 				///
		time2(df(1) start(agec) logtoff) 	///
		time3(df(1) start(yearc) logtoff) 	///
		nohr noorthog
timer off 4

timer list
