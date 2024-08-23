//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

clear all
tr:do ./build/buildmlib.do
mata mata clear

use "`drive'/multistate/data/multistate_example",clear
rename pid id

msset, id(id) states(rfi osi) times(rf os) cr

gen stime = _stop

bys id: gen d1 = _status[1]
bys id: gen d2 = _status[2]

bys id: drop if _n==2

set seed 35867
// keep if runiform()<0.3

timer clear
timer on 1
merlin 	(stime chemo , family(rp, df(3) failure(d1)) )		///
	(stime chemo , family(rp, df(3) failure(d2)) )		///
	, evaltype(gf2)
timer off 1

timer on 2
uhtred 	(stime chemo , family(rp, df(3) failure(d1)) )		///
	(stime chemo , family(rp, df(3) failure(d2)) )		///
	, evaltype(gf2)
timer off 2

timer list
