//survsim competing risks
clear 
set seed 4541
set obs 5000
gen trt=runiform()>0.5
gen age = rnormal(50,5)
survsim time state event, hazard1(distribution(weib) lambda(0.1) gamma(0.8) ) ///
	hazard2(distribution(exp) lambda(0.2)  covariates(trt -0.5)) maxtime(10)	

gen id=_n	
expand 2	
bysort id: gen cause=_n		
order id


gen state=state1-1
gen event=(cause==state)  

gen cancer=(cause==1)
gen other=(cause==2)


keep id trt age time0 time1 cancer other cause state event

stset time1, failure(event) 


uhtred 	(time1  age	///
		if cause==1								///
		, family(rp, failure(event) df(1)))		///
		(time1 	trt	///
		if cause==2								///
		, family(rp, failure(event) df(1)))		///
		,

assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap2)'"'              == `"0"'
assert `"`e(ndistap2)'"'          == `"0"'
assert `"`e(constant2)'"'         == `"1"'
assert `"`e(orthog2)'"'           == `"orthog"'
assert `"`e(knots2)'"'            == `"-7.07174 2.302453"'
assert `"`e(timevar2)'"'          == `"time1"'
assert `"`e(failure2)'"'          == `"event"'
assert `"`e(response2)'"'         == `"time1 event"'
assert `"`e(family2)'"'           == `"rp"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-7.175559 2.301569"'
assert `"`e(timevar1)'"'          == `"time1"'
assert `"`e(failure1)'"'          == `"event"'
assert `"`e(response1)'"'         == `"time1 event"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0 0"'
assert `"`e(hastb)'"'             == `"1 1"'
assert `"`e(hasxb)'"'             == `"1 1"'
assert `"`e(allvars)'"'           == `"age event time1 trt"'
assert `"`e(title)'"'             == `"Fixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"0"'
assert `"`e(from)'"'              == `"0"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"uhtred_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf2"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 6
assert         e(N)          == 10000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -13628.62172669903) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.0010769303914529
mat T_b[1,2] =  -1.59358639066376
mat T_b[1,3] =   1.08171467238767
mat T_b[1,4] = -.5148021552991292
mat T_b[1,5] = -.8157532864747798
mat T_b[1,6] =  1.356674188472324
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:age xb1:_cons tb1:_rcs1_1 xb2:trt xb2:_cons tb2:_rcs2_1"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0000275562525869
mat T_V[1,2] =  -.001378649387389
mat T_V[1,3] =  6.67779523831e-07
mat T_V[2,1] =  -.001378649387389
mat T_V[2,2] =  .0699001101165307
mat T_V[2,3] = -.0004006743285855
mat T_V[3,1] =  6.67779523831e-07
mat T_V[3,2] = -.0004006743285855
mat T_V[3,3] =  .0005930327165903
mat T_V[4,4] =  .0013463123026857
mat T_V[4,5] = -.0005635144673324
mat T_V[4,6] = -.0000681530695734
mat T_V[5,4] = -.0005635144673324
mat T_V[5,5] =  .0007545075005375
mat T_V[5,6] = -.0002487074283526
mat T_V[6,4] = -.0000681530695734
mat T_V[6,5] = -.0002487074283526
mat T_V[6,6] =   .000412609602891
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:age xb1:_cons tb1:_rcs1_1 xb2:trt xb2:_cons tb2:_rcs2_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:age xb1:_cons tb1:_rcs1_1 xb2:trt xb2:_cons tb2:_rcs2_1"'
mat drop C_V T_V

qui {
mat T_rcsrmat_2 = J(2,2,0)
mat T_rcsrmat_2[1,1] =  1.346942846937757
mat T_rcsrmat_2[2,1] =  .7779056998937681
mat T_rcsrmat_2[2,2] =                  1
}
matrix C_rcsrmat_2 = e(rcsrmat_2)
assert mreldif( C_rcsrmat_2 , T_rcsrmat_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_2'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_2'"' `"c1 c2"'
mat drop C_rcsrmat_2 T_rcsrmat_2

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  1.346942846937757
mat T_rcsrmat_1[2,1] =  .7779056998937681
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -6.73377730628e-06
mat T_gradient[1,2] = -1.31076996002e-07
mat T_gradient[1,3] = -7.62147954566e-08
mat T_gradient[1,4] = -2.09666307016e-06
mat T_gradient[1,5] = -.0001026950832327
mat T_gradient[1,6] = -.0000582478845717
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:age xb1:_cons tb1:_rcs1_1 xb2:trt xb2:_cons tb2:_rcs2_1"'
mat drop C_gradient T_gradient
