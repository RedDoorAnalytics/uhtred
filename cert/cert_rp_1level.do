
clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen bmi = rnormal(30,3)

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) //ltruncated(t0)

stset stime, f(died) //enter(t0)

//============================================================================//
//RP PH model- simple weibull

uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, df(1) failure(died))) 	///
                , 


_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3                                        , family(rp, df(1) failure(died)))                      ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4.415933 2.302107"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"_t died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"_t bmi died trt x1 x2 x3"'
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

assert         e(rank)       == 7
assert         e(N)          == 5000
assert         e(k)          == 7
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -7490.386596280393) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,7,0)
mat T_b[1,1] = -.5246209510357646
mat T_b[1,2] = -.0520931035801331
mat T_b[1,3] =  .0917190702873983
mat T_b[1,4] = -.3982079775890079
mat T_b[1,5] =  .5540648299391859
mat T_b[1,6] =  .5482119599812287
mat T_b[1,7] =  .9177061878953923
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(7,7,0)
mat T_V[1,1] =  .0020538130399537
mat T_V[1,2] = -8.91247804861e-06
mat T_V[1,3] =  4.43368367942e-07
mat T_V[1,4] =  .0000309741980021
mat T_V[1,5] = -.0000311667945861
mat T_V[1,6] = -.0005461491671776
mat T_V[1,7] = -.0000420438405964
mat T_V[2,1] = -8.91247804861e-06
mat T_V[2,2] =  .0000542688440129
mat T_V[2,3] =  3.75603585991e-06
mat T_V[2,4] =  4.12136647695e-07
mat T_V[2,5] = -4.26062483341e-06
mat T_V[2,6] = -.0016030827335666
mat T_V[2,7] = -3.95790696029e-06
mat T_V[3,1] =  4.43368367942e-07
mat T_V[3,2] =  3.75603585991e-06
mat T_V[3,3] =  .0005064012913387
mat T_V[3,4] =  9.28344046167e-06
mat T_V[3,5] =  7.88227885011e-06
mat T_V[3,6] = -.0001515780992091
mat T_V[3,7] =  6.38379606266e-06
mat T_V[4,1] =  .0000309741980021
mat T_V[4,2] =  4.12136647695e-07
mat T_V[4,3] =  9.28344046167e-06
mat T_V[4,4] =  .0004989955270776
mat T_V[4,5] = -.0000275860599974
mat T_V[4,6] =  .0001358698858308
mat T_V[4,7] = -.0000323049650197
mat T_V[5,1] = -.0000311667945861
mat T_V[5,2] = -4.26062483341e-06
mat T_V[5,3] =  7.88227885011e-06
mat T_V[5,4] = -.0000275860599974
mat T_V[5,5] =  .0005346231648273
mat T_V[5,6] =   -.00008866497689
mat T_V[5,7] =  .0000478925567356
mat T_V[6,1] = -.0005461491671776
mat T_V[6,2] = -.0016030827335666
mat T_V[6,3] = -.0001515780992091
mat T_V[6,4] =  .0001358698858308
mat T_V[6,5] =   -.00008866497689
mat T_V[6,6] =  .0483418807390042
mat T_V[6,7] =  5.01111915374e-06
mat T_V[7,1] = -.0000420438405964
mat T_V[7,2] = -3.95790696029e-06
mat T_V[7,3] =  6.38379606266e-06
mat T_V[7,4] = -.0000323049650197
mat T_V[7,5] =  .0000478925567356
mat T_V[7,6] =  5.01111915374e-06
mat T_V[7,7] =  .0003472640599698
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7971729577304844
mat T_rcsrmat_1[2,1] =    1.8714958263367
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,7,0)
mat T_gradient[1,1] = -1.79367982273e-10
mat T_gradient[1,2] = -1.77666105161e-08
mat T_gradient[1,3] = -1.42126899909e-10
mat T_gradient[1,4] =  4.51078270222e-10
mat T_gradient[1,5] = -5.11701273523e-10
mat T_gradient[1,6] = -6.01650600868e-10
mat T_gradient[1,7] = -5.99408092605e-11
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1"'
mat drop C_gradient T_gradient



//============================================================================//
//RP PH model-3 df

uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, df(3) failure(died))) 	///
                , 




_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3                                        , family(rp, df(3) failure(died)))                      ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4.415933 1.043513 1.805705 2.302107"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"_t died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"_t bmi died trt x1 x2 x3"'
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

assert         e(rank)       == 9
assert         e(N)          == 5000
assert         e(k)          == 9
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -7490.039919617051) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5236311568622735
mat T_b[1,2] = -.0520184063514713
mat T_b[1,3] =  .0915971770101783
mat T_b[1,4] = -.3974864349729968
mat T_b[1,5] =  .5529681576746152
mat T_b[1,6] =  .5445505321248796
mat T_b[1,7] =  .9253189352630145
mat T_b[1,8] =  .0118160786649667
mat T_b[1,9] =  .0015132906401527
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =  .0020552254995932
mat T_V[1,2] = -8.78660638854e-06
mat T_V[1,3] =  2.28081994814e-07
mat T_V[1,4] =  .0000319901288127
mat T_V[1,5] = -.0000325159548646
mat T_V[1,6] = -.0005522161988129
mat T_V[1,7] = -.0000318845477472
mat T_V[1,8] =   .000013826034635
mat T_V[1,9] =  4.11010088246e-06
mat T_V[2,1] = -8.78660638854e-06
mat T_V[2,2] =  .0000542777097081
mat T_V[2,3] =  3.72842728224e-06
mat T_V[2,4] =  4.87169247655e-07
mat T_V[2,5] = -4.36184312899e-06
mat T_V[2,6] = -.0016035311938525
mat T_V[2,7] = -3.18062326894e-06
mat T_V[2,8] =  1.05708383684e-06
mat T_V[2,9] =  2.95105240791e-07
mat T_V[3,1] =  2.28081994814e-07
mat T_V[3,2] =  3.72842728224e-06
mat T_V[3,3] =  .0005064686461931
mat T_V[3,4] =  9.16308877118e-06
mat T_V[3,5] =  8.07081049676e-06
mat T_V[3,6] = -.0001504719036784
mat T_V[3,7] =  5.07301123621e-06
mat T_V[3,8] = -1.90597338507e-06
mat T_V[3,9] = -3.52351126012e-07
mat T_V[4,1] =  .0000319901288127
mat T_V[4,2] =  4.87169247655e-07
mat T_V[4,3] =  9.16308877118e-06
mat T_V[4,4] =  .0004996418296812
mat T_V[4,5] =  -.000028620720164
mat T_V[4,6] =  .0001319462743491
mat T_V[4,7] = -.0000248081113608
mat T_V[4,8] =  .0000103649379176
mat T_V[4,9] =  2.78112042326e-06
mat T_V[5,1] = -.0000325159548646
mat T_V[5,2] = -4.36184312899e-06
mat T_V[5,3] =  8.07081049676e-06
mat T_V[5,4] =  -.000028620720164
mat T_V[5,5] =  .0005362114897186
mat T_V[5,6] = -.0000831742852139
mat T_V[5,7] =  .0000365526016021
mat T_V[5,8] = -.0000155449443249
mat T_V[5,9] = -4.38482716948e-06
mat T_V[6,1] = -.0005522161988129
mat T_V[6,2] = -.0016035311938525
mat T_V[6,3] = -.0001504719036784
mat T_V[6,4] =  .0001319462743491
mat T_V[6,5] = -.0000831742852139
mat T_V[6,6] =  .0483641007509306
mat T_V[6,7] = -.0000347294569263
mat T_V[6,8] = -.0000549776114735
mat T_V[6,9] = -.0000141940957856
mat T_V[7,1] = -.0000318845477472
mat T_V[7,2] = -3.18062326894e-06
mat T_V[7,3] =  5.07301123621e-06
mat T_V[7,4] = -.0000248081113608
mat T_V[7,5] =  .0000365526016021
mat T_V[7,6] = -.0000347294569263
mat T_V[7,7] =  .0004455021735947
mat T_V[7,8] =  .0001614271630293
mat T_V[7,9] = -1.50797239462e-06
mat T_V[8,1] =   .000013826034635
mat T_V[8,2] =  1.05708383684e-06
mat T_V[8,3] = -1.90597338507e-06
mat T_V[8,4] =  .0000103649379176
mat T_V[8,5] = -.0000155449443249
mat T_V[8,6] = -.0000549776114735
mat T_V[8,7] =  .0001614271630293
mat T_V[8,8] =  .0002934403674717
mat T_V[8,9] = -.0000392463960342
mat T_V[9,1] =  4.11010088246e-06
mat T_V[9,2] =  2.95105240791e-07
mat T_V[9,3] = -3.52351126012e-07
mat T_V[9,4] =  2.78112042326e-06
mat T_V[9,5] = -4.38482716948e-06
mat T_V[9,6] = -.0000141940957856
mat T_V[9,7] = -1.50797239462e-06
mat T_V[9,8] = -.0000392463960342
mat T_V[9,9] =  .0000507009984826
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .7971729577304844
mat T_rcsrmat_1[1,2] = -12.32939739719427
mat T_rcsrmat_1[1,3] = -5.092327053364209
mat T_rcsrmat_1[2,2] =  2.651955395098022
mat T_rcsrmat_1[2,3] =  1.181226618289126
mat T_rcsrmat_1[3,3] =  .0917647908310291
mat T_rcsrmat_1[4,1] =    1.8714958263367
mat T_rcsrmat_1[4,2] = -47.20905966714736
mat T_rcsrmat_1[4,3] = -19.07834032616201
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] = -1.72233595441e-10
mat T_gradient[1,2] = -1.69431210173e-08
mat T_gradient[1,3] = -1.35350944489e-10
mat T_gradient[1,4] =  4.26876293482e-10
mat T_gradient[1,5] = -4.81572941707e-10
mat T_gradient[1,6] = -5.73569754192e-10
mat T_gradient[1,7] = -4.11831541092e-11
mat T_gradient[1,8] = -1.08123173609e-10
mat T_gradient[1,9] = -2.28014176984e-10
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_gradient T_gradient


//============================================================================//
//RP PH model factor

uhtred (_t	trt bmi x1 x2 x3 i.agecat			///
		, family(rp, df(1) failure(died))) 	///
		,

_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3 i.agecat                                       , family(rp, df(1) failure(died)))                      ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4.415933 2.302107"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"_t died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"_t bmi died i.agecat trt x1 x2 x3"'
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

assert         e(rank)       == 11
assert         e(N)          == 5000
assert         e(k)          == 12
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 1
assert reldif( e(ll)          , -7487.502134773174) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,12,0)
mat T_b[1,1] = -.5267118791024842
mat T_b[1,2] = -.0519575049346563
mat T_b[1,3] =  .0924302483881762
mat T_b[1,4] = -.3990042160316856
mat T_b[1,5] =   .551867942946791
mat T_b[1,7] =  .0563507562402874
mat T_b[1,8] =  .0743899195269524
mat T_b[1,9] =  .1849886808304233
mat T_b[1,10] =  .1248491870843496
mat T_b[1,11] =  .4716397651316285
mat T_b[1,12] =  .9181368092280155
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_rcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(12,12,0)
mat T_V[1,1] =  .0020555196323709
mat T_V[1,2] = -8.82479297774e-06
mat T_V[1,3] = -1.98546398925e-06
mat T_V[1,4] =  .0000292154216045
mat T_V[1,5] = -.0000316408783072
mat T_V[1,7] =  .0000278024619091
mat T_V[1,8] =   .000026584028736
mat T_V[1,9] = -.0000492257045225
mat T_V[1,10] = -.0000662813623974
mat T_V[1,11] =  -.000558566468899
mat T_V[1,12] = -.0000423740013561
mat T_V[2,1] = -8.82479297774e-06
mat T_V[2,2] =  .0000543578445681
mat T_V[2,3] =  3.70468266221e-06
mat T_V[2,4] =  5.62544265431e-07
mat T_V[2,5] = -4.14462141876e-06
mat T_V[2,7] = -2.46378657495e-06
mat T_V[2,8] =  5.94106395807e-06
mat T_V[2,9] =  3.25311604916e-06
mat T_V[2,10] = -1.56555459171e-06
mat T_V[2,11] = -.0016074151921426
mat T_V[2,12] = -3.95547726036e-06
mat T_V[3,1] = -1.98546398925e-06
mat T_V[3,2] =  3.70468266221e-06
mat T_V[3,3] =  .0005082805078332
mat T_V[3,4] =  9.20676493105e-06
mat T_V[3,5] =  7.17508783854e-06
mat T_V[3,7] =  2.88747296676e-06
mat T_V[3,8] = -.0000377185676147
mat T_V[3,9] =  .0000200021033231
mat T_V[3,10] =  .0000478420560563
mat T_V[3,11] = -.0001414117187263
mat T_V[3,12] =  6.56038385383e-06
mat T_V[4,1] =  .0000292154216045
mat T_V[4,2] =  5.62544265431e-07
mat T_V[4,3] =  9.20676493105e-06
mat T_V[4,4] =  .0004999597905542
mat T_V[4,5] = -.0000266589132079
mat T_V[4,7] = -.0000288158124338
mat T_V[4,8] = -.0000105366467145
mat T_V[4,9] =  -.000032382155629
mat T_V[4,10] = -.0000309881673088
mat T_V[4,11] =   .000150904107812
mat T_V[4,12] = -.0000324311913663
mat T_V[5,1] = -.0000316408783072
mat T_V[5,2] = -4.14462141876e-06
mat T_V[5,3] =  7.17508783854e-06
mat T_V[5,4] = -.0000266589132079
mat T_V[5,5] =  .0005350260315053
mat T_V[5,7] = -.0000359389474495
mat T_V[5,8] = -.0000168084180521
mat T_V[5,9] = -.0000817744963088
mat T_V[5,10] =  .0000317193740995
mat T_V[5,11] = -.0000621259336114
mat T_V[5,12] =  .0000476531022021
mat T_V[7,1] =  .0000278024619091
mat T_V[7,2] = -2.46378657495e-06
mat T_V[7,3] =  2.88747296676e-06
mat T_V[7,4] = -.0000288158124338
mat T_V[7,5] = -.0000359389474495
mat T_V[7,7] =  .0047232880814423
mat T_V[7,8] =  .0032275775784718
mat T_V[7,9] =  .0032332525805377
mat T_V[7,10] =  .0032252306719719
mat T_V[7,11] = -.0031635641144773
mat T_V[7,12] =  4.96877022447e-06
mat T_V[8,1] =   .000026584028736
mat T_V[8,2] =  5.94106395807e-06
mat T_V[8,3] = -.0000377185676147
mat T_V[8,4] = -.0000105366467145
mat T_V[8,5] = -.0000168084180521
mat T_V[8,7] =  .0032275775784718
mat T_V[8,8] =  .0046969204185773
mat T_V[8,9] =  .0032277156002027
mat T_V[8,10] =  .0032211168725273
mat T_V[8,11] = -.0034112138526474
mat T_V[8,12] =  7.45454322637e-06
mat T_V[9,1] = -.0000492257045225
mat T_V[9,2] =  3.25311604916e-06
mat T_V[9,3] =  .0000200021033231
mat T_V[9,4] =  -.000032382155629
mat T_V[9,5] = -.0000817744963088
mat T_V[9,7] =  .0032332525805377
mat T_V[9,8] =  .0032277156002027
mat T_V[9,9] =  .0064391517780374
mat T_V[9,10] =  .0032269924374884
mat T_V[9,11] = -.0032913939320573
mat T_V[9,12] =  .0000139664301752
mat T_V[10,1] = -.0000662813623974
mat T_V[10,2] = -1.56555459171e-06
mat T_V[10,3] =  .0000478420560563
mat T_V[10,4] = -.0000309881673088
mat T_V[10,5] =  .0000317193740995
mat T_V[10,7] =  .0032252306719719
mat T_V[10,8] =  .0032211168725273
mat T_V[10,9] =  .0032269924374884
mat T_V[10,10] =  .0249751691509285
mat T_V[10,11] = -.0031790126634683
mat T_V[10,12] =  .0000125826299809
mat T_V[11,1] =  -.000558566468899
mat T_V[11,2] = -.0016074151921426
mat T_V[11,3] = -.0001414117187263
mat T_V[11,4] =   .000150904107812
mat T_V[11,5] = -.0000621259336114
mat T_V[11,7] = -.0031635641144773
mat T_V[11,8] = -.0034112138526474
mat T_V[11,9] = -.0032913939320573
mat T_V[11,10] = -.0031790126634683
mat T_V[11,11] =   .051252669361297
mat T_V[11,12] = -1.39060223768e-06
mat T_V[12,1] = -.0000423740013561
mat T_V[12,2] = -3.95547726036e-06
mat T_V[12,3] =  6.56038385383e-06
mat T_V[12,4] = -.0000324311913663
mat T_V[12,5] =  .0000476531022021
mat T_V[12,7] =  4.96877022447e-06
mat T_V[12,8] =  7.45454322637e-06
mat T_V[12,9] =  .0000139664301752
mat T_V[12,10] =  .0000125826299809
mat T_V[12,11] = -1.39060223768e-06
mat T_V[12,12] =  .0003475158456889
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_rcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_rcs1_1"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7971729577304844
mat T_rcsrmat_1[2,1] =    1.8714958263367
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,12,0)
mat T_gradient[1,1] = -1.99722446914e-10
mat T_gradient[1,2] = -1.92159709256e-08
mat T_gradient[1,3] = -1.49191985785e-10
mat T_gradient[1,4] =  4.78327614957e-10
mat T_gradient[1,5] = -5.53774551368e-10
mat T_gradient[1,7] = -1.55000373647e-10
mat T_gradient[1,8] = -2.77938914939e-10
mat T_gradient[1,9] = -1.45493943282e-10
mat T_gradient[1,10] = -6.25630103279e-12
mat T_gradient[1,11] = -6.50224027399e-10
mat T_gradient[1,12] = -6.70012031967e-11
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_rcs1_1"'
mat drop C_gradient T_gradient




//============================================================================//
//RP nonPH model


//mkassert
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(3) failure(died))) 	///
		,
		
_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3                                        c.trt#rcs(_t, df(1) log orthog)                         , family(rp, df(3) failure(died)))                      ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-4.415933 2.302585"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4.415933 1.043513 1.805705 2.302107"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"_t died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"_t bmi c.trt died trt x1 x2 x3"'
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

assert         e(rank)       == 10
assert         e(N)          == 5000
assert         e(k)          == 10
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -7490.02474791766 ) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,10,0)
mat T_b[1,1] = -.5256870771463352
mat T_b[1,2] = -.0520225080665863
mat T_b[1,3] =  .0915997016900165
mat T_b[1,4] =   -.39745034178866
mat T_b[1,5] =   .552862792345569
mat T_b[1,6] =  .5453353726723186
mat T_b[1,7] =  .0066054127474835
mat T_b[1,8] =  .9227870207174348
mat T_b[1,9] =  .0118956440620156
mat T_b[1,10] =  .0015506229453361
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(10,10,0)
mat T_V[1,1] =  .0021949998597706
mat T_V[1,2] = -8.53664978969e-06
mat T_V[1,3] =  9.97851571885e-08
mat T_V[1,4] =  .0000293401959497
mat T_V[1,5] = -.0000250870616693
mat T_V[1,6] = -.0006044725068366
mat T_V[1,7] = -.0004484120974972
mat T_V[1,8] =  .0001388242651215
mat T_V[1,9] =  8.53717148454e-06
mat T_V[1,10] =  1.56222518637e-06
mat T_V[2,1] = -8.53664978969e-06
mat T_V[2,2] =   .000054277115965
mat T_V[2,3] =  3.72969291149e-06
mat T_V[2,4] =  4.83605494585e-07
mat T_V[2,5] = -4.35523365089e-06
mat T_V[2,6] = -.0016035882014122
mat T_V[2,7] = -9.04349776858e-07
mat T_V[2,8] = -2.83555132823e-06
mat T_V[2,9] =  1.04407886506e-06
mat T_V[2,10] =  2.89990496410e-07
mat T_V[3,1] =  9.97851571885e-08
mat T_V[3,2] =  3.72969291149e-06
mat T_V[3,3] =  .0005064782536227
mat T_V[3,4] =  9.15902995379e-06
mat T_V[3,5] =  8.06777614065e-06
mat T_V[3,6] = -.0001504775972805
mat T_V[3,7] =  5.69588086947e-07
mat T_V[3,8] =  4.85756396234e-06
mat T_V[3,9] = -1.89985123852e-06
mat T_V[3,10] = -3.48232394658e-07
mat T_V[4,1] =  .0000293401959497
mat T_V[4,2] =  4.83605494585e-07
mat T_V[4,3] =  9.15902995379e-06
mat T_V[4,4] =   .000499697144827
mat T_V[4,5] = -.0000287456950827
mat T_V[4,6] =  .0001329269883443
mat T_V[4,7] =  7.74853893081e-06
mat T_V[4,8] = -.0000277788308492
mat T_V[4,9] =  .0000104584575841
mat T_V[4,10] =  2.82182152057e-06
mat T_V[5,1] = -.0000250870616693
mat T_V[5,2] = -4.35523365089e-06
mat T_V[5,3] =  8.06777614065e-06
mat T_V[5,4] = -.0000287456950827
mat T_V[5,5] =  .0005365414545369
mat T_V[5,6] = -.0000857688049006
mat T_V[5,7] = -.0000227769152717
mat T_V[5,8] =  .0000452645274177
mat T_V[5,9] = -.0000158278591408
mat T_V[5,10] = -4.50926518935e-06
mat T_V[6,1] = -.0006044725068366
mat T_V[6,2] = -.0016035882014122
mat T_V[6,3] = -.0001504775972805
mat T_V[6,4] =  .0001329269883443
mat T_V[6,5] = -.0000857688049006
mat T_V[6,6] =  .0483823771020776
mat T_V[6,7] =  .0001704356356081
mat T_V[6,8] = -.0000995725012408
mat T_V[6,9] = -.0000529124102401
mat T_V[6,10] = -.0000132251966369
mat T_V[7,1] = -.0004484120974972
mat T_V[7,2] = -9.04349776858e-07
mat T_V[7,3] =  5.69588086947e-07
mat T_V[7,4] =  7.74853893081e-06
mat T_V[7,5] = -.0000227769152717
mat T_V[7,6] =  .0001704356356081
mat T_V[7,7] =  .0014393827344227
mat T_V[7,8] = -.0005493611381038
mat T_V[7,9] =  .0000173710276111
mat T_V[7,10] =  8.13871452146e-06
mat T_V[8,1] =  .0001388242651215
mat T_V[8,2] = -2.83555132823e-06
mat T_V[8,3] =  4.85756396234e-06
mat T_V[8,4] = -.0000277788308492
mat T_V[8,5] =  .0000452645274177
mat T_V[8,6] = -.0000995725012408
mat T_V[8,7] = -.0005493611381038
mat T_V[8,8] =  .0006551960075844
mat T_V[8,9] =  .0001547951297958
mat T_V[8,10] = -4.61522374494e-06
mat T_V[9,1] =  8.53717148454e-06
mat T_V[9,2] =  1.04407886506e-06
mat T_V[9,3] = -1.89985123852e-06
mat T_V[9,4] =  .0000104584575841
mat T_V[9,5] = -.0000158278591408
mat T_V[9,6] = -.0000529124102401
mat T_V[9,7] =  .0000173710276111
mat T_V[9,8] =  .0001547951297958
mat T_V[9,9] =  .0002936884486765
mat T_V[9,10] =  -.000039159066631
mat T_V[10,1] =  1.56222518637e-06
mat T_V[10,2] =  2.89990496410e-07
mat T_V[10,3] = -3.48232394658e-07
mat T_V[10,4] =  2.82182152057e-06
mat T_V[10,5] = -4.50926518935e-06
mat T_V[10,6] = -.0000132251966369
mat T_V[10,7] =  8.13871452146e-06
mat T_V[10,8] = -4.61522374494e-06
mat T_V[10,9] =  -.000039159066631
mat T_V[10,10] =  .0000507610923185
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_V T_V

qui {
mat T_rmat_1_6_2 = J(2,2,0)
mat T_rmat_1_6_2[1,1] =  .7971729577304844
mat T_rmat_1_6_2[2,1] =    1.8714958263367
mat T_rmat_1_6_2[2,2] =                  1
}
matrix C_rmat_1_6_2 = e(rmat_1_6_2)
assert mreldif( C_rmat_1_6_2 , T_rmat_1_6_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_6_2'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rmat_1_6_2'"' `"c1 c2"'
mat drop C_rmat_1_6_2 T_rmat_1_6_2

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .7971729577304844
mat T_rcsrmat_1[1,2] = -12.32939739719427
mat T_rcsrmat_1[1,3] = -5.092327053364209
mat T_rcsrmat_1[2,2] =  2.651955395098022
mat T_rcsrmat_1[2,3] =  1.181226618289126
mat T_rcsrmat_1[3,3] =  .0917647908310291
mat T_rcsrmat_1[4,1] =    1.8714958263367
mat T_rcsrmat_1[4,2] = -47.20905966714736
mat T_rcsrmat_1[4,3] = -19.07834032616201
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,10,0)
mat T_gradient[1,1] = -1.56476807206e-10
mat T_gradient[1,2] = -1.57146944346e-08
mat T_gradient[1,3] = -1.26829418393e-10
mat T_gradient[1,4] =  3.93535033743e-10
mat T_gradient[1,5] = -4.39309819127e-10
mat T_gradient[1,6] = -5.32234759737e-10
mat T_gradient[1,7] = -3.42032097089e-11
mat T_gradient[1,8] =  5.39703420843e-11
mat T_gradient[1,9] = -1.82828264680e-10
mat T_gradient[1,10] = -2.91886520162e-10
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_gradient T_gradient

//============================================================================//
//user-defined knots (PH RP)

uhtred (_t 	trt bmi x1 x2 x3 			///
		, family(rp, knots(-4 1 1.5 2) failure(died))) 	///
,


_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3                                        , family(rp, knots(-4 1 1.5 2) failure(died)))  ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4 1 1.5 2"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"_t died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"_t bmi died trt x1 x2 x3"'
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

assert         e(rank)       == 9
assert         e(N)          == 5000
assert         e(k)          == 9
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -7489.57641500612 ) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5232515334671862
mat T_b[1,2] = -.0519960287671335
mat T_b[1,3] =  .0915590243850024
mat T_b[1,4] = -.3971858141172385
mat T_b[1,5] =  .5525142126090495
mat T_b[1,6] =  .5433832146558263
mat T_b[1,7] =  .9253127291277446
mat T_b[1,8] =  .0077929225644105
mat T_b[1,9] =  .0077803938494558
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =  .0020550265685279
mat T_V[1,2] = -8.79697383885e-06
mat T_V[1,3] =  2.00763762245e-07
mat T_V[1,4] =  .0000319118275019
mat T_V[1,5] = -.0000322775221337
mat T_V[1,6] = -.0005517926555807
mat T_V[1,7] = -.0000320781053172
mat T_V[1,8] =  .0000135254446508
mat T_V[1,9] =  4.10490333887e-06
mat T_V[2,1] = -8.79697383885e-06
mat T_V[2,2] =  .0000542815027505
mat T_V[2,3] =  3.71759518641e-06
mat T_V[2,4] =  4.76351377250e-07
mat T_V[2,5] = -4.33918586896e-06
mat T_V[2,6] = -.0016036355508163
mat T_V[2,7] = -3.19359786747e-06
mat T_V[2,8] =  1.05877102614e-06
mat T_V[2,9] =  2.52700172253e-07
mat T_V[3,1] =  2.00763762245e-07
mat T_V[3,2] =  3.71759518641e-06
mat T_V[3,3] =  .0005065083124033
mat T_V[3,4] =  9.19381493957e-06
mat T_V[3,5] =  8.08518250859e-06
mat T_V[3,6] = -.0001501521204417
mat T_V[3,7] =  5.09956250236e-06
mat T_V[3,8] = -1.84708634841e-06
mat T_V[3,9] = -3.87720546739e-07
mat T_V[4,1] =  .0000319118275019
mat T_V[4,2] =  4.76351377250e-07
mat T_V[4,3] =  9.19381493957e-06
mat T_V[4,4] =  .0004995137205075
mat T_V[4,5] = -.0000284777304384
mat T_V[4,6] =   .000132266500082
mat T_V[4,7] = -.0000249195703296
mat T_V[4,8] =  9.95785599210e-06
mat T_V[4,9] =  3.10394407821e-06
mat T_V[5,1] = -.0000322775221337
mat T_V[5,2] = -4.33918586896e-06
mat T_V[5,3] =  8.08518250859e-06
mat T_V[5,4] = -.0000284777304384
mat T_V[5,5] =  .0005360154471145
mat T_V[5,6] = -.0000839412798216
mat T_V[5,7] =   .000036736115735
mat T_V[5,8] = -.0000149989169264
mat T_V[5,9] = -4.74818698159e-06
mat T_V[6,1] = -.0005517926555807
mat T_V[6,2] = -.0016036355508163
mat T_V[6,3] = -.0001501521204417
mat T_V[6,4] =   .000132266500082
mat T_V[6,5] = -.0000839412798216
mat T_V[6,6] =  .0483668898148939
mat T_V[6,7] = -.0000340131539025
mat T_V[6,8] = -.0000536965659333
mat T_V[6,9] =  -.000014314515558
mat T_V[7,1] = -.0000320781053172
mat T_V[7,2] = -3.19359786747e-06
mat T_V[7,3] =  5.09956250236e-06
mat T_V[7,4] = -.0000249195703296
mat T_V[7,5] =   .000036736115735
mat T_V[7,6] = -.0000340131539025
mat T_V[7,7] =  .0004424868684348
mat T_V[7,8] =   .000154932623717
mat T_V[7,9] =  4.60097201699e-06
mat T_V[8,1] =  .0000135254446508
mat T_V[8,2] =  1.05877102614e-06
mat T_V[8,3] = -1.84708634841e-06
mat T_V[8,4] =  9.95785599210e-06
mat T_V[8,5] = -.0000149989169264
mat T_V[8,6] = -.0000536965659333
mat T_V[8,7] =   .000154932623717
mat T_V[8,8] =  .0002888657661029
mat T_V[8,9] = -.0000382904452303
mat T_V[9,1] =  4.10490333887e-06
mat T_V[9,2] =  2.52700172253e-07
mat T_V[9,3] = -3.87720546739e-07
mat T_V[9,4] =  3.10394407821e-06
mat T_V[9,5] = -4.74818698159e-06
mat T_V[9,6] =  -.000014314515558
mat T_V[9,7] =  4.60097201699e-06
mat T_V[9,8] = -.0000382904452303
mat T_V[9,9] =  .0000602153371009
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .7971729577304844
mat T_rcsrmat_1[1,2] = -9.148641337221349
mat T_rcsrmat_1[1,3] = -4.782100700190186
mat T_rcsrmat_1[2,2] =  2.034361945785094
mat T_rcsrmat_1[2,3] =  1.133127065236349
mat T_rcsrmat_1[3,3] =  .0704895934254864
mat T_rcsrmat_1[4,1] =    1.8714958263367
mat T_rcsrmat_1[4,2] = -33.88946390180183
mat T_rcsrmat_1[4,3] = -17.37104851006879
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] = -1.71896327034e-10
mat T_gradient[1,2] = -1.68999721994e-08
mat T_gradient[1,3] = -1.35325339596e-10
mat T_gradient[1,4] =  4.25705319411e-10
mat T_gradient[1,5] = -4.80654573624e-10
mat T_gradient[1,6] = -5.72150198747e-10
mat T_gradient[1,7] = -1.04392917921e-11
mat T_gradient[1,8] = -1.40471972596e-10
mat T_gradient[1,9] = -2.68727403177e-10
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_gradient T_gradient

