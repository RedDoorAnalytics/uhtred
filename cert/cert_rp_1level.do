
clear 
set seed 725
set obs 5000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
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
assert reldif( e(ll)          , -7490.039919617052) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5236311568622737
mat T_b[1,2] = -.0520184063514714
mat T_b[1,3] =  .0915971770101783
mat T_b[1,4] =  -.397486434972997
mat T_b[1,5] =  .5529681576746154
mat T_b[1,6] =  .5445505321248828
mat T_b[1,7] =  .9253189352630145
mat T_b[1,8] =  .0118160786649674
mat T_b[1,9] =   .001513290640151
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =  .0020552254995932
mat T_V[1,2] = -8.78660638853e-06
mat T_V[1,3] =  2.28081994815e-07
mat T_V[1,4] =  .0000319901288127
mat T_V[1,5] = -.0000325159548646
mat T_V[1,6] = -.0005522161988129
mat T_V[1,7] = -.0000318845477472
mat T_V[1,8] =   .000013826034635
mat T_V[1,9] =  4.11010088246e-06
mat T_V[2,1] = -8.78660638853e-06
mat T_V[2,2] =  .0000542777097081
mat T_V[2,3] =  3.72842728224e-06
mat T_V[2,4] =  4.87169247655e-07
mat T_V[2,5] = -4.36184312899e-06
mat T_V[2,6] = -.0016035311938525
mat T_V[2,7] = -3.18062326894e-06
mat T_V[2,8] =  1.05708383684e-06
mat T_V[2,9] =  2.95105240791e-07
mat T_V[3,1] =  2.28081994815e-07
mat T_V[3,2] =  3.72842728224e-06
mat T_V[3,3] =  .0005064686461931
mat T_V[3,4] =  9.16308877118e-06
mat T_V[3,5] =  8.07081049676e-06
mat T_V[3,6] = -.0001504719036784
mat T_V[3,7] =  5.07301123621e-06
mat T_V[3,8] = -1.90597338507e-06
mat T_V[3,9] = -3.52351126013e-07
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
mat T_V[6,6] =  .0483641007509301
mat T_V[6,7] = -.0000347294569263
mat T_V[6,8] = -.0000549776114735
mat T_V[6,9] = -.0000141940957856
mat T_V[7,1] = -.0000318845477472
mat T_V[7,2] = -3.18062326894e-06
mat T_V[7,3] =  5.07301123621e-06
mat T_V[7,4] = -.0000248081113608
mat T_V[7,5] =  .0000365526016021
mat T_V[7,6] = -.0000347294569263
mat T_V[7,7] =  .0004455021735946
mat T_V[7,8] =  .0001614271630292
mat T_V[7,9] = -1.50797239454e-06
mat T_V[8,1] =   .000013826034635
mat T_V[8,2] =  1.05708383684e-06
mat T_V[8,3] = -1.90597338507e-06
mat T_V[8,4] =  .0000103649379176
mat T_V[8,5] = -.0000155449443249
mat T_V[8,6] = -.0000549776114735
mat T_V[8,7] =  .0001614271630292
mat T_V[8,8] =  .0002934403674718
mat T_V[8,9] = -.0000392463960342
mat T_V[9,1] =  4.11010088246e-06
mat T_V[9,2] =  2.95105240791e-07
mat T_V[9,3] = -3.52351126013e-07
mat T_V[9,4] =  2.78112042326e-06
mat T_V[9,5] = -4.38482716948e-06
mat T_V[9,6] = -.0000141940957856
mat T_V[9,7] = -1.50797239454e-06
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
mat T_rcsrmat_1[1,2] = -12.32939739719497
mat T_rcsrmat_1[1,3] = -5.092327053364387
mat T_rcsrmat_1[2,2] =  2.651955395098098
mat T_rcsrmat_1[2,3] =  1.181226618289067
mat T_rcsrmat_1[3,3] =  .0917647908310305
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
mat T_gradient[1,1] = -1.72065420939e-10
mat T_gradient[1,2] = -1.69325902188e-08
mat T_gradient[1,3] = -1.35318197889e-10
mat T_gradient[1,4] =  4.26769544319e-10
mat T_gradient[1,5] = -4.81534958717e-10
mat T_gradient[1,6] = -5.73222792150e-10
mat T_gradient[1,7] = -4.07976569194e-11
mat T_gradient[1,8] = -1.08515156930e-10
mat T_gradient[1,9] = -2.28649767522e-10
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_gradient T_gradient


		
uhtred (_t 	trt bmi x1 x2 x3 			///
		c.trt#rcs(_t, df(1) log orthog) 	///
		, family(rp, df(1) failure(died))) 	///
                , 

_assert_streq `"`e(cmdline)'"' `"uhtred (_t      trt bmi x1 x2 x3                                        c.trt#rcs(_t, df(1) log orthog)                         , family(rp, df(1) failure(died)))                      ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-4.415933 2.302585"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4.415933 2.302107"'
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
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"uhtred_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf2"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 8
assert         e(N)          == 5000
assert         e(k)          == 8
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -7490.377223570948) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5262461531745911
mat T_b[1,2] = -.0520968439076628
mat T_b[1,3] =  .0917219706118695
mat T_b[1,4] = -.3981852864678235
mat T_b[1,5] =   .553991106313936
mat T_b[1,6] =  .5488564690490514
mat T_b[1,7] =  .0051823869147301
mat T_b[1,8] =   .915671519440337
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0021950862205026
mat T_V[1,2] = -8.61041054867e-06
mat T_V[1,3] =  2.26371995434e-07
mat T_V[1,4] =  .0000288445200233
mat T_V[1,5] = -.0000245458944991
mat T_V[1,6] = -.0006012872415715
mat T_V[1,7] = -.0004499898124209
mat T_V[1,8] =  .0001336554556258
mat T_V[2,1] = -8.61041054867e-06
mat T_V[2,2] =   .000054268675059
mat T_V[2,3] =  3.75685652131e-06
mat T_V[2,4] =  4.08779813611e-07
mat T_V[2,5] = -4.25216205059e-06
mat T_V[2,6] = -.0016031722639078
mat T_V[2,7] = -1.04307923213e-06
mat T_V[2,8] = -3.54883762204e-06
mat T_V[3,1] =  2.26371995434e-07
mat T_V[3,2] =  3.75685652131e-06
mat T_V[3,3] =  .0005064086156895
mat T_V[3,4] =  9.28119109496e-06
mat T_V[3,5] =  7.87532347601e-06
mat T_V[3,6] = -.0001515385921176
mat T_V[3,7] =  8.17644046249e-07
mat T_V[3,8] =  6.06526056261e-06
mat T_V[4,1] =  .0000288445200233
mat T_V[4,2] =  4.08779813611e-07
mat T_V[4,3] =  9.28119109496e-06
mat T_V[4,4] =  .0004990343341869
mat T_V[4,5] = -.0000276765357845
mat T_V[4,6] =  .0001366849371142
mat T_V[4,7] =  6.18663847585e-06
mat T_V[4,8] =     -.000034734704
mat T_V[5,1] = -.0000245458944991
mat T_V[5,2] = -4.25216205059e-06
mat T_V[5,3] =  7.87532347601e-06
mat T_V[5,4] = -.0000276765357845
mat T_V[5,5] =  .0005348867051181
mat T_V[5,6] = -.0000910923606914
mat T_V[5,7] = -.0000202533550802
mat T_V[5,8] =   .000055834101125
mat T_V[6,1] = -.0006012872415715
mat T_V[6,2] = -.0016031722639078
mat T_V[6,3] = -.0001515385921176
mat T_V[6,4] =  .0001366849371142
mat T_V[6,5] = -.0000910923606914
mat T_V[6,6] =  .0483624226690458
mat T_V[6,7] =  .0001778641950534
mat T_V[6,8] = -.0000644363983593
mat T_V[7,1] = -.0004499898124209
mat T_V[7,2] = -1.04307923213e-06
mat T_V[7,3] =  8.17644046249e-07
mat T_V[7,4] =  6.18663847585e-06
mat T_V[7,5] = -.0000202533550802
mat T_V[7,6] =  .0001778641950534
mat T_V[7,7] =  .0014338671059938
mat T_V[7,8] = -.0005611128188747
mat T_V[8,1] =  .0001336554556258
mat T_V[8,2] = -3.54883762204e-06
mat T_V[8,3] =  6.06526056261e-06
mat T_V[8,4] =     -.000034734704
mat T_V[8,5] =   .000055834101125
mat T_V[8,6] = -.0000644363983593
mat T_V[8,7] = -.0005611128188747
mat T_V[8,8] =  .0005668676059698
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1"'
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
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] = -1.65007243980e-10
mat T_gradient[1,2] = -1.65462937796e-08
mat T_gradient[1,3] = -1.33790760263e-10
mat T_gradient[1,4] =  4.19194233280e-10
mat T_gradient[1,5] = -4.70951295723e-10
mat T_gradient[1,6] = -5.60642386682e-10
mat T_gradient[1,7] = -4.21006077223e-11
mat T_gradient[1,8] =  2.41453246286e-11
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_rcs1_1"'
mat drop C_gradient T_gradient
