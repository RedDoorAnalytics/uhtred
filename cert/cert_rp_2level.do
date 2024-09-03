set seed 98798
clear
set obs 1000
gen id 	= _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.1))
gen u1 	= rnormal(0,sd1)
gen age = rnormal()
expand 10
sort id 
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.1 u1 1) maxt(5) 
stset stime, f(died)

uhtred (stime trt age M1[id]@1, family(rp, df(1) failure(died)))	///
	, 
	
_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id]@1, family(rp, df(1) failure(died)))                ,"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd(M1)"'
assert `"`e(re_ivscale1)'"'       == `"exp"'
assert `"`e(re_eqns1)'"'          == `"lns1_1"'
assert `"`e(Nreparams1)'"'        == `"1"'
assert `"`e(Nres1)'"'             == `"1"'
assert `"`e(latents1)'"'          == `"M1"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-5.24262 1.609222"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age died id stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"0"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"uhtred_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 5
assert         e(N)          == 10000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -13338.61300872465) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] =  -.509013901253118
mat T_b[1,2] =  .1075778917214951
mat T_b[1,3] = -.8548997014008141
mat T_b[1,4] =  .8918369629326414
mat T_b[1,5] =                  1
mat T_b[1,6] = -2.398681885782334
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0010224654131747
mat T_V[1,2] = -.0000148505006577
mat T_V[1,3] = -.0004183027229316
mat T_V[1,4] = -.0000154516506638
mat T_V[1,6] = -.0006718711581376
mat T_V[2,1] = -.0000148505006577
mat T_V[2,2] =  .0002519890363482
mat T_V[2,3] = -.0000141735940032
mat T_V[2,4] =  3.66900112425e-06
mat T_V[2,6] =  .0000822842268142
mat T_V[3,1] = -.0004183027229316
mat T_V[3,2] = -.0000141735940032
mat T_V[3,3] =  .0004548558602652
mat T_V[3,4] = -.0000601132492111
mat T_V[3,6] = -.0021270697784458
mat T_V[4,1] = -.0000154516506638
mat T_V[4,2] =  3.66900112425e-06
mat T_V[4,3] = -.0000601132492111
mat T_V[4,4] =  .0001688021560809
mat T_V[4,6] =  .0010032030448185
mat T_V[6,1] = -.0006718711581376
mat T_V[6,2] =  .0000822842268142
mat T_V[6,3] = -.0021270697784458
mat T_V[6,4] =  .0010032030448185
mat T_V[6,6] =  .5081109589346494
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7258941009043415
mat T_rcsrmat_1[2,1] =  1.217614314980523
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -.0009563866775216
mat T_gradient[1,2] =  .0000280214283434
mat T_gradient[1,3] = -.0019662461439669
mat T_gradient[1,4] =  .0015224727192591
mat T_gradient[1,6] =  .0028729909633992
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  1.184565247291186
mat T_ml_scale[1,2] =  4.407351521420006
mat T_ml_scale[1,3] =  .3847001785034195
mat T_ml_scale[1,4] =  .3052708474003165
mat T_ml_scale[1,5] =  6.303021725562637
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale




uhtred (stime trt age M1[id]@1, family(rp, df(1) failure(died)))	///
	, intmethod(gh)

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id]@1, family(rp, df(1) failure(died)))                , intmethod(gh)"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"ghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd(M1)"'
assert `"`e(re_ivscale1)'"'       == `"exp"'
assert `"`e(re_eqns1)'"'          == `"lns1_1"'
assert `"`e(Nreparams1)'"'        == `"1"'
assert `"`e(Nres1)'"'             == `"1"'
assert `"`e(latents1)'"'          == `"M1"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-5.24262 1.609222"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age died id stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"uhtred_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 5
assert         e(N)          == 10000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -13338.61300964758) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.5090135381739073
mat T_b[1,2] =  .1075778475352739
mat T_b[1,3] = -.8548985836193904
mat T_b[1,4] =   .891836431782292
mat T_b[1,5] =                  1
mat T_b[1,6] =  -2.39893900327269
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0010224484517523
mat T_V[1,2] = -.0000148505147674
mat T_V[1,3] = -.0004182924726809
mat T_V[1,4] = -.0000154522945516
mat T_V[1,6] = -.0006725717932497
mat T_V[2,1] = -.0000148505147674
mat T_V[2,2] =   .000251984587817
mat T_V[2,3] = -.0000141738057893
mat T_V[2,4] =  3.66908346179e-06
mat T_V[2,6] =  .0000823769755409
mat T_V[3,1] = -.0004182924726809
mat T_V[3,2] = -.0000141738057893
mat T_V[3,3] =  .0004548508829956
mat T_V[3,4] = -.0000601153854098
mat T_V[3,6] = -.0021291501072252
mat T_V[4,1] = -.0000154522945516
mat T_V[4,2] =  3.66908346179e-06
mat T_V[4,3] = -.0000601153854098
mat T_V[4,4] =  .0001688030887502
mat T_V[4,6] =   .001004239423796
mat T_V[6,1] = -.0006725717932497
mat T_V[6,2] =  .0000823769755409
mat T_V[6,3] = -.0021291501072252
mat T_V[6,4] =   .001004239423796
mat T_V[6,6] =   .508888932214358
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7258941009043415
mat T_rcsrmat_1[2,1] =  1.217614314980523
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -.0010490813898815
mat T_gradient[1,2] =  .0000186970843794
mat T_gradient[1,3] = -.0021371724695113
mat T_gradient[1,4] =  .0015973664924546
mat T_gradient[1,6] =  .0034266098702524
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  1.185592618417029
mat T_ml_scale[1,2] =  4.665021874177106
mat T_ml_scale[1,3] =  .3823129419630767
mat T_ml_scale[1,4] =  .3034110717937598
mat T_ml_scale[1,5] =  5.952471511850828
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale
