set seed 725488
clear
set obs 1000
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,0.25\0.25,1)
drawnorm u1 u2, means(0 0) sds(1 0.5) corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt
gen age = rnormal() + u2

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(age 0.01 trtui 1 u1 1) maxt(5) 
stset stime, f(dead)

//random intercept
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(unstr) 




_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead))),               cov(unstr)"'
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
assert `"`e(knots1)'"'            == `"-6.056272 1.608266"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"dead"'
assert `"`e(response1)'"'         == `"stime dead"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id1"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age dead id1 stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
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
assert         e(N)          == 5000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -6777.468734858823) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4313294896398842
mat T_b[1,2] =  .0737414484770948
mat T_b[1,3] =  -.999213524547046
mat T_b[1,4] =  1.065943738267433
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0203946488854133
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0021520822321426
mat T_V[1,2] =  .0000192282475498
mat T_V[1,3] =  -.000918172031606
mat T_V[1,4] = -.0000617493829473
mat T_V[1,6] = -.0000600833757343
mat T_V[2,1] =  .0000192282475498
mat T_V[2,2] =  .0004884427145159
mat T_V[2,3] =  -.000023521986417
mat T_V[2,4] = -1.04971810800e-06
mat T_V[2,6] = -.0000689062303826
mat T_V[3,1] =  -.000918172031606
mat T_V[3,2] =  -.000023521986417
mat T_V[3,3] =  .0021451675935357
mat T_V[3,4] = -.0001605526928973
mat T_V[3,6] = -.0004982015037962
mat T_V[4,1] = -.0000617493829473
mat T_V[4,2] = -1.04971810800e-06
mat T_V[4,3] = -.0001605526928973
mat T_V[4,4] =  .0004089132632236
mat T_V[4,6] =   .000211349769869
mat T_V[6,1] = -.0000600833757343
mat T_V[6,2] = -.0000689062303826
mat T_V[6,3] = -.0004982015037962
mat T_V[6,4] =   .000211349769869
mat T_V[6,6] =  .0015313280175845
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .9075200192305893
mat T_rcsrmat_1[2,1] =  1.049222979798123
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] =  .0000623383947618
mat T_gradient[1,2] =  .0000194370517971
mat T_gradient[1,3] =  .0001432282908904
mat T_gradient[1,4] =  .0000128371226488
mat T_gradient[1,6] =  5.00086507045e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  2.177544289258381
mat T_ml_scale[1,2] =  6.139947824211223
mat T_ml_scale[1,3] =  .5431294461628475
mat T_ml_scale[1,4] =  .2118345580003686
mat T_ml_scale[1,5] =   27.7669200951121
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale



//diagonal covariance matrix
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(diag) 


_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead))),               cov(diag)"'
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
assert `"`e(knots1)'"'            == `"-6.056272 1.608266"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"dead"'
assert `"`e(response1)'"'         == `"stime dead"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id1"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age dead id1 stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
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
assert         e(N)          == 5000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -6777.468734858823) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4313294896398842
mat T_b[1,2] =  .0737414484770948
mat T_b[1,3] =  -.999213524547046
mat T_b[1,4] =  1.065943738267433
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0203946488854133
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0021520822321426
mat T_V[1,2] =  .0000192282475498
mat T_V[1,3] =  -.000918172031606
mat T_V[1,4] = -.0000617493829473
mat T_V[1,6] = -.0000600833757343
mat T_V[2,1] =  .0000192282475498
mat T_V[2,2] =  .0004884427145159
mat T_V[2,3] =  -.000023521986417
mat T_V[2,4] = -1.04971810800e-06
mat T_V[2,6] = -.0000689062303826
mat T_V[3,1] =  -.000918172031606
mat T_V[3,2] =  -.000023521986417
mat T_V[3,3] =  .0021451675935357
mat T_V[3,4] = -.0001605526928973
mat T_V[3,6] = -.0004982015037962
mat T_V[4,1] = -.0000617493829473
mat T_V[4,2] = -1.04971810800e-06
mat T_V[4,3] = -.0001605526928973
mat T_V[4,4] =  .0004089132632236
mat T_V[4,6] =   .000211349769869
mat T_V[6,1] = -.0000600833757343
mat T_V[6,2] = -.0000689062303826
mat T_V[6,3] = -.0004982015037962
mat T_V[6,4] =   .000211349769869
mat T_V[6,6] =  .0015313280175845
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .9075200192305893
mat T_rcsrmat_1[2,1] =  1.049222979798123
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] =  .0000623383947618
mat T_gradient[1,2] =  .0000194370517971
mat T_gradient[1,3] =  .0001432282908904
mat T_gradient[1,4] =  .0000128371226488
mat T_gradient[1,6] =  5.00086507045e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  2.177544289258381
mat T_ml_scale[1,2] =  6.139947824211223
mat T_ml_scale[1,3] =  .5431294461628475
mat T_ml_scale[1,4] =  .2118345580003686
mat T_ml_scale[1,5] =   27.7669200951121
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale


//identiy covariance matrix

uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	cov(iden) 
	


_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead))),               cov(iden)"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd"'
assert `"`e(re_ivscale1)'"'       == `"exp"'
assert `"`e(re_eqns1)'"'          == `"lns1_1"'
assert `"`e(Nreparams1)'"'        == `"1"'
assert `"`e(Nres1)'"'             == `"1"'
assert `"`e(latents1)'"'          == `"M1"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-6.056272 1.608266"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"dead"'
assert `"`e(response1)'"'         == `"stime dead"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id1"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age dead id1 stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
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
assert         e(N)          == 5000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -6777.468734858823) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4313294896398842
mat T_b[1,2] =  .0737414484770948
mat T_b[1,3] =  -.999213524547046
mat T_b[1,4] =  1.065943738267433
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0203946488854133
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0021520822321426
mat T_V[1,2] =  .0000192282475498
mat T_V[1,3] =  -.000918172031606
mat T_V[1,4] = -.0000617493829473
mat T_V[1,6] = -.0000600833757343
mat T_V[2,1] =  .0000192282475498
mat T_V[2,2] =  .0004884427145159
mat T_V[2,3] =  -.000023521986417
mat T_V[2,4] = -1.04971810800e-06
mat T_V[2,6] = -.0000689062303826
mat T_V[3,1] =  -.000918172031606
mat T_V[3,2] =  -.000023521986417
mat T_V[3,3] =  .0021451675935357
mat T_V[3,4] = -.0001605526928973
mat T_V[3,6] = -.0004982015037962
mat T_V[4,1] = -.0000617493829473
mat T_V[4,2] = -1.04971810800e-06
mat T_V[4,3] = -.0001605526928973
mat T_V[4,4] =  .0004089132632236
mat T_V[4,6] =   .000211349769869
mat T_V[6,1] = -.0000600833757343
mat T_V[6,2] = -.0000689062303826
mat T_V[6,3] = -.0004982015037962
mat T_V[6,4] =   .000211349769869
mat T_V[6,6] =  .0015313280175845
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .9075200192305893
mat T_rcsrmat_1[2,1] =  1.049222979798123
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] =  .0000623383947618
mat T_gradient[1,2] =  .0000194370517971
mat T_gradient[1,3] =  .0001432282908904
mat T_gradient[1,4] =  .0000128371226488
mat T_gradient[1,6] =  5.00086507045e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  2.177544289258381
mat T_ml_scale[1,2] =  6.139947824211223
mat T_ml_scale[1,3] =  .5431294461628475
mat T_ml_scale[1,4] =  .2118345580003686
mat T_ml_scale[1,5] =   27.7669200951121
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale


//intmethod - ghermite
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))), 	///
	intmethod(ghermite)

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead))),               intmethod(ghermite)"'
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
assert `"`e(knots1)'"'            == `"-6.056272 1.608266"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"dead"'
assert `"`e(response1)'"'         == `"stime dead"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id1"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age dead id1 stime trt"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"uhtred"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"uhtred_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
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
assert         e(N)          == 5000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -6777.715872128853) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4288531503102935
mat T_b[1,2] =  .0759828027094849
mat T_b[1,3] = -.9961364945163345
mat T_b[1,4] =  1.063901997473898
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0101973074565366
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0021547581526675
mat T_V[1,2] =  .0000113485028967
mat T_V[1,3] = -.0008998466316245
mat T_V[1,4] = -.0000610450239076
mat T_V[1,6] =  5.28117697684e-07
mat T_V[2,1] =  .0000113485028967
mat T_V[2,2] =  .0004895657880288
mat T_V[2,3] = -.0000248525446527
mat T_V[2,4] =  4.66351688681e-08
mat T_V[2,6] = -.0000450736535445
mat T_V[3,1] = -.0008998466316245
mat T_V[3,2] = -.0000248525446527
mat T_V[3,3] =  .0020760020011815
mat T_V[3,4] = -.0001474690280692
mat T_V[3,6] = -.0004705633458991
mat T_V[4,1] = -.0000610450239076
mat T_V[4,2] =  4.66351688681e-08
mat T_V[4,3] = -.0001474690280692
mat T_V[4,4] =  .0004044973538514
mat T_V[4,6] =   .000190569871074
mat T_V[6,1] =  5.28117697684e-07
mat T_V[6,2] = -.0000450736535445
mat T_V[6,3] = -.0004705633458991
mat T_V[6,4] =   .000190569871074
mat T_V[6,6] =  .0013614384183595
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .9075200192305893
mat T_rcsrmat_1[2,1] =  1.049222979798123
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] =  -.000134000759785
mat T_gradient[1,2] = -6.81766376460e-06
mat T_gradient[1,3] = -.0003029311737403
mat T_gradient[1,4] = -.0001290999217295
mat T_gradient[1,6] = -.0001002355568886
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,5,0)
mat T_ml_scale[1,1] =  1.754900305256777
mat T_ml_scale[1,2] =  4.802873579846938
mat T_ml_scale[1,3] =  .9205772382461238
mat T_ml_scale[1,4] =  .2581279331648406
mat T_ml_scale[1,5] =  89.95532708421773
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5"'
mat drop C_ml_scale T_ml_scale



	
	