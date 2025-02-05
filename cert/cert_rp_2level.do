set seed 98798
clear
set obs 1000
gen id1 = _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.1))
gen u1 	= rnormal(0,sd1)
gen age = rnormal()
expand 100
sort id1 
survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.1 u1 1) maxt(5) 


uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))) 	///
	,

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead)))                ,"'
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
assert `"`e(knots1)'"'            == `"-9.125254 1.609348"'
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
assert `"`e(hasopts)'"'           == `"0"'
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
assert         e(N)          == 100000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -134291.9624137046) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4921848428391861
mat T_b[1,2] =  .0990174471193304
mat T_b[1,3] = -.8553735012897704
mat T_b[1,4] =  .9004400477682779
mat T_b[1,5] =                  1
mat T_b[1,6] = -2.403852348528028
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-5
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0001305355523448
mat T_V[1,2] = -1.58152133276e-06
mat T_V[1,3] = -.0000569707843571
mat T_V[1,4] = -1.36065118746e-06
mat T_V[1,6] = -.0000147055433378
mat T_V[2,1] = -1.58152133276e-06
mat T_V[2,2] =  .0000324270684163
mat T_V[2,3] = -1.24285255279e-06
mat T_V[2,4] =  2.82982280814e-07
mat T_V[2,6] =  5.18085838397e-06
mat T_V[3,1] = -.0000569707843571
mat T_V[3,2] = -1.24285255279e-06
mat T_V[3,3] =  .0000596463190987
mat T_V[3,4] = -5.75165969224e-06
mat T_V[3,6] = -.0000319316780844
mat T_V[4,1] = -1.36065118746e-06
mat T_V[4,2] =  2.82982280814e-07
mat T_V[4,3] = -5.75165969224e-06
mat T_V[4,4] =  .0000168649223197
mat T_V[4,6] =  .0000110230343581
mat T_V[6,1] = -.0000147055433378
mat T_V[6,2] =  5.18085838397e-06
mat T_V[6,3] = -.0000319316780844
mat T_V[6,4] =  .0000110230343581
mat T_V[6,6] =  .0075509125694396
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-6
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7481013813490104
mat T_rcsrmat_1[2,1] =  1.205948296908631
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -9.69092831434e-06
mat T_gradient[1,2] =  .0000265881332365
mat T_gradient[1,3] = -.0007384310012594
mat T_gradient[1,4] = -.0000136022756857
mat T_gradient[1,6] =  1.41846863523e-06
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-3
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient


// intmethod ghermite
uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))) 	///
	, intmethod(ghermite)


_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt age M1[id1]@1,         family(rp, df(1) failure(dead)))                , intmethod(ghermite)"'
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
assert `"`e(knots1)'"'            == `"-9.125254 1.609348"'
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
assert         e(N)          == 100000
assert         e(k)          == 6
assert         e(k_eq)       == 4
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -134291.9624188172) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] = -.4921850458143355
mat T_b[1,2] =  .0990174675848569
mat T_b[1,3] = -.8553732327752838
mat T_b[1,4] =  .9004400423362672
mat T_b[1,5] =                  1
mat T_b[1,6] = -2.403861656043603
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-4
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =  .0001305351120025
mat T_V[1,2] = -1.58168634287e-06
mat T_V[1,3] = -.0000569705175292
mat T_V[1,4] = -1.36069106796e-06
mat T_V[1,6] = -.0000147270796171
mat T_V[2,1] = -1.58168634287e-06
mat T_V[2,2] =  .0000324270826439
mat T_V[2,3] = -1.24270262296e-06
mat T_V[2,4] =  2.82986790589e-07
mat T_V[2,6] =  5.18145598297e-06
mat T_V[3,1] = -.0000569705175292
mat T_V[3,2] = -1.24270262296e-06
mat T_V[3,3] =  .0000596458925368
mat T_V[3,4] = -5.75161768943e-06
mat T_V[3,6] = -.0000319093738354
mat T_V[4,1] = -1.36069106796e-06
mat T_V[4,2] =  2.82986790589e-07
mat T_V[4,3] = -5.75161768943e-06
mat T_V[4,4] =  .0000168649215695
mat T_V[4,6] =  .0000110230785301
mat T_V[6,1] = -.0000147270796171
mat T_V[6,2] =  5.18145598297e-06
mat T_V[6,3] = -.0000319093738354
mat T_V[6,4] =  .0000110230785301
mat T_V[6,6] =  .0075511411958584
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-6
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .7481013813490104
mat T_rcsrmat_1[2,1] =  1.205948296908631
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -.0001613780398727
mat T_gradient[1,2] = -3.41148490390e-06
mat T_gradient[1,3] = -.0002625436828968
mat T_gradient[1,4] =  .0000301797368251
mat T_gradient[1,6] = -.0009169437781076
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-3
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 lns1_1:_cons"'
mat drop C_gradient T_gradient
