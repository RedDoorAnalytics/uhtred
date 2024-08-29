
set seed 7254
clear 
set obs 100
gen id1 = _n
gen age = runiform()
gen sd1 = exp(log(1))
gen u1 = rnormal(0,sd1)
expand 10
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(1))
gen u2 = rnormal(0,sd2)
expand 5
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.02 u1 1 u2 1) maxt(10)
stset stime1, f(dead1)
replace id2 = _n

uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 , family(rp, df(1) failure(dead1))) ///
			, evaltype(gf0) //intmethod(gh) //intpoints(15) 


_assert_streq `"`e(cmdline)'"' `"uhtred (stime1 trt age M2[id1>id2]@1 M1[id1]@1 , family(rp, df(1) failure(dead1)))                         , evaltype(gf0)"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod2)'"'        == `"mvaghermite"'
assert `"`e(intpoints2)'"'        == `"7"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist2)'"'          == `"normal"'
assert `"`e(re_label2)'"'         == `"sd(M2)"'
assert `"`e(re_ivscale2)'"'       == `"exp"'
assert `"`e(re_eqns2)'"'          == `"lns2_1"'
assert `"`e(Nreparams2)'"'        == `"1"'
assert `"`e(Nres2)'"'             == `"1"'
assert `"`e(latents2)'"'          == `"M2"'
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
assert `"`e(knots1)'"'            == `"-6.34343 2.29853"'
assert `"`e(timevar1)'"'          == `"stime1"'
assert `"`e(failure1)'"'          == `"dead1"'
assert `"`e(response1)'"'         == `"stime1 dead1"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(levelvars)'"'         == `"id1 id2"'
assert `"`e(haszb)'"'             == `"1"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age dead1 id1 id2 stime1 trt"'
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

assert         e(rank)       == 6
assert         e(N)          == 5000
assert         e(k)          == 8
assert         e(k_eq)       == 6
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -9585.643594362182) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 3

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5491320542814528
mat T_b[1,2] = -.2336421414039778
mat T_b[1,3] = -.5269678611249945
mat T_b[1,4] =  1.531185159748634
mat T_b[1,5] =                  1
mat T_b[1,6] =                  1
mat T_b[1,7] =  .0541256906912157
mat T_b[1,8] = -.0095884008530593
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0031701264152182
mat T_V[1,2] =  .0002984245396177
mat T_V[1,3] =  -.001094187152354
mat T_V[1,4] = -.0016090275086677
mat T_V[1,7] = -.0011910195751481
mat T_V[1,8] = -.0026640127443669
mat T_V[2,1] =  .0002984245396177
mat T_V[2,2] =   .144717978845459
mat T_V[2,3] = -.0741557994086873
mat T_V[2,4] = -.0006839046467606
mat T_V[2,7] = -.0003181964435937
mat T_V[2,8] = -.0011311114894928
mat T_V[3,1] =  -.001094187152354
mat T_V[3,2] = -.0741557994086873
mat T_V[3,3] =  .0503760324643051
mat T_V[3,4] = -.0002373611422725
mat T_V[3,7] = -.0003520385658214
mat T_V[3,8] = -.0003270774586831
mat T_V[4,1] = -.0016090275086677
mat T_V[4,2] = -.0006839046467606
mat T_V[4,3] = -.0002373611422725
mat T_V[4,4] =  .0039735950789341
mat T_V[4,7] =  .0027072019765935
mat T_V[4,8] =   .006224768083078
mat T_V[7,1] = -.0011910195751481
mat T_V[7,2] = -.0003181964435937
mat T_V[7,3] = -.0003520385658214
mat T_V[7,4] =  .0027072019765935
mat T_V[7,7] =  .0076234323539407
mat T_V[7,8] =  .0044688105615897
mat T_V[8,1] = -.0026640127443669
mat T_V[8,2] = -.0011311114894928
mat T_V[8,3] = -.0003270774586831
mat T_V[8,4] =   .006224768083078
mat T_V[8,7] =  .0044688105615897
mat T_V[8,8] =   .011013187820039
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  1.309925998200917
mat T_rcsrmat_1[2,1] =  1.113771202883215
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] = -1.62074233638e-07
mat T_gradient[1,2] =  1.53596960230e-08
mat T_gradient[1,3] =  3.77986646102e-08
mat T_gradient[1,4] =  4.25507025540e-06
mat T_gradient[1,7] =  2.81754926565e-07
mat T_gradient[1,8] =  7.73956053510e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,6,0)
mat T_ml_scale[1,1] =  2.251210976186122
mat T_ml_scale[1,2] =  16.70076034796118
mat T_ml_scale[1,3] =  5.554319597255222
mat T_ml_scale[1,4] =  .3617906385379899
mat T_ml_scale[1,5] =  37.17990566517763
mat T_ml_scale[1,6] =  74.87616561807033
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5 c6"'
mat drop C_ml_scale T_ml_scale
	