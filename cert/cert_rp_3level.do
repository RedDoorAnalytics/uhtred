clear 
set seed 7254

set obs 100
gen id1 = _n
gen age = runiform()
gen u1 = rnormal(0,0.5)
expand 30
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen u2 = rnormal(0,0.5)
expand 10

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.02 u1 1 u2 1) maxt(10)

uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 ///
	, family(rp, df(1) failure(dead1)))
_assert_streq `"`e(cmdline)'"' `"uhtred (stime1 trt age M2[id1>id2]@1 M1[id1]@1         , family(rp, df(1) failure(dead1)))"'
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
assert `"`e(knots1)'"'            == `"-10.19983 2.302559"'
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

assert         e(rank)       == 6
assert         e(N)          == 30000
assert         e(k)          == 8
assert         e(k_eq)       == 6
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -62035.20258137111) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 3

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5266248309557517
mat T_b[1,2] = -.0126131958863629
mat T_b[1,3] = -.5344043680778759
mat T_b[1,4] =  1.250338988613769
mat T_b[1,5] =                  1
mat T_b[1,6] =                  1
mat T_b[1,7] = -.6301954196905883
mat T_b[1,8] = -.6625471841659873
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0005802724805071
mat T_V[1,2] =  .0000190567546559
mat T_V[1,3] = -.0002780645449972
mat T_V[1,4] = -.0000114406897714
mat T_V[1,7] = -8.04525482318e-06
mat T_V[1,8] = -.0000229714036142
mat T_V[2,1] =  .0000190567546559
mat T_V[2,2] =    .03670823087781
mat T_V[2,3] = -.0188387465675493
mat T_V[2,4] = -1.26242895676e-06
mat T_V[2,7] =  .0000239913628653
mat T_V[2,8] =  7.11606759352e-07
mat T_V[3,1] = -.0002780645449972
mat T_V[3,2] = -.0188387465675493
mat T_V[3,3] =   .012784512782954
mat T_V[3,4] = -.0000228212650069
mat T_V[3,7] = -.0000470439676194
mat T_V[3,8] = -.0000315314900431
mat T_V[4,1] = -.0000114406897714
mat T_V[4,2] = -1.26242895676e-06
mat T_V[4,3] = -.0000228212650069
mat T_V[4,4] =  .0000566675354212
mat T_V[4,7] =  .0000220231305721
mat T_V[4,8] =  .0000349586630569
mat T_V[7,1] = -8.04525482318e-06
mat T_V[7,2] =  .0000239913628653
mat T_V[7,3] = -.0000470439676194
mat T_V[7,4] =  .0000220231305721
mat T_V[7,7] =  .0055597862615126
mat T_V[7,8] =  .0000295152340157
mat T_V[8,1] = -.0000229714036142
mat T_V[8,2] =  7.11606759352e-07
mat T_V[8,3] = -.0000315314900431
mat T_V[8,4] =  .0000349586630569
mat T_V[8,7] =  .0000295152340157
mat T_V[8,8] =  .0004852772962281
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_V T_V

mat T_VCV_2=  .2649973911943641
matrix C_VCV_2 = e(VCV_2)
assert mreldif( C_VCV_2 , T_VCV_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_VCV_2'"' `"r1"'
_assert_streq `"`: colfullnames C_VCV_2'"' `"c1"'
mat drop C_VCV_2 T_VCV_2

mat T_VCV_1=  .2813739094685654
matrix C_VCV_1 = e(VCV_1)
assert mreldif( C_VCV_1 , T_VCV_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_VCV_1'"' `"r1"'
_assert_streq `"`: colfullnames C_VCV_1'"' `"c1"'
mat drop C_VCV_1 T_VCV_1

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  1.044445335194847
mat T_rcsrmat_1[2,1] =    1.3918523564047
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] = -2.98585765501e-07
mat T_gradient[1,2] = -2.27578111697e-07
mat T_gradient[1,3] = -4.44040648362e-07
mat T_gradient[1,4] = -2.39608692998e-07
mat T_gradient[1,7] =  5.63992212307e-09
mat T_gradient[1,8] = -8.67212537122e-06
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,6,0)
mat T_ml_scale[1,1] =  2.158221177195881
mat T_ml_scale[1,2] =  390.3911933953053
mat T_ml_scale[1,3] =   6.97173256613322
mat T_ml_scale[1,4] =  .3909230917934726
mat T_ml_scale[1,5] =  6.083706978491136
mat T_ml_scale[1,6] =  2.216066468505855
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5 c6"'
mat drop C_ml_scale T_ml_scale
