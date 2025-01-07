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
assert reldif( e(ll)          , -62035.27099400835) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 3

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5275958779155659
mat T_b[1,2] = -.0133006804956048
mat T_b[1,3] =  -.533877227357005
mat T_b[1,4] =  1.250286569093285
mat T_b[1,5] =                  1
mat T_b[1,6] =                  1
mat T_b[1,7] = -.6308212558685061
mat T_b[1,8] = -.6642988116574527
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-6
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =   .000578238521227
mat T_V[1,2] =  .0000167702571783
mat T_V[1,3] =  -.000275017278571
mat T_V[1,4] = -.0000115731830805
mat T_V[1,7] = -8.29013118604e-06
mat T_V[1,8] = -.0000250634538564
mat T_V[2,1] =  .0000167702571783
mat T_V[2,2] =  .0366474304602872
mat T_V[2,3] = -.0188043397826051
mat T_V[2,4] = -1.10282104241e-06
mat T_V[2,7] =  .0000338142635785
mat T_V[2,8] = -4.14981000004e-06
mat T_V[3,1] =  -.000275017278571
mat T_V[3,2] = -.0188043397826051
mat T_V[3,3] =   .012760076022562
mat T_V[3,4] = -.0000229600537787
mat T_V[3,7] = -.0000512715135316
mat T_V[3,8] = -.0000287834206469
mat T_V[4,1] = -.0000115731830805
mat T_V[4,2] = -1.10282104241e-06
mat T_V[4,3] = -.0000229600537787
mat T_V[4,4] =  .0000566343208537
mat T_V[4,7] =  .0000218224407806
mat T_V[4,8] =  .0000344584877208
mat T_V[7,1] = -8.29013118604e-06
mat T_V[7,2] =  .0000338142635785
mat T_V[7,3] = -.0000512715135316
mat T_V[7,4] =  .0000218224407806
mat T_V[7,7] =  .0055602329201207
mat T_V[7,8] =  .0000283438786424
mat T_V[8,1] = -.0000250634538564
mat T_V[8,2] = -4.14981000004e-06
mat T_V[8,3] = -.0000287834206469
mat T_V[8,4] =  .0000344584877208
mat T_V[8,7] =  .0000283438786424
mat T_V[8,8] =  .0004758276088331
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-6
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_V T_V

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
mat T_gradient[1,1] =   .001082696402701
mat T_gradient[1,2] =  .0007923357168968
mat T_gradient[1,3] =  .0022716928726266
mat T_gradient[1,4] =  .0007726865567516
mat T_gradient[1,7] =  .0017767524425885
mat T_gradient[1,8] =  .0001820719241663
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-5
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_gradient T_gradient
