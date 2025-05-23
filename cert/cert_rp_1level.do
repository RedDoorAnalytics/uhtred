
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
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		maxt(10) 


//============================================================================//
// rp ph model - df1

uhtred (stime trt bmi x1 x2 x3, family(rp, df(1) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3, family(rp, df(1) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.988795 2.302106"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5773.855544144517) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,7,0)
mat T_b[1,1] = -.5165850072254397
mat T_b[1,2] = -.0536506242228911
mat T_b[1,3] =   .089041697808268
mat T_b[1,4] =  -.413644694215166
mat T_b[1,5] =  .5372759072971561
mat T_b[1,6] =  .2441742218503316
mat T_b[1,7] =  .7518285636230789
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(7,7,0)
mat T_V[1,1] =  .0030145006933191
mat T_V[1,2] = -.0000126756772554
mat T_V[1,3] = -1.11665596459e-06
mat T_V[1,4] =   .000037503043091
mat T_V[1,5] = -.0000201625291158
mat T_V[1,6] = -.0007976501800245
mat T_V[1,7] = -.0000323949728692
mat T_V[2,1] = -.0000126756772554
mat T_V[2,2] =  .0000790203453125
mat T_V[2,3] =  5.15071356371e-06
mat T_V[2,4] =  2.24157685316e-07
mat T_V[2,5] = -4.53027950291e-06
mat T_V[2,6] = -.0023309437739427
mat T_V[2,7] = -3.33946982304e-06
mat T_V[3,1] = -1.11665596459e-06
mat T_V[3,2] =  5.15071356371e-06
mat T_V[3,3] =   .000744333706991
mat T_V[3,4] =  .0000137453122933
mat T_V[3,5] =  9.22703592544e-06
mat T_V[3,6] = -.0002131006546416
mat T_V[3,7] =  5.64230701876e-06
mat T_V[4,1] =   .000037503043091
mat T_V[4,2] =  2.24157685316e-07
mat T_V[4,3] =  .0000137453122933
mat T_V[4,4] =  .0007151959219153
mat T_V[4,5] = -.0000263065396501
mat T_V[4,6] =  .0002403789029903
mat T_V[4,7] = -.0000254427058897
mat T_V[5,1] = -.0000201625291158
mat T_V[5,2] = -4.53027950291e-06
mat T_V[5,3] =  9.22703592544e-06
mat T_V[5,4] = -.0000263065396501
mat T_V[5,5] =   .000756729999554
mat T_V[5,6] = -.0002026337664198
mat T_V[5,7] =  .0000369537907309
mat T_V[6,1] = -.0007976501800245
mat T_V[6,2] = -.0023309437739427
mat T_V[6,3] = -.0002131006546416
mat T_V[6,4] =  .0002403789029903
mat T_V[6,5] = -.0002026337664198
mat T_V[6,6] =  .0702030657604157
mat T_V[6,7] =  3.76370988110e-06
mat T_V[7,1] = -.0000323949728692
mat T_V[7,2] = -3.33946982304e-06
mat T_V[7,3] =  5.64230701876e-06
mat T_V[7,4] = -.0000254427058897
mat T_V[7,5] =  .0000369537907309
mat T_V[7,6] =  3.76370988110e-06
mat T_V[7,7] =  .0003598549070492
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[2,1] =  2.026577982062658
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,7,0)
mat T_gradient[1,1] = -1.06020658631e-07
mat T_gradient[1,2] = -6.88013358885e-06
mat T_gradient[1,3] = -9.96469736826e-09
mat T_gradient[1,4] =  3.57457493663e-08
mat T_gradient[1,5] = -3.30037669777e-08
mat T_gradient[1,6] = -2.29013358877e-07
mat T_gradient[1,7] = -3.65018775197e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1"'
mat drop C_gradient T_gradient


//============================================================================//
// rp ph model - 3 df

uhtred (stime trt bmi x1 x2 x3, family(rp, df(3) failure(died)))
		
_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3, family(rp, df(3) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.988795 1.180037 1.83553 2.302106"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5771.510731746175) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5153995342969059
mat T_b[1,2] = -.0535922914694009
mat T_b[1,3] =  .0888076242418419
mat T_b[1,4] =   -.41259930828416
mat T_b[1,5] =  .5357406858418392
mat T_b[1,6] =  .2404661331914562
mat T_b[1,7] =  .7581216834622467
mat T_b[1,8] = -.0001661999583854
mat T_b[1,9] =  .0138938827563727
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =  .0030150973575314
mat T_V[1,2] = -.0000124942593313
mat T_V[1,3] = -1.24733504632e-06
mat T_V[1,4] =  .0000381589707626
mat T_V[1,5] = -.0000205617823404
mat T_V[1,6] =  -.000804420227059
mat T_V[1,7] = -.0000245130930593
mat T_V[1,8] =  .0000104753772738
mat T_V[1,9] =  1.79671182040e-06
mat T_V[2,1] = -.0000124942593313
mat T_V[2,2] =  .0000790410110184
mat T_V[2,3] =  5.10001886569e-06
mat T_V[2,4] =  2.97155357224e-07
mat T_V[2,5] = -4.56423056819e-06
mat T_V[2,6] = -.0023316911930153
mat T_V[2,7] = -2.73896825730e-06
mat T_V[2,8] =  8.71475640722e-07
mat T_V[2,9] =  2.06451770979e-08
mat T_V[3,1] = -1.24733504632e-06
mat T_V[3,2] =  5.10001886569e-06
mat T_V[3,3] =  .0007442741860809
mat T_V[3,4] =  .0000137162658281
mat T_V[3,5] =  9.43677039895e-06
mat T_V[3,6] =  -.000211379442549
mat T_V[3,7] =  4.41074266657e-06
mat T_V[3,8] = -1.47462713097e-06
mat T_V[3,9] = -4.34778507512e-07
mat T_V[4,1] =  .0000381589707626
mat T_V[4,2] =  2.97155357224e-07
mat T_V[4,3] =  .0000137162658281
mat T_V[4,4] =  .0007153967565176
mat T_V[4,5] = -.0000266476550588
mat T_V[4,6] =   .000236902356508
mat T_V[4,7] = -.0000193326320682
mat T_V[4,8] =  7.60786493674e-06
mat T_V[4,9] =  1.82043175771e-06
mat T_V[5,1] = -.0000205617823404
mat T_V[5,2] = -4.56423056819e-06
mat T_V[5,3] =  9.43677039895e-06
mat T_V[5,4] = -.0000266476550588
mat T_V[5,5] =   .000757442943514
mat T_V[5,6] = -.0002000301639045
mat T_V[5,7] =  .0000278998654524
mat T_V[5,8] = -.0000113317375079
mat T_V[5,9] = -2.66020510717e-06
mat T_V[6,1] =  -.000804420227059
mat T_V[6,2] = -.0023316911930153
mat T_V[6,3] =  -.000211379442549
mat T_V[6,4] =   .000236902356508
mat T_V[6,5] = -.0002000301639045
mat T_V[6,6] =  .0702316269602105
mat T_V[6,7] = -.0000277418773817
mat T_V[6,8] = -.0000421808970207
mat T_V[6,9] = -6.46074378931e-06
mat T_V[7,1] = -.0000245130930593
mat T_V[7,2] = -2.73896825730e-06
mat T_V[7,3] =  4.41074266657e-06
mat T_V[7,4] = -.0000193326320682
mat T_V[7,5] =  .0000278998654524
mat T_V[7,6] = -.0000277418773817
mat T_V[7,7] =  .0004511324474014
mat T_V[7,8] =  .0001357674727107
mat T_V[7,9] =  9.88719018117e-06
mat T_V[8,1] =  .0000104753772738
mat T_V[8,2] =  8.71475640722e-07
mat T_V[8,3] = -1.47462713097e-06
mat T_V[8,4] =  7.60786493674e-06
mat T_V[8,5] = -.0000113317375079
mat T_V[8,6] = -.0000421808970207
mat T_V[8,7] =  .0001357674727107
mat T_V[8,8] =  .0002428672306644
mat T_V[8,9] = -.0000274336469453
mat T_V[9,1] =  1.79671182040e-06
mat T_V[9,2] =  2.06451770979e-08
mat T_V[9,3] = -4.34778507512e-07
mat T_V[9,4] =  1.82043175771e-06
mat T_V[9,5] = -2.66020510717e-06
mat T_V[9,6] = -6.46074378931e-06
mat T_V[9,7] =  9.88719018117e-06
mat T_V[9,8] = -.0000274336469453
mat T_V[9,9] =  .0000437562124964
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[1,2] = -8.435695793572606
mat T_rcsrmat_1[1,3] = -3.658962257730614
mat T_rcsrmat_1[2,2] =  1.843410155446673
mat T_rcsrmat_1[2,3] =  .8567385999018583
mat T_rcsrmat_1[3,3] =  .0586748184202542
mat T_rcsrmat_1[4,1] =  2.026577982062658
mat T_rcsrmat_1[4,2] = -38.92980822103527
mat T_rcsrmat_1[4,3] = -16.56962891904038
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] = -4.31813265700e-11
mat T_gradient[1,2] = -4.00527078170e-09
mat T_gradient[1,3] = -3.87018547070e-11
mat T_gradient[1,4] =  1.24819546788e-10
mat T_gradient[1,5] = -1.45581349201e-10
mat T_gradient[1,6] = -1.40575849572e-10
mat T_gradient[1,7] =  2.27253352724e-08
mat T_gradient[1,8] =  5.92493904494e-08
mat T_gradient[1,9] = -3.56981144723e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient

//orthog matrices check
qui {
	mat A = J(4,4,0)
	mat A[1,1] = .79717296
	mat A[1,2] = -12.329397
	mat A[1,3] = -5.0923271
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = 2.6519554
	mat A[2,3] = 1.1812266
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .09176479 
	mat A[3,4] =  0

	mat A[4,1] =  1.8714958
	mat A[4,2] = -47.20906
	mat A[4,3] =  -19.07834
	mat A[4,4] =  1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-0
matrix drop A


//============================================================================//
// rp ph model - factor variable

uhtred (stime trt bmi x1 x2 x3 i.agecat, family(rp, df(1) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3 i.agecat, family(rp, df(1) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.988795 2.302106"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died i.agecat stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5772.343694828301) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,12,0)
mat T_b[1,1] = -.5180644548945068
mat T_b[1,2] =  -.053773239619956
mat T_b[1,3] =  .0904214727189516
mat T_b[1,4] = -.4142380762020072
mat T_b[1,5] =  .5364411107983733
mat T_b[1,7] = -.0451782949293751
mat T_b[1,8] = -.0991370775358976
mat T_b[1,9] =  .0279503297215027
mat T_b[1,10] = -.1107925166469488
mat T_b[1,11] =  .2945961761497692
mat T_b[1,12] =  .7518785689821531
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_brcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(12,12,0)
mat T_V[1,1] =  .0030171271029654
mat T_V[1,2] = -.0000126519467403
mat T_V[1,3] = -4.31444621278e-06
mat T_V[1,4] =  .0000374385649067
mat T_V[1,5] = -.0000201913142852
mat T_V[1,7] =    .00007006524133
mat T_V[1,8] =  .0000442009137439
mat T_V[1,9] = -.0000416487286064
mat T_V[1,10] = -.0000965629083946
mat T_V[1,11] = -.0008284195040046
mat T_V[1,12] = -.0000327033877174
mat T_V[2,1] = -.0000126519467403
mat T_V[2,2] =   .000079206903353
mat T_V[2,3] =  4.82061276228e-06
mat T_V[2,4] =  4.64103870022e-07
mat T_V[2,5] = -4.50484443274e-06
mat T_V[2,7] = -3.68646073918e-06
mat T_V[2,8] =  6.95004975422e-06
mat T_V[2,9] =  3.91625633058e-06
mat T_V[2,10] =  6.77604987470e-06
mat T_V[2,11] = -.0023381158432773
mat T_V[2,12] = -3.33649397277e-06
mat T_V[3,1] = -4.31444621278e-06
mat T_V[3,2] =  4.82061276228e-06
mat T_V[3,3] =  .0007461960905315
mat T_V[3,4] =  .0000131732556022
mat T_V[3,5] =  8.74979594274e-06
mat T_V[3,7] =  6.71934881638e-06
mat T_V[3,8] = -.0000566018126755
mat T_V[3,9] =  .0000141394157482
mat T_V[3,10] =  .0001058328647808
mat T_V[3,11] = -.0001906684373963
mat T_V[3,12] =  5.76520505718e-06
mat T_V[4,1] =  .0000374385649067
mat T_V[4,2] =  4.64103870022e-07
mat T_V[4,3] =  .0000131732556022
mat T_V[4,4] =  .0007164154237498
mat T_V[4,5] = -.0000266638079506
mat T_V[4,7] = -.0000364141439575
mat T_V[4,8] =  3.18871688792e-06
mat T_V[4,9] = -.0000358553267062
mat T_V[4,10] = -.0000229014733413
mat T_V[4,11] =  .0002510938294609
mat T_V[4,12] = -.0000253471110321
mat T_V[5,1] = -.0000201913142852
mat T_V[5,2] = -4.50484443274e-06
mat T_V[5,3] =  8.74979594274e-06
mat T_V[5,4] = -.0000266638079506
mat T_V[5,5] =  .0007588485978328
mat T_V[5,7] = -.0000659324848737
mat T_V[5,8] = -.0000445270391569
mat T_V[5,9] =  -.000127217247061
mat T_V[5,10] = -.0000249606572409
mat T_V[5,11] = -.0001477972630655
mat T_V[5,12] =   .000036762352694
mat T_V[7,1] =    .00007006524133
mat T_V[7,2] = -3.68646073918e-06
mat T_V[7,3] =  6.71934881638e-06
mat T_V[7,4] = -.0000364141439575
mat T_V[7,5] = -.0000659324848737
mat T_V[7,7] =  .0064475117459744
mat T_V[7,8] =  .0042956458568713
mat T_V[7,9] =  .0043041077972424
mat T_V[7,10] =  .0042940084977788
mat T_V[7,11] = -.0042009954155355
mat T_V[7,12] = -3.23073432842e-06
mat T_V[8,1] =  .0000442009137439
mat T_V[8,2] =  6.95004975422e-06
mat T_V[8,3] = -.0000566018126755
mat T_V[8,4] =  3.18871688792e-06
mat T_V[8,5] = -.0000445270391569
mat T_V[8,7] =  .0042956458568713
mat T_V[8,8] =  .0065421254026263
mat T_V[8,9] =  .0042976978413333
mat T_V[8,10] =  .0042843045402803
mat T_V[8,11] = -.0044947818541108
mat T_V[8,12] = -3.23673945130e-07
mat T_V[9,1] = -.0000416487286064
mat T_V[9,2] =  3.91625633058e-06
mat T_V[9,3] =  .0000141394157482
mat T_V[9,4] = -.0000358553267062
mat T_V[9,5] =  -.000127217247061
mat T_V[9,7] =  .0043041077972424
mat T_V[9,8] =  .0042976978413333
mat T_V[9,9] =  .0090114572349993
mat T_V[9,10] =  .0043010819436228
mat T_V[9,11] = -.0043601118799485
mat T_V[9,12] =  1.66246939235e-06
mat T_V[10,1] = -.0000965629083946
mat T_V[10,2] =  6.77604987470e-06
mat T_V[10,3] =  .0001058328647808
mat T_V[10,4] = -.0000229014733413
mat T_V[10,5] = -.0000249606572409
mat T_V[10,7] =  .0042940084977788
mat T_V[10,8] =  .0042843045402803
mat T_V[10,9] =  .0043010819436228
mat T_V[10,10] =  .0387951779639044
mat T_V[10,11] = -.0044574872062519
mat T_V[10,12] = -9.75712895456e-06
mat T_V[11,1] = -.0008284195040046
mat T_V[11,2] = -.0023381158432773
mat T_V[11,3] = -.0001906684373963
mat T_V[11,4] =  .0002510938294609
mat T_V[11,5] = -.0001477972630655
mat T_V[11,7] = -.0042009954155355
mat T_V[11,8] = -.0044947818541108
mat T_V[11,9] = -.0043601118799485
mat T_V[11,10] = -.0044574872062519
mat T_V[11,11] =  .0740250187587558
mat T_V[11,12] =  5.06340979681e-06
mat T_V[12,1] = -.0000327033877174
mat T_V[12,2] = -3.33649397277e-06
mat T_V[12,3] =  5.76520505718e-06
mat T_V[12,4] = -.0000253471110321
mat T_V[12,5] =   .000036762352694
mat T_V[12,7] = -3.23073432842e-06
mat T_V[12,8] = -3.23673945130e-07
mat T_V[12,9] =  1.66246939235e-06
mat T_V[12,10] = -9.75712895456e-06
mat T_V[12,11] =  5.06340979681e-06
mat T_V[12,12] =  .0003599014850002
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_brcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_brcs1_1"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[2,1] =  2.026577982062658
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,12,0)
mat T_gradient[1,1] = -1.07003580122e-07
mat T_gradient[1,2] = -6.93463021539e-06
mat T_gradient[1,3] = -9.44598904083e-09
mat T_gradient[1,4] =  3.44995891856e-08
mat T_gradient[1,5] = -3.30384605245e-08
mat T_gradient[1,7] = -7.40298805486e-08
mat T_gradient[1,8] = -8.21853810701e-08
mat T_gradient[1,9] = -3.62032437431e-08
mat T_gradient[1,10] = -4.96870770916e-09
mat T_gradient[1,11] = -2.30740340310e-07
mat T_gradient[1,12] = -3.65599845563e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:0b.agecat xb1:1.agecat xb1:2.agecat xb1:3.agecat xb1:4.agecat xb1:_cons tb1:_brcs1_1"'
mat drop C_gradient T_gradient


//============================================================================//
// user-defined knots PH model
uhtred (stime trt bmi x1 x2 x3, family(rp, knots(-3.5 1 1.5 2) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3, family(rp, knots(-3.5 1 1.5 2) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.5 1 1.5 2"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5770.540651590882) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5151777682780363
mat T_b[1,2] = -.0535720215855634
mat T_b[1,3] =  .0887722448408774
mat T_b[1,4] = -.4125777134629683
mat T_b[1,5] =  .5356764428624713
mat T_b[1,6] =  .2396267866949445
mat T_b[1,7] =  .7582992465605269
mat T_b[1,8] = -.0032400707879088
mat T_b[1,9] =  .0196690009653912
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =   .003015112729449
mat T_V[1,2] = -.0000124855307004
mat T_V[1,3] = -1.18030742386e-06
mat T_V[1,4] =  .0000381415678833
mat T_V[1,5] = -.0000203979974521
mat T_V[1,6] = -.0008047995001103
mat T_V[1,7] = -.0000245192887961
mat T_V[1,8] =  9.95239885772e-06
mat T_V[1,9] =  3.00953179885e-06
mat T_V[2,1] = -.0000124855307004
mat T_V[2,2] =  .0000790287397036
mat T_V[2,3] =  5.11275693772e-06
mat T_V[2,4] =  2.91573861909e-07
mat T_V[2,5] = -4.56597553587e-06
mat T_V[2,6] = -.0023313413570737
mat T_V[2,7] = -2.73206862120e-06
mat T_V[2,8] =  8.25076198477e-07
mat T_V[2,9] =  1.22734398975e-07
mat T_V[3,1] = -1.18030742386e-06
mat T_V[3,2] =  5.11275693772e-06
mat T_V[3,3] =  .0007442340932892
mat T_V[3,4] =  .0000137231446756
mat T_V[3,5] =  9.41131495358e-06
mat T_V[3,6] =  -.000211761724297
mat T_V[3,7] =  4.39668885315e-06
mat T_V[3,8] = -1.41023160456e-06
mat T_V[3,9] = -6.19129559506e-07
mat T_V[4,1] =  .0000381415678833
mat T_V[4,2] =  2.91573861909e-07
mat T_V[4,3] =  .0000137231446756
mat T_V[4,4] =  .0007154740851298
mat T_V[4,5] = -.0000266673888608
mat T_V[4,6] =  .0002371612447248
mat T_V[4,7] = -.0000194736018315
mat T_V[4,8] =  7.40741242409e-06
mat T_V[4,9] =  2.26746196107e-06
mat T_V[5,1] = -.0000203979974521
mat T_V[5,2] = -4.56597553587e-06
mat T_V[5,3] =  9.41131495358e-06
mat T_V[5,4] = -.0000266673888608
mat T_V[5,5] =  .0007574866416573
mat T_V[5,6] = -.0002001287615774
mat T_V[5,7] =  .0000280770067926
mat T_V[5,8] =  -.000010994331113
mat T_V[5,9] = -3.43127415157e-06
mat T_V[6,1] = -.0008047995001103
mat T_V[6,2] = -.0023313413570737
mat T_V[6,3] =  -.000211761724297
mat T_V[6,4] =  .0002371612447248
mat T_V[6,5] = -.0002001287615774
mat T_V[6,6] =  .0702218138493397
mat T_V[6,7] = -.0000278725462354
mat T_V[6,8] = -.0000394949761767
mat T_V[6,9] =   -.00001216567775
mat T_V[7,1] = -.0000245192887961
mat T_V[7,2] = -2.73206862120e-06
mat T_V[7,3] =  4.39668885315e-06
mat T_V[7,4] = -.0000194736018315
mat T_V[7,5] =  .0000280770067926
mat T_V[7,6] = -.0000278725462354
mat T_V[7,7] =   .000448500045115
mat T_V[7,8] =  .0001269967380345
mat T_V[7,9] =  .0000222420175756
mat T_V[8,1] =  9.95239885772e-06
mat T_V[8,2] =  8.25076198477e-07
mat T_V[8,3] = -1.41023160456e-06
mat T_V[8,4] =  7.40741242409e-06
mat T_V[8,5] =  -.000010994331113
mat T_V[8,6] = -.0000394949761767
mat T_V[8,7] =  .0001269967380345
mat T_V[8,8] =  .0002327314964656
mat T_V[8,9] = -.0000213071042735
mat T_V[9,1] =  3.00953179885e-06
mat T_V[9,2] =  1.22734398975e-07
mat T_V[9,3] = -6.19129559506e-07
mat T_V[9,4] =  2.26746196107e-06
mat T_V[9,5] = -3.43127415157e-06
mat T_V[9,6] =   -.00001216567775
mat T_V[9,7] =  .0000222420175756
mat T_V[9,8] = -.0000213071042735
mat T_V[9,9] =  .0000599670800285
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[1,2] = -6.742606420838285
mat T_rcsrmat_1[1,3] =  -3.55544856779441
mat T_rcsrmat_1[2,2] =  1.510744869937544
mat T_rcsrmat_1[2,3] =  .8562440553834105
mat T_rcsrmat_1[3,3] =  .0564600884291515
mat T_rcsrmat_1[4,1] =  2.026577982062658
mat T_rcsrmat_1[4,2] = -30.08658057770258
mat T_rcsrmat_1[4,3] =  -15.5275355165639
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] = -2.92930584379e-06
mat T_gradient[1,2] =  -.000337473021886
mat T_gradient[1,3] = -2.51304785585e-06
mat T_gradient[1,4] =  .0000107194804764
mat T_gradient[1,5] = -.0000129092827629
mat T_gradient[1,6] = -.0000116401620042
mat T_gradient[1,7] =  .0001903263379313
mat T_gradient[1,8] =  .0005300765013116
mat T_gradient[1,9] = -.0004017122443233
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient

//============================================================================//
// user-defined BOUNDARY knots PH model

uhtred (stime trt bmi x1 x2 x3, family(rp, knots(-3.5 1 1.5 2) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3, family(rp, knots(-3.5 1 1.5 2) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.5 1 1.5 2"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5770.540651590882) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] = -.5151777682780363
mat T_b[1,2] = -.0535720215855634
mat T_b[1,3] =  .0887722448408774
mat T_b[1,4] = -.4125777134629683
mat T_b[1,5] =  .5356764428624713
mat T_b[1,6] =  .2396267866949445
mat T_b[1,7] =  .7582992465605269
mat T_b[1,8] = -.0032400707879088
mat T_b[1,9] =  .0196690009653912
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =   .003015112729449
mat T_V[1,2] = -.0000124855307004
mat T_V[1,3] = -1.18030742386e-06
mat T_V[1,4] =  .0000381415678833
mat T_V[1,5] = -.0000203979974521
mat T_V[1,6] = -.0008047995001103
mat T_V[1,7] = -.0000245192887961
mat T_V[1,8] =  9.95239885772e-06
mat T_V[1,9] =  3.00953179885e-06
mat T_V[2,1] = -.0000124855307004
mat T_V[2,2] =  .0000790287397036
mat T_V[2,3] =  5.11275693772e-06
mat T_V[2,4] =  2.91573861909e-07
mat T_V[2,5] = -4.56597553587e-06
mat T_V[2,6] = -.0023313413570737
mat T_V[2,7] = -2.73206862120e-06
mat T_V[2,8] =  8.25076198477e-07
mat T_V[2,9] =  1.22734398975e-07
mat T_V[3,1] = -1.18030742386e-06
mat T_V[3,2] =  5.11275693772e-06
mat T_V[3,3] =  .0007442340932892
mat T_V[3,4] =  .0000137231446756
mat T_V[3,5] =  9.41131495358e-06
mat T_V[3,6] =  -.000211761724297
mat T_V[3,7] =  4.39668885315e-06
mat T_V[3,8] = -1.41023160456e-06
mat T_V[3,9] = -6.19129559506e-07
mat T_V[4,1] =  .0000381415678833
mat T_V[4,2] =  2.91573861909e-07
mat T_V[4,3] =  .0000137231446756
mat T_V[4,4] =  .0007154740851298
mat T_V[4,5] = -.0000266673888608
mat T_V[4,6] =  .0002371612447248
mat T_V[4,7] = -.0000194736018315
mat T_V[4,8] =  7.40741242409e-06
mat T_V[4,9] =  2.26746196107e-06
mat T_V[5,1] = -.0000203979974521
mat T_V[5,2] = -4.56597553587e-06
mat T_V[5,3] =  9.41131495358e-06
mat T_V[5,4] = -.0000266673888608
mat T_V[5,5] =  .0007574866416573
mat T_V[5,6] = -.0002001287615774
mat T_V[5,7] =  .0000280770067926
mat T_V[5,8] =  -.000010994331113
mat T_V[5,9] = -3.43127415157e-06
mat T_V[6,1] = -.0008047995001103
mat T_V[6,2] = -.0023313413570737
mat T_V[6,3] =  -.000211761724297
mat T_V[6,4] =  .0002371612447248
mat T_V[6,5] = -.0002001287615774
mat T_V[6,6] =  .0702218138493397
mat T_V[6,7] = -.0000278725462354
mat T_V[6,8] = -.0000394949761767
mat T_V[6,9] =   -.00001216567775
mat T_V[7,1] = -.0000245192887961
mat T_V[7,2] = -2.73206862120e-06
mat T_V[7,3] =  4.39668885315e-06
mat T_V[7,4] = -.0000194736018315
mat T_V[7,5] =  .0000280770067926
mat T_V[7,6] = -.0000278725462354
mat T_V[7,7] =   .000448500045115
mat T_V[7,8] =  .0001269967380345
mat T_V[7,9] =  .0000222420175756
mat T_V[8,1] =  9.95239885772e-06
mat T_V[8,2] =  8.25076198477e-07
mat T_V[8,3] = -1.41023160456e-06
mat T_V[8,4] =  7.40741242409e-06
mat T_V[8,5] =  -.000010994331113
mat T_V[8,6] = -.0000394949761767
mat T_V[8,7] =  .0001269967380345
mat T_V[8,8] =  .0002327314964656
mat T_V[8,9] = -.0000213071042735
mat T_V[9,1] =  3.00953179885e-06
mat T_V[9,2] =  1.22734398975e-07
mat T_V[9,3] = -6.19129559506e-07
mat T_V[9,4] =  2.26746196107e-06
mat T_V[9,5] = -3.43127415157e-06
mat T_V[9,6] =   -.00001216567775
mat T_V[9,7] =  .0000222420175756
mat T_V[9,8] = -.0000213071042735
mat T_V[9,9] =  .0000599670800285
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[1,2] = -6.742606420838285
mat T_rcsrmat_1[1,3] =  -3.55544856779441
mat T_rcsrmat_1[2,2] =  1.510744869937544
mat T_rcsrmat_1[2,3] =  .8562440553834105
mat T_rcsrmat_1[3,3] =  .0564600884291515
mat T_rcsrmat_1[4,1] =  2.026577982062658
mat T_rcsrmat_1[4,2] = -30.08658057770258
mat T_rcsrmat_1[4,3] =  -15.5275355165639
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] = -2.92930584379e-06
mat T_gradient[1,2] =  -.000337473021886
mat T_gradient[1,3] = -2.51304785585e-06
mat T_gradient[1,4] =  .0000107194804764
mat T_gradient[1,5] = -.0000129092827629
mat T_gradient[1,6] = -.0000116401620042
mat T_gradient[1,7] =  .0001903263379313
mat T_gradient[1,8] =  .0005300765013116
mat T_gradient[1,9] = -.0004017122443233
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient



//============================================================================//
//interaction term in PH - cat#cat
uhtred (stime i.trt bmi x1 x2 x3 i.trt#i.agecat, family(rp, df(3) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime i.trt bmi x1 x2 x3 i.trt#i.agecat, family(rp, df(3) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.988795 1.180037 1.83553 2.302106"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi died i.agecat i.trt stime x1 x2 x3"'
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

assert         e(rank)       == 17
assert         e(N)          == 5000
assert         e(k)          == 20
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 3
assert reldif( e(ll)          , -5767.921852194574) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,20,0)
mat T_b[1,2] = -.3956418143707023
mat T_b[1,3] =  -.053763346868748
mat T_b[1,4] =  .0929874847388738
mat T_b[1,5] = -.4105455114311095
mat T_b[1,6] =  .5351831359735395
mat T_b[1,8] =  .0408815748552386
mat T_b[1,9] = -.0914495742683547
mat T_b[1,10] =   .125518600490574
mat T_b[1,11] =    .04316551333753
mat T_b[1,13] = -.1764056753351244
mat T_b[1,14] = -.1068873814712446
mat T_b[1,15] = -.1122447642506615
mat T_b[1,16] = -.3308562547961565
mat T_b[1,17] =  .2423434373286982
mat T_b[1,18] =  .7583280723451924
mat T_b[1,19] = -.0002308868199841
mat T_b[1,20] =  .0138515964498071
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
mat drop C_b T_b

qui {
mat T_V = J(20,20,0)
mat T_V[2,2] =  .0174991401683504
mat T_V[2,3] = -.0000114844454355
mat T_V[2,4] =  .0000813716186937
mat T_V[2,5] =  .0001883049179755
mat T_V[2,6] = -5.25434259465e-06
mat T_V[2,8] =  .0076524571911134
mat T_V[2,9] =  .0076407420439943
mat T_V[2,10] =  .0076589002582105
mat T_V[2,11] =   .007671253134403
mat T_V[2,13] = -.0098523228149457
mat T_V[2,14] = -.0098295207268675
mat T_V[2,15] = -.0098621726314457
mat T_V[2,16] = -.0098443900940558
mat T_V[2,17] = -.0072477249706814
mat T_V[2,18] = -.0000209773678072
mat T_V[2,19] =  .0000129908942358
mat T_V[2,20] =  8.50722928302e-07
mat T_V[3,2] = -.0000114844454355
mat T_V[3,3] =  .0000793502512754
mat T_V[3,4] =  4.76873129876e-06
mat T_V[3,5] =  5.17347120868e-07
mat T_V[3,6] = -4.56739388714e-06
mat T_V[3,8] = -2.41404702276e-06
mat T_V[3,9] =  8.48151989528e-06
mat T_V[3,10] =  3.02707351260e-07
mat T_V[3,11] =  .0000100901123649
mat T_V[3,13] = -5.41653570810e-06
mat T_V[3,14] =  4.78402383146e-06
mat T_V[3,15] =  8.87045553977e-06
mat T_V[3,16] =  2.65237376398e-06
mat T_V[3,17] = -.0023429092574887
mat T_V[3,18] = -2.74835223530e-06
mat T_V[3,19] =  8.65594082841e-07
mat T_V[3,20] =  1.20488758535e-08
mat T_V[4,2] =  .0000813716186937
mat T_V[4,3] =  4.76873129876e-06
mat T_V[4,4] =  .0007493714334375
mat T_V[4,5] =  .0000140250620044
mat T_V[4,6] =  .0000109497299603
mat T_V[4,8] =  .0000697470019531
mat T_V[4,9] = -.0000545763046329
mat T_V[4,10] =  .0001016035775154
mat T_V[4,11] =  .0000851909927575
mat T_V[4,13] =  -.000088803975597
mat T_V[4,14] =  -.000057186970738
mat T_V[4,15] = -.0001107862706458
mat T_V[4,16] =  .0001231585301238
mat T_V[4,17] = -.0002265533289198
mat T_V[4,18] =  4.60670206809e-06
mat T_V[4,19] = -1.55856198343e-06
mat T_V[4,20] = -4.91588363519e-07
mat T_V[5,2] =  .0001883049179755
mat T_V[5,3] =  5.17347120868e-07
mat T_V[5,4] =  .0000140250620044
mat T_V[5,5] =  .0007167458996986
mat T_V[5,6] = -.0000260012878418
mat T_V[5,8] =  .0000443067204204
mat T_V[5,9] =   .000054916258798
mat T_V[5,10] =  .0000616508516482
mat T_V[5,11] =  .0001140175618738
mat T_V[5,13] = -.0001552929604081
mat T_V[5,14] = -.0000705453756511
mat T_V[5,15] = -.0001741289345878
mat T_V[5,16] = -.0002064233888734
mat T_V[5,17] =  .0001842690315893
mat T_V[5,18] = -.0000191073855683
mat T_V[5,19] =  7.55905619895e-06
mat T_V[5,20] =  1.87313961099e-06
mat T_V[6,2] = -5.25434259465e-06
mat T_V[6,3] = -4.56739388714e-06
mat T_V[6,4] =  .0000109497299603
mat T_V[6,5] = -.0000260012878418
mat T_V[6,6] =    .00075878578015
mat T_V[6,8] = -.0000745088935235
mat T_V[6,9] = -.0000361037384793
mat T_V[6,10] =  -.000091248266086
mat T_V[6,11] =  .0000317120148642
mat T_V[6,13] = -.0000519589658094
mat T_V[6,14] = -.0000576446427742
mat T_V[6,15] = -.0001750881726765
mat T_V[6,16] = -.0000981665299393
mat T_V[6,17] = -.0001502545433779
mat T_V[6,18] =  .0000279464776124
mat T_V[6,19] = -.0000112491979797
mat T_V[6,20] = -2.58131420285e-06
mat T_V[8,2] =  .0076524571911134
mat T_V[8,3] = -2.41404702276e-06
mat T_V[8,4] =  .0000697470019531
mat T_V[8,5] =  .0000443067204204
mat T_V[8,6] = -.0000745088935235
mat T_V[8,8] =  .0110175975062059
mat T_V[8,9] =  .0076349209467269
mat T_V[8,10] =  .0076560625631301
mat T_V[8,11] =  .0076446896416344
mat T_V[8,13] =  -.000011765331673
mat T_V[8,14] = -3.70219176244e-06
mat T_V[8,15] = -3.03272726673e-06
mat T_V[8,16] =  9.82345503895e-06
mat T_V[8,17] = -.0075311300079388
mat T_V[8,18] = -2.03182361248e-06
mat T_V[8,19] =  2.42252691602e-06
mat T_V[8,20] = -4.47399004783e-06
mat T_V[9,2] =  .0076407420439943
mat T_V[9,3] =  8.48151989528e-06
mat T_V[9,4] = -.0000545763046329
mat T_V[9,5] =   .000054916258798
mat T_V[9,6] = -.0000361037384793
mat T_V[9,8] =  .0076349209467269
mat T_V[9,9] =  .0115053781672752
mat T_V[9,10] =  .0076349203287407
mat T_V[9,11] =  .0076357345272178
mat T_V[9,13] = -3.63070240006e-06
mat T_V[9,14] =  1.93472830826e-06
mat T_V[9,15] =  3.56100992587e-06
mat T_V[9,16] = -.0000203151632291
mat T_V[9,17] = -.0078502745030134
mat T_V[9,18] = -1.87496851373e-06
mat T_V[9,19] =  1.67787787822e-06
mat T_V[9,20] = -1.04866577293e-06
mat T_V[10,2] =  .0076589002582105
mat T_V[10,3] =  3.02707351260e-07
mat T_V[10,4] =  .0001016035775154
mat T_V[10,5] =  .0000616508516482
mat T_V[10,6] =  -.000091248266086
mat T_V[10,8] =  .0076560625631301
mat T_V[10,9] =  .0076349203287407
mat T_V[10,10] =   .015416174131799
mat T_V[10,11] =  .0076501048195739
mat T_V[10,13] = -.0000183308354648
mat T_V[10,14] = -6.12342587294e-06
mat T_V[10,15] = -7.73864249475e-06
mat T_V[10,16] =  .0000119925334288
mat T_V[10,17] = -.0076060665566235
mat T_V[10,18] =  4.90152262965e-06
mat T_V[10,19] = -3.62856727324e-06
mat T_V[10,20] = -2.79376249517e-06
mat T_V[11,2] =   .007671253134403
mat T_V[11,3] =  .0000100901123649
mat T_V[11,4] =  .0000851909927575
mat T_V[11,5] =  .0001140175618738
mat T_V[11,6] =  .0000317120148642
mat T_V[11,8] =  .0076446896416344
mat T_V[11,9] =  .0076357345272178
mat T_V[11,10] =  .0076501048195739
mat T_V[11,11] =  .0632195594112993
mat T_V[11,13] =   -.00003719912952
mat T_V[11,14] = -.0000195851600714
mat T_V[11,15] = -.0000464886144121
mat T_V[11,16] = -.0000229458233224
mat T_V[11,17] = -.0079212312477821
mat T_V[11,18] = -6.48702227247e-06
mat T_V[11,19] =  5.53497291183e-06
mat T_V[11,20] = -2.60643059173e-06
mat T_V[13,2] = -.0098523228149457
mat T_V[13,3] = -5.41653570810e-06
mat T_V[13,4] =  -.000088803975597
mat T_V[13,5] = -.0001552929604081
mat T_V[13,6] = -.0000519589658094
mat T_V[13,8] =  -.000011765331673
mat T_V[13,9] = -3.63070240006e-06
mat T_V[13,10] = -.0000183308354648
mat T_V[13,11] =   -.00003719912952
mat T_V[13,13] =  .0157692822083538
mat T_V[13,14] =  .0098296333118654
mat T_V[13,15] =  .0098669291070668
mat T_V[13,16] =  .0098424392049923
mat T_V[13,17] =   .000149542796929
mat T_V[13,18] = -7.53708536134e-06
mat T_V[13,19] = -8.67809482637e-07
mat T_V[13,20] = -9.99545258217e-07
mat T_V[14,2] = -.0098295207268675
mat T_V[14,3] =  4.78402383146e-06
mat T_V[14,4] =  -.000057186970738
mat T_V[14,5] = -.0000705453756511
mat T_V[14,6] = -.0000576446427742
mat T_V[14,8] = -3.70219176244e-06
mat T_V[14,9] =  1.93472830826e-06
mat T_V[14,10] = -6.12342587294e-06
mat T_V[14,11] = -.0000195851600714
mat T_V[14,13] =  .0098296333118654
mat T_V[14,14] =  .0151677201043228
mat T_V[14,15] =  .0098437965824987
mat T_V[14,16] =  .0098232817407615
mat T_V[14,17] = -.0001327920802304
mat T_V[14,18] =  1.17997611220e-07
mat T_V[14,19] = -3.46670354095e-06
mat T_V[14,20] = -1.06858721235e-06
mat T_V[15,2] = -.0098621726314457
mat T_V[15,3] =  8.87045553977e-06
mat T_V[15,4] = -.0001107862706458
mat T_V[15,5] = -.0001741289345878
mat T_V[15,6] = -.0001750881726765
mat T_V[15,8] = -3.03272726673e-06
mat T_V[15,9] =  3.56100992587e-06
mat T_V[15,10] = -7.73864249475e-06
mat T_V[15,11] = -.0000464886144121
mat T_V[15,13] =  .0098669291070668
mat T_V[15,14] =  .0098437965824987
mat T_V[15,15] =   .021810420977936
mat T_V[15,16] =  .0098621673107138
mat T_V[15,17] = -.0002303077529231
mat T_V[15,18] = -.0000102735321997
mat T_V[15,19] = -4.42701651442e-07
mat T_V[15,20] = -3.67891952912e-06
mat T_V[16,2] = -.0098443900940558
mat T_V[16,3] =  2.65237376398e-06
mat T_V[16,4] =  .0001231585301238
mat T_V[16,5] = -.0002064233888734
mat T_V[16,6] = -.0000981665299393
mat T_V[16,8] =  9.82345503895e-06
mat T_V[16,9] = -.0000203151632291
mat T_V[16,10] =  .0000119925334288
mat T_V[16,11] = -.0000229458233224
mat T_V[16,13] =  .0098424392049923
mat T_V[16,14] =  .0098232817407615
mat T_V[16,15] =  .0098621673107138
mat T_V[16,16] =  .1008105653519338
mat T_V[16,17] = -.0001080021740539
mat T_V[16,18] = -.0000141293723342
mat T_V[16,19] =  5.97544706334e-07
mat T_V[16,20] = -5.46919850615e-06
mat T_V[17,2] = -.0072477249706814
mat T_V[17,3] = -.0023429092574887
mat T_V[17,4] = -.0002265533289198
mat T_V[17,5] =  .0001842690315893
mat T_V[17,6] = -.0001502545433779
mat T_V[17,8] = -.0075311300079388
mat T_V[17,9] = -.0078502745030134
mat T_V[17,10] = -.0076060665566235
mat T_V[17,11] = -.0079212312477821
mat T_V[17,13] =   .000149542796929
mat T_V[17,14] = -.0001327920802304
mat T_V[17,15] = -.0002303077529231
mat T_V[17,16] = -.0001080021740539
mat T_V[17,17] =  .0769970038108876
mat T_V[17,18] = -.0000266620005464
mat T_V[17,19] = -.0000430592888625
mat T_V[17,20] = -3.88688838372e-06
mat T_V[18,2] = -.0000209773678072
mat T_V[18,3] = -2.74835223530e-06
mat T_V[18,4] =  4.60670206809e-06
mat T_V[18,5] = -.0000191073855683
mat T_V[18,6] =  .0000279464776124
mat T_V[18,8] = -2.03182361248e-06
mat T_V[18,9] = -1.87496851373e-06
mat T_V[18,10] =  4.90152262965e-06
mat T_V[18,11] = -6.48702227247e-06
mat T_V[18,13] = -7.53708536134e-06
mat T_V[18,14] =  1.17997611220e-07
mat T_V[18,15] = -.0000102735321997
mat T_V[18,16] = -.0000141293723342
mat T_V[18,17] = -.0000266620005464
mat T_V[18,18] =  .0004511792854032
mat T_V[18,19] =  .0001357453514948
mat T_V[18,20] =  9.85798621555e-06
mat T_V[19,2] =  .0000129908942358
mat T_V[19,3] =  8.65594082841e-07
mat T_V[19,4] = -1.55856198343e-06
mat T_V[19,5] =  7.55905619895e-06
mat T_V[19,6] = -.0000112491979797
mat T_V[19,8] =  2.42252691602e-06
mat T_V[19,9] =  1.67787787822e-06
mat T_V[19,10] = -3.62856727324e-06
mat T_V[19,11] =  5.53497291183e-06
mat T_V[19,13] = -8.67809482637e-07
mat T_V[19,14] = -3.46670354095e-06
mat T_V[19,15] = -4.42701651442e-07
mat T_V[19,16] =  5.97544706334e-07
mat T_V[19,17] = -.0000430592888625
mat T_V[19,18] =  .0001357453514948
mat T_V[19,19] =  .0002429329248651
mat T_V[19,20] = -.0000274277547687
mat T_V[20,2] =  8.50722928302e-07
mat T_V[20,3] =  1.20488758535e-08
mat T_V[20,4] = -4.91588363519e-07
mat T_V[20,5] =  1.87313961099e-06
mat T_V[20,6] = -2.58131420285e-06
mat T_V[20,8] = -4.47399004783e-06
mat T_V[20,9] = -1.04866577293e-06
mat T_V[20,10] = -2.79376249517e-06
mat T_V[20,11] = -2.60643059173e-06
mat T_V[20,13] = -9.99545258217e-07
mat T_V[20,14] = -1.06858721235e-06
mat T_V[20,15] = -3.67891952912e-06
mat T_V[20,16] = -5.46919850615e-06
mat T_V[20,17] = -3.88688838372e-06
mat T_V[20,18] =  9.85798621555e-06
mat T_V[20,19] = -.0000274277547687
mat T_V[20,20] =  .0000437888398354
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8

mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[1,2] = -8.435695793572606
mat T_rcsrmat_1[1,3] = -3.658962257730614
mat T_rcsrmat_1[2,2] =  1.843410155446673
mat T_rcsrmat_1[2,3] =  .8567385999018583
mat T_rcsrmat_1[3,3] =  .0586748184202542
mat T_rcsrmat_1[4,1] =  2.026577982062658
mat T_rcsrmat_1[4,2] = -38.92980822103527
mat T_rcsrmat_1[4,3] = -16.56962891904038
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,20,0)
mat T_gradient[1,2] = -4.44886072426e-11
mat T_gradient[1,3] = -4.07134237435e-09
mat T_gradient[1,4] = -3.75621184769e-11
mat T_gradient[1,5] =  1.25761241285e-10
mat T_gradient[1,6] = -1.47028118636e-10
mat T_gradient[1,8] = -4.57583450297e-11
mat T_gradient[1,9] = -2.67394960340e-11
mat T_gradient[1,10] = -8.78697808959e-12
mat T_gradient[1,11] = -6.11428963015e-12
mat T_gradient[1,13] = -1.27065788447e-11
mat T_gradient[1,14] = -1.25619410707e-11
mat T_gradient[1,15] = -4.10594128142e-12
mat T_gradient[1,16] = -5.03624919546e-14
mat T_gradient[1,17] = -1.42856657059e-10
mat T_gradient[1,18] =  2.29208527098e-08
mat T_gradient[1,19] =  5.97503193606e-08
mat T_gradient[1,20] = -3.59957102371e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
mat drop C_gradient T_gradient

//============================================================================//
// interaction term in PH - cts cts

uhtred (stime trt bmi c.x1#c.x2 x3, family(rp, df(3) failure(died)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi c.x1#c.x2 x3, family(rp, df(3) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.988795 1.180037 1.83553 2.302106"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi c.x1 c.x2 died stime trt x3"'
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
assert reldif( e(ll)          , -5896.743918118266) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.4924789946927949
mat T_b[1,2] = -.0538330979952552
mat T_b[1,3] =  .0008855040311123
mat T_b[1,4] =  .5205139911584153
mat T_b[1,5] =  .3067393877773399
mat T_b[1,6] =  .7468687450482938
mat T_b[1,7] =  .0042627980130086
mat T_b[1,8] =  .0148616038477548
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:c.x1#c.x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0030119953307758
mat T_V[1,2] = -9.78036160519e-06
mat T_V[1,3] = -.0000190994098222
mat T_V[1,4] =  -.000015118033791
mat T_V[1,5] = -.0008996523848148
mat T_V[1,6] = -.0000202258433948
mat T_V[1,7] =  8.99677473944e-06
mat T_V[1,8] =  1.76580837330e-06
mat T_V[2,1] = -9.78036160519e-06
mat T_V[2,2] =  .0000781110306601
mat T_V[2,3] = -3.83315217802e-06
mat T_V[2,4] = -1.52901017138e-06
mat T_V[2,5] = -.0023063163992274
mat T_V[2,6] = -2.44399925009e-06
mat T_V[2,7] =  8.61134241872e-07
mat T_V[2,8] =  6.90277869444e-08
mat T_V[3,1] = -.0000190994098222
mat T_V[3,2] = -3.83315217802e-06
mat T_V[3,3] =  .0007682756043243
mat T_V[3,4] =  3.24822216792e-06
mat T_V[3,5] =  .0001370988899443
mat T_V[3,6] =  1.42004233860e-06
mat T_V[3,7] = -5.43850301277e-07
mat T_V[3,8] = -5.52266184805e-08
mat T_V[4,1] =  -.000015118033791
mat T_V[4,2] = -1.52901017138e-06
mat T_V[4,3] =  3.24822216792e-06
mat T_V[4,4] =  .0007564947812323
mat T_V[4,5] =  -.000281570046167
mat T_V[4,6] =  .0000234393078203
mat T_V[4,7] = -9.61278165788e-06
mat T_V[4,8] = -2.16009314249e-06
mat T_V[5,1] = -.0008996523848148
mat T_V[5,2] = -.0023063163992274
mat T_V[5,3] =  .0001370988899443
mat T_V[5,4] =  -.000281570046167
mat T_V[5,5] =   .069459539384129
mat T_V[5,6] = -.0000353802217926
mat T_V[5,7] = -.0000425161843822
mat T_V[5,8] = -8.07744043539e-06
mat T_V[6,1] = -.0000202258433948
mat T_V[6,2] = -2.44399925009e-06
mat T_V[6,3] =  1.42004233860e-06
mat T_V[6,4] =  .0000234393078203
mat T_V[6,5] = -.0000353802217926
mat T_V[6,6] =  .0004477235557392
mat T_V[6,7] =  .0001380895742111
mat T_V[6,8] =  .0000102416254696
mat T_V[7,1] =  8.99677473944e-06
mat T_V[7,2] =  8.61134241872e-07
mat T_V[7,3] = -5.43850301277e-07
mat T_V[7,4] = -9.61278165788e-06
mat T_V[7,5] = -.0000425161843822
mat T_V[7,6] =  .0001380895742111
mat T_V[7,7] =  .0002417868805974
mat T_V[7,8] = -.0000280221876481
mat T_V[8,1] =  1.76580837330e-06
mat T_V[8,2] =  6.90277869444e-08
mat T_V[8,3] = -5.52266184805e-08
mat T_V[8,4] = -2.16009314249e-06
mat T_V[8,5] = -8.07744043539e-06
mat T_V[8,6] =  .0000102416254696
mat T_V[8,7] = -.0000280221876481
mat T_V[8,8] =  .0000422858857359
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:c.x1#c.x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:c.x1#c.x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6536567793092846
mat T_rcsrmat_1[1,2] = -8.435695793572606
mat T_rcsrmat_1[1,3] = -3.658962257730614
mat T_rcsrmat_1[2,2] =  1.843410155446673
mat T_rcsrmat_1[2,3] =  .8567385999018583
mat T_rcsrmat_1[3,3] =  .0586748184202542
mat T_rcsrmat_1[4,1] =  2.026577982062658
mat T_rcsrmat_1[4,2] = -38.92980822103527
mat T_rcsrmat_1[4,3] = -16.56962891904038
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] = -3.50799216065e-11
mat T_gradient[1,2] = -3.32807748027e-09
mat T_gradient[1,3] =  7.01977275827e-12
mat T_gradient[1,4] = -1.18770542473e-10
mat T_gradient[1,5] = -1.15975454934e-10
mat T_gradient[1,6] =  2.30987047874e-08
mat T_gradient[1,7] =  5.99913448709e-08
mat T_gradient[1,8] = -3.59491867975e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:c.x1#c.x2 xb1:x3 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient
		
		
//============================================================================//
// rp nonPH model - df1

//sim nonPH data
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
egen agecat=cut(age), at(0,50,55,60,65,100) icodes

gen t0 = runiform() * 2

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
		cov(trt -0.5 bmi -0.05 x1 0.1 x2 -0.4 x3 0.5) ///
		tde(trt 0.01) tdefunc(log({t}))	///
		maxt(10) 

//model
uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(1) log orthog), family(rp, df(1) failure(died))) 	

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(1) log orthog), family(rp, df(1) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-3.948971 2.302585"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.948971 2.301309"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi c.trt died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5786.684804478459) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5110810019875132
mat T_b[1,2] =  -.053126377722847
mat T_b[1,3] =  .0870718097259109
mat T_b[1,4] = -.4156085485858795
mat T_b[1,5] =  .5368070381086772
mat T_b[1,6] =  .2285884237108362
mat T_b[1,7] =  .0127008723190856
mat T_b[1,8] =  .7473732667807018
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0030953016685342
mat T_V[1,2] = -.0000124398076723
mat T_V[1,3] = -5.52700541709e-07
mat T_V[1,4] =  .0000350939666563
mat T_V[1,5] =  -.000015133452605
mat T_V[1,6] = -.0008354593006785
mat T_V[1,7] = -.0003768404343652
mat T_V[1,8] =  .0001107563231214
mat T_V[2,1] = -.0000124398076723
mat T_V[2,2] =  .0000788096467773
mat T_V[2,3] =  5.12976874850e-06
mat T_V[2,4] =  2.30673509505e-07
mat T_V[2,5] = -4.58690096755e-06
mat T_V[2,6] = -.0023251624749239
mat T_V[2,7] = -1.51352409499e-06
mat T_V[2,8] = -2.78081512358e-06
mat T_V[3,1] = -5.52700541709e-07
mat T_V[3,2] =  5.12976874850e-06
mat T_V[3,3] =   .000742144637881
mat T_V[3,4] =  .0000137200763818
mat T_V[3,5] =  9.28438046888e-06
mat T_V[3,6] = -.0002110700634643
mat T_V[3,7] = -1.07885241693e-06
mat T_V[3,8] =  6.06200190622e-06
mat T_V[4,1] =  .0000350939666563
mat T_V[4,2] =  2.30673509505e-07
mat T_V[4,3] =  .0000137200763818
mat T_V[4,4] =  .0007133427175258
mat T_V[4,5] = -.0000263600058623
mat T_V[4,6] =  .0002416071801877
mat T_V[4,7] =  5.16896767134e-06
mat T_V[4,8] = -.0000275147809765
mat T_V[5,1] =  -.000015133452605
mat T_V[5,2] = -4.58690096755e-06
mat T_V[5,3] =  9.28438046888e-06
mat T_V[5,4] = -.0000263600058623
mat T_V[5,5] =  .0007546612856878
mat T_V[5,6] = -.0002011986272068
mat T_V[5,7] = -.0000145423668311
mat T_V[5,8] =  .0000425414591869
mat T_V[6,1] = -.0008354593006785
mat T_V[6,2] = -.0023251624749239
mat T_V[6,3] = -.0002110700634643
mat T_V[6,4] =  .0002416071801877
mat T_V[6,5] = -.0002011986272068
mat T_V[6,6] =  .0700551596455408
mat T_V[6,7] =  .0001652747421438
mat T_V[6,8] = -.0000584844455981
mat T_V[7,1] = -.0003768404343652
mat T_V[7,2] = -1.51352409499e-06
mat T_V[7,3] = -1.07885241693e-06
mat T_V[7,4] =  5.16896767134e-06
mat T_V[7,5] = -.0000145423668311
mat T_V[7,6] =  .0001652747421438
mat T_V[7,7] =  .0015130587425177
mat T_V[7,8] = -.0005755730912234
mat T_V[8,1] =  .0001107563231214
mat T_V[8,2] = -2.78081512358e-06
mat T_V[8,3] =  6.06200190622e-06
mat T_V[8,4] = -.0000275147809765
mat T_V[8,5] =  .0000425414591869
mat T_V[8,6] = -.0000584844455981
mat T_V[8,7] = -.0005755730912234
mat T_V[8,8] =  .0005781086756115
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1"'
mat drop C_V T_V

qui {
mat T_rmat_1_6_2 = J(2,2,0)
mat T_rmat_1_6_2[1,1] =  .6530109699495861
mat T_rmat_1_6_2[2,1] =  2.026098456985472
mat T_rmat_1_6_2[2,2] =                  1
}
matrix C_rmat_1_6_2 = e(rmat_1_6_2)
assert mreldif( C_rmat_1_6_2 , T_rmat_1_6_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_6_2'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rmat_1_6_2'"' `"c1 c2"'
mat drop C_rmat_1_6_2 T_rmat_1_6_2

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .6530109699495861
mat T_rcsrmat_1[2,1] =  2.026098456985472
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] = -1.43187987196e-08
mat T_gradient[1,2] = -1.58633856895e-06
mat T_gradient[1,3] = -1.28272872010e-08
mat T_gradient[1,4] =  5.05825903198e-08
mat T_gradient[1,5] = -5.96113961319e-08
mat T_gradient[1,6] = -5.49978114890e-08
mat T_gradient[1,7] = -1.56280010287e-09
mat T_gradient[1,8] =  3.51212383950e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1"'
mat drop C_gradient T_gradient


//============================================================================//
//RP nonPH model - 3 df
uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(1) log orthog), family(rp, df(3) failure(died))) 

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(1) log orthog), family(rp, df(3) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-3.948971 2.302585"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.948971 1.184131 1.832773 2.301309"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi c.trt died stime trt x1 x2 x3"'
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
assert reldif( e(ll)          , -5784.111202172127) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,10,0)
mat T_b[1,1] = -.5108055399694648
mat T_b[1,2] = -.0530650298966401
mat T_b[1,3] =   .086817423516179
mat T_b[1,4] = -.4144731659313165
mat T_b[1,5] =  .5351276258204217
mat T_b[1,6] =  .2249816301410384
mat T_b[1,7] =  .0164615402329426
mat T_b[1,8] =  .7528560022899016
mat T_b[1,9] =  .0004319194503738
mat T_b[1,10] =   .014511808684217
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(10,10,0)
mat T_V[1,1] =  .0030952942063676
mat T_V[1,2] = -.0000122864560126
mat T_V[1,3] = -5.85160725453e-07
mat T_V[1,4] =  .0000353638994943
mat T_V[1,5] = -.0000148892385709
mat T_V[1,6] = -.0008406335069469
mat T_V[1,7] =  -.000375886712645
mat T_V[1,8] =  .0001145822818169
mat T_V[1,9] =  7.33827959835e-06
mat T_V[1,10] = -7.24097544690e-07
mat T_V[2,1] = -.0000122864560126
mat T_V[2,2] =  .0000788301630957
mat T_V[2,3] =  5.07757005165e-06
mat T_V[2,4] =  3.08332633591e-07
mat T_V[2,5] = -4.62234601488e-06
mat T_V[2,6] = -.0023258920412754
mat T_V[2,7] = -1.46014612243e-06
mat T_V[2,8] = -2.20969605647e-06
mat T_V[2,9] =  8.58571600203e-07
mat T_V[2,10] =  1.33996803467e-08
mat T_V[3,1] = -5.85160725453e-07
mat T_V[3,2] =  5.07757005165e-06
mat T_V[3,3] =  .0007420876573538
mat T_V[3,4] =  .0000136922088864
mat T_V[3,5] =  9.49445010025e-06
mat T_V[3,6] =  -.000209330100484
mat T_V[3,7] = -1.34371650504e-06
mat T_V[3,8] =  4.91235244749e-06
mat T_V[3,9] = -1.48689146856e-06
mat T_V[3,10] = -4.45357741607e-07
mat T_V[4,1] =  .0000353638994943
mat T_V[4,2] =  3.08332633591e-07
mat T_V[4,3] =  .0000136922088864
mat T_V[4,4] =  .0007135471975347
mat T_V[4,5] = -.0000267024576777
mat T_V[4,6] =  .0002380963188367
mat T_V[4,7] =  6.15577828848e-06
mat T_V[4,8] = -.0000216987604027
mat T_V[4,9] =  7.67410550138e-06
mat T_V[4,10] =  1.85293808873e-06
mat T_V[5,1] = -.0000148892385709
mat T_V[5,2] = -4.62234601488e-06
mat T_V[5,3] =  9.49445010025e-06
mat T_V[5,4] = -.0000267024576777
mat T_V[5,5] =  .0007553681694915
mat T_V[5,6] = -.0001987124865245
mat T_V[5,7] = -.0000162102159998
mat T_V[5,8] =  .0000339584140971
mat T_V[5,9] = -.0000114695865121
mat T_V[5,10] = -2.74075485800e-06
mat T_V[6,1] = -.0008406335069469
mat T_V[6,2] = -.0023258920412754
mat T_V[6,3] =  -.000209330100484
mat T_V[6,4] =  .0002380963188367
mat T_V[6,5] = -.0001987124865245
mat T_V[6,6] =  .0700823777885477
mat T_V[6,7] =  .0001608131766956
mat T_V[6,8] = -.0000869288622969
mat T_V[6,9] = -.0000409029146607
mat T_V[6,10] = -5.50510820957e-06
mat T_V[7,1] =  -.000375886712645
mat T_V[7,2] = -1.46014612243e-06
mat T_V[7,3] = -1.34371650504e-06
mat T_V[7,4] =  6.15577828848e-06
mat T_V[7,5] = -.0000162102159998
mat T_V[7,6] =  .0001608131766956
mat T_V[7,7] =  .0015151307873975
mat T_V[7,8] = -.0005629755813107
mat T_V[7,9] =  .0000125128501111
mat T_V[7,10] =  9.38488797486e-06
mat T_V[8,1] =  .0001145822818169
mat T_V[8,2] = -2.20969605647e-06
mat T_V[8,3] =  4.91235244749e-06
mat T_V[8,4] = -.0000216987604027
mat T_V[8,5] =  .0000339584140971
mat T_V[8,6] = -.0000869288622969
mat T_V[8,7] = -.0005629755813107
mat T_V[8,8] =  .0006599567939796
mat T_V[8,9] =   .000131130242408
mat T_V[8,10] =  6.54121977806e-06
mat T_V[9,1] =  7.33827959835e-06
mat T_V[9,2] =  8.58571600203e-07
mat T_V[9,3] = -1.48689146856e-06
mat T_V[9,4] =  7.67410550138e-06
mat T_V[9,5] = -.0000114695865121
mat T_V[9,6] = -.0000409029146607
mat T_V[9,7] =  .0000125128501111
mat T_V[9,8] =   .000131130242408
mat T_V[9,9] =  .0002427201037381
mat T_V[9,10] = -.0000274753663788
mat T_V[10,1] = -7.24097544690e-07
mat T_V[10,2] =  1.33996803467e-08
mat T_V[10,3] = -4.45357741607e-07
mat T_V[10,4] =  1.85293808873e-06
mat T_V[10,5] = -2.74075485800e-06
mat T_V[10,6] = -5.50510820957e-06
mat T_V[10,7] =  9.38488797486e-06
mat T_V[10,8] =  6.54121977806e-06
mat T_V[10,9] = -.0000274753663788
mat T_V[10,10] =  .0000438843829879
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rmat_1_6_2 = J(2,2,0)
mat T_rmat_1_6_2[1,1] =  .6530109699495861
mat T_rmat_1_6_2[2,1] =  2.026098456985472
mat T_rmat_1_6_2[2,2] =                  1
}
matrix C_rmat_1_6_2 = e(rmat_1_6_2)
assert mreldif( C_rmat_1_6_2 , T_rmat_1_6_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_6_2'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rmat_1_6_2'"' `"c1 c2"'
mat drop C_rmat_1_6_2 T_rmat_1_6_2

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6530109699495861
mat T_rcsrmat_1[1,2] = -8.323745356845835
mat T_rcsrmat_1[1,3] = -3.641373298096312
mat T_rcsrmat_1[2,2] =  1.825614318150785
mat T_rcsrmat_1[2,3] =  .8554267137907268
mat T_rcsrmat_1[3,3] =  .0583013797329871
mat T_rcsrmat_1[4,1] =  2.026098456985472
mat T_rcsrmat_1[4,2] = -38.24004846936445
mat T_rcsrmat_1[4,3] = -16.41651361041069
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,10,0)
mat T_gradient[1,1] = -4.42879725948e-11
mat T_gradient[1,2] = -6.39255703927e-09
mat T_gradient[1,3] = -5.72379279119e-11
mat T_gradient[1,4] =  1.93247642132e-10
mat T_gradient[1,5] = -2.21298069086e-10
mat T_gradient[1,6] = -2.24207187538e-10
mat T_gradient[1,7] =  8.70347900347e-09
mat T_gradient[1,8] =  2.78514287397e-08
mat T_gradient[1,9] =  6.83932134884e-08
mat T_gradient[1,10] = -4.67634382736e-08
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient


// orthog matrix certs
qui {
	mat A = J(4,4,0)
	mat A[1,1] = .79717296
	mat A[1,2] = -12.329397
	mat A[1,3] = -5.0923271
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = 2.6519554
	mat A[2,3] = 1.1812266
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .09176479 
	mat A[3,4] = 0

	mat A[4,1] = 1.8714958
	mat A[4,2] = -47.20906
	mat A[4,3] = -19.07834
	mat A[4,4] = 1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-0

qui {
	mat B = J(2,2,0)
	mat B[1,1] = .79717296
	mat B[1,2] = 0

	mat B[2,1] = 1.8714958 
	mat B[2,2] = 1
}
assert mreldif( e(rmat_1_6_2) , B ) < 1E-0
matrix drop A B


//============================================================================//
// non-PH model >1 df in tde

uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(2) log orthog event), family(rp, df(3) failure(died))) 	

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, df(2) log orthog event), family(rp, df(3) failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-3.948971 1.533659 2.301309"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-3.948971 1.184131 1.832773 2.301309"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi c.trt died stime trt x1 x2 x3"'
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
assert         e(k)          == 11
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -5783.967569217161) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,11,0)
mat T_b[1,1] = -.5090297044896325
mat T_b[1,2] = -.0530846675500467
mat T_b[1,3] =  .0868092949001664
mat T_b[1,4] = -.4144549646834603
mat T_b[1,5] =  .5350088746506344
mat T_b[1,6] =  .2247842995056785
mat T_b[1,7] =  .0057063197834344
mat T_b[1,8] = -.0160753750929419
mat T_b[1,9] =  .7567963839078035
mat T_b[1,10] =  .0065754799062569
mat T_b[1,11] =  .0147884328822219
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(11,11,0)
mat T_V[1,1] =  .0031061281602259
mat T_V[1,2] = -.0000124471359292
mat T_V[1,3] = -5.78902058642e-07
mat T_V[1,4] =  .0000351720980054
mat T_V[1,5] = -.0000152211391529
mat T_V[1,6] = -.0008415347937852
mat T_V[1,7] = -.0004400602259618
mat T_V[1,8] = -.0000964080041412
mat T_V[1,9] =  .0001420819520509
mat T_V[1,10] =  .0000496877571728
mat T_V[1,11] =  8.59128883838e-07
mat T_V[2,1] = -.0000124471359292
mat T_V[2,2] =  .0000788329468562
mat T_V[2,3] =  5.07543336430e-06
mat T_V[2,4] =  3.24187000198e-07
mat T_V[2,5] = -4.63515358818e-06
mat T_V[2,6] = -.0023258891101427
mat T_V[2,7] = -7.68782036573e-07
mat T_V[2,8] =  1.11860838000e-06
mat T_V[2,9] = -2.46074017241e-06
mat T_V[2,10] =  4.33638145949e-07
mat T_V[2,11] = -2.93978712356e-09
mat T_V[3,1] = -5.78902058642e-07
mat T_V[3,2] =  5.07543336430e-06
mat T_V[3,3] =  .0007420803720971
mat T_V[3,4] =  .0000137002469157
mat T_V[3,5] =  9.51051179032e-06
mat T_V[3,6] = -.0002092730222164
mat T_V[3,7] = -9.38989781704e-07
mat T_V[3,8] =  4.12073512600e-07
mat T_V[3,9] =  4.76878453451e-06
mat T_V[3,10] = -1.63674733551e-06
mat T_V[3,11] = -4.50261054659e-07
mat T_V[4,1] =  .0000351720980054
mat T_V[4,2] =  3.24187000198e-07
mat T_V[4,3] =  .0000137002469157
mat T_V[4,4] =  .0007135768234068
mat T_V[4,5] = -.0000267274843358
mat T_V[4,6] =  .0002377309966909
mat T_V[4,7] =  5.03204929610e-06
mat T_V[4,8] = -8.60819445874e-07
mat T_V[4,9] = -.0000213086770212
mat T_V[4,10] =  7.95329874996e-06
mat T_V[4,11] =  1.85846453035e-06
mat T_V[5,1] = -.0000152211391529
mat T_V[5,2] = -4.63515358818e-06
mat T_V[5,3] =  9.51051179032e-06
mat T_V[5,4] = -.0000267274843358
mat T_V[5,5] =  .0007553523249043
mat T_V[5,6] = -.0001981715132763
mat T_V[5,7] = -.0000111183944105
mat T_V[5,8] =  6.41919609689e-06
mat T_V[5,9] =  .0000321281482583
mat T_V[5,10] = -.0000138360422636
mat T_V[5,11] = -2.82415401396e-06
mat T_V[6,1] = -.0008415347937852
mat T_V[6,2] = -.0023258891101427
mat T_V[6,3] = -.0002092730222164
mat T_V[6,4] =  .0002377309966909
mat T_V[6,5] = -.0001981715132763
mat T_V[6,6] =  .0700825591075442
mat T_V[6,7] =  .0001715081039564
mat T_V[6,8] =  .0000122775727261
mat T_V[6,9] = -.0000923808866074
mat T_V[6,10] = -.0000476394528346
mat T_V[6,11] = -5.74530147128e-06
mat T_V[7,1] = -.0004400602259618
mat T_V[7,2] = -7.68782036573e-07
mat T_V[7,3] = -9.38989781704e-07
mat T_V[7,4] =  5.03204929610e-06
mat T_V[7,5] = -.0000111183944105
mat T_V[7,6] =  .0001715081039564
mat T_V[7,7] =  .0019027695336156
mat T_V[7,8] =  .0005898710044534
mat T_V[7,9] = -.0007245948529272
mat T_V[7,10] = -.0002365236374354
mat T_V[7,11] =  2.97267935992e-07
mat T_V[8,1] = -.0000964080041412
mat T_V[8,2] =  1.11860838000e-06
mat T_V[8,3] =  4.12073512600e-07
mat T_V[8,4] = -8.60819445874e-07
mat T_V[8,5] =  6.41919609689e-06
mat T_V[8,6] =  .0000122775727261
mat T_V[8,7] =  .0005898710044534
mat T_V[8,8] =  .0008929297878447
mat T_V[8,9] = -.0002298025445344
mat T_V[8,10] = -.0003532550407536
mat T_V[8,11] = -.0000128679958288
mat T_V[9,1] =  .0001420819520509
mat T_V[9,2] = -2.46074017241e-06
mat T_V[9,3] =  4.76878453451e-06
mat T_V[9,4] = -.0000213086770212
mat T_V[9,5] =  .0000321281482583
mat T_V[9,6] = -.0000923808866074
mat T_V[9,7] = -.0007245948529272
mat T_V[9,8] = -.0002298025445344
mat T_V[9,9] =  .0007266608647824
mat T_V[9,10] =  .0002277259270977
mat T_V[9,11] =  .0000100576174553
mat T_V[10,1] =  .0000496877571728
mat T_V[10,2] =  4.33638145949e-07
mat T_V[10,3] = -1.63674733551e-06
mat T_V[10,4] =  7.95329874996e-06
mat T_V[10,5] = -.0000138360422636
mat T_V[10,6] = -.0000476394528346
mat T_V[10,7] = -.0002365236374354
mat T_V[10,8] = -.0003532550407536
mat T_V[10,9] =  .0002277259270977
mat T_V[10,10] =  .0003821219448046
mat T_V[10,11] = -.0000224170060769
mat T_V[11,1] =  8.59128883838e-07
mat T_V[11,2] = -2.93978712356e-09
mat T_V[11,3] = -4.50261054659e-07
mat T_V[11,4] =  1.85846453035e-06
mat T_V[11,5] = -2.82415401396e-06
mat T_V[11,6] = -5.74530147128e-06
mat T_V[11,7] =  2.97267935992e-07
mat T_V[11,8] = -.0000128679958288
mat T_V[11,9] =  .0000100576174553
mat T_V[11,10] = -.0000224170060769
mat T_V[11,11] =   .000044040028345
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rmat_1_6_2 = J(3,3,0)
mat T_rmat_1_6_2[1,1] =  .6530109699495861
mat T_rmat_1_6_2[1,2] = -5.876622008643532
mat T_rmat_1_6_2[2,2] =  1.347170758724998
mat T_rmat_1_6_2[3,1] =  2.026098456985472
mat T_rmat_1_6_2[3,2] = -26.67741778371033
mat T_rmat_1_6_2[3,3] =                  1
}
matrix C_rmat_1_6_2 = e(rmat_1_6_2)
assert mreldif( C_rmat_1_6_2 , T_rmat_1_6_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_6_2'"' `"r1 r2 r3"'
_assert_streq `"`: colfullnames C_rmat_1_6_2'"' `"c1 c2 c3"'
mat drop C_rmat_1_6_2 T_rmat_1_6_2

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6530109699495861
mat T_rcsrmat_1[1,2] = -8.323745356845835
mat T_rcsrmat_1[1,3] = -3.641373298096312
mat T_rcsrmat_1[2,2] =  1.825614318150785
mat T_rcsrmat_1[2,3] =  .8554267137907268
mat T_rcsrmat_1[3,3] =  .0583013797329871
mat T_rcsrmat_1[4,1] =  2.026098456985472
mat T_rcsrmat_1[4,2] = -38.24004846936445
mat T_rcsrmat_1[4,3] = -16.41651361041069
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,11,0)
mat T_gradient[1,1] = -2.30345089093e-07
mat T_gradient[1,2] = -.0000118676801543
mat T_gradient[1,3] =  1.30126447188e-08
mat T_gradient[1,4] = -6.41832438802e-08
mat T_gradient[1,5] =  1.06498117534e-07
mat T_gradient[1,6] = -3.89109463910e-07
mat T_gradient[1,7] =  3.41810816925e-08
mat T_gradient[1,8] =  3.23362582834e-07
mat T_gradient[1,9] = -4.26881568809e-09
mat T_gradient[1,10] =  3.40912919028e-07
mat T_gradient[1,11] = -1.45150248278e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient
	


//============================================================================//
// user-defined knots (nonPH RP)

uhtred (stime trt bmi x1 x2 x3 c.trt#rcs(stime, knots(-4 1.5 2.3) log orthog), family(rp, knots(-4 1 1.5 2) failure(died))) 

assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_6_2)'"'       == `"-4 1.5 2.3"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-4 1 1.5 2"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"bmi c.trt died stime trt x1 x2 x3"'
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
assert         e(k)          == 11
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -5783.156786091448) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,11,0)
mat T_b[1,1] = -.5088128235502523
mat T_b[1,2] = -.0530633072373552
mat T_b[1,3] =  .0867777351540062
mat T_b[1,4] = -.4144468961214485
mat T_b[1,5] =  .5349608014305245
mat T_b[1,6] =  .2239229419479673
mat T_b[1,7] =  .0057428551224929
mat T_b[1,8] = -.0161740479473846
mat T_b[1,9] =  .7572573383190734
mat T_b[1,10] =  .0041483604972014
mat T_b[1,11] =  .0205506025553737
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(11,11,0)
mat T_V[1,1] =  .0031054783995551
mat T_V[1,2] = -.0000124353279362
mat T_V[1,3] = -5.03611537484e-07
mat T_V[1,4] =  .0000351454607059
mat T_V[1,5] = -.0000150198927439
mat T_V[1,6] = -.0008417441908869
mat T_V[1,7] = -.0004369548517222
mat T_V[1,8] = -.0000931596149362
mat T_V[1,9] =  .0001410131431981
mat T_V[1,10] =  .0000481408726667
mat T_V[1,11] =  3.81274743210e-06
mat T_V[2,1] = -.0000124353279362
mat T_V[2,2] =  .0000788190039908
mat T_V[2,3] =  5.09017297518e-06
mat T_V[2,4] =  3.17624756324e-07
mat T_V[2,5] = -4.63824593684e-06
mat T_V[2,6] = -.0023254917581469
mat T_V[2,7] = -8.15621385342e-07
mat T_V[2,8] =  1.06884263961e-06
mat T_V[2,9] = -2.43809616082e-06
mat T_V[2,10] =  4.03889020049e-07
mat T_V[2,11] =  8.73356941993e-08
mat T_V[3,1] = -5.03611537484e-07
mat T_V[3,2] =  5.09017297518e-06
mat T_V[3,3] =  .0007420390085219
mat T_V[3,4] =  .0000136992192317
mat T_V[3,5] =  9.47782043873e-06
mat T_V[3,6] = -.0002097177078766
mat T_V[3,7] = -9.64080339245e-07
mat T_V[3,8] =  3.80994247528e-07
mat T_V[3,9] =  4.76421840887e-06
mat T_V[3,10] = -1.55701315895e-06
mat T_V[3,11] = -6.59796665456e-07
mat T_V[4,1] =  .0000351454607059
mat T_V[4,2] =  3.17624756324e-07
mat T_V[4,3] =  .0000136992192317
mat T_V[4,4] =  .0007136630274617
mat T_V[4,5] =  -.000026750876469
mat T_V[4,6] =  .0002380361111035
mat T_V[4,7] =  5.09969080692e-06
mat T_V[4,8] = -5.53773751020e-07
mat T_V[4,9] = -.0000214964328904
mat T_V[4,10] =  7.59838833193e-06
mat T_V[4,11] =  2.40656354668e-06
mat T_V[5,1] = -.0000150198927439
mat T_V[5,2] = -4.63824593684e-06
mat T_V[5,3] =  9.47782043873e-06
mat T_V[5,4] =  -.000026750876469
mat T_V[5,5] =   .000755386433733
mat T_V[5,6] = -.0001982544717354
mat T_V[5,7] = -.0000114096599233
mat T_V[5,8] =  5.85132203181e-06
mat T_V[5,9] =  .0000324440370147
mat T_V[5,10] = -.0000132268202247
mat T_V[5,11] = -3.84765903278e-06
mat T_V[6,1] = -.0008417441908869
mat T_V[6,2] = -.0023254917581469
mat T_V[6,3] = -.0002097177078766
mat T_V[6,4] =  .0002380361111035
mat T_V[6,5] = -.0001982544717354
mat T_V[6,6] =  .0700712804388588
mat T_V[6,7] =  .0001714934307254
mat T_V[6,8] =  .0000120995418471
mat T_V[6,9] = -.0000925303511562
mat T_V[6,10] = -.0000449564628886
mat T_V[6,11] = -.0000120160257195
mat T_V[7,1] = -.0004369548517222
mat T_V[7,2] = -8.15621385342e-07
mat T_V[7,3] = -9.64080339245e-07
mat T_V[7,4] =  5.09969080692e-06
mat T_V[7,5] = -.0000114096599233
mat T_V[7,6] =  .0001714934307254
mat T_V[7,7] =  .0018869706424658
mat T_V[7,8] =  .0005699966900075
mat T_V[7,9] = -.0007191908815213
mat T_V[7,10] =  -.000230323547784
mat T_V[7,11] = -9.57972749419e-06
mat T_V[8,1] = -.0000931596149362
mat T_V[8,2] =  1.06884263961e-06
mat T_V[8,3] =  3.80994247528e-07
mat T_V[8,4] = -5.53773751020e-07
mat T_V[8,5] =  5.85132203181e-06
mat T_V[8,6] =  .0000120995418471
mat T_V[8,7] =  .0005699966900075
mat T_V[8,8] =  .0008671701241979
mat T_V[8,9] = -.0002224980578225
mat T_V[8,10] = -.0003444484351099
mat T_V[8,11] = -.0000289588393175
mat T_V[9,1] =  .0001410131431981
mat T_V[9,2] = -2.43809616082e-06
mat T_V[9,3] =  4.76421840887e-06
mat T_V[9,4] = -.0000214964328904
mat T_V[9,5] =  .0000324440370147
mat T_V[9,6] = -.0000925303511562
mat T_V[9,7] = -.0007191908815213
mat T_V[9,8] = -.0002224980578225
mat T_V[9,9] =  .0007226337808184
mat T_V[9,10] =  .0002172147686427
mat T_V[9,11] =  .0000272269348046
mat T_V[10,1] =  .0000481408726667
mat T_V[10,2] =  4.03889020049e-07
mat T_V[10,3] = -1.55701315895e-06
mat T_V[10,4] =  7.59838833193e-06
mat T_V[10,5] = -.0000132268202247
mat T_V[10,6] = -.0000449564628886
mat T_V[10,7] =  -.000230323547784
mat T_V[10,8] = -.0003444484351099
mat T_V[10,9] =  .0002172147686427
mat T_V[10,10] =  .0003708965611793
mat T_V[10,11] = -9.29195609783e-06
mat T_V[11,1] =  3.81274743210e-06
mat T_V[11,2] =  8.73356941993e-08
mat T_V[11,3] = -6.59796665456e-07
mat T_V[11,4] =  2.40656354668e-06
mat T_V[11,5] = -3.84765903278e-06
mat T_V[11,6] = -.0000120160257195
mat T_V[11,7] = -9.57972749419e-06
mat T_V[11,8] = -.0000289588393175
mat T_V[11,9] =  .0000272269348046
mat T_V[11,10] = -9.29195609783e-06
mat T_V[11,11] =   .000061855423384
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rmat_1_6_2 = J(3,3,0)
mat T_rmat_1_6_2[1,1] =  .6530109699495861
mat T_rmat_1_6_2[1,2] = -6.185840605415935
mat T_rmat_1_6_2[2,2] =  1.401310623310694
mat T_rmat_1_6_2[3,1] =  2.026098456985472
mat T_rmat_1_6_2[3,2] = -28.25392226794522
mat T_rmat_1_6_2[3,3] =                  1
}
matrix C_rmat_1_6_2 = e(rmat_1_6_2)
assert mreldif( C_rmat_1_6_2 , T_rmat_1_6_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_6_2'"' `"r1 r2 r3"'
_assert_streq `"`: colfullnames C_rmat_1_6_2'"' `"c1 c2 c3"'
mat drop C_rmat_1_6_2 T_rmat_1_6_2

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .6530109699495861
mat T_rcsrmat_1[1,2] = -7.654538235848746
mat T_rcsrmat_1[1,3] = -4.011651034585714
mat T_rcsrmat_1[2,2] =   1.57434089988249
mat T_rcsrmat_1[2,3] =  .8874627806223095
mat T_rcsrmat_1[3,3] =  .0575000626812369
mat T_rcsrmat_1[4,1] =  2.026098456985472
mat T_rcsrmat_1[4,2] = -35.88175669347928
mat T_rcsrmat_1[4,3] = -18.42471283735619
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,11,0)
mat T_gradient[1,1] = -1.51809950073e-07
mat T_gradient[1,2] = -7.56465908192e-06
mat T_gradient[1,3] = -1.02097382736e-08
mat T_gradient[1,4] =  3.38519620715e-08
mat T_gradient[1,5] = -4.02604464928e-08
mat T_gradient[1,6] = -2.51523815289e-07
mat T_gradient[1,7] =  3.88436426283e-06
mat T_gradient[1,8] =  9.14169020718e-06
mat T_gradient[1,9] =  3.89954809104e-06
mat T_gradient[1,10] =  9.57961996552e-06
mat T_gradient[1,11] = -4.85144329558e-06
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:x1 xb1:x2 xb1:x3 xb1:_cons tb1:c.trt#c._rcs1_6_2_1 tb1:c.trt#c._rcs1_6_2_2 tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient
