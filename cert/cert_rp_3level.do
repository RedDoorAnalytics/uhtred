clear 
set seed 7254

set obs 100
gen id1 = _n
gen age = runiform()
gen sd1 = exp(log(0.1))
gen u1 = rnormal(0,sd1)
expand 100
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(0.1))
gen u2 = rnormal(0,sd2)
expand 10
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.02 u1 1 u2 1) maxt(10)

replace id2 = _n

uhtred 	(stime1 trt age M2[id1>id2]@1 M1[id1]@1 ///
	, family(rp, df(1) failure(dead1))) ///
			, 

_assert_streq `"`e(cmdline)'"' `"uhtred (stime1 trt age M2[id1>id2]@1 M1[id1]@1         , family(rp, df(1) failure(dead1)))                         ,"'
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
assert `"`e(knots1)'"'            == `"-6.64526 2.302561"'
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
assert         e(N)          == 100000
assert         e(k)          == 8
assert         e(k_eq)       == 6
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -219867.4829783649) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 3

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.5025782698792253
mat T_b[1,2] =  .0231598212948998
mat T_b[1,3] = -.4947432452337845
mat T_b[1,4] =   1.14860319107674
mat T_b[1,5] =                  1
mat T_b[1,6] =                  1
mat T_b[1,7] = -2.300749422222957
mat T_b[1,8] = -1.756934058281157
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-6
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0000768590456688
mat T_V[1,2] = -5.64569464016e-08
mat T_V[1,3] = -.0000218517963494
mat T_V[1,4] = -.0000263432992817
mat T_V[1,7] = -.0000421157120797
mat T_V[1,8] =  -.001281522512179
mat T_V[2,1] = -5.64569464016e-08
mat T_V[2,2] =  .0014159878843411
mat T_V[2,3] = -.0007266265359429
mat T_V[2,4] =  3.75708379570e-07
mat T_V[2,7] =  1.77139826637e-06
mat T_V[2,8] =  .0000151700397751
mat T_V[3,1] = -.0000218517963494
mat T_V[3,2] = -.0007266265359429
mat T_V[3,3] =  .0005030284565406
mat T_V[3,4] =  -.000011295949492
mat T_V[3,7] = -.0000138813098987
mat T_V[3,8] = -.0002535592683424
mat T_V[4,1] = -.0000263432992817
mat T_V[4,2] =  3.75708379570e-07
mat T_V[4,3] =  -.000011295949492
mat T_V[4,4] =   .000047152954011
mat T_V[4,7] =  .0000539265572805
mat T_V[4,8] =  .0017545216557454
mat T_V[7,1] = -.0000421157120797
mat T_V[7,2] =  1.77139826637e-06
mat T_V[7,3] = -.0000138813098987
mat T_V[7,4] =  .0000539265572805
mat T_V[7,7] =   .006631424160519
mat T_V[7,8] =  .0025874211043052
mat T_V[8,1] =  -.001281522512179
mat T_V[8,2] =  .0000151700397751
mat T_V[8,3] = -.0002535592683424
mat T_V[8,4] =  .0017545216557454
mat T_V[8,7] =  .0025874211043052
mat T_V[8,8] =  .0932614352010597
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-6
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(2,2,0)
mat T_rcsrmat_1[1,1] =  .9552445360205324
mat T_rcsrmat_1[2,1] =  1.489067572024473
mat T_rcsrmat_1[2,2] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-6
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] =  .0000877645011674
mat T_gradient[1,2] = -.0005136161027755
mat T_gradient[1,3] =  .0008011548224963
mat T_gradient[1,4] =  .0149135507323603
mat T_gradient[1,7] =  .0004072959220272
mat T_gradient[1,8] = -.0102175185052282
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-4
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:age xb1:_cons tb1:_rcs1_1 zb1_1:_re_M1 zb1_2:_re_M2 lns1_1:_cons lns2_1:_cons"'
mat drop C_gradient T_gradient

