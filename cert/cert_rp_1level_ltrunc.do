set seed 72549

clear
set obs 500
gen id 	= _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen bmi = rnormal(30,3)

survsim stime died, dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1) maxt(10)

gen t0 = 0
replace t0 = runiform() * 2 //if _n>200

drop if stime<t0

//PH model
uhtred (stime trt bmi age x1 ///
		, family(rp, df(3) failure(died) ltruncated(t0)))

_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi age x1                 , family(rp, df(3) failure(died) ltruncated(t0)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-1.126265 1.352654 1.920142 2.296919"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(ltruncated1)'"'       == `"t0"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died t0"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age bmi died stime t0 trt x1"'
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
assert         e(N)          == 483
assert         e(k)          == 8
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -719.1719298640731) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,8,0)
mat T_b[1,1] = -.4950876558110083
mat T_b[1,2] = -.0820916113445683
mat T_b[1,3] =  .0053576960194915
mat T_b[1,4] = -.0605198694899665
mat T_b[1,5] =  1.393239089930795
mat T_b[1,6] =  .6250287401096962
mat T_b[1,7] = -.0359662562368731
mat T_b[1,8] =  .0054913619625753
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_b T_b

qui {
mat T_V = J(8,8,0)
mat T_V[1,1] =  .0224425423499754
mat T_V[1,2] =   .000031388744325
mat T_V[1,3] = -.0002097024888498
mat T_V[1,4] = -.0002830615447071
mat T_V[1,5] =  .0014014137772944
mat T_V[1,6] = -.0001426195145807
mat T_V[1,7] =  .0001437342833462
mat T_V[1,8] =  4.73869947113e-06
mat T_V[2,1] =   .000031388744325
mat T_V[2,2] =  .0005958087417274
mat T_V[2,3] =  .0000106682740195
mat T_V[2,4] =  .0002905960668185
mat T_V[2,5] = -.0181494524318568
mat T_V[2,6] =  .0000203664407633
mat T_V[2,7] =  .0000340599071961
mat T_V[2,8] =  5.48971481482e-06
mat T_V[3,1] = -.0002097024888498
mat T_V[3,2] =  .0000106682740195
mat T_V[3,3] =  .0001959851158537
mat T_V[3,4] = -9.94632689721e-06
mat T_V[3,5] = -.0110275766134051
mat T_V[3,6] =  .0000118155822653
mat T_V[3,7] =  2.85535287105e-06
mat T_V[3,8] =  1.34772536814e-06
mat T_V[4,1] = -.0002830615447071
mat T_V[4,2] =  .0002905960668185
mat T_V[4,3] = -9.94632689721e-06
mat T_V[4,4] =  .0053750154913617
mat T_V[4,5] = -.0078547954144265
mat T_V[4,6] =  .0002541240482058
mat T_V[4,7] =  .0001046996209033
mat T_V[4,8] =  8.39862891781e-06
mat T_V[5,1] =  .0014014137772944
mat T_V[5,2] = -.0181494524318568
mat T_V[5,3] = -.0110275766134051
mat T_V[5,4] = -.0078547954144265
mat T_V[5,5] =  1.156867751289587
mat T_V[5,6] =  -.013916569286144
mat T_V[5,7] = -.0068370516492207
mat T_V[5,8] =  .0000196773492065
mat T_V[6,1] = -.0001426195145807
mat T_V[6,2] =  .0000203664407633
mat T_V[6,3] =  .0000118155822653
mat T_V[6,4] =  .0002541240482058
mat T_V[6,5] =  -.013916569286144
mat T_V[6,6] =  .0165153704085045
mat T_V[6,7] =  .0070723775092447
mat T_V[6,8] = -.0003227570086795
mat T_V[7,1] =  .0001437342833462
mat T_V[7,2] =  .0000340599071961
mat T_V[7,3] =  2.85535287105e-06
mat T_V[7,4] =  .0001046996209033
mat T_V[7,5] = -.0068370516492207
mat T_V[7,6] =  .0070723775092447
mat T_V[7,7] =  .0043757663577268
mat T_V[7,8] = -.0003381784215482
mat T_V[8,1] =  4.73869947113e-06
mat T_V[8,2] =  5.48971481482e-06
mat T_V[8,3] =  1.34772536814e-06
mat T_V[8,4] =  8.39862891781e-06
mat T_V[8,5] =  .0000196773492065
mat T_V[8,6] = -.0003227570086795
mat T_V[8,7] = -.0003381784215482
mat T_V[8,8] =  .0002919684829957
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_V T_V

qui {
mat T_rcsrmat_1 = J(4,4,0)
mat T_rcsrmat_1[1,1] =  .5338731215410862
mat T_rcsrmat_1[1,2] = -2.695161102636694
mat T_rcsrmat_1[1,3] = -1.176727157036661
mat T_rcsrmat_1[2,2] =  .6773003379716674
mat T_rcsrmat_1[2,3] =  .3248116475136131
mat T_rcsrmat_1[3,3] =  .0323402935563342
mat T_rcsrmat_1[4,1] =  2.012986809123404
mat T_rcsrmat_1[4,2] = -8.575807511173112
mat T_rcsrmat_1[4,3] = -3.624400613598444
mat T_rcsrmat_1[4,4] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3 r4"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3 c4"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,8,0)
mat T_gradient[1,1] =  .0000397357126856
mat T_gradient[1,2] =  .0030157194461089
mat T_gradient[1,3] =  .0056622200466103
mat T_gradient[1,4] = -6.58084179691e-07
mat T_gradient[1,5] =  .0001028975284411
mat T_gradient[1,6] = -.0001490017698843
mat T_gradient[1,7] = -.0000177541689498
mat T_gradient[1,8] = -.0001623811374384
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:_brcs1_1 tb1:_brcs1_2 tb1:_brcs1_3"'
mat drop C_gradient T_gradient


//nonPH model
set seed 72549

clear
set obs 500
gen id 	= _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen x1 = rnormal()
gen bmi = rnormal(30,3)


survsim stime died, dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.01 bmi -0.05 x1 0.1) ///
	tde(trt 0.03) tdefunction(log({t})) maxt(10)

gen t0 = 0
replace t0 = runiform() * 2 //if _n>200

drop if stime<t0

uhtred (stime trt bmi age x1 c.trt#rcs(stime, df(2) log orthog event) ///
	, family(rp, df(2) failure(died) ltruncated(t0)))
		
_assert_streq `"`e(cmdline)'"' `"uhtred (stime trt bmi age x1 c.trt#rcs(stime, df(2) log orthog event)         , family(rp, df(2) failure(died) ltruncated(t0)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(knots_1_5_2)'"'       == `"-1.126265 1.700885 2.296919"'
assert `"`e(orthog1)'"'           == `"orthog"'
assert `"`e(knots1)'"'            == `"-1.126265 1.700885 2.296919"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(ltruncated1)'"'       == `"t0"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died t0"'
assert `"`e(family1)'"'           == `"rp"'
assert `"`e(haszb)'"'             == `"0"'
assert `"`e(hastb)'"'             == `"1"'
assert `"`e(hasxb)'"'             == `"1"'
assert `"`e(allvars)'"'           == `"age bmi c.trt died stime t0 trt x1"'
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
assert         e(N)          == 483
assert         e(k)          == 9
assert         e(k_eq)       == 2
assert         e(noconstant) == 0
assert         e(consonly)   == 0
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -731.1511284717557) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,9,0)
mat T_b[1,1] =   .027392591636159
mat T_b[1,2] = -.0788900677858013
mat T_b[1,3] =  .0050175454770941
mat T_b[1,4] = -.0531196075039316
mat T_b[1,5] =  1.316616522073824
mat T_b[1,6] = -.2976930045773955
mat T_b[1,7] =  -.065349013140771
mat T_b[1,8] =  .6184681106138655
mat T_b[1,9] = -.0231208424790223
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:c.trt#c._rcs1_5_2_1 tb1:c.trt#c._rcs1_5_2_2 tb1:_brcs1_1 tb1:_brcs1_2"'
mat drop C_b T_b

qui {
mat T_V = J(9,9,0)
mat T_V[1,1] =  .1738451306239696
mat T_V[1,2] =   .000161248761339
mat T_V[1,3] = -.0001281710671851
mat T_V[1,4] =  .0004718079997034
mat T_V[1,5] =  -.019599903671928
mat T_V[1,6] = -.0749731316119751
mat T_V[1,7] =   .005202476800208
mat T_V[1,8] =  .0165788128210868
mat T_V[1,9] =   .007395068186089
mat T_V[2,1] =   .000161248761339
mat T_V[2,2] =  .0005833527494787
mat T_V[2,3] =  .0000105920727484
mat T_V[2,4] =  .0002784622319359
mat T_V[2,5] = -.0178181156709266
mat T_V[2,6] = -.0000868332264093
mat T_V[2,7] = -.0000574627713388
mat T_V[2,8] =  .0000474684817221
mat T_V[2,9] =  .0000663215660604
mat T_V[3,1] = -.0001281710671851
mat T_V[3,2] =  .0000105920727484
mat T_V[3,3] =  .0001916043088975
mat T_V[3,4] = -8.29099568505e-06
mat T_V[3,5] =  -.010781472950983
mat T_V[3,6] = -.0000347854515708
mat T_V[3,7] =  2.63003164094e-06
mat T_V[3,8] =  7.42065488693e-06
mat T_V[3,9] =  1.96685983232e-06
mat T_V[4,1] =  .0004718079997034
mat T_V[4,2] =  .0002784622319359
mat T_V[4,3] = -8.29099568505e-06
mat T_V[4,4] =  .0052628874436704
mat T_V[4,5] = -.0077243943104274
mat T_V[4,6] = -.0005667747216009
mat T_V[4,7] = -.0001290851515348
mat T_V[4,8] =  .0003774698396953
mat T_V[4,9] =  .0001650682315334
mat T_V[5,1] =  -.019599903671928
mat T_V[5,2] = -.0178181156709266
mat T_V[5,3] =  -.010781472950983
mat T_V[5,4] = -.0077243943104274
mat T_V[5,5] =   1.13746219727877
mat T_V[5,6] =  .0210047889746963
mat T_V[5,7] =   .008920660103391
mat T_V[5,8] = -.0183481147589076
mat T_V[5,9] = -.0094332886419427
mat T_V[6,1] = -.0749731316119751
mat T_V[6,2] = -.0000868332264093
mat T_V[6,3] = -.0000347854515708
mat T_V[6,4] = -.0005667747216009
mat T_V[6,5] =  .0210047889746963
mat T_V[6,6] =  .0463616622984609
mat T_V[6,7] =  .0037931667011061
mat T_V[6,8] = -.0217812096878099
mat T_V[6,9] = -.0091898379911611
mat T_V[7,1] =   .005202476800208
mat T_V[7,2] = -.0000574627713388
mat T_V[7,3] =  2.63003164094e-06
mat T_V[7,4] = -.0001290851515348
mat T_V[7,5] =   .008920660103391
mat T_V[7,6] =  .0037931667011061
mat T_V[7,7] =  .0071523829571082
mat T_V[7,8] = -.0091776773171354
mat T_V[7,9] =  -.005890256211081
mat T_V[8,1] =  .0165788128210868
mat T_V[8,2] =  .0000474684817221
mat T_V[8,3] =  7.42065488693e-06
mat T_V[8,4] =  .0003774698396953
mat T_V[8,5] = -.0183481147589076
mat T_V[8,6] = -.0217812096878099
mat T_V[8,7] = -.0091776773171354
mat T_V[8,8] =  .0217651537661072
mat T_V[8,9] =  .0091807711012099
mat T_V[9,1] =   .007395068186089
mat T_V[9,2] =  .0000663215660604
mat T_V[9,3] =  1.96685983232e-06
mat T_V[9,4] =  .0001650682315334
mat T_V[9,5] = -.0094332886419427
mat T_V[9,6] = -.0091898379911611
mat T_V[9,7] =  -.005890256211081
mat T_V[9,8] =  .0091807711012099
mat T_V[9,9] =  .0058921153761378
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:c.trt#c._rcs1_5_2_1 tb1:c.trt#c._rcs1_5_2_2 tb1:_brcs1_1 tb1:_brcs1_2"'
_assert_streq `"`: colfullnames C_V'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:c.trt#c._rcs1_5_2_1 tb1:c.trt#c._rcs1_5_2_2 tb1:_brcs1_1 tb1:_brcs1_2"'
mat drop C_V T_V

qui {
mat T_rmat_1_5_2 = J(3,3,0)
mat T_rmat_1_5_2[1,1] =  .5335825466231103
mat T_rmat_1_5_2[1,2] = -1.818182862046044
mat T_rmat_1_5_2[2,2] =  .4882334283001356
mat T_rmat_1_5_2[3,1] =  2.009764734973028
mat T_rmat_1_5_2[3,2] = -5.629430830503892
mat T_rmat_1_5_2[3,3] =                  1
}
matrix C_rmat_1_5_2 = e(rmat_1_5_2)
assert mreldif( C_rmat_1_5_2 , T_rmat_1_5_2 ) < 1E-8
_assert_streq `"`: rowfullnames C_rmat_1_5_2'"' `"r1 r2 r3"'
_assert_streq `"`: colfullnames C_rmat_1_5_2'"' `"c1 c2 c3"'
mat drop C_rmat_1_5_2 T_rmat_1_5_2

qui {
mat T_rcsrmat_1 = J(3,3,0)
mat T_rcsrmat_1[1,1] =  .5335825466231103
mat T_rcsrmat_1[1,2] = -1.818182862046044
mat T_rcsrmat_1[2,2] =  .4882334283001356
mat T_rcsrmat_1[3,1] =  2.009764734973028
mat T_rcsrmat_1[3,2] = -5.629430830503892
mat T_rcsrmat_1[3,3] =                  1
}
matrix C_rcsrmat_1 = e(rcsrmat_1)
assert mreldif( C_rcsrmat_1 , T_rcsrmat_1 ) < 1E-8
_assert_streq `"`: rowfullnames C_rcsrmat_1'"' `"r1 r2 r3"'
_assert_streq `"`: colfullnames C_rcsrmat_1'"' `"c1 c2 c3"'
mat drop C_rcsrmat_1 T_rcsrmat_1

qui {
mat T_gradient = J(1,9,0)
mat T_gradient[1,1] =  1.05284858193e-07
mat T_gradient[1,2] =  6.97118445991e-06
mat T_gradient[1,3] =  .0000131089264084
mat T_gradient[1,4] = -1.30048686671e-09
mat T_gradient[1,5] =  2.37478889914e-07
mat T_gradient[1,6] =  3.70498951063e-06
mat T_gradient[1,7] =  .0000110123955905
mat T_gradient[1,8] =  3.51936952705e-06
mat T_gradient[1,9] =  .0000109429094997
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:c.trt#c._rcs1_5_2_1 tb1:c.trt#c._rcs1_5_2_2 tb1:_brcs1_1 tb1:_brcs1_2"'
mat drop C_gradient T_gradient
