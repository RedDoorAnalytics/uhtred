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

qui {
	mat A = J(4,4,0)
	mat A[1,1] = .53387312
	mat A[1,2] = -2.6951611
	mat A[1,3] = -1.1767272
	mat A[1,4] = 0

	mat A[2,1] = 0
	mat A[2,2] = .67730034
	mat A[2,3] = .32481165
	mat A[2,4] = 0

	mat A[3,1] = 0
	mat A[3,2] = 0
	mat A[3,3] = .03234029
	mat A[3,4] =  0

	mat A[4,1] =  2.0129868
	mat A[4,2] = -8.5758075
	mat A[4,3] =  -3.6244006
	mat A[4,4] =  1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-5
matrix drop A






assert `"`r(PT_rseps)'"'            == `"`""' `""' `""' `""' `""' `""' `"sep"' `""' `""' `""'"'
assert `"`r(PT_rnotes)'"'           == `"`""' `""' `""' `""' `""' `""' `""' `""' `""' `""'"'
_assert_streq `"`r(PT_raligns)'"' `"`"left"' `"right"' `"right"' `"right"' `"right"' `"right"' `"left"' `"right"' `"right"' `"right"'"'

_assert_streq `"`r(PT_rtitles)'"' `"`"xb1"' `"trt"' `"bmi"' `"age"' `"x1"' `"_cons"' `"tb1"' `"_rcs1_1"' `"_rcs1_2"' `"_rcs1_3"'"'
assert `"`r(PT_cformats)'"'         == `"`"%9.0g"' `"%9.0g"' `"%8.2f"' `"%5.3f"' `"%9.0g"' `"%9.0g"'"'
assert `"`r(PT_cspans1)'"'          == `"`"1"' `"1"' `"1"' `"1"' `"2"' `"0"'"'
assert `"`r(PT_ctitles1)'"'         == `"`"Coefficient"' `"Std. err."' `"z"' `"P>|z|"' `"[95% conf. interval]"' `""'"'
assert `"`r(put_tables)'"'          == `"PT"'
assert `"`r(citype)'"'              == `"normal"'
assert `"`r(_collect_prefix_get)'"' == `"ignore"'

assert         r(PT_has_legend) == 0
assert         r(PT_has_cnotes) == 0
assert         r(PT_k_ctitles)  == 1
assert         r(level)         == 95
assert         r(k_eform)       == 1

qui {
mat T_PT = J(10,6,0)
mat T_PT[1,1] =                 .b
mat T_PT[1,2] =                 .b
mat T_PT[1,3] =                 .b
mat T_PT[1,4] =                 .b
mat T_PT[1,5] =                 .b
mat T_PT[1,6] =                 .b
mat T_PT[2,1] = -.4950876558110081
mat T_PT[2,2] =  .1498083520701547
mat T_PT[2,3] =  -3.30480676791078
mat T_PT[2,4] =  .0009504190178285
mat T_PT[2,5] = -.7887066304518077
mat T_PT[2,6] = -.2014686811702086
mat T_PT[3,1] = -.0820916113445683
mat T_PT[3,2] =  .0244091937951136
mat T_PT[3,3] = -3.363143085905685
mat T_PT[3,4] =  .0007706041097999
mat T_PT[3,5] = -.1299327520746495
mat T_PT[3,6] = -.0342504706144871
mat T_PT[4,1] =  .0053576960194916
mat T_PT[4,2] =  .0139994684132558
mat T_PT[4,3] =   .382707104393941
mat T_PT[4,4] =  .7019369462429386
mat T_PT[4,5] = -.0220807578731958
mat T_PT[4,6] =  .0327961499121789
mat T_PT[5,1] = -.0605198694899665
mat T_PT[5,2] =  .0733144971432098
mat T_PT[5,3] = -.8254829787858909
mat T_PT[5,4] =  .4090974422802846
mat T_PT[5,5] = -.2042136434353225
mat T_PT[5,6] =  .0831739044553894
mat T_PT[6,1] =  1.393239089930791
mat T_PT[6,2] =   1.07557786853839
mat T_PT[6,3] =  1.295340049925044
mat T_PT[6,4] =  .1952029489945581
mat T_PT[6,5] = -.7148547949728108
mat T_PT[6,6] =  3.501332974834392
mat T_PT[7,1] =                 .b
mat T_PT[7,2] =                 .b
mat T_PT[7,3] =                 .b
mat T_PT[7,4] =                 .b
mat T_PT[7,5] =                 .b
mat T_PT[7,6] =                 .b
mat T_PT[8,1] =  .6250287401096971
mat T_PT[8,2] =   .128512141093767
mat T_PT[8,3] =  4.863577361563481
mat T_PT[8,4] =  1.15282937459e-06
mat T_PT[8,5] =  .3731495719897839
mat T_PT[8,6] =  .8769079082296103
mat T_PT[9,1] = -.0359662562368732
mat T_PT[9,2] =  .0661495756428317
mat T_PT[9,3] =  -.543711065223904
mat T_PT[9,4] =  .5866403167684658
mat T_PT[9,5] = -.1656170420894313
mat T_PT[9,6] =  .0936845296156849
mat T_PT[10,1] =  .0054913619625757
mat T_PT[10,2] =  .0170870852691649
mat T_PT[10,3] =  .3213749961490113
mat T_PT[10,4] =  .7479262291326123
mat T_PT[10,5] = -.0279987097657525
mat T_PT[10,6] =  .0389814336909039
}
matrix C_PT = r(PT)
assert mreldif( C_PT , T_PT ) < 1E-8
_assert_streq `"`: rowfullnames C_PT'"' `"r1 r2 r3 r4 r5 r6 r7 r8 r9 r10"'
_assert_streq `"`: colfullnames C_PT'"' `"c1 c2 c3 c4 c5 c6"'
mat drop C_PT T_PT

qui {
mat T_table = J(9,8,0)
mat T_table[1,1] = -.4950876558110081
mat T_table[1,2] = -.0820916113445683
mat T_table[1,3] =  .0053576960194916
mat T_table[1,4] = -.0605198694899665
mat T_table[1,5] =  1.393239089930791
mat T_table[1,6] =  .6250287401096971
mat T_table[1,7] = -.0359662562368732
mat T_table[1,8] =  .0054913619625757
mat T_table[2,1] =  .1498083520701547
mat T_table[2,2] =  .0244091937951136
mat T_table[2,3] =  .0139994684132558
mat T_table[2,4] =  .0733144971432098
mat T_table[2,5] =   1.07557786853839
mat T_table[2,6] =   .128512141093767
mat T_table[2,7] =  .0661495756428317
mat T_table[2,8] =  .0170870852691649
mat T_table[3,1] =  -3.30480676791078
mat T_table[3,2] = -3.363143085905685
mat T_table[3,3] =   .382707104393941
mat T_table[3,4] = -.8254829787858909
mat T_table[3,5] =  1.295340049925044
mat T_table[3,6] =  4.863577361563481
mat T_table[3,7] =  -.543711065223904
mat T_table[3,8] =  .3213749961490113
mat T_table[4,1] =  .0009504190178285
mat T_table[4,2] =  .0007706041097999
mat T_table[4,3] =  .7019369462429386
mat T_table[4,4] =  .4090974422802846
mat T_table[4,5] =  .1952029489945581
mat T_table[4,6] =  1.15282937459e-06
mat T_table[4,7] =  .5866403167684658
mat T_table[4,8] =  .7479262291326123
mat T_table[5,1] = -.7887066304518077
mat T_table[5,2] = -.1299327520746495
mat T_table[5,3] = -.0220807578731958
mat T_table[5,4] = -.2042136434353225
mat T_table[5,5] = -.7148547949728108
mat T_table[5,6] =  .3731495719897839
mat T_table[5,7] = -.1656170420894313
mat T_table[5,8] = -.0279987097657525
mat T_table[6,1] = -.2014686811702086
mat T_table[6,2] = -.0342504706144871
mat T_table[6,3] =  .0327961499121789
mat T_table[6,4] =  .0831739044553894
mat T_table[6,5] =  3.501332974834392
mat T_table[6,6] =  .8769079082296103
mat T_table[6,7] =  .0936845296156849
mat T_table[6,8] =  .0389814336909039
mat T_table[7,1] =                  .
mat T_table[7,2] =                  .
mat T_table[7,3] =                  .
mat T_table[7,4] =                  .
mat T_table[7,5] =                  .
mat T_table[7,6] =                  .
mat T_table[7,7] =                  .
mat T_table[7,8] =                  .
mat T_table[8,1] =  1.959963984540054
mat T_table[8,2] =  1.959963984540054
mat T_table[8,3] =  1.959963984540054
mat T_table[8,4] =  1.959963984540054
mat T_table[8,5] =  1.959963984540054
mat T_table[8,6] =  1.959963984540054
mat T_table[8,7] =  1.959963984540054
mat T_table[8,8] =  1.959963984540054
}
matrix C_table = r(table)
assert mreldif( C_table , T_table ) < 1E-8
_assert_streq `"`: rowfullnames C_table'"' `"b se z pvalue ll ul df crit eform"'
_assert_streq `"`: colfullnames C_table'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:_rcs1_1 tb1:_rcs1_2 tb1:_rcs1_3"'
mat drop C_table T_table


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
	tde(trt 0.03) tdefunction(log({t}))maxt(10)

gen t0 = 0
replace t0 = runiform() * 2 //if _n>200

drop if stime<t0





uhtred (stime trt bmi age x1 c.trt#rcs(stime, df(2) log orthog event) ///
	, family(rp, df(2) failure(died) ltruncated(t0)))
		
qui {
	mat A = J(3,3,0)
	mat A[1,1] = .53358255
	mat A[1,2] = -1.8181829
	mat A[1,3] = 0

	mat A[2,1] = 0
	mat A[2,2] = .48823343
	mat A[2,3] = 0

	mat A[3,1] =  2.0097647
	mat A[3,2] = -5.6294308 
	mat A[3,3] = 1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-5
matrix drop A

qui {
	mat B = J(3,3,0)
	mat B[1,1] = .53358255
	mat B[1,2] = -1.8181829
	mat B[1,3] = 0

	mat B[2,1] = 0
	mat B[2,2] = .48823343 
	mat B[2,3] = 0

	mat B[3,1] = 2.0097647
	mat B[3,2] = -5.6294308  
	mat B[3,3] =  1
}
assert mreldif( e(rmat_1_5_2) , B ) < 1E-5
matrix drop B
	




assert `"`r(PT_cformats)'"'         == `"`"%9.0g"' `"%9.0g"' `"%8.2f"' `"%5.3f"' `"%9.0g"' `"%9.0g"'"'
assert `"`r(PT_cspans1)'"'          == `"`"1"' `"1"' `"1"' `"1"' `"2"' `"0"'"'
assert `"`r(PT_ctitles1)'"'         == `"`"Coefficient"' `"Std. err."' `"z"' `"P>|z|"' `"[95% conf. interval]"' `""'"'
assert `"`r(put_tables)'"'          == `"PT"'
assert `"`r(citype)'"'              == `"normal"'
assert `"`r(_collect_prefix_get)'"' == `"ignore"'

assert         r(PT_has_legend) == 0
assert         r(PT_has_cnotes) == 0
assert         r(PT_k_ctitles)  == 1
assert         r(level)         == 95
assert         r(k_eform)       == 1

qui {
mat T_PT = J(17,6,0)
mat T_PT[1,1] =                 .b
mat T_PT[1,2] =                 .b
mat T_PT[1,3] =                 .b
mat T_PT[1,4] =                 .b
mat T_PT[1,5] =                 .b
mat T_PT[1,6] =                 .b
mat T_PT[2,1] =  .0273925916361476
mat T_PT[2,2] =  .4169473955116808
mat T_PT[2,3] =  .0656979559796296
mat T_PT[2,4] =  .9476182997912777
mat T_PT[2,5] = -.7898092870145239
mat T_PT[2,6] =  .8445944702868191
mat T_PT[3,1] = -.0788900677858013
mat T_PT[3,2] =  .0241526965260348
mat T_PT[3,3] = -3.266304766458013
mat T_PT[3,4] =  .0010896091256243
mat T_PT[3,5] = -.1262284831063552
mat T_PT[3,6] = -.0315516524652475
mat T_PT[4,1] =  .0050175454770941
mat T_PT[4,2] =  .0138421208236846
mat T_PT[4,3] =  .3624838665263512
mat T_PT[4,4] =  .7169904776276192
mat T_PT[4,5] = -.0221125128069796
mat T_PT[4,6] =  .0321476037611678
mat T_PT[5,1] = -.0531196075039317
mat T_PT[5,2] =  .0725457610317132
mat T_PT[5,3] = -.7322220726406126
mat T_PT[5,4] =  .4640330341303085
mat T_PT[5,5] = -.1953066863571389
mat T_PT[5,6] =  .0890674713492756
mat T_PT[6,1] =  1.316616522073825
mat T_PT[6,2] =   1.06651872804875
mat T_PT[6,3] =  1.234499205168804
mat T_PT[6,4] =  .2170169427206989
mat T_PT[6,5] = -.7737217737391928
mat T_PT[6,6] =  3.406954817886843
mat T_PT[7,1] =                 .b
mat T_PT[7,2] =                 .b
mat T_PT[7,3] =                 .b
mat T_PT[7,4] =                 .b
mat T_PT[7,5] =                 .b
mat T_PT[7,6] =                 .b
mat T_PT[8,1] =                 .b
mat T_PT[8,2] =                 .b
mat T_PT[8,3] =                 .b
mat T_PT[8,4] =                 .b
mat T_PT[8,5] =                 .b
mat T_PT[8,6] =                 .b
mat T_PT[9,1] =                 .b
mat T_PT[9,2] =                 .b
mat T_PT[9,3] =                 .b
mat T_PT[9,4] =                 .b
mat T_PT[9,5] =                 .b
mat T_PT[9,6] =                 .b
mat T_PT[10,1] = -.2976930045773885
mat T_PT[10,2] =  .2153175847404536
mat T_PT[10,3] = -1.382576369395148
mat T_PT[10,4] =  .1667947988925097
mat T_PT[10,5] = -.7197077159068286
mat T_PT[10,6] =  .1243217067520516
mat T_PT[11,1] =                 .b
mat T_PT[11,2] =                 .b
mat T_PT[11,3] =                 .b
mat T_PT[11,4] =                 .b
mat T_PT[11,5] =                 .b
mat T_PT[11,6] =                 .b
mat T_PT[12,1] =                 .b
mat T_PT[12,2] =                 .b
mat T_PT[12,3] =                 .b
mat T_PT[12,4] =                 .b
mat T_PT[12,5] =                 .b
mat T_PT[12,6] =                 .b
mat T_PT[13,1] =                 .b
mat T_PT[13,2] =                 .b
mat T_PT[13,3] =                 .b
mat T_PT[13,4] =                 .b
mat T_PT[13,5] =                 .b
mat T_PT[13,6] =                 .b
mat T_PT[14,1] = -.0653490131407713
mat T_PT[14,2] =  .0845717621733655
mat T_PT[14,3] = -.7727048776258315
mat T_PT[14,4] =  .4396970593419642
mat T_PT[14,5] = -.2311066211096545
mat T_PT[14,6] =  .1004085948281119
mat T_PT[15,1] =                 .b
mat T_PT[15,2] =                 .b
mat T_PT[15,3] =                 .b
mat T_PT[15,4] =                 .b
mat T_PT[15,5] =                 .b
mat T_PT[15,6] =                 .b
mat T_PT[16,1] =  .6184681106138655
mat T_PT[16,2] =  .1475301791705923
mat T_PT[16,3] =  4.192146407540911
mat T_PT[16,4] =  .0000276327565074
mat T_PT[16,5] =  .3293142728067633
mat T_PT[16,6] =  .9076219484209678
mat T_PT[17,1] = -.0231208424790228
mat T_PT[17,2] =  .0767601157902843
mat T_PT[17,3] = -.3012090620367362
mat T_PT[17,4] =  .7632550800259732
mat T_PT[17,5] = -.1735679048771044
mat T_PT[17,6] =  .1273262199190588
}
matrix C_PT = r(PT)
assert mreldif( C_PT , T_PT ) < 1E-8
_assert_streq `"`: rowfullnames C_PT'"' `"r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15 r16 r17"'
_assert_streq `"`: colfullnames C_PT'"' `"c1 c2 c3 c4 c5 c6"'
mat drop C_PT T_PT

qui {
mat T_table = J(9,9,0)
mat T_table[1,1] =  .0273925916361476
mat T_table[1,2] = -.0788900677858013
mat T_table[1,3] =  .0050175454770941
mat T_table[1,4] = -.0531196075039317
mat T_table[1,5] =  1.316616522073825
mat T_table[1,6] = -.2976930045773885
mat T_table[1,7] = -.0653490131407713
mat T_table[1,8] =  .6184681106138655
mat T_table[1,9] = -.0231208424790228
mat T_table[2,1] =  .4169473955116808
mat T_table[2,2] =  .0241526965260348
mat T_table[2,3] =  .0138421208236846
mat T_table[2,4] =  .0725457610317132
mat T_table[2,5] =   1.06651872804875
mat T_table[2,6] =  .2153175847404536
mat T_table[2,7] =  .0845717621733655
mat T_table[2,8] =  .1475301791705923
mat T_table[2,9] =  .0767601157902843
mat T_table[3,1] =  .0656979559796296
mat T_table[3,2] = -3.266304766458013
mat T_table[3,3] =  .3624838665263512
mat T_table[3,4] = -.7322220726406126
mat T_table[3,5] =  1.234499205168804
mat T_table[3,6] = -1.382576369395148
mat T_table[3,7] = -.7727048776258315
mat T_table[3,8] =  4.192146407540911
mat T_table[3,9] = -.3012090620367362
mat T_table[4,1] =  .9476182997912777
mat T_table[4,2] =  .0010896091256243
mat T_table[4,3] =  .7169904776276192
mat T_table[4,4] =  .4640330341303085
mat T_table[4,5] =  .2170169427206989
mat T_table[4,6] =  .1667947988925097
mat T_table[4,7] =  .4396970593419642
mat T_table[4,8] =  .0000276327565074
mat T_table[4,9] =  .7632550800259732
mat T_table[5,1] = -.7898092870145239
mat T_table[5,2] = -.1262284831063552
mat T_table[5,3] = -.0221125128069796
mat T_table[5,4] = -.1953066863571389
mat T_table[5,5] = -.7737217737391928
mat T_table[5,6] = -.7197077159068286
mat T_table[5,7] = -.2311066211096545
mat T_table[5,8] =  .3293142728067633
mat T_table[5,9] = -.1735679048771044
mat T_table[6,1] =  .8445944702868191
mat T_table[6,2] = -.0315516524652475
mat T_table[6,3] =  .0321476037611678
mat T_table[6,4] =  .0890674713492756
mat T_table[6,5] =  3.406954817886843
mat T_table[6,6] =  .1243217067520516
mat T_table[6,7] =  .1004085948281119
mat T_table[6,8] =  .9076219484209678
mat T_table[6,9] =  .1273262199190588
mat T_table[7,1] =                  .
mat T_table[7,2] =                  .
mat T_table[7,3] =                  .
mat T_table[7,4] =                  .
mat T_table[7,5] =                  .
mat T_table[7,6] =                  .
mat T_table[7,7] =                  .
mat T_table[7,8] =                  .
mat T_table[7,9] =                  .
mat T_table[8,1] =  1.959963984540054
mat T_table[8,2] =  1.959963984540054
mat T_table[8,3] =  1.959963984540054
mat T_table[8,4] =  1.959963984540054
mat T_table[8,5] =  1.959963984540054
mat T_table[8,6] =  1.959963984540054
mat T_table[8,7] =  1.959963984540054
mat T_table[8,8] =  1.959963984540054
mat T_table[8,9] =  1.959963984540054
}
matrix C_table = r(table)
assert mreldif( C_table , T_table ) < 1E-8
_assert_streq `"`: rowfullnames C_table'"' `"b se z pvalue ll ul df crit eform"'
_assert_streq `"`: colfullnames C_table'"' `"xb1:trt xb1:bmi xb1:age xb1:x1 xb1:_cons tb1:c.trt#c._rcs1_5_2_1 tb1:c.trt#c._rcs1_5_2_2 tb1:_rcs1_1 tb1:_rcs1_2"'
mat drop C_table T_table
