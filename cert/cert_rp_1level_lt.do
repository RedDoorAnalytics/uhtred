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


stset stime, enter(t0) f(died)

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
uhtred (stime trt bmi age x1 ///
		c.trt#rcs(stime, df(2) log orthog event) ///
		, family(rp, df(2) failure(died) ltruncated(t0)))
		
		
qui {
	mat A = J(3,3,0)
	mat A[1,1] = .53387312
	mat A[1,2] = -1.8618457
	mat A[1,3] = 0

	mat A[2,1] = 0
	mat A[2,2] = .49964344
	mat A[2,3] = 0

	mat A[3,1] = 2.0129868
	mat A[3,2] = -5.7941996 
	mat A[3,3] = 1
}
assert mreldif( e(rcsrmat_1) , A ) < 1E-5
matrix drop A

qui {
	mat B = J(3,3,0)
	mat B[1,1] = .53387312
	mat B[1,2] = -1.8618457
	mat B[1,3] = 0

	mat B[2,1] = 0
	mat B[2,2] = .49964344 
	mat B[2,3] = 0

	mat B[3,1] = 2.0129868
	mat B[3,2] = -5.7941996  
	mat B[3,3] =  1
}
assert mreldif( e(rmat_1_5_2) , B ) < 1E-5
matrix drop B	
	

assert `"`r(PT_rseps)'"'            == `"`""' `""' `""' `""' `""' `""' `"sep"' `""' `""' `""' `""' `""' `""'"'
assert `"`r(PT_rnotes)'"'           == `"`""' `""' `""' `""' `""' `""' `""' `""' `""' `""' `""' `""' `""'"'
_assert_streq `"`r(PT_raligns)'"' `"`"left"' `"right"' `"right"' `"right"' `"right"' `"right"' `"left"' `"right"' `"right"' `"right"' `"right"' `"right"' `"right"'"'
_assert_streq `"`r(PT_rtitles)'"' `"`"xb1"' `"trt"' `"bmi"' `"age"' `"x1"' `"_cons"' `"tb1"' `"c.trt#c._rcs1_5_2_1"' `""' `"c.trt#c._rcs1_5_2_2"' `""' `"_rcs1_1"' `"_rcs1_2"'"'


qui {
mat T_PT = J(13,6,0)
mat T_PT[1,1] =                 .b
mat T_PT[1,2] =                 .b
mat T_PT[1,3] =                 .b
mat T_PT[1,4] =                 .b
mat T_PT[1,5] =                 .b
mat T_PT[1,6] =                 .b
mat T_PT[2,1] = -.0118461796404981
mat T_PT[2,2] =  .4568296001397835
mat T_PT[2,3] = -.0259312873703309
mat T_PT[2,4] =   .979312144719611
mat T_PT[2,5] = -.9072157429863077
mat T_PT[2,6] =  .8835233837053116
mat T_PT[3,1] = -.0819249136206437
mat T_PT[3,2] =  .0244165098837342
mat T_PT[3,3] = -3.355308109584508
mat T_PT[3,4] =  .0007927663461861
mat T_PT[3,5] =  -.129780393620929
mat T_PT[3,6] = -.0340694336203584
mat T_PT[4,1] =   .005456519145374
mat T_PT[4,2] =  .0139953849861674
mat T_PT[4,3] =  .3898798890324932
mat T_PT[4,4] =  .6966253657077921
mat T_PT[4,5] = -.0219739313772868
mat T_PT[4,6] =  .0328869696680348
mat T_PT[5,1] = -.0602156334114575
mat T_PT[5,2] =  .0733416863701172
mat T_PT[5,3] = -.8210287544736922
mat T_PT[5,4] =  .4116298897121654
mat T_PT[5,5] = -.2039626972623193
mat T_PT[5,6] =  .0835314304394043
mat T_PT[6,1] =  1.387150784895123
mat T_PT[6,2] =  1.077983740047763
mat T_PT[6,3] =  1.286801213563445
mat T_PT[6,4] =  .1981635838581495
mat T_PT[6,5] = -.7256585215182807
mat T_PT[6,6] =  3.499960091308526
mat T_PT[7,1] =                 .b
mat T_PT[7,2] =                 .b
mat T_PT[7,3] =                 .b
mat T_PT[7,4] =                 .b
mat T_PT[7,5] =                 .b
mat T_PT[7,6] =                 .b
mat T_PT[8,1] = -.3090360177999907
mat T_PT[8,2] =  .2246716003172449
mat T_PT[8,3] = -1.375501030675974
mat T_PT[8,4] =  .1689761664019615
mat T_PT[8,5] = -.7493842627707684
mat T_PT[8,6] =  .1313122271707869
mat T_PT[9,1] =                 .b
mat T_PT[9,2] =                 .b
mat T_PT[9,3] =                 .b
mat T_PT[9,4] =                 .b
mat T_PT[9,5] =                 .b
mat T_PT[9,6] =                 .b
mat T_PT[10,1] = -.0583891939377452
mat T_PT[10,2] =  .0856723962345034
mat T_PT[10,3] = -.6815403385930953
mat T_PT[10,4] =  .4955296513562518
mat T_PT[10,5] = -.2263040050266168
mat T_PT[10,6] =  .1095256171511263
mat T_PT[11,1] =                 .b
mat T_PT[11,2] =                 .b
mat T_PT[11,3] =                 .b
mat T_PT[11,4] =                 .b
mat T_PT[11,5] =                 .b
mat T_PT[11,6] =                 .b
mat T_PT[12,1] =  .6167617037517593
mat T_PT[12,2] =  .1506737419521023
mat T_PT[12,3] =  4.093358907538262
mat T_PT[12,4] =  .0000425168831717
mat T_PT[12,5] =  .3214465961097571
mat T_PT[12,6] =  .9120768113937616
mat T_PT[13,1] = -.0244152477892564
mat T_PT[13,2] =  .0777824948357272
mat T_PT[13,3] = -.3138912918751218
mat T_PT[13,4] =  .7536036044578849
mat T_PT[13,5] = -.1768661362949545
mat T_PT[13,6] =  .1280356407164417
}
matrix C_PT = r(PT)
assert mreldif( C_PT , T_PT ) < 1E-8
_assert_streq `"`: rowfullnames C_PT'"' `"r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13"'
_assert_streq `"`: colfullnames C_PT'"' `"c1 c2 c3 c4 c5 c6"'
mat drop C_PT T_PT

qui {
mat T_table = J(9,9,0)
mat T_table[1,1] = -.0118461796404981
mat T_table[1,2] = -.0819249136206437
mat T_table[1,3] =   .005456519145374
mat T_table[1,4] = -.0602156334114575
mat T_table[1,5] =  1.387150784895123
mat T_table[1,6] = -.3090360177999907
mat T_table[1,7] = -.0583891939377452
mat T_table[1,8] =  .6167617037517593
mat T_table[1,9] = -.0244152477892564
mat T_table[2,1] =  .4568296001397835
mat T_table[2,2] =  .0244165098837342
mat T_table[2,3] =  .0139953849861674
mat T_table[2,4] =  .0733416863701172
mat T_table[2,5] =  1.077983740047763
mat T_table[2,6] =  .2246716003172449
mat T_table[2,7] =  .0856723962345034
mat T_table[2,8] =  .1506737419521023
mat T_table[2,9] =  .0777824948357272
mat T_table[3,1] = -.0259312873703309
mat T_table[3,2] = -3.355308109584508
mat T_table[3,3] =  .3898798890324932
mat T_table[3,4] = -.8210287544736922
mat T_table[3,5] =  1.286801213563445
mat T_table[3,6] = -1.375501030675974
mat T_table[3,7] = -.6815403385930953
mat T_table[3,8] =  4.093358907538262
mat T_table[3,9] = -.3138912918751218
mat T_table[4,1] =   .979312144719611
mat T_table[4,2] =  .0007927663461861
mat T_table[4,3] =  .6966253657077921
mat T_table[4,4] =  .4116298897121654
mat T_table[4,5] =  .1981635838581495
mat T_table[4,6] =  .1689761664019615
mat T_table[4,7] =  .4955296513562518
mat T_table[4,8] =  .0000425168831717
mat T_table[4,9] =  .7536036044578849
mat T_table[5,1] = -.9072157429863077
mat T_table[5,2] =  -.129780393620929
mat T_table[5,3] = -.0219739313772868
mat T_table[5,4] = -.2039626972623193
mat T_table[5,5] = -.7256585215182807
mat T_table[5,6] = -.7493842627707684
mat T_table[5,7] = -.2263040050266168
mat T_table[5,8] =  .3214465961097571
mat T_table[5,9] = -.1768661362949545
mat T_table[6,1] =  .8835233837053116
mat T_table[6,2] = -.0340694336203584
mat T_table[6,3] =  .0328869696680348
mat T_table[6,4] =  .0835314304394043
mat T_table[6,5] =  3.499960091308526
mat T_table[6,6] =  .1313122271707869
mat T_table[6,7] =  .1095256171511263
mat T_table[6,8] =  .9120768113937616
mat T_table[6,9] =  .1280356407164417
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

		
