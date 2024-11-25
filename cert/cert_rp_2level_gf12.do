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
	family(rp, df(1) failure(dead))) ///
	, evaltype(gf0)
	
matrix C_b0 = e(b)

uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))) ///
	, evaltype(gf1)
	
matrix C_b1 = e(b)

uhtred (stime trt age M1[id1]@1, ///
	family(rp, df(1) failure(dead))) ///
	, evaltype(gf2)
	
matrix C_b2 = e(b)

assert mreldif( C_b0 , C_b1 ) < 1E-4
assert mreldif( C_b0 , C_b2 ) < 1E-4
