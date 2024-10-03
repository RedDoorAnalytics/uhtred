set seed 725488
clear
set obs 1000
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,0.25\0.25,1)
drawnorm u1 u2, means(0 0) sds(1 0.5) corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt
gen age = rnormal() + u2

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(age 0.01 trtui 1 u1 1) maxt(5) 
stset stime, f(dead)


//============================================================================//
//error checks

//2 level random intercept
rcof "uhtred (stime trt age M1[id1]@1, family(rp, df(1) failure(dead))), cov(bla) " == 198
rcof "uhtred (stime trt age M1[id1]@1, family(rp, df(1) failure(dead))), intmethod(bla) " == 198
rcof "uhtred (stime trt age M[id1]@1, family(rp, df(1) failure(dead))) " == 1986
rcof "uhtred (stime trt age M1[id1]@trt, family(rp, df(1) failure(dead))) " == 1986
rcof "uhtred (stime trt age M1[bla]@1, family(rp, df(1) failure(dead))) " == 3598

