*! version 1.0.0

local gml 	struct uhtred_struct scalar
local pgml	pointer(struct uhtred_struct scalar) scalar
local Egml	struct uhtred_ereturn_struct scalar
local RS 	real scalar
local SS 	string scalar
local PS 	pointer scalar
local RM 	real matrix
local SM	string matrix
local PC 	pointer colvector
local SC 	string colvector
local SR	string rowvector
local TR	transmorphic
local RC	real colvector
local RR	real rowvector
local PM	pointer matrix

version 17

mata:

`RM' _logl_weibull(`TR' M, `RR' b, `gml' gml, `RM' G, `RM' H)
{	
	// Setup //
	
	//model info
	model 	= 1
	hasxb 	= gml.hasxb[model]
	if (hasxb) xeqn = gml.xeqn[model]
	hastb 	= gml.hastb[model]
	if (hastb) teqn = gml.teqn[model]
	haszb 	= gml.haszb[model]
	if (haszb) zeqn = gml.zeqn[model,.]
	haslt	= gml.hasltrunc[model]
	hasic	= gml.haslint[model]
	hasbh	= gml.hasbh[model,]
	Nrelevs = gml.Nrelevels
	
	//data
	y 	= uhtred_util_depvar(gml)
	Nobs	= uhtred_util_nobs(gml)

	//logl
	logl 	= (Nrelevs ? J(Nobs,gml.ndim[Nrelevs],0) : J(Nobs,1,0))

	//linear predictors
	if (hastb) {
		xzb   = uhtred_util_xzb(M,b,gml)
		xtzb  = xzb + uhtred_util_tb(M,b,gml)
	}
	else 	xtzb  = uhtred_util_xtzb(M,b,gml)
	gam     = exp(uhtred_util_dap(gml,1))

	// log likelihood //
	
	//exactly observed events -> log hazard function
	gml.survind 	= 1
	Nobs1 		= uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 		= uhtred_get_surv_index(gml)
		logl[index1,] 	= xtzb[index1,] :+ log(gam) :+ 
			(gam - 1) :* log(y[index1,1])
		if (hasbh[1]) {
			logl[index1,] = log(exp(logl[index1,]) :+ 
				uhtred_util_bhazard(gml))
		}
		if (hasbh[2]) {
			logl[index1,] = log(exp(logl[index1,]) :* 
				uhtred_util_bhazard(gml))
		}
	}

	//exactly observed events and/or right censoring -> survival function
	gml.survind 	= 2
	Nobs2 		= uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (hastb | hasbh[2]) {
			Ngq 	= gml.chip
			chq2 	= J(Nobs2,Ngq,0)
			gq 	= uhtred_glegendre(Ngq)
			y2	= y[index2,1] :/ 2
			if (haslt) {
				y02   = y[index2,3] :/ 2
				y2m02 = (y2 :- y02)
				qp2 = y2m02 :* J(Nobs2,1,gq[,1]') :+ y2 :+ y02
			}
			else {
				qp2 = y2 :* J(Nobs2,1,gq[,1]') :+ y2
			}

			for (q=1;q<=Ngq;q++) {
				chq2[,q] = uhtred_util_tb(M,b,gml,qp2[,q]) 
			}
			if (hasbh[2]) {
				chq2 = chq2 :+ log(uhtred_util_bhazard(gml))
			}
			chq2 = exp(chq2 :+ (gam-1) :* log(qp2))
			if (haslt) 	chq = y2m02 :* (chq2 * gq[,2]) 
			else 		chq = y2 :* (chq2 * gq[,2]) 
			expxzb1   	= exp(xzb[index2,])
			expxzb1gam 	= expxzb1 :* gam
			logl[index2,] 	= logl[index2,] :- chq :* expxzb1gam
		}
		else {
			logl[index2,] = logl[index2,] :- exp(xtzb[index2,]) :* 
						y[index2,1] :^ gam
			if (haslt) {
				gml.survind 	= 4
				index4 		= uhtred_get_surv_index(gml)
				logl[index4,] 	= logl[index4,] :+ 
					exp(xtzb[index4,]) :* y[index4,3] :^ gam
			}
		}
	}

	//interval censoring -> cdf function
	gml.survind = 3
	Nobs3 = uhtred_util_nobs(gml)
	if (Nobs3) {
		//exit times
		index3 = uhtred_get_surv_index(gml)
		/*if (hast) {
			Ngq 	= gml.chip
			chq3 	= J(Nobs3,Ngq,.)
			gq 		= uhtred_glegendre(Ngq)
			qp3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,1]') :+ y[index3,1]:/2
			for (q=1;q<=Ngq;q++) {
				chq3[,q] = exp(uhtred_util_xzb(gml,qp3[,q])) :* gam :* qp3[,q] :^(gam-1)
			}
			ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
			logl[index3,] = -exp(-ch3)
		}
		else */
		logl[index3,] = -exp(-exp(xtzb[index3,]) :* y[index3,1] :^gam)
		
		//entry times
		gml.survind = 5
		if (gml.hasltrunc[model]) tind = 4
		else 			  tind = 3
		index5 = uhtred_get_surv_index(gml)
		/*if (hast) {
			Nobs5 	= uhtred_get_nobs(gml)
			chq5 	= J(Nobs5,Ngq,.)
			qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
			for (q=1;q<=Ngq;q++) {
				chq5[,q] = exp(uhtred_util_xzb(gml,qp5[,q])) :* gam :* qp5[,q] :^(gam-1)
			}
			ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
			logl[index5,] = logl[index5,] :+ exp(-ch5)
		}
		else*/
		logl[index5,] = logl[index5,] :+ exp(-exp(xtzb[index5,]) :* y[index5,tind] :^gam)
		//logL
		logl[index3,] = log(logl[index3,])
	}	
	
	if (gml.todo==0) return(logl)

	// score //
	
	gml.survind 	= 0

	//get coefficient indices
	if (hasxb) {
		sindex1 = merlin_util_score_indices(M,xeqn)
		Nxbs 	= cols(sindex1)
		if (hastb) geqn = 3
		else geqn = 2
	}
	else geqn = 2
	
	if (hastb) {
		sindex2 = merlin_util_score_indices(M,teqn)
		Ntbs 	= cols(sindex2)
	}
	sindex4 = merlin_util_score_indices(M,geqn)
	
	
	//exactly observed events -> hazard function
	if (Nobs1) {
		if (hasxb) G[index1,sindex1] = gml.X[index1,]
		if (hastb) G[index1,sindex2] = gml.XT[index1,]
		G[index1,sindex4] = 1 :+ gam :* log(y[index1,1])
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		if (hastb) {
			dch2 = dch22 = J(Nobs2,1,0)
			if (haslt) {
				qw2  = (y[index2,1]:-y[index2,3]):/2 :* 
					J(Nobs2,1,gq[,2]')
			}
			else qw2  = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
			dchq2temp = chq2 :* qw2
			for (q=1;q<=Ngq;q++) {	
				dch2 = dch2 :+ dchq2temp[,q] :* 
					uhtred_util_t(gml,qp2[,q]) 
				dch22 = dch22 :+ dchq2temp[,q] :* 
					(gam:*log(qp2[,q]) :+ 1)
			}
			if (hasxb) {
				Gtemp = -chq :* expxzb1gam

				G[index2,sindex1] = G[index2,sindex1] :+ 
					Gtemp :* gml.X[index2,]
			}
			G[index2,sindex2] = G[index2,sindex2] :- 
						dch2 :* expxzb1gam
			G[index2,sindex4] = G[index2,sindex4] :- dch22 :* 
						expxzb1gam
		}
		else {
			Gtemp = exp(xtzb[index2,]) :* y[index2,1] :^ gam
			G[index2,sindex1] = G[index2,sindex1] :- 
				Gtemp :* gml.X[index2,]
			
			G[index2,sindex4] = G[index2,sindex4] :- 
				Gtemp :* gam :* log(y[index2,1])
			if (haslt) {
				Gtemp0 = exp(xtzb[index4,]) :* y[index4,3] :^ gam
				G[index4,sindex1] = G[index4,sindex1] :+ 
					Gtemp0 :* gml.X[index4,]
				G[index4,sindex4] = G[index4,sindex4] :+ 
					Gtemp0 :* gam :* log(y[index4,3])
			}
		}
	}

	if (gml.todo==1) return(logl)
	
	// Hessian //
	
	if (hasxb) {
		Hxb2 = Hxbg = J(gml.Nobs,1,0)
		if (hastb) {
			xtbind  = _uhtred_get_hessian_indices(Nxbs,Ntbs)
			Hxbtb  	= J(Nobs,cols(xtbind),0)
			tbind   = _uhtred_get_hessian_indices(Ntbs)
			Htb2    = J(Nobs,cols(tbind),0)
			Htbg 	= J(gml.Nobs,Ntbs,0)
		}
	}
	else {
		if (hastb) {
			Htb2 = Htbg = J(gml.Nobs,1,0)
		}
	}
	Hg2 = J(gml.Nobs,1,0)
	
	if (Nobs1) {
		 Hg2[index1] = gam :* log(y[index1,1])
	}  

	gml.survind = 2
	if (Nobs2) {
		if (hastb) {

			if (hasxb) {
				Hxb2[index2] = Gtemp
				Hxbg[index2] = -dch22 :* expxzb1gam
			}
			
			logqp2	= log(qp2)
			if (!haslt) {
				hazq = y[index2,1] :/ 2 :* 
					J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* 
					gam :* (gam :* logqp2 :+ 1)
			}
			else {
				hazq = (y[index2,1] :- y[index2,3]) :/ 2 :* 
					J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* 
					gam :* (gam :* logqp2 :+ 1)
			}
			hazq2 = logqp2 :* gam :* qp2:^(gam:-1) :* 
				((gam:^2 :* logqp2 :+ 2:*gam) :+ gam :+ 
				1:/logqp2)

			dchtg = J(Nobs2,1,0)
			chdxdt = chdt2 = J(Nobs2,1,0)
			for (q=1;q<=Ngq;q++) {	
				tq = uhtred_util_t(gml,qp2[,q]) 
				exptqb = exp(uhtred_util_tb(M,b,gml,qp2[,q]))
				chdxdt = chdxdt :+ dchq2temp[,q] :* 
					gml.X[index2,xtbind[1,]] :* 
					tq[index2,xtbind[2,]]
				chdt2 = chdt2 :+ dchq2temp[,q] :* 
					tq[index2,tbind[1,]] :* 
					tq[index2,tbind[2,]]
				dchtg = dchtg :+ exptqb :* hazq[,q] :* 
					tq
				hazq2[,q] = hazq2[,q] :* exptqb

			}
			Hxbtb[index2,] 	= - chdxdt :* expxzb1gam
			Htb2[index2,] 	= - chdt2 :* expxzb1gam
			Htbg[index2,] 	= - dchtg :* exp(xzb[index2])
			
			if (!haslt) {
				Hg2[index2] = Hg2[index2] :- 
					y[index2,1]:/2 :* (hazq2 * gq[,2]) :* 
					exp(xzb)
			}
			else {
				Hg2[index2] = Hg2[index2] :- 
					(y[index2,1] :- y[index2,3]) :/ 2 :* 
					(hazq2 * gq[,2]) :* exp(xzb)
			}
		}
		else {
			Hxb2[index2,] = -Gtemp
			Hxbg[index2,] = -Gtemp :* gam :* log(y[index2,1]) 
			Hg2[index2] = Hg2[index2] :- log(y[index2,1]) :* gam :* Gtemp :* (gam :* log(y[index2,1]) :+ 1)
			if (haslt) {
				Hxb2[index4,] = Hxb2[index4,] :+ Gtemp0
				Hxbg[index4,] = Hxbg[index4,] :+ Gtemp0 :* gam :* log(y[index4,3])
				Hg2[index4]   = Hg2[index4] :+ log(y[index4,3]) :* gam :* Gtemp0 :* (gam :* log(y[index4,3]) :+ 1)
			}
		}
	}
	
	if (hasxb) {
		//dxb dxb
		H[|moptimize_util_eq_indices(M,xeqn,xeqn)|] = 
			moptimize_util_matsum(M,xeqn,xeqn,Hxb2, 
				moptimize_util_sum(M,logl))

		//dxb dg
		H[|moptimize_util_eq_indices(M,xeqn,geqn)|] = 
			moptimize_util_matsum(M,xeqn,geqn,Hxbg,
				moptimize_util_sum(M,logl))
		H[|moptimize_util_eq_indices(M,geqn,xeqn)|] = 
			H[|moptimize_util_eq_indices(M,xeqn,geqn)|]'

		if (hastb) {
			//dxb dtb
			Hxt 	= J(Nxbs,Ntbs,.)
			Hxtsum 	= quadcolsum(Hxbtb,1)
			el = 1				
			for (e1=1;e1<=Nxbs;e1++) {
				for (e2=1;e2<=Ntbs;e2++) {
					Hxt[e1,e2] = Hxtsum[el++]
				}					
			}
			H[|moptimize_util_eq_indices(M,xeqn,teqn)|] = Hxt
			H[|moptimize_util_eq_indices(M,teqn,xeqn)|] = 
				H[|moptimize_util_eq_indices(M,xeqn,teqn)|]'
			
			//dtb dtb
			Ht	= J(Ntbs,Ntbs,.)
			Htsum 	= quadcolsum(Htb2,1)
			el 	= 1
			for (e1=1;e1<=Ntbs;e1++) {
				e2 = 1
				while (e2<=e1) {
					if (e1==e2) Ht[e1,e1] = Htsum[el++]
					else 	Ht[e2,e1] = Ht[e1,e2] = Htsum[el++]
					e2++
				}
			}
			H[|moptimize_util_eq_indices(M,teqn,teqn)|] = Ht
			
			//dt/dg
			H[|moptimize_util_eq_indices(M,geqn,teqn)|] = quadcolsum(Htbg,1)
			H[|moptimize_util_eq_indices(M,teqn,geqn)|] = 
				H[|moptimize_util_eq_indices(M,geqn,teqn)|]'
		}
		else {
			
		}
		
	}
	else {
		if (hastb) {
			//dtb dtb
			Ht	= J(Ntbs,Ntbs,.)
			Htsum 	= quadcolsum(Htb2,1)
			el 	= 1
			for (e1=1;e1<=Ntbs;e1++) {
				e2 = 1
				while (e2<=e1) {
					if (e1==e2) Ht[e1,e1] = Htsum[el++]
					else 	Ht[e2,e1] = Ht[e1,e2] = Htsum[el++]
					e2++
				}
			}
			H[|moptimize_util_eq_indices(M,teqn,teqn)|] = Ht
		}
		
	}
	
	//dg dg
	H[|moptimize_util_eq_indices(M,geqn,geqn)|] = moptimize_util_matsum(M,geqn,geqn, Hg2, moptimize_util_sum(M,logl))
	
// 	H[|moptimize_util_eq_indices(M, 1, 2)|] = moptimize_util_matsum(M, 1, 2, Hxbg, moptimize_util_sum(M,logl))
// 	H[|moptimize_util_eq_indices(M, 2, 1)|] = H[|moptimize_util_eq_indices(M, 1, 2)|]'
// 	H[|moptimize_util_eq_indices(M, 2, 2)|] = moptimize_util_matsum(M, 2, 2, Hg2, moptimize_util_sum(M,logl))
	return(logl)
}


`RM' _logl_weibull2(`TR' M, `RR' b, `gml' gml, `RM' G, `RM' H)
{	
	// Setup //
	
	//model info
	model 	= 1
	hasxb 	= gml.hasxb[model]
	if (hasxb) xeqn = gml.xeqn[model]
	hastb 	= gml.hastb[model]
	if (hastb) teqn = gml.teqn[model]
	haszb 	= gml.haszb[model]
	if (haszb) zeqn = gml.zeqn[model,.]
	haslt	= gml.hasltrunc[model]
	hasic	= gml.haslint[model]
	hasbh	= gml.hasbh[model,]
	Nrelevs = gml.Nrelevels
	
	//data
	y 	= uhtred_util_depvar(gml)
	Nobs	= uhtred_util_nobs(gml)

	//logl
	logl 	= (Nrelevs ? J(Nobs,gml.ndim[Nrelevs],0) : J(Nobs,1,0))

	//linear predictors
	if (hastb) {
		xzb   = uhtred_util_xzb(M,b,gml)
		xtzb  = xzb + uhtred_util_tb(M,b,gml)
	}
	else 	xtzb  = uhtred_util_xtzb(M,b,gml)
	gam     = exp(uhtred_util_dap(gml,1))

	// log likelihood //
	
	//exactly observed events -> log hazard function
	gml.survind 	= 1
	Nobs1 		= uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 		= uhtred_get_surv_index(gml)
		logl[index1,] 	= xtzb[index1,] :+ log(gam) :+ 
			(gam - 1) :* log(y[index1,1])
		if (hasbh[1]) {
			logl[index1,] = log(exp(logl[index1,]) :+ 
				uhtred_util_bhazard(gml))
		}
		if (hasbh[2]) {
			logl[index1,] = log(exp(logl[index1,]) :* 
				uhtred_util_bhazard(gml))
		}
	}

	//exactly observed events and/or right censoring -> survival function
	gml.survind 	= 2
	Nobs2 		= uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (hastb | hasbh[2]) {
			Ngq 	= gml.chip
			chq2 	= J(Nobs2,Ngq,0)
// 			gq 	= uhtred_glegendre(Ngq)
// 			y2	= y[index2,1] :/ 2
// 			if (haslt) {
// 				y02   = y[index2,3] :/ 2
// 				y2m02 = (y2 :- y02)
// 				qp2 = y2m02 :* J(Nobs2,1,gq[,1]') :+ y2 :+ y02
// 			}
// 			else {
// 				qp2 = y2 :* J(Nobs2,1,gq[,1]') :+ y2
// 			}

			for (q=1;q<=Ngq;q++) {
// 				chq2[,q] = uhtred_util_tb(M,b,gml,qp2[,q]) 
				chq2[,q] = uhtred_util_tb(M,b,gml,gml.chqp[,q]) 
			}
			if (hasbh[2]) {
				chq2 = chq2 :+ log(uhtred_util_bhazard(gml))
			}
// 			chq2 = exp(chq2 :+ (gam-1) :* log(qp2))
			chq2 = exp(chq2 :+ (gam-1) :* log(gml.chqp))
// 			if (haslt) 	chq = y2m02 :* (chq2 * gq[,2]) 
// 			else 		chq = y2 :* (chq2 * gq[,2]) 

			chq = gml.y2 :* (chq2 * gml.chqw)

			expxzb1   	= exp(xzb[index2,])
			expxzb1gam 	= expxzb1 :* gam
			logl[index2,] 	= logl[index2,] :- chq :* expxzb1gam
		}
		else {
			logl[index2,] = logl[index2,] :- exp(xtzb[index2,]) :* 
						y[index2,1] :^ gam
			if (haslt) {
				gml.survind 	= 4
				index4 		= uhtred_get_surv_index(gml)
				logl[index4,] 	= logl[index4,] :+ 
					exp(xtzb[index4,]) :* y[index4,3] :^ gam
			}
		}
	}

	//interval censoring -> cdf function
	gml.survind = 3
	Nobs3 = uhtred_util_nobs(gml)
	if (Nobs3) {
		//exit times
		index3 = uhtred_get_surv_index(gml)
		/*if (hast) {
			Ngq 	= gml.chip
			chq3 	= J(Nobs3,Ngq,.)
			gq 		= uhtred_glegendre(Ngq)
			qp3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,1]') :+ y[index3,1]:/2
			for (q=1;q<=Ngq;q++) {
				chq3[,q] = exp(uhtred_util_xzb(gml,qp3[,q])) :* gam :* qp3[,q] :^(gam-1)
			}
			ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
			logl[index3,] = -exp(-ch3)
		}
		else */
		logl[index3,] = -exp(-exp(xtzb[index3,]) :* y[index3,1] :^gam)
		
		//entry times
		gml.survind = 5
		if (gml.hasltrunc[model]) tind = 4
		else 			  tind = 3
		index5 = uhtred_get_surv_index(gml)
		/*if (hast) {
			Nobs5 	= uhtred_get_nobs(gml)
			chq5 	= J(Nobs5,Ngq,.)
			qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
			for (q=1;q<=Ngq;q++) {
				chq5[,q] = exp(uhtred_util_xzb(gml,qp5[,q])) :* gam :* qp5[,q] :^(gam-1)
			}
			ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
			logl[index5,] = logl[index5,] :+ exp(-ch5)
		}
		else*/
		logl[index5,] = logl[index5,] :+ exp(-exp(xtzb[index5,]) :* y[index5,tind] :^gam)
		//logL
		logl[index3,] = log(logl[index3,])
	}	
	
	if (gml.todo==0) return(logl)

	// score //
	
	gml.survind 	= 0

	//get coefficient indices
	if (hasxb) {
		sindex1 = merlin_util_score_indices(M,xeqn)
		Nxbs 	= cols(sindex1)
		if (hastb) geqn = 3
		else geqn = 2
	}
	else geqn = 2
	
	if (hastb) {
		sindex2 = merlin_util_score_indices(M,teqn)
		Ntbs 	= cols(sindex2)
	}
	sindex4 = merlin_util_score_indices(M,geqn)
	
	
	//exactly observed events -> hazard function
	if (Nobs1) {
		if (hasxb) G[index1,sindex1] = gml.X[index1,]
		if (hastb) G[index1,sindex2] = gml.XT[index1,]
		G[index1,sindex4] = 1 :+ gam :* log(y[index1,1])
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		if (hastb) {
			dch2 = dch22 = J(Nobs2,1,0)
			if (haslt) {
				qw2  = (y[index2,1]:-y[index2,3]):/2 :* 
					J(Nobs2,1,gq[,2]')
			}
			else qw2  = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
			dchq2temp = chq2 :* qw2
			for (q=1;q<=Ngq;q++) {	
				dch2 = dch2 :+ dchq2temp[,q] :* 
					uhtred_util_t(gml,qp2[,q]) 
				dch22 = dch22 :+ dchq2temp[,q] :* 
					(gam:*log(qp2[,q]) :+ 1)
			}
			if (hasxb) {
				Gtemp = -chq :* expxzb1gam

				G[index2,sindex1] = G[index2,sindex1] :+ 
					Gtemp :* gml.X[index2,]
			}
			G[index2,sindex2] = G[index2,sindex2] :- 
						dch2 :* expxzb1gam
			G[index2,sindex4] = G[index2,sindex4] :- dch22 :* 
						expxzb1gam
		}
		else {
			Gtemp = exp(xtzb[index2,]) :* y[index2,1] :^ gam
			G[index2,sindex1] = G[index2,sindex1] :- 
				Gtemp :* gml.X[index2,]
			
			G[index2,sindex4] = G[index2,sindex4] :- 
				Gtemp :* gam :* log(y[index2,1])
			if (haslt) {
				Gtemp0 = exp(xtzb[index4,]) :* y[index4,3] :^ gam
				G[index4,sindex1] = G[index4,sindex1] :+ 
					Gtemp0 :* gml.X[index4,]
				G[index4,sindex4] = G[index4,sindex4] :+ 
					Gtemp0 :* gam :* log(y[index4,3])
			}
		}
	}

	if (gml.todo==1) return(logl)
	
	// Hessian //
	
	if (hasxb) {
		Hxb2 = Hxbg = J(gml.Nobs,1,0)
		if (hastb) {
			xtbind  = _uhtred_get_hessian_indices(Nxbs,Ntbs)
			Hxbtb  	= J(Nobs,cols(xtbind),0)
			tbind   = _uhtred_get_hessian_indices(Ntbs)
			Htb2    = J(Nobs,cols(tbind),0)
			Htbg 	= J(gml.Nobs,Ntbs,0)
		}
	}
	else {
		if (hastb) {
			Htb2 = Htbg = J(gml.Nobs,1,0)
		}
	}
	Hg2 = J(gml.Nobs,1,0)
	
	if (Nobs1) {
		 Hg2[index1] = gam :* log(y[index1,1])
	}  

	gml.survind = 2
	if (Nobs2) {
		if (hastb) {

			if (hasxb) {
				Hxb2[index2] = Gtemp
				Hxbg[index2] = -dch22 :* expxzb1gam
			}
			
			logqp2	= log(qp2)
			if (!haslt) {
				hazq = y[index2,1] :/ 2 :* 
					J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* 
					gam :* (gam :* logqp2 :+ 1)
			}
			else {
				hazq = (y[index2,1] :- y[index2,3]) :/ 2 :* 
					J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* 
					gam :* (gam :* logqp2 :+ 1)
			}
			hazq2 = logqp2 :* gam :* qp2:^(gam:-1) :* 
				((gam:^2 :* logqp2 :+ 2:*gam) :+ gam :+ 
				1:/logqp2)

			dchtg = J(Nobs2,1,0)
			chdxdt = chdt2 = J(Nobs2,1,0)
			for (q=1;q<=Ngq;q++) {	
				tq = uhtred_util_t(gml,qp2[,q]) 
				exptqb = exp(uhtred_util_tb(M,b,gml,qp2[,q]))
				chdxdt = chdxdt :+ dchq2temp[,q] :* 
					gml.X[index2,xtbind[1,]] :* 
					tq[index2,xtbind[2,]]
				chdt2 = chdt2 :+ dchq2temp[,q] :* 
					tq[index2,tbind[1,]] :* 
					tq[index2,tbind[2,]]
				dchtg = dchtg :+ exptqb :* hazq[,q] :* 
					tq
				hazq2[,q] = hazq2[,q] :* exptqb

			}
			Hxbtb[index2,] 	= - chdxdt :* expxzb1gam
			Htb2[index2,] 	= - chdt2 :* expxzb1gam
			Htbg[index2,] 	= - dchtg :* exp(xzb[index2])
			
			if (!haslt) {
				Hg2[index2] = Hg2[index2] :- 
					y[index2,1]:/2 :* (hazq2 * gq[,2]) :* 
					exp(xzb)
			}
			else {
				Hg2[index2] = Hg2[index2] :- 
					(y[index2,1] :- y[index2,3]) :/ 2 :* 
					(hazq2 * gq[,2]) :* exp(xzb)
			}
		}
		else {
			Hxb2[index2,] = -Gtemp
			Hxbg[index2,] = -Gtemp :* gam :* log(y[index2,1]) 
			Hg2[index2] = Hg2[index2] :- log(y[index2,1]) :* gam :* Gtemp :* (gam :* log(y[index2,1]) :+ 1)
			if (haslt) {
				Hxb2[index4,] = Hxb2[index4,] :+ Gtemp0
				Hxbg[index4,] = Hxbg[index4,] :+ Gtemp0 :* gam :* log(y[index4,3])
				Hg2[index4]   = Hg2[index4] :+ log(y[index4,3]) :* gam :* Gtemp0 :* (gam :* log(y[index4,3]) :+ 1)
			}
		}
	}
	
	if (hasxb) {
		//dxb dxb
		H[|moptimize_util_eq_indices(M,xeqn,xeqn)|] = 
			moptimize_util_matsum(M,xeqn,xeqn,Hxb2, 
				moptimize_util_sum(M,logl))

		//dxb dg
		H[|moptimize_util_eq_indices(M,xeqn,geqn)|] = 
			moptimize_util_matsum(M,xeqn,geqn,Hxbg,
				moptimize_util_sum(M,logl))
		H[|moptimize_util_eq_indices(M,geqn,xeqn)|] = 
			H[|moptimize_util_eq_indices(M,xeqn,geqn)|]'

		if (hastb) {
			//dxb dtb
			Hxt 	= J(Nxbs,Ntbs,.)
			Hxtsum 	= quadcolsum(Hxbtb,1)
			el = 1				
			for (e1=1;e1<=Nxbs;e1++) {
				for (e2=1;e2<=Ntbs;e2++) {
					Hxt[e1,e2] = Hxtsum[el++]
				}					
			}
			H[|moptimize_util_eq_indices(M,xeqn,teqn)|] = Hxt
			H[|moptimize_util_eq_indices(M,teqn,xeqn)|] = 
				H[|moptimize_util_eq_indices(M,xeqn,teqn)|]'
			
			//dtb dtb
			Ht	= J(Ntbs,Ntbs,.)
			Htsum 	= quadcolsum(Htb2,1)
			el 	= 1
			for (e1=1;e1<=Ntbs;e1++) {
				e2 = 1
				while (e2<=e1) {
					if (e1==e2) Ht[e1,e1] = Htsum[el++]
					else 	Ht[e2,e1] = Ht[e1,e2] = Htsum[el++]
					e2++
				}
			}
			H[|moptimize_util_eq_indices(M,teqn,teqn)|] = Ht
			
			//dt/dg
			H[|moptimize_util_eq_indices(M,geqn,teqn)|] = quadcolsum(Htbg,1)
			H[|moptimize_util_eq_indices(M,teqn,geqn)|] = 
				H[|moptimize_util_eq_indices(M,geqn,teqn)|]'
		}
		else {
			
		}
		
	}
	else {
		if (hastb) {
			//dtb dtb
			Ht	= J(Ntbs,Ntbs,.)
			Htsum 	= quadcolsum(Htb2,1)
			el 	= 1
			for (e1=1;e1<=Ntbs;e1++) {
				e2 = 1
				while (e2<=e1) {
					if (e1==e2) Ht[e1,e1] = Htsum[el++]
					else 	Ht[e2,e1] = Ht[e1,e2] = Htsum[el++]
					e2++
				}
			}
			H[|moptimize_util_eq_indices(M,teqn,teqn)|] = Ht
		}
		
	}
	
	//dg dg
	H[|moptimize_util_eq_indices(M,geqn,geqn)|] = moptimize_util_matsum(M,geqn,geqn, Hg2, moptimize_util_sum(M,logl))
	
// 	H[|moptimize_util_eq_indices(M, 1, 2)|] = moptimize_util_matsum(M, 1, 2, Hxbg, moptimize_util_sum(M,logl))
// 	H[|moptimize_util_eq_indices(M, 2, 1)|] = H[|moptimize_util_eq_indices(M, 1, 2)|]'
// 	H[|moptimize_util_eq_indices(M, 2, 2)|] = moptimize_util_matsum(M, 2, 2, Hg2, moptimize_util_sum(M,logl))
	return(logl)
}

// `RM' uhtred_weibull_logh(`gml' gml,`RC' t)
// {
// 	gam 	= asarray(gml.distancb,(gml.model,1))
// 	logh 	= uhtred_util_xzb(gml,t) :+ log(gam) :+ (gam - 1) :* log(t)
// 	if (gml.hasbh[gml.model,1]) {
// 		logh = log(exp(logh) :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))])
// 	}
// 	return(logh)
// }
//
// `RM' uhtred_weibull_ch(`gml' gml,`RC' t, | `RC' t0)
// {
// 	gam 	= asarray(gml.distancb,(gml.model,1))
// 	hast0 	= args()==3
//	
// 	if (gml.NI[gml.model]) {
// 		nobs 	= rows(t)
// 		ch 	= J(nobs,1,0)
// 		Ngq 	= gml.chip
// 		gq 	= uhtred_glegendre(Ngq)
// 		qp	= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
// 		loghazq = log(gam) :+ (gam-1) :* log(qp) :+ log(t:/2 :* J(nobs,1,gq[,2]'))
//		
// 		if (!hast0) {
// 			for (q=1;q<=Ngq;q++) {
// 				ch = ch :+ exp(uhtred_util_xzb(gml,qp[,q]) :+ loghazq[,q])
// 			}
// 		}
// 		else {
// 			for (q=1;q<=Ngq;q++) {
// 				ch = ch :+ exp(uhtred_util_xzb(gml,qp[,q],t0) :+ loghazq[,q])
// 			}
// 		}
// 		return(ch)
// 	}
// 	else {
// 		if (!hast0)	return(exp(uhtred_util_xzb(gml,t)) :* t :^ gam)
// 		else 		return(exp(uhtred_util_xzb(gml,t,t0)) :* t :^ gam)
// 	}
// }
//
// `RM' uhtred_weibull_cdf(`gml' gml,`RC' t)
// {
// 	return(1:-uhtred_weibull_s(gml,t))
// }
//
// `RM' uhtred_weibull_s(`gml' gml, `RC' t)
// {
// 	return(exp(-uhtred_weibull_ch(gml,t))) 	
// }
//

end
