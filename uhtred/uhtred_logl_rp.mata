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


/*
-> functions for family(rp, ...)
*/

mata:

`RM' _logl_rp(`TR' M, `RR' b, `gml' gml, | `RM' G, `RM' H)
{	
	// Setup //
	
	//model info
	model 	= gml.model
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
		xtzb  = xzb :+ uhtred_util_tb(M,b,gml)
		brcs  = b[|uhtred_util_bindices(gml,teqn)|]'
	}
	else 	xtzb  = uhtred_util_xtzb(M,b,gml)
	expxtzb	= exp(xtzb)

	// log likelihood //
		
	//exactly observed events -> hazard function
	gml.survind = 1
	Nobs1 = uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 		= uhtred_get_surv_index(gml)
		dXT		= asarray(gml.dXT,model)[index1,]
		dxb 		= dXT * brcs
		logl[index1,] 	= xtzb[index1,] :+ log(dxb) :- log(y[index1,1])
		if (hasbh[1]) {
			totalh1		= exp(logl[index1,]) :+ uhtred_util_bhazard(gml)
			logl[index1,] 	= log(totalh1)
		}
		if (hasbh[2]) {
			logl[index1,] = log(expxtzb[index1,] :* uhtred_util_bhazard(gml) :+ 
				exp(logl[index1]) :* uhtred_util_bHazard(gml))
		}
	}

	//exactly observed events and/or right censoring -> CH function
	gml.survind = 2
	Nobs2 = uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (hasbh[2])	{
			logl[index2,] = logl[index2,] :- 
				expxtzb[index2,] :* 
				uhtred_util_bHazard(gml)
		}
		else	logl[index2,] = logl[index2,] :- expxtzb[index2,]
		if (haslt) {
			gml.survind = 4
			index4 = uhtred_get_surv_index(gml)
			XT0    = asarray(gml.XT0,model)[index4,]
			expxtzb0 = exp(xzb[index4] + 
				XT0 * brcs)
			if (hasbh[2]) {
				logl[index4,] = logl[index4,] :+ 
					expxtzb0 :* uhtred_util_bHazard(gml)
			}
			else	logl[index4,] = logl[index4,] :+ expxtzb0
		}
	}

		//interval censoring -> cdf function
		gml.survind = 3
		Nobs3 = uhtred_util_nobs(gml)
		if (Nobs3) {
			index3 	= uhtred_get_surv_index(gml)
			//exit times
			logl[index3,] 	= 1:-exp(-expxb[index3])
			//entry times
			gml.survind = 5
			index5 	= uhtred_get_surv_index(gml)
			expxbl0	= exp(xb1[index5] :+ gml.XTL[index5,] * brcs)
			logl[index5,] = logl[index5,] :- (1:-exp(-expxbl0))
			//logL
			logl[index3,] = log(logl[index3,])
		}

	if (gml.todo==0) return(logl)

	//===================================================================//
	// score

	if (hasxb) {
		sindex1 = uhtred_util_score_indices(M,xeqn)
		X = asarray(gml.X,model)
	}
	if (hastb) {
		sindex2 = uhtred_util_score_indices(M,teqn)
		XT = asarray(gml.XT,model)
	}
	
	//exactly observed events -> hazard function
	if (Nobs1) {		
		if (hasbh[1]) {
			c1 = expxtzb[index1] :/ totalh1 :/ y[index1,1]
			G[index1,sindex1] = G[index1,sindex1] :+ 
						c1 :* X[index1,] :* dxb
			G[index1,sindex2] = G[index1,sindex2] :+ 
						c1 :* (dXT :+ 
						XT[index1,] :* dxb)
		}
		else {
			G[index1,sindex1] = G[index1,sindex1] :+ X[index1,]
			G[index1,sindex2] = G[index1,sindex2] :+ 
						XT[index1,] :+ 
						dXT :/ dxb
		}
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		G[index2,sindex1] = G[index2,sindex1] :- 
					X[index2,] :* expxtzb[index2]
		G[index2,sindex2] = G[index2,sindex2] :- 
					XT[index2,] :* expxtzb[index2]
		if (haslt) {
			G[index4,sindex1] = G[index4,sindex1] :+ 
						X[index4,] :* expxtzb0
			G[index4,sindex2] = G[index4,sindex2] :+ 
						XT0 :* expxtzb0
		} 
	}

	//interval censoring -> cdf function
	if (Nobs3) {
		//exit times
		expexpxb3 		= exp(-expxb[index3]) :* expxb[index3]
		G[index3,sindex1] 	= G[index3,sindex1] :+ 
						expexpxb3 :* X[index3,]
		G[index3,sindex2] 	= expexpxb3 :* gml.XT[index3,]
		//entry times
		expexpxbl05 		= exp(-expxbl0) :* expxbl0
		G[index5,sindex1] 	= G[index5,sindex1] :- 
						expexpxbl05 :* X[index5,]
		G[index5,sindex2] 	= G[index5,sindex2] :- expexpxbl05 :* gml.XTL[index5,]
		//finish up
		explogl3 		= exp(logl[index3,])
		G[index3,] 		= G[index3,] :/ explogl3
	}

	if (gml.todo==1) return(logl)	
	
	Nxb = cols(sindex1)
	xindex1 = xindex2 = J(1,0,.)
	for (i=1; i<=Nxb; i++) {
		refind = 1
		while (refind<=i) {
			xindex1 = xindex1,i
			xindex2 = xindex2,refind
			refind++
		}
	}
	Hxb2 = J(Nobs,cols(xindex2),0)

	//xbrcs indices
	xindex3 = xindex4 = J(1,0,.)
	Nrcsb = cols(sindex2)
	for (i=1; i<=Nxb; i++) {
		for (j=1;j<=Nrcsb;j++) {
			xindex3 = xindex3,i
			xindex4 = xindex4,j
		}
	}
	Hxbrcs = J(Nobs,cols(xindex3),0)
	
	//rcs indices
	xindex5 = xindex6 = J(1,0,.)
	
	for (i=1; i<=Nrcsb; i++) {
		refind = 1
		while (refind<=i) {
			xindex5 = xindex5,i
			xindex6 = xindex6,refind
			refind++
		}
	}
	Hrcs2 = J(Nobs,cols(xindex6),0)
		
	//exactly observed events -> hazard function
	if (Nobs1) {
		if (hasbh[1]) {
			Hxb2[index1] = c1 :* dxb :* 
				(1 :- (expxtzb[index1] :* dxb :/ 
				y[index1,1]):/totalh1) 
			Hxbrcs[index1,] = c1 :* X[index1,xindex3] :* 
				(dXT[,xindex4] :+ dxb :* XT[index1,xindex4] :-
				dxb :/ totalh1 :* (expxtzb[index1] :* dxb :/ y[index1,1] :* 
				(XT[index1,xindex4] + dXT[,xindex4]:/dxb)))
			Hrcs2[index1,] = c1 :* 
				((XT[index1,xindex5] :* dXT[,xindex6]) :+
				(dXT[,xindex5] :+ XT[index1,xindex5] :* dxb) :* XT[index1,xindex6] :- 
				(dXT[,xindex5] :+ XT[index1,xindex5] :* dxb) :/ totalh1 :* 
				(expxtzb[index1] :* dxb :/ y[index1,1] :* 
				(XT[index1,xindex6] + dXT[index1,xindex6]:/dxb)))
		}
		else {
			Hrcs2[index1,] = - dXT[,xindex5] :* 
				dXT[,xindex6] :/ (dxb:^2)
		}
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		//xb
		Hxb2[index2,] 	= Hxb2[index2,] :- expxtzb[index2]
		Hxbrcs[index2,] = Hxbrcs[index2,] :- 
					X[index2,xindex3] :* 
					XT[index2,xindex4] :* expxtzb[index2]
		Hrcs2[index2,] 	= Hrcs2[index2,] :- XT[index2,xindex5] :* 
					XT[index2,xindex6] :* expxtzb[index2]
		if (haslt) {
			Hxb2[index4,] 	= Hxb2[index4,] :+ expxtzb0
			Hxbrcs[index4,] = Hxbrcs[index4,] :+ 
				X[index4,xindex3] :* 
				XT0[index4,xindex4] :* expxtzb0[index4]
			Hrcs2[index4,] 	= Hrcs2[index4,] :+ 
				XT0[index4,xindex5] :* 
				XT0[index4,xindex6] :* expxtzb0
		}
	}

	if (hasxb) {
		Hxbf = J(Nxb,Nxb,.)
		Hxbsum = quadcolsum(Hxb2 :* X[,xindex1] :* X[,xindex2],1)
		el 	= 1
		for (e1=1;e1<=Nxb;e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) Hxbf[e1,e1] = Hxbsum[el++]
				else 	Hxbf[e2,e1] = Hxbf[e1,e2] = Hxbsum[el++]
				e2++
			}
		}
		H[|uhtred_util_bindices(gml,xeqn,xeqn)|] = Hxbf

		Hoff = J(Nxb,Nrcsb,.)
		Hxbrcssum = quadcolsum(Hxbrcs,1)
		el = 1				
		for (e1=1;e1<=Nxb;e1++) {
			for (e2=1;e2<=Nrcsb;e2++) {
				Hoff[e1,e2] = Hxbrcssum[el++]
			}					
		}
		H[|uhtred_util_bindices(gml,xeqn,teqn)|] = Hoff
		H[|uhtred_util_bindices(gml,teqn,xeqn)|] = H[|uhtred_util_bindices(gml,xeqn,teqn)|]'
	}

	H2	= J(Nrcsb,Nrcsb,.)
	Hrcssum = quadcolsum(Hrcs2,1)
	el 	= 1
	for (e1=1;e1<=Nrcsb;e1++) {
		e2 = 1
		while (e2<=e1) {
			if (e1==e2) H2[e1,e1] = Hrcssum[el++]
			else 	H2[e2,e1] = H2[e1,e2] = Hrcssum[el++]
			e2++
		}
	}
	H[|uhtred_util_bindices(gml,teqn,teqn)|] = H2
	return(logl)
}

end
