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

`RM' _score_rp_xb(`gml' gml, `RS' xindex)
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

	//linear predictors
	if (hastb) {
		xzb   = uhtred_util_p_xzb(gml)
		xtzb  = xzb :+ uhtred_util_p_tb(gml)
		brcs  = gml.myb[|uhtred_util_bindices(gml,teqn)|]'
	}
	else 	xtzb  = uhtred_util_p_xtzb(gml)
	expxtzb	= exp(xtzb)

	// log likelihood //
		
	//exactly observed events -> hazard function
	gml.survind = 1
	Nobs1 = uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 	= uhtred_get_surv_index(gml)
	}

	//exactly observed events and/or right censoring -> CH function
	gml.survind = 2
	Nobs2 = uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (haslt) {
			gml.survind = 4
			index4 = uhtred_get_surv_index(gml)
			expxtzb0 = exp(xzb[index4,] :+ gml.XT0[index4,] * brcs)
		}
	}

	//===================================================================//
	// score

	res = J(gml.Nobs[gml.Nlevels,1],gml.ndim[gml.Nrelevels],0)

	//exactly observed events -> hazard function
	if (Nobs1) {		
		res[index1,] = res[index1,] :+ gml.X[index1,xindex]
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		res[index2,] = res[index2,] :- gml.X[index2,xindex] :* expxtzb[index2,]
		if (haslt) {
			res[index4,] = res[index4,] :+ gml.X[index4,xindex] :* expxtzb0
		} 
	}
	return(res)
}

`RM' _score_rp_tb(`gml' gml, `RS' xindex)
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

	//linear predictors
	xzb   = uhtred_util_p_xzb(gml)
	xtzb  = xzb :+ uhtred_util_p_tb(gml)
	brcs  = gml.myb[|uhtred_util_bindices(gml,teqn)|]'
	expxtzb	= exp(xtzb)

	// log likelihood //
		
	//exactly observed events -> hazard function
	gml.survind = 1
	Nobs1 = uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 		= uhtred_get_surv_index(gml)
		dXT		= gml.dXT[index1,]
		dxb 		= dXT * brcs
	}

	//exactly observed events and/or right censoring -> CH function
	gml.survind = 2
	Nobs2 = uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (haslt) {
			gml.survind = 4
			index4 = uhtred_get_surv_index(gml)
			XT0    = gml.XT0[index4,]
			expxtzb0 = exp(xzb[index4,] :+ XT0 * brcs)
		}
	}

	//===================================================================//
	// score

	res 	= J(gml.Nobs[gml.Nlevels,1],gml.ndim[gml.Nrelevels],0)
	XT 	= gml.XT[,xindex]
	
	//exactly observed events -> hazard function
	if (Nobs1) {		
		res[index1,] = res[index1,] :+ XT[index1,] :+ 
						dXT[,xindex] :/ dxb
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		res[index2,] = res[index2,] :- 
					XT[index2,] :* expxtzb[index2,]
		if (haslt) {
			res[index4,] = res[index4,] :+ 
						XT0[,xindex] :* expxtzb0
		} 
	}
	
	return(res)
}

`RM' _score_rp_vcv(`gml' gml,`RS' xindex)
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

	//linear predictors
	if (hastb) {
		xzb   = uhtred_util_p_xzb(gml)
		xtzb  = xzb :+ uhtred_util_p_tb(gml)
		brcs  = gml.myb[|uhtred_util_bindices(gml,teqn)|]'
	}
	else 	xtzb  = uhtred_util_p_xtzb(gml)
	expxtzb	= exp(xtzb)
	
	zb = uhtred_util_p_zb_lev(gml,xindex)

	// log likelihood //
		
	//exactly observed events -> hazard function
	gml.survind = 1
	Nobs1 = uhtred_util_nobs(gml)
	if (Nobs1) {
		index1 = uhtred_get_surv_index(gml)
	}

	//exactly observed events and/or right censoring -> CH function
	gml.survind = 2
	Nobs2 = uhtred_util_nobs(gml)
	if (Nobs2) {
		index2 	= uhtred_get_surv_index(gml)
		if (haslt) {
			gml.survind = 4
			index4 = uhtred_get_surv_index(gml)
			XT0    = gml.XT0[index4,]
			expxtzb0 = exp(xzb[index4,] :+ XT0 * brcs)
		}
	}

	//===================================================================//
	// score
	
	res 	= J(gml.Nobs[gml.Nlevels,1],gml.ndim[gml.Nrelevels],0)
	
	//exactly observed events -> hazard function
	if (Nobs1) {		
		res[index1,] = res[index1,] :+ zb[index1,]
	}  

	//exactly observed events and/or right censoring -> survival function
	if (Nobs2) {
		res[index2,] = res[index2,] :- zb[index2,] :* expxtzb[index2,]
		if (haslt) {
			res[index4,] = res[index4,] :+ 
						zb[index4,] :* expxtzb0

		} 
	}
	
	return(res)
}

end

