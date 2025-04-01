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

// prolog
void uhtred_prolog(	`RR' b,		///
                        `TS' M, 	///
                        `RS' lnl) 
{
	`pgml' p
	p = &moptimize_util_userinfo(M,1)

	if ((*p).gridsearch) uhtred_vcv_gridsearch(M,*p,b)

	//post starting values
	if ((*p).iter==0) {
		p->svs = b
	}
	
	if (sum((*p).adapt) & (*p).iter<501 & (*p).iter<=(*p).adaptit) {

		`gml' 	gml
		`RS' 	oldlnl, newlnl , it
		
		gml = moptimize_util_userinfo(M,1)
		
		// need to fill up linear predictors, first time through
		gml.myb = b
		uhtred_xb(M,gml,b)
		
		gml.tofix = 0
		gml.fxls  = J(1,0,.)
		
		for (i=1; i<=gml.Nrelevels; i++) {
			
			if (i>1) {
				gml.tofix = 1
				gml.fxls = gml.fxls,(i-1)
			}
			
			oldlnl 	= quadsum(uhtred_logl_panels(i,M,b,gml),1)
			
			it = 0	
			if (gml.showadapt) uhtred_di_adaptiter(it,oldlnl)
			
			it = 1
			(*gml.Pupdateip)(gml,i)
			
			newlnl = quadsum(uhtred_logl_panels(i,M,b,gml),1)
			if (gml.showadapt) uhtred_di_adaptiter(it,newlnl)
			
			while (reldif(oldlnl,newlnl)>gml.atol & it<=1000) {
				it++
				swap(oldlnl,newlnl)
				(*gml.Pupdateip)(gml,i)
				newlnl = quadsum(uhtred_logl_panels(i,M,b,gml),1)
				if (gml.showadapt) uhtred_di_adaptiter(it,newlnl)
			}
		
			//update stored nodes
			uhtred_gh_post_ip(gml,i)
		}		
		gml.tofix = 0
	}

	p->iter = (*p).iter :+ 1	
}

/*
	display adapted log-likelihood
*/
void uhtred_di_adaptiter(`RS' iter, `RS' val)
{
	st_numscalar("mliter",iter)
	st_numscalar("ll",val)		
	stata(`"di as txt "-- Iteration " _col(14) mliter as txt ":" _col(19) "Adapted log likelihood = " as res %10.0g ll"')
}

void uhtred_gh_post_ip(`gml' gml, `RS' i)
{
	asarray(gml.Pgml->stackednodes,i,asarray(gml.stackednodes,i))
	for (j=1; j<=gml.Nobs[i,1]; j++) {
		asarray(gml.Pgml->aghip,(i,j),asarray(gml.aghip,(i,j)))
	}
	for (r=1; r<=gml.Nres[i]; r++) {
		asarray(gml.Pgml->aghip2,(i,r),asarray(gml.aghip2,(i,r)))
	}
	asarray(gml.Pgml->aghlogl,i,asarray(gml.aghlogl,i))

	//now store them so we can condition on them for lower levels
	blups 		= J(gml.Nobs[i,1],gml.Nres[i],.)
	baseweights 	= asarray(gml.baseGHweights,i)
	L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
	numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
	baseweights2 	= baseweights'	
	cholV		= cholesky(asarray(gml.vcvs,i))
	for (j=1; j<=gml.Nobs[i,1]; j++) {
		ipij 		= cholV * asarray(gml.aghip,(i,j))
		blups[j,] 	= numer[j,] * (ipij :* baseweights2)'
	}
	asarray(gml.blups,i,blups)	
}	


// Update adaptive quadrature locations and scales
void uhtred_gh_update_ip(`gml' gml,i)
{
	ndim = gml.ndim[i]
	ind1 = 1; ind2 = ndim
	stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
	detvcv 		= J(gml.Nobs[i,1],1,.)
	basenodes 	= asarray(gml.baseGHnodes,i)
	baseweights 	= asarray(gml.baseGHweights,i)
	L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
	numer 		= exp(log(asarray(gml.Li_ip,gml.qind)) :- log(L_i))	//c cancels out
	baseweights2 	= baseweights'
	
	for (j=1; j<=gml.Nobs[i,1]; j++) {
		ipij 	 = asarray(gml.aghip,(i,j))
		newblups = numer[j,] * (ipij :* baseweights2)'
		vcv_new  = (numer[j,] * 
				(uhtred_outerprod_by_col(ipij') :* 
				baseweights)) :- 
				newblups # newblups
		vcv_new  = rowshape(vcv_new,gml.Nres[i])
		nodes 	 = newblups' :+ cholesky(vcv_new) * basenodes
		asarray(gml.aghip,(i,j),nodes)
		detvcv[j] = det(vcv_new)
		stackednodes[|ind1,.\ind2,.|] = nodes'
		ind1 = ind1 + ndim
		ind2 = ind2 + ndim
	}
	
	//update logl extra contribution
	asarray(gml.aghlogl,
		i,
		rowshape(uhtred_lndmvnorm(stackednodes,I(gml.Nres[i])),
			 gml.Nobs[i,1]) 
		:+ log(detvcv):/2)
	//update stacked nodes, and RE specific stacked nodes
	asarray(gml.stackednodes,i,stackednodes')
	res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
	for (r=1;r<=gml.Nres[i];r++) {
		asarray(gml.aghip2,
			(i,r),
			rowshape(res[r,],gml.Nobs[i,1]))
	}	
}

end
