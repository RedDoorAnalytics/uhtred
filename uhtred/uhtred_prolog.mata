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

	if (sum((*p).adapt) & (*p).iter<501 & (*p).iter<=(*p).adaptit) {

		`gml' 	gml
		`RS' 	oldlnl, newlnl , it
		
		gml = moptimize_util_userinfo(M,1)
		
		// need to fill up linear predictors, first time through
		uhtred_xb(M,gml,b)
		gml.myb = b
		oldlnl 	= quadsum(uhtred_logl_panels(1,M,b,gml),1)
		
		it = 0	
		if (gml.showadapt) uhtred_di_adaptiter(it,oldlnl)
		
		it = 1
		(*gml.Pupdateip)(gml)
		newlnl = quadsum(uhtred_logl_panels(1,M,b,gml),1)
		if (gml.showadapt) uhtred_di_adaptiter(it,newlnl)
		
		while (reldif(oldlnl,newlnl)>gml.atol & it<=1000) {
			it++
			swap(oldlnl,newlnl)
			(*gml.Pupdateip)(gml)
			newlnl = quadsum(uhtred_logl_panels(1,M,b,gml),1)
			if (gml.showadapt) uhtred_di_adaptiter(it,newlnl)
		}
		
		//update stored nodes -> needs work!
		uhtred_gh_post_ip(gml)
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

void uhtred_mc_update(	`gml' gml, 	///
                        `RS' lev, 	///
                        `RR' b, 	///
                        `TS' M, 	///
                        `RS' lnl) 
{
	//increase in steps of ?
	gml.Pgml->ndim[lev] = gml.ndim[lev] = gml.ndim[lev] + 2			//must update external and passed struct
	uhtred_update_Zs_bs(gml,lev)	//!!re-write to add 3 new draws to existing
}

// Update adaptive quadrature locations and scales
void uhtred_gh_update_ip(`gml' gml)
{
	for (i=1; i<=1; i++) {
		ndim = gml.ndim[i]
		ind1 = 1; ind2 = ndim
		stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		detvcv 		= J(gml.Nobs[i,1],1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
		numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
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
}

void uhtred_gh_post_ip(`gml' gml)
{
	for (i=1; i<=1; i++) {
		asarray(gml.Pgml->stackednodes,i,asarray(gml.stackednodes,i))
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			asarray(gml.Pgml->aghip,(i,j),asarray(gml.aghip,(i,j)))
		}
		for (r=1; r<=gml.Nres[i]; r++) {
			asarray(gml.Pgml->aghip2,(i,r),asarray(gml.aghip2,(i,r)))
		}
		asarray(gml.Pgml->aghlogl,i,asarray(gml.aghlogl,i))
	}
}	

void uhtred_mc_update_ip(`gml' gml) 
{

	ip = gml.ip
	ndim = ip

	for (i=1; i<=1; i++) {
		
		ind1 = 1; ind2 = ndim
		stackednodes = J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		stackednodes2 = J(gml.Nobs[i,1],ndim,.)
		L_i = quadrowsum(asarray(gml.Li_ip,gml.qind))
		numer = asarray(gml.Li_ip,gml.qind) :/ L_i
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			ipij 		= asarray(gml.aghip,(i,j))
			newblups 	=  numer[j,] * (ipij')
			vcv_new 	= (numer[j,] * (uhtred_outerprod_by_col(ipij'))) :- rowshape((newblups')*newblups,1)
			vcv_new 	= rowshape(vcv_new,gml.Nres[i])
			
			/*rseed(gml.seed)
			nodes = uhtred_drawnorm(newblups',vcv_new,ip)
			asarray(gml.aghip,(i,j),nodes)
			nodes2 = cholesky(asarray(gml.vcvs,i)) * nodes
			asarray(gml.aghip2,(i,j),nodes2)*/

			nodes = newblups' :+ cholesky(vcv_new) * asarray(gml.bdraws,i)
			asarray(gml.aghip,(i,j),nodes)
			nodes2 = cholesky(asarray(gml.vcvs,i)) * nodes
			asarray(gml.aghip2,(i,j),nodes2)
			for (k=1;k<=ndim;k++) stackednodes2[j,k] = lnmvnormalden(newblups', vcv_new, nodes[,k])
			stackednodes[|ind1,.\ind2,.|] = nodes'
			ind1 = ind1 + ndim
			ind2 = ind2 + ndim
		}
		//update logl extra contribution
		//asarray(gml.aghlogl,i,rowshape(uhtred_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :+ logdetChol)
		
		//update stacked nodes, and RE specific stacked nodes
		asarray(gml.stackednodes,i,stackednodes')
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
		}
		
		asarray(gml.aghlogl,i,rowshape(uhtred_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :- stackednodes2)

	}	
}

void uhtred_update_ip_newNip(`gml' gml, `RS' i)
{
	//update nodes
	if (gml.adapt[i]) {
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)	
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
		}	
	}
	else {
		if (gml.usegh[i]) {
// 			if (gml.hasRExb[i]) {
// 				vcv = asarray(gml.vcvs,i)
// 				ndim = gml.ndim[i]
// 				ind1 = 1; ind2 = ndim
// 				stackednodes = J(gml.ndim[i]*gml.Nobs[i,1],gml.Nres[i],.)
// 				for (j=1;j<=gml.Nobs[i,1];j++) {
// 					nodes = cholesky(invvech(vcv[j,])) * asarray(gml.baseGHnodes,i)
// 					stackednodes[|ind1,.\ind2,.|] = nodes'
// 					ind1 = ind1 + ndim
// 					ind2 = ind2 + ndim
// 				}
// 				stackednodes = stackednodes'			//!!fix
// 				for (r=1;r<=gml.Nres[i];r++) {
// 					asarray(gml.aghip2,(i,r),rowshape(stackednodes[r,],gml.Nobs[i,1]))
// 				}				
// 			}
// 			else {
				bnodes = cholesky(asarray(gml.vcvs,i)) * asarray(gml.baseGHnodes,i)
				asarray(gml.b,i,bnodes)	
// 			}
		}
		else {
			asarray(gml.b,i,cholesky(asarray(gml.vcvs,i)) * asarray(gml.bdraws,i))
		}
			
	}	
}

// Update adaptive quadrature locations and scales
void uhtred_gh_update_ip_alllevs(`gml' gml)
{
	for (i=1; i<=gml.Nrelevels; i++) {

		ndim = gml.ndim[i]
		ind1 = 1; ind2 = ndim
		stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		detvcv 		= J(gml.Nobs[i,1],1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
		numer 		= exp(log(asarray(gml.Li_ip,gml.qind)) :- log(L_i))	//c cancels out
		baseweights2 	= baseweights'
//gml.qind
// st_view(res1=.,.,"res1","flag")
		
		
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			ipij 	 = asarray(gml.aghip,(i,j))

			newblups = numer[j,] * (ipij :* baseweights2)'
			res1[j,] = newblups
			res1[j,] = newblups
// 			if (missing(newblups)) -999
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
}

end
