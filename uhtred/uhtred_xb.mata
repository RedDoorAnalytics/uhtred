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

//Evaluate linear predictors
void uhtred_xb(`TR' M, `gml' gml, `RR' b)
{
	eqnind = 1
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.hasxb[i]) eqnind++
		if (gml.hastb[i]) eqnind++
		if (gml.haszb[i]) eqnind = eqnind + gml.Nrelevels
// 		if (gml.Ndap[i]) {
// 			for (k=1;k<=gml.Ndap[i];k++) {
// 				eqnind2 = moptimize_util_eq_indices(M,eqnind++)[1,2]
// 				asarray(gml.distancb,(i,k),b[eqnind2])
// 			}
// 		}
// 		if (gml.Nap[i]) {
// 			for (k=1;k<=gml.Nap[i];k++) {
// 				eqnind2 = moptimize_util_eq_indices(M,eqnind++)[1,2]
// 				asarray(gml.apxb,(i,k),b[eqnind2])
// 				eqnind++
// 			}
// 		}
	}
	if (gml.Nrelevels) {
		eqnind = moptimize_util_eq_indices(M,eqnind)[1,2]
		uhtred_fillvcv(gml,b,eqnind)
	}
}

//Evaluate linear predictors for predict
void uhtred_p_xb(`gml' gml, `RR' b)
{
	eqnind = 1
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.hasxb[i]) eqnind++
		if (gml.hastb[i]) eqnind++
		if (gml.haszb[i]) eqnind = eqnind + gml.Nrelevels
// 		if (gml.Ndap[i]) {
// 			for (k=1;k<=gml.Ndap[i];k++) {
// 				eqnind2 = uhtred_util_bindices(gml,eqnind++)[1,2]
// 				asarray(gml.distancb,(i,k),b[eqnind2])
// 			}
// 		}
// 		if (gml.Nap[i]) {
// 			for (k=1;k<=gml.Nap[i];k++) {
// 				eqnind2 = uhtred_util_bindices(gml,eqnind++)[1,2]
// 				asarray(gml.apxb,(i,k),b[eqnind2])
// 				eqnind++
// 			}
// 		}
	}
	if (gml.Nrelevels) {
		eqnind = uhtred_util_bindices(gml,eqnind-1)[1,2] + 1
		uhtred_fillvcv(gml,b,eqnind)
	}
}

void uhtred_fillvcv(`gml' gml, 
		`RR' b,
		`RS' eqnind)
{
	`RM' covariances, vcv, sdcor
	`RC' Nres
	`RS' adapt

	covariances = gml.covariances	//indep \ exch \ unstr
	Nres 		= gml.Nres
	adapt 		= gml.adapt

	for (i=1 ; i < gml.Nlevels ; i++) {
		if (Nres[i]) {
			sdcor = vcv = asarray(gml.vcvs,i)
			if (covariances[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					sdcor[j,j] = exp(b[eqnind++])
					vcv[j,j] = sdcor[j,j] :^2
				}
			}
			else if (covariances[2,i]) {
				var_xb = exp(b[eqnind++]):^2
				cov_xb = var_xb :* tanh(b[eqnind++])
				for (j=1;j<=Nres[i];j++) {
					ind = 1
					while (ind<j) {
						vcv[ind,j] = vcv[j,ind] = cov_xb
						ind++  
					}	
					vcv[j,j] = var_xb
				}
			}
			else if (covariances[3,i]) {
				for (j=1;j<=Nres[i];j++) {
					sdcor[j,j] = exp(b[eqnind++])
					vcv[j,j] = sdcor[j,j] :^2
				}
				for (j=1;j<=Nres[i];j++) {
					ind = 1
					while (ind<j) {
						sdcor[ind,j] = sdcor[j,ind] = tanh(b[eqnind++])
						ind++
					}
					ind = 1
					while (ind<j) {
						vcv[ind,j] = vcv[j,ind] = sdcor[ind,ind]:*sdcor[j,j]:*sdcor[ind,j]
						ind++
					}
				}
			
			}
			else {
				var_xb = exp(b[eqnind++]):^2
				for (j=1;j<=Nres[i];j++) vcv[j,j] = var_xb
			}			
			asarray(gml.vcvs,i,vcv)
			if (!gml.predict) asarray(gml.Pgml->vcvs,i,vcv)
			uhtred_update_ip(gml,i)
		}
	}
}


void uhtred_update_ip(`gml' gml, `RS' i)
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

end
