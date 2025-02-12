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

void uhtred_p_ebmeans(`gml' gml, | `RS' getses)
{


	getses 		= args()==2
	gml.blups 	= asarray_create("real",1)
	gml.fixedonly 	= 0

	for (i=1; i<=gml.Nrelevels; i++) {

		if ((i-1)>0) {
			//fix the previous levels posterior means
			gml.fixedlevels = 1::(i-1)
		}
		
		//adaptive updates for the current level, conditioning on 
		//any higher levels ebmeans (not nodes)
		
		oldlnl 	= quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
// asarray(gml.vcvs,i)
// oldlnl
		_ebmeans_update(gml,i)
		exit(1986)
		newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
// newlnl
		while (reldif(oldlnl,newlnl)>gml.atol) {
			swap(oldlnl,newlnl)
			_ebmeans_update(gml,i)
			newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		}	
		
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
	
	
}

// Update adaptive quadrature locations and scales
void _ebmeans_update(`gml' gml,`RS' i)
{
// 	for (i=1; i<=gml.Nrelevels; i++) {
i
		ndim = gml.ndim[i]
		ind1 = 1; ind2 = ndim
		stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		detvcv 		= J(gml.Nobs[i,1],1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
		numer 		= exp(log(asarray(gml.Li_ip,gml.qind)) :- log(L_i))	//c cancels out
		baseweights2 	= baseweights'
gml.qind
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			ipij 	 = asarray(gml.aghip,(i,j))

			newblups = numer[j,] * (ipij :* baseweights2)'
			newblups
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
// 	}
}


end

