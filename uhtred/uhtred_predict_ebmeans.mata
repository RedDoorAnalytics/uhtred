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

void uhtred_p_ebmeans(`gml' gml, `RS' getses)
{
	gml.blups 	= asarray_create("real",1)
	gml.fixedonly 	= 0
	
	
	gml.tofix = 0
	gml.fxls  = J(1,0,.)

	//start with initial starting values
	
	gml.myb = st_matrix("e(b_svs)")
	uhtred_p_xb(gml,gml.myb)
	gml.Pgml = &gml
		
	for (i=1; i<=gml.Nrelevels; i++) {

		if (i>1) {
			gml.tofix = 1
			gml.fxls = gml.fxls,(i-1)
		}
			
		//adaptive updates for the current level, conditioning on 
		//any higher levels ebmeans (not nodes)
		oldlnl 	= quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		_ebmeans_update(gml,i)
		
		newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		
		while (reldif(oldlnl,newlnl)>gml.atol) {
			swap(oldlnl,newlnl)
			_ebmeans_update(gml,i)
			newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		}	
		
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
	
	//now move on to estimated coefficients
	
	gml.myb = st_matrix("e(b)")
	uhtred_p_xb(gml,gml.myb)
	
	gml.tofix = 0
	gml.fxls  = J(1,0,.)

	for (i=1; i<=gml.Nrelevels; i++) {

		if (i>1) {
			gml.tofix = 1
			gml.fxls = gml.fxls,(i-1)
		}
			
		//adaptive updates for the current level, conditioning on 
		//any higher levels ebmeans (not nodes)
		oldlnl 	= quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		_ebmeans_update(gml,i)

		newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		
		while (reldif(oldlnl,newlnl)>gml.atol) {
			swap(oldlnl,newlnl)
			_ebmeans_update(gml,i)
			newlnl = quadsum(uhtred_logl_panels(i,M=.,gml.myb,gml),1)
		}	
		
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
	gml.tofix = 0
	
	if (getses) {
		for (i=1; i<=gml.Nrelevels; i++) {
			seblups 	= J(gml.Nobs[i,1],gml.Nres[i],.)
			baseweights 	= asarray(gml.baseGHweights,i)
			L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
			numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
			baseweights2 	= baseweights'	
			cholV		= cholesky(asarray(gml.vcvs,i))
			for (j=1; j<=gml.Nobs[i,1]; j++) {
				ipij 		= cholV * asarray(gml.aghip,(i,j))
				blups 		= numer[j,] * (ipij :* baseweights2)'
				vcv_new 	= (numer[j,] * (uhtred_outerprod_by_col(ipij') :* baseweights)) :- blups # blups
				seblups[j,] = sqrt(diagonal(rowshape(vcv_new,gml.Nres[i])))'
			}
			asarray(gml.blups,i,seblups)
		}
	}
	
}

// Update adaptive quadrature locations and scales
void _ebmeans_update(`gml' gml,`RS' i)
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


void uhtred_p_init_ip(`gml' gml)
{
	
	for (i=1;i<gml.Nlevels;i++) {
		
		if (gml.usegh[i]) {
						
			if (gml.Nres[i]) {
			
				qmat = _gauss_hermite_nodes(gml.ip[i])
				gml.ndim[i] = gml.ip[i]:^gml.Nres[i]
				
				x = J(gml.Nres[i],1,qmat[1,]):*sqrt(2)
				asarray(gml.baseGHnodes,i,uhtred_expand_matrix(x))
				
				w = J(gml.Nres[i],1,qmat[2,]):/sqrt(pi())
				asarray(gml.baseGHweights,i,uhtred_expand_matrix(w,1)')
				
				if (gml.adapt[i]) {
                                        //!! check cholesky() needed in below
					asarray(gml.stackednodes,i,J(1,gml.Nobs[i,1],asarray(gml.vcvs,i) * asarray(gml.baseGHnodes,i)))
					for (r=1;r<=gml.Nres[i];r++) {
						asarray(gml.aghip2,(i,r),rowshape(asarray(gml.stackednodes,i)[r,],gml.Nobs[i,1]))
					}	
					for (j=1; j<=gml.Nobs[i,1]; j++) {
						asarray(gml.aghip,(i,j),asarray(gml.vcvs,i) * asarray(gml.baseGHnodes,i))
					}
					asarray(gml.baseGHweights,i,((2:*pi()):^(gml.Nres[i]:/2):*exp(quadcolsum(asarray(gml.baseGHnodes,i):^2):/2) :* asarray(gml.baseGHweights,i)')')
					asarray(gml.aghlogl,i,J(gml.Nobs[i,1],1,uhtred_lndmvnorm(asarray(gml.baseGHnodes,i)',I(gml.Nres[i]))') :+ log(sqrt(det(asarray(gml.vcvs,i)))))

				}
				else {
					asarray(gml.b,i,asarray(gml.baseGHnodes,i))				
				}
			}
		}
		else {
			gml.ndim[i] = gml.ip[i]
			if (gml.adapt[i]) {
				rseed(gml.seed)			//reset seed
				baseip = invnormal(halton(gml.ndim[i],gml.Nres[i])')
				asarray(gml.bdraws,i,baseip)
				asarray(gml.stackednodes,i,J(1,gml.Nobs[i,1],baseip))
				for (r=1;r<=gml.Nres[i];r++) {
					asarray(gml.aghip2,(i,r),rowshape(asarray(gml.stackednodes,i)[r,],gml.Nobs[i,1]))
				}
				for (j=1; j<=gml.Nobs[i,]; j++) {
					asarray(gml.aghip,(i,j),baseip)
				}
				asarray(gml.aghlogl,i,J(gml.Nobs[i,1],gml.ndim[i],0))
			}
			else {
				rseed(gml.seed)			//reset seed
				asarray(gml.bdraws,i,invnormal(halton(gml.ndim[i],gml.Nres[i])'))
				asarray(gml.b,i,asarray(gml.bdraws,i))
			}
		
		}
	}
	//note; Npanels is indexed by level and model -> they are the same across all models (until ob level which isn't used here)
}

end

