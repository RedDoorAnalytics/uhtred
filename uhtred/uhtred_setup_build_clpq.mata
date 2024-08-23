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

version 18

mata:

void uhtred_build_clpq(`gml' gml)
{
	if (!gml.NI) return
	
	gml.chq = asarray_create("real",2)
	Ngq 	= gml.chip
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		gml.model 	= mod
		gml.survind 	= 2
		Nobs2 		= uhtred_util_nobs(gml)
		index2 		= uhtred_get_surv_index(gml)
		haslt		= gml.hasltrunc[mod]
		gq 		= uhtred_glegendre(Ngq)	
		y  		= uhtred_util_depvar(gml)
		y2		= y[index2,1] :/ 2
		if (haslt) {
			y02   = y[index2,3] :/ 2
			y2m02 = (y2 :- y02)
			qp2 = y2m02 :* J(Nobs2,1,gq[,1]') :+ y2 :+ y02
		}
		else 	qp2 = y2 :* J(Nobs2,1,gq[,1]') :+ y2
		
		for (q=1;q<=Ngq;q++) {
			asarray(gml.chq,(mod,q),uhtred_util_t(gml,qp2[,q]))
		}	
	}
}

end
