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

//gf core function
void uhtred_gf(	`TS' M,
                `RS' todo,
                `RR' b,
                `RC' lnfi,
                `RM' G,
                `RM' H)
{
	`gml' gml

	gml  		= moptimize_util_userinfo(M,1)
	gml.todo 	= todo
	gml.myb		= b
	gml.survind     = 0
	
	uhtred_xb(M,gml,b)
	
	if (gml.Nrelevels) {	
		if (gml.todo) {
			nb = cols(b)
			G = J(gml.Nobs[1],nb,0)
			if (gml.todo==2) {
				H = J(nb,nb,0)
			}
		}
		`pgml' Pgml
		Pgml = findexternal(gml.GML)
		Pgml->lnfi1 = lnfi = uhtred_logl_panels(1,M,b,gml)
			
		if (todo==0) return
                        
			gml.survind = 0
			uhtred_score(M,b,G,lnfi,gml)
			if (todo==1) return
			
				gml.survind = 0
				uhtred_hessian_panels(M,b,G,H,lnfi,gml)
				return
	}
	
	if (gml.todo) {
		nb = cols(b)
		G = J(gml.N,nb,0)
		if (gml.todo==2) {
			H = J(nb,nb,0)
		}
	}

	lnfi = uhtred_logl_ob(M,b,gml,G,H)
}

`RC' uhtred_logl_ob(`TR' M, `RR' b, `gml' gml, `RM' G, `RM' H)
{
	`RC' logl, index
	logl = J(gml.N,1,0)

	for (j=1;j<=gml.Nmodels;j++) {
		gml.model   = gml.modtoind = j
		gml.survind = 0
		index = uhtred_util_index(gml)
		logl[index] = logl[index] :+ (*gml.Plnl[j])(M,b,gml,G,H)
	}
	return(logl)
}

`RC' uhtred_logl_panels(`RS' index,	///	-level-
			`TR' M, `RR' b, `gml' gml)	//	
{
	`RS' index2
	`RM' res, panelindex

	index2 = index+1
	res = J(gml.Nobs[index,1],gml.ndim[index],0)

	if (index<gml.Nrelevels) {
		panelindex = asarray(gml.panelindexes,(index,1))
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			res[,q] = panelsum(uhtred_logl_panels(index2,M,b,gml),panelindex)
		}
	}
	else {
		j = 1
		gml.survind = 0
		res2 = (*gml.Plnl[j])(M,b,gml)
		res = res :+ panelsum(res2,asarray(gml.panelindexes,(index,j)))
	}
	
	if (gml.adapt[index]) 	res = res :+ asarray(gml.aghlogl,index)
	if (gml.usegh[index]) 	c   = rowmax(res :+ log(asarray(gml.baseGHweights,index)'))
	else			c   = rowmax(res)
	expres 	= exp(res :- c)

	if (gml.adapt[index] | gml.todo) {
		if (index==1) gml.qind[index2] = 0
		asarray(gml.Li_ip,gml.qind,expres) 
		if (gml.todo) asarray(gml.Li_ip2,gml.qind,exp(res)) 
	}

	if (gml.usegh[index]) {			//GHQ
		return(c :+ log(expres * asarray(gml.baseGHweights,index)))
	}
	else {					//MCI
		return(c :+ log(quadrowsum(expres):/gml.ndim[index]))	
	}
}

end
