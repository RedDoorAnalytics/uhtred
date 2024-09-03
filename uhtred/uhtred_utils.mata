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

`RC' uhtred_util_index(`gml' gml)
{
	if (gml.survind==-99) 	return(1::gml.N)
	else 			index = asarray(gml.xbindex,gml.modtoind)
	if (!gml.survind) 	return(index)
	else 		        return(index[uhtred_get_surv_index(gml)])
}

//get model specific index for survival outcomes
//can't if else, as it could be 0
`RC' uhtred_get_surv_index(`gml' gml)
{
	mod = gml.model
	if 	(gml.survind==1) return(asarray(gml.surv_index,(mod,1)))
	else if (gml.survind==2) return(asarray(gml.surv_index,(mod,2)))
	else if (gml.survind==3) return(asarray(gml.surv_index,(mod,3)))
	else if (gml.survind==5) return(asarray(gml.surv_index,(mod,5)))
	else if (gml.survind==4) return(asarray(gml.surv_index,(mod,4)))
	else if (gml.survind==6) return(asarray(gml.surv_index,(mod,6)))
	else if (gml.survind==7) return(asarray(gml.surv_index,(mod,7)))
}

`RS' uhtred_util_nobs(`gml' gml)
{
	mod = gml.model
	if 	(gml.survind==0) 	return(gml.Nobs[gml.Nlevels,mod])
	else if (gml.survind==1) 	return(gml.Nsurv[mod,1])
	else if (gml.survind==2) 	return(gml.Nsurv[mod,2])
	else if (gml.survind==3) 	return(gml.Nsurv[mod,3])
	else if (gml.survind==5) 	return(gml.Nsurv[mod,5])
	else if (gml.survind==4) 	return(gml.Nsurv[mod,4])
	else if (gml.survind==6) 	return(gml.Nsurv[mod,6])
	else if (gml.survind==7) 	return(gml.Nsurv[mod,7])
}

`RM' uhtred_util_xtzb(`TR' M, `RR' b, `gml' gml, | `RC' t)
{
	`RM' x
	`RS' hast, mod
	
	hast = args()==4
	mod  = gml.model
	xb   = 0
	
	if (gml.hasxb[mod]) {
		xb = moptimize_util_xb(M,b,gml.xeqn[mod])
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	if (gml.hastb[mod]) {
		if (hast)	xb = xb :+ uhtred_util_tb(M,b,gml,t)
		else 		xb = xb :+ uhtred_util_tb(M,b,gml)
	}
	if (gml.haszb[mod]) xb = xb :+ uhtred_util_zb(M,b,gml)
	return(xb)
}

`RM' uhtred_util_p_xtzb(`gml' gml, | `RC' t)
{
	`RM' x
	`RS' hast, mod
	
	hast = args()==4
	mod  = gml.model
	xb   = 0

	if (gml.hasxb[mod]) {
		xb = asarray(gml.X,mod) * 
			gml.myb[|uhtred_util_bindices(gml,gml.xeqn[mod])|]'
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	if (gml.hastb[mod]) {
		if (hast)	xb = xb :+ uhtred_util_p_tb(gml,t)
		else 		xb = xb :+ uhtred_util_p_tb(gml)
	}
	if (gml.haszb[mod]) xb = xb :+ uhtred_util_p_zb(gml)
	return(xb)
}

`RM' uhtred_util_xzb(`TR' M, `RR' b, `gml' gml)
{
	`RM' x
	`RS' mod
	
	mod  = gml.model
	xb   = 0
	
	if (gml.hasxb[mod]) {
		xb = moptimize_util_xb(M,b,gml.xeqn[mod])
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	if (gml.haszb[mod]) xb = xb :+ uhtred_util_zb(M,b,gml)
	return(xb)
}

`RM' uhtred_util_p_xzb(`gml' gml)
{
	`RM' x
	`RS' mod
	
	mod  = gml.model
	xb   = 0
	
	if (gml.hasxb[mod]) {
		xb = asarray(gml.X,mod) * 
			gml.myb[|uhtred_util_bindices(gml,gml.xeqn[mod])|]'
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	if (gml.haszb[mod]) xb = xb :+ uhtred_util_zb(gml)
	return(xb)
}

`RM' uhtred_util_zb(`TR' M, `RR' b, `gml' gml)
{
	mod = gml.model
	zb = 0
	for (lev=1;lev<=gml.Nrelevels;lev++) {
		bre 	= b[|moptimize_util_eq_indices(M,gml.zeqn[mod,lev])|]
		reindex = asarray(gml.Zbindex,(mod,lev))
		if (gml.adapt[1]) {
			Zbeta 	= asarray(gml.Z,(mod,lev)) :* bre	
			indexr  = uhtred_get_adpanelindex(gml,lev)
			indexc  = .
			if (lev!=gml.Nrelevels) indexc = gml.qind[,lev+1]

			for (r=1;r<=gml.Nres[lev];r++)  {
				zb = zb :+ Zbeta[,r] :* 
					asarray(gml.aghip2,
						(lev,reindex[r]))[indexr,indexc]
			}			
		}
		else {
			if (lev==gml.Nrelevels) {
				zb = zb :+ (asarray(gml.Z,(mod,lev)) :* bre) * 
					asarray(gml.b,lev)[reindex,]
			}
			else {
				zb = zb :+ (asarray(gml.Z,(mod,lev)) :* bre) * 
					asarray(gml.b,lev)[reindex,
						gml.qind[,lev+1]]
			}
		}
	}
	return(zb)
}

`RM' uhtred_util_p_zb(`gml' gml)
{
	mod = gml.model
	zb = 0
	for (lev=1;lev<=gml.Nrelevels;lev++) {
		bre 	= gml.myb[|uhtred_util_bindices(gml,gml.zeqn[mod,lev])|]
		reindex = asarray(gml.Zbindex,(mod,lev))
		if (gml.adapt[1]) {
			Zbeta 	= asarray(gml.Z,(mod,lev)) :* bre	
			indexr  = uhtred_get_adpanelindex(gml,lev)
			indexc  = .
			if (lev!=gml.Nrelevels) indexc = gml.qind[,lev+1]

			for (r=1;r<=gml.Nres[lev];r++)  {
				zb = zb :+ Zbeta[,r] :* 
					asarray(gml.aghip2,
						(lev,reindex[r]))[indexr,indexc]
			}			
		}
		else {
			if (lev==gml.Nrelevels) {
				zb = zb :+ (asarray(gml.Z,(mod,lev)) :* bre) * 
					asarray(gml.b,lev)[reindex,]
			}
			else {
				zb = zb :+ (asarray(gml.Z,(mod,lev)) :* bre) * 
					asarray(gml.b,lev)[reindex,
						gml.qind[,lev+1]]
			}
		}
	}
	return(zb)
}

`RM' uhtred_util_tb(`TR' M, `RR' b, `gml' gml, | `RC' t)
{
	`RM' xb
	`RS' hast
	
	xb   = 0
	hast = args()==4
	teqn = gml.teqn[gml.model]
	
	if (hast) {
		bt = b[|moptimize_util_eq_indices(M,teqn)|]'
		xb = _uhtred_util_tb_update(bt,gml,t)
	}
	else {
		xb = moptimize_util_xb(M,b,teqn)
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	
	return(xb)
}

`RM' uhtred_util_p_tb(`gml' gml, | `RC' t)
{
	`RM' xb
	`RS' hast
	
	xb   = 0
	hast = args()==4
	teqn = gml.teqn[gml.model]
	
	if (hast) {
		bt = gml.myb[|uhtred_util_bindices(gml,teqn)|]'
		xb = _uhtred_util_tb_update(bt,gml,t)
	}
	else {
		xb = asarray(gml.XT,gml.model) * 
			gml.myb[|uhtred_util_bindices(gml,teqn)|]'
		if (gml.Nmodels>1) {
			xb = xb[uhtred_util_index(gml),]
		}
	}
	
	return(xb)
}

`RM' _uhtred_util_tb_update(`RC' bt, `gml' gml, `RC' t)
{
	return(_uhtred_util_t_update(gml,t) * bt)
}

`RM' uhtred_util_t(`gml' gml, | `RC' t)
{
	if (args()==2) 	return(_uhtred_util_t_update(gml,t))
	else 		return(asarray(gml.XT,gml.model))
}

`RM' _uhtred_util_t_update(`gml' gml, `RC' t)
{
	/*
	- loop over each component in the t equation
	  -> loop over each element
	     -> if time then replace
	- read in new T design matrix with st_data
	*/
	
	mod 	= gml.model
	Ncmps 	= gml.Ncmps[mod]		        //# of components
	Nels 	= asarray(gml.Nels,mod)			//# els per component
	Nobs	= uhtred_util_nobs(gml)
	newtvarlist = gml.tvarlist[mod]

	//now rebuild
	for (i=1;i<=Ncmps;i++) {
		istdep 	= sum(asarray(gml.eltvar,(mod,i))[,1])
		if (istdep) {
			eltypes = asarray(gml.elindex,(mod,i))
			
			//need comp varlist entry for st_data/view
			for (j=1;j<=Nels[i];j++) {
				if (eltypes[j]==8) {
					elvars = uhtred_xz_rcs(gml,i,j,0,t)
					//now post to Stata
					Nvars = cols(elvars)
					names = J(1,0,"")
					stub  = "_temp_"+strofreal(mod)+
						"_"+strofreal(i)+
						"_"+strofreal(j)
					for (k=1;k<=Nvars;k++) {
						names = names,(stub+"_"+strofreal(k))
					}

					stata("cap drop "+invtokens(names))
					id = st_addvar("double",names)
					st_store(.,id,gml.touse,elvars)
					oldstub = "_rcs" + strofreal(mod) + 
							"_"+strofreal(i) + 
							"_" + strofreal(j) + 
							"_*"
					newtvarlist = subinstr(newtvarlist,oldstub,stub+"*",1)
				}
			}

		}
	}	

	T = st_data(.,newtvarlist,gml.modeltouses[gml.modtoind])
	stata("cap drop "+invtokens(names))
	return(T)
}

`RM' uhtred_util_tb_qs(`TR' M, `RR' b, `gml' gml,`RS' Nobs2)
{
	`RM' xb
	
	chq2 = J(Nobs2,gml.chip,0)
	teqn = gml.teqn[gml.model]
	bt   = b[|moptimize_util_eq_indices(M,teqn)|]'
	
	for (q=1;q<=gml.chip;q++) {
		chq2[,q] = asarray(gml.chq,(gml.model,q)) * bt
	}	
	
	return(chq2)
}

`RM' uhtred_util_p_tb_qs(`gml' gml,`RS' Nobs2)
{
	`RM' xb
	
	chq2 = J(Nobs2,gml.chip,0)
	teqn = gml.teqn[gml.model]
	bt   = gml.myb[|uhtred_util_bindices(gml,teqn)|]'
	
	for (q=1;q<=gml.chip;q++) {
		chq2[,q] = asarray(gml.chq,(gml.model,q)) * bt
	}	
	
	return(chq2)
}
/*
	uhtred_util_depvar()
	-> extract dependent variable(s) for current model
	-> for survival outcomes, first column in time, second is event indicator
*/

`RM' uhtred_util_depvar(`gml' gml)
{
	if (!gml.survind) return(asarray(gml.y,gml.model))
	else return(asarray(gml.y,gml.model)[uhtred_util_index(gml),])
}

/*
	uhtred_util_ap()
	-> extract ancillary parameters for current model
*/

`RC' uhtred_util_ap(`gml' gml, `RR' b, `RS' index)
{
// 	return(asarray(gml.apxb,(1,index)))
}

/*
	uhtred_util_dap()
	-> extract distributional ancillary parameters for current model
*/

`RC' uhtred_util_dap(`gml' gml, `RS' index)
{
	return(asarray(gml.distancb,(gml.model,index)))
}

/*
	uhtred_util_expval()
	-> extract expected value for current model
	-> possibly dependent on new timevar
*/


/*
	uhtred_util_xzb_deriv()
	-> extract d/dt of complex linear predictor for the current model
	-> possibly dependent on new timevar
*/

`RC' uhtred_util_timevar(`gml' gml)
{
	if (!gml.survind) return(asarray(gml.timevar,gml.model))
	else return(asarray(gml.timevar,gml.model)[uhtred_util_index(gml)])
}


/*
	updates particular element values -> called from overoutcome() in predict
*/

`RR' uhtred_util_istimevar(`gml' gml, `RS' c, `RS' el)
{
	return(asarray(gml.eltvar,(gml.model,c))[el,])
}

`RM' uhtred_util_bhazard(`gml' gml)
{
	if 	(gml.hasbh[gml.model,1]) return(gml.bhazard[uhtred_util_index(gml),])
	else if (gml.hasbh[gml.model,2]) return(asarray(gml.bhazard,(gml.model,gml.survind)))
}

`RM' uhtred_util_bHazard(`gml' gml)
{
	return(asarray(gml.bhazard,(1,gml.survind+2)))
}

`RR' uhtred_util_score_indices(`TR' M, `RS' eqn)
{
	sindex1 = moptimize_util_eq_indices(M,eqn)
	if (sindex1[1,2]==sindex1[2,2]) sindex1	= sindex1[2,2]
	else sindex1 = (sindex1[1,2]..sindex1[2,2])
	return(sindex1)
}

`RM' uhtred_util_bindices(`gml' gml, `RS' eqn, | `RS' eqn2)
{
	if (args()==2) 	return(asarray(gml.bindices,(gml.model,eqn)))
	else 		return(asarray(gml.Hindices,(gml.model,eqn,eqn2)))
}

end
