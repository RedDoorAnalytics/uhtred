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

void uhtred_predict(`SS' object, `SS' newvar, `SS' touse, `SS' stat, `SS' predtype)
{
	`gml' gml
	
	swap(gml,*findexternal(object))
	gml.Pgml = &gml

	uhtred_predict_setup(gml,stat,touse)
	stand	= st_local("standardise")!=""

	if 	(stat=="reffects" | stat=="reses") {
		pred = uhtred_predict_blups(gml,stat=="reses")
	}
	else {
		pf 	= uhtred_p_getpf(gml,stat)
		pred 	= uhtred_predict_core(gml,pf,predtype,stand)
	}
	id = st_addvar("double",tokens(newvar))
	if (stand & !gml.issurv[gml.model])	{
		st_store(1,id,pred)
	}
	else	st_store(.,id,touse,pred)
}

void uhtred_predict_setup(`gml' gml, `SS' stat, `SS' touse)
{

	gml.myb = st_matrix(st_local("best"))	//get parameter estimates
	uhtred_p_xb(gml,gml.myb)		//and fill up (also updates NI points)
	
	k = strtoreal(st_local("outcome"))
	gml.model = gml.modtoind = k		//keep here uhtred_xb() changes it

	if (stat!="tprob" & stat!="tlos") {
		uhtred_predict_error_check(gml,stat)
	}
	
	if (st_global("e(cmd2)")=="stexcess") {
		gml.Puserst[1,3] = &uhtred_p_stexcess_ch()
		gml.Puserst[1,1] = &uhtred_p_stexcess_h()
	}
	
	//update N
	if (st_local("npredict")=="") {
		gml.N = gml.Nobs[gml.Nlevels,k]
	}

	//fill in timevar if specified
	if (st_local("timevar")!="") {
		asarray(gml.timevar,gml.model,st_data(.,st_local("timevar"),touse))
	}	
}

void uhtred_predict_error_check(`gml' gml, `SS' stat)
{
	f = gml.familys[gml.model]
	
	if (stat=="mu") {
		if (f=="rp" | f=="weibull" | f=="loghazard") {
			uhtred_error("mu not supported with family("+f+")")
		}
	}
	
}

`PS' uhtred_p_getpf(`gml' gml, `SS' stat)
{

	if 	(stat=="eta")		return(&uhtred_p_eta())
	else if (stat=="cif") 		return(&uhtred_p_cif())
	else if (stat=="totalsurvival") return(&uhtred_p_totalsurvival())
	else if (stat=="rmst")		return(&uhtred_p_rmst())
	else if (stat=="timelost")	return(&uhtred_p_timelost())
	else if (stat=="totaltimelost")	return(&uhtred_p_totaltimelost())
	else if (stat=="tprob")		return(&uhtred_p_transprob())
	else if (stat=="tlos")		return(&uhtred_p_los())
	else if (stat=="median")	return(&uhtred_p_med())
	else {
		f = gml.familys[gml.model]
		/*
		if (f=="weibull") {
			if 	(stat=="hazard")	return(&uhtred_p_weibull_h())
			else if (stat=="survival")	return(&uhtred_p_weibull_s())
			else if (stat=="chazard")	return(&uhtred_p_weibull_ch())
			else if (stat=="logchazard")	return(&uhtred_p_weibull_logch())
			else if (stat=="mu") 		return(&uhtred_p_weibull_mu())
                        else if (stat=="density") 	return(&uhtred_p_weibull_dens())
		}
		else if (f=="loghazard") {
			if 	(stat=="hazard")	return(&uhtred_p_loghazard_h())
			else if (stat=="survival")	return(&uhtred_p_loghazard_s())
			else if (stat=="chazard")	return(&uhtred_p_loghazard_ch())
			else if (stat=="logchazard")	return(&uhtred_p_loghazard_logch())
                        else if (stat=="density")	return(&uhtred_p_loghazard_dens())
		}
		else*/
		if (f=="rp") {
			if 	(stat=="hazard")	return(&uhtred_p_rp_h())
			else if (stat=="survival")	return(&uhtred_p_rp_s())
			else if (stat=="chazard")	return(&uhtred_p_rp_ch())
			else if (stat=="logchazard")	return(&uhtred_p_rp_logch())
                        else if (stat=="density")	return(&uhtred_p_rp_dens())
		}
	}
}

`RC' uhtred_predict_core(`gml' gml, `PS' pf, `SS' predtype, `RS' stand)
{
	gml.fixedonly = 1
	if (predtype=="fitted") {
		if (st_local("panel")!="") {
			uhtred_predict_getspecblups(gml)
			gml.fixedonly = 3
		}
		else {
			(void) uhtred_predict_blups(gml,0)	//fills up all blups 
			gml.fixedonly = 2
		}
		//reset
		gml.model = gml.modtoind = strtoreal(st_local("outcome"))
	}

	gml.survind = 0

	if (gml.issurv[gml.model] & stand) {
		t 	= asarray(gml.timevar,gml.model)
		Nt 	= rows(t)
		Nobs	= gml.Nobs[gml.Nlevels,gml.model]
		pred 	= J(Nt,1,.)
		for (i=1;i<=Nt;i++) {
			asarray(gml.timevar,gml.model,J(Nobs,1,t[i]))
			pred[i] = mean((*pf)(gml))
		}
		
	}
	else {
		if (predtype=="marginal") {
			pred = uhtred_predict_marginal(1,gml,pf)
		}
		else    pred = (*pf)(gml)
		if (stand) pred = mean(pred)
	}
		
	return(pred)
}

//integrate over random effects
`RC' uhtred_predict_marginal(	`RS' index,	///	-level-
                                `gml' gml,	///
                                `PS' pf)	//	
{
	gml.fixedonly = 0
        index2 = index + 1
	if (index<gml.Nrelevels) {
		result = J(gml.Nobs[gml.Nlevels,1],gml.ndim[index],0)
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			result[,q] = uhtred_predict_marginal(index2,gml,pf)
		}
	}
	else result = (*pf)(gml)

        //GHQ or MCI
	if (gml.usegh[index]) {
                return(result * asarray(gml.baseGHweights,index))
        }
	else 	return(quadrowsum(result):/gml.ndim[index])
}

`RM' uhtred_p_eta(`gml' gml)
{
	if (gml.hastb[gml.model]) {
		t = uhtred_util_timevar(gml)
		return(uhtred_util_p_xtzb(gml,t))
	}
	else return(uhtred_util_p_xtzb(gml))
}

`RM' uhtred_p_cif(`gml' gml, | `RC' t)
{
	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	
	refmod 	= gml.model

	//get pointers to cause-specific h and ch functions -> gml.model indexes this
	hf 	= uhtred_p_getpf(gml, "hazard")
	chfs 	= J(1,0,NULL)
		
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			chfs 	= chfs,uhtred_p_getpf(gml, "chazard")
			modind 	= modind,gml.model
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			chfs 	= chfs,uhtred_p_getpf(gml, "chazard")
		}
	}
	
	if (Nsurvmodels==1) {
		
		gml.model = modind
		pf = uhtred_p_getpf(gml, "survival")
		return(1:-(*pf)(gml,t))
		
	}
	else {
	
		gml.model 	= refmod
		Ngq 		= gml.chip
		gq 		= uhtred_gq(Ngq,"legendre")
		result 		= J(gml.N,1,0)
		qp		= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2

		for (q=1; q<=Ngq; q++) {						//cif integral
			gml.model 	= refmod
			haz		= (*hf)(gml,qp[,q])
			ochres 		= J(gml.N,1,0)

			for (k=1; k<=Nsurvmodels; k++) {			//overall cumulative hazard integral
				gml.model = modind[k]
				ochres = ochres :+ (*chfs[k])(gml,qp[,q])
			}
			result = result :+ haz :* exp(-ochres) :* gq[q,2] :* t :/ 2
		}

		gml.model 	= refmod
		return(result)	
	
	}	
	
}

`RM' uhtred_p_totalsurvival(`gml' gml, | `RC' t)
{
	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	
	refmod 	= gml.model

	//get pointers to cause-specific h and ch functions -> gml.model indexes this
	chfs 	= J(1,0,NULL)
		
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			chfs 	= chfs,uhtred_p_getpf(gml, "chazard")
			modind 	= modind,gml.model
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			f = gml.familys[modm]
			chfs 	= chfs,uhtred_p_getpf(gml, "chazard")
		}
	}
	
	if (Nsurvmodels==1) {
		gml.model = modind
		pf = uhtred_p_getpf(gml, "survival")
		return(1:-(*pf)(gml,t))
	}
	else {
		gml.model 	= refmod
		result 		= J(gml.N,1,0)
		ochres 		= J(gml.N,1,0)
		for (k=1; k<=Nsurvmodels; k++) {			//overall cumulative hazard integral
			gml.model = modind[k]
			result = result :+ (*chfs[k])(gml,t)
		}
		gml.model = refmod
		return(exp(-result))	
	}	
	
}

`RM' uhtred_p_timelost(`gml' gml, | `RC' t)
{
	if (gml.familys[gml.model]=="cox") {
		return(uhtred_p_cox_timelost(gml))
	}
	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)

	Ngq 	= gml.chip
	gq 	= uhtred_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	
	//integrate cif
	result 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		result = result :+ uhtred_p_cif(gml,qp[,q]) :* gq[q,2] :* t :/ 2
	}
	return(result)	
}

`RM' uhtred_p_totaltimelost(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	refmod 	= gml.model
	
	if (st_local("causes")=="") {
		modind 	= J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			modind 	= modind,gml.model
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			f = gml.familys[modm]
			modind 	= modind,modm
		}
		
	}
	
	//total time lost
	result 	= J(gml.N,1,0)
	for (q=1; q<=Nsurvmodels; q++) {		
		gml.model 	= modind[q]
		result 		= result :+ uhtred_p_timelost(gml,t)
	}
	return(result)	
}

`RM' uhtred_p_rmst(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)

	return(t:-uhtred_p_totaltimelost(gml,t))	
}

/*
- calculate blups for all panels and store within gml.blups array
*/

`RM' uhtred_predict_blups(`gml' gml, `RS' getses)
{
	uhtred_p_ebmeans(gml,getses)
	
	gml.survind = 0
	gml.model = gml.modtoind = strtoreal(st_local("outcome"))	//getblups changes it
	
	//gml.blups will contain blups or seblups
	res = asarray(gml.blups,1)[uhtred_get_adpanelindex(gml,1),]
	if (gml.Nrelevels>1) {
		for (i=2;i<=gml.Nrelevels;i++) { 
			res = res,asarray(gml.blups,i)[uhtred_get_adpanelindex(gml,i),]
		}
	}
	
	return(res)
}

`RM' uhtred_predict_getspecblups(`gml' gml)
{
	gml.blups = asarray_create("real",1)
	
	stata("tempvar specb spectouse")
	if (st_local("blupif")!="") {
		stata("qui gen byte "+st_local("spectouse")+ "= "+st_global("e(levelvars)")+"=="+st_local("panel")+" & "+st_local("touse")+" & "+st_local("blupif"))
	}
	else {
		stata("qui gen byte "+st_local("spectouse")+ "= "+st_global("e(levelvars)")+"=="+st_local("panel")+" & "+st_local("touse"))
	}
	stata("qui predict "+st_local("specb")+"* if "+st_local("spectouse")+", reffects")
	specblups = st_data(.,st_local("specb")+"*",st_local("spectouse"))
	asarray(gml.blups,1,specblups[1,])
}

end
