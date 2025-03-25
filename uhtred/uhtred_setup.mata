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

struct uhtred_ereturn_struct {
	`SR' levelvars		//level vars 
	`RC' Nreparams		//# of VCV parameters at each level
	`SC' reeqns		//eqn names for re parameters
	`SC' reivscale		//inverse scale of re parameters
	`SC' relabel		//labels for re parameters
	`SM' userfunctions			
	`SC' Nvars
	`RS' uhtred
}

struct uhtred_struct {

// [op] - denotes an option which is not always completely filled in

	`pgml' Pgml		//self-pointer
	
	//model definition
	`TR' y			//response variables
	`SC' familys		//family of each model
	`SC' responses		//response var names 
                                // -> includes event and ltruncated if there
	
	`RC' issurv
	`SC' failures		//[op]	//failure variable name of each model 			-> for eret list
	`SC' ltruncated		//[op]	//ltruncated variable name of each model		-> for eret list
	`SC' linterval		//[op]	//linterval variable name of each model			-> for eret list
	`TR' surv_index
	`RM' Nsurv
	`RS' survind
	
	`RC' hasltrunc		//has left truncation/delayed entry -> handle with counting process notation
	`RC' haslint		//has interval censoring
	`RS' hasanyltrunc
	`RS' hasanylint
	
	`SS' touse		//global touse variable name
	`SR' levelvars		//level varables from highest to lowest
	`RM' hasbh		//presence of bhazard or bhfile for each model
	`TR' bhazard		//stored bhazard variable
	`SC' bhvarnames		//bhazard() varnames
	`RC' timevar		//stored timevar()
	`SC' tvarnames		//for eret list
	`SR' modeltouses
	`TR' xbindex
	
	`RS' Nmodels		//# of models
	`RS' Nlevels		//# of levels (includes ob. level)
	`RS' Nrelevels		//# of levels with random effects
	`RS' N			//# of observations (ob. level, main touse)
	`RM' Nobs		//# of observations for each model = level by models
	
	`RC' Ndap		//# of distributional ancillary parameters for each model
	`RC' Nap		//[op]	//# of anciliary parameters in a user model
	`TR' distancb
	`TR' apxb
	
	`PC' Plnl, Pxb		//pointers to logL and linear predictor functions for each model
	`RS' todo		//gf0,1,2
	`PC' Ph			//[op]	//pointer to hazard functions
	`PC' Pch		//[op]	//pointer to cumulative hazard functions
	`PC' Plogh		//[op]	//pointer to log hazard functions
	`PC' Pcdf
	`PM' Puserst
	
	//complex syntax
	`RR' myb
	`RC' hascons		//whether each complex predictor has a constant term
	`RC' Ncmps		//number of components in each complex syntax
	`TR' Nels		//# of elements per component per model
	`TR' elindex		//index of each element within each component, for each model
	`TR' elvarname		//varname for each element (used in matching when at() from predict mediation etc.)
	`TR' elinfo
	`TR' hasconstraint	//each model clp, whether each parameter is constraned 0/1
	`TR' eltvar		//RR for each mod,comp - flag for whether input var matches timevar(), ltruncated()
	
	//integration and random effects
	`RS' gridsearch
	`RM' covariances
	`TR' panelindexes
	`TR' vcvs
	`TR' baseGHnodes
	`TR' baseGHweights
	`RC' ndim
	`RS' adapt
	`RS' iter
	`RS' atol
	`RS' showadapt
	`RS' df
	`TR' b, bdraws
	`TR' Li_ip, Li_ip2
	`TR' c
	`TR' aghip
	`TR' aghip2
	`TR' aghlogl
	`TR' stackednodes
	`TR' adpanelindexes
	`PS' Pupdateip
	`SS' seed
	`RC' usegh
	`RC' ip
	`RC' Nres
	
	`RC' NI, y2
	`RS' chip
	`RM' chqp, chqw
	`TR' chq
	
	`TR' latlevs		//unique REs at each level for eret list
	`RS' adaptit		//number of iterations to allow adapting nodes
	
	`RS' tofix
	`RR' fxls
	
	//delayed entry
	`RS' ltflag
	
	//gets updated
	`RS' lnfi1
	`RR' qind
	`RS' model, modtoind
	
	//for predictions
	`RS' predict				
	`RS' fixedonly		//predictions based only on fixed effects - so skip res in utils
	`TR' blups
	
	//ereturn list
	`Egml' E
	`SS' allvars
	
	//morgana
	`RS' morgana
	
	//design matrices
	`TR' X			//design matrix for fixed non-timevar covariates
	`TR' XT			//design matrix for fixed timevar covariates
	`TR' dXT		//design matrix for d/dt of fixed timevar covariates
	`TR' XT0		//design matrix for fixed timevar covariates at 
				//  ltruncated() times
	`TR' XTL		//design matrix for fixed timevar covariates at 
				//  linterval() times
	`TR' Z			//Z design matrix (at each level) for simple models
	`TR' Zbindex		//indices to extract quad. points for Z * u[indices]
	`RC' hasxb		//has xb linear predictor
	`RC' hastb		//has tb linear predictor
	`RC' haszb		//has zb linear predictor
	`RC' xeqn		//b equation index for each x equation
	`RC' teqn		//b equation index for each t equation
	`RM' zeqn		//b equation index for each z equation
	`SC' tvarlist		//time equation parsed varlist
	`TR' bindices		//model,equation coefficient indices
	`TR' Hindices		//model,equation,equation coefficient indices
	
	`TR' Sb
	
	//===================================================================//
	//development
	
        `SS' indicatorvar       //indicator() varname from stexcess
        `RC' indicator          //indicator view for stexcess
	
	`RR' svs
	
}

void uhtred_setup(`SS' GML,`SS' touse)
{
	//declarations
	`pgml' 	pGML
	`gml' 	gml
	
	//initialise
	gml.Pgml = pGML = crexternal(GML)
	
	//core setup
	uhtred_setup_core(gml,touse)
	
	uhtred_setup_family_general(gml)
	uhtred_setup_check_clp(gml)
	uhtred_setup_levels(gml)
	uhtred_setup_error_checks(gml)
	uhtred_setup_touses(gml)
	uhtred_setup_levelvars(gml)
	
	uhtred_get_ys(gml)
       
//         gml.indicatorvar = st_local("indicator")
//         if (gml.indicatorvar!="") {
//                 st_view(gml.indicator=.,.,gml.indicatorvar,gml.touse)
//         }
	
	//clp
	uhtred_get_noconstants(gml)
        uhtred_get_timevars(gml)
	uhtred_get_latents(gml)
	uhtred_setup_evaltype(gml)
        uhtred_build_clp(gml,coef=.)

	//survival extras
	uhtred_setup_survival(gml)
	uhtred_build_clpq(gml)

	//latents
	uhtred_setup_vcv(gml,coef)
	

	//pointers
	uhtred_get_logl_p(gml)

 	//ml stuff
	uhtred_setup_mleqns(gml)
	uhtred_setup_wrappers(gml)
	uhtred_starting_values(gml)


// if (gml.predict) {
// 	asarray(gml.vcvs,1)
// 	asarray(gml.vcvs,2)
// 	gml.myb = st_matrix(st_local("from"))
// 	gml.myb
// // 	uhtred_p_xb(gml,gml.myb)
// // 	uhtred_p_init_ip(gml)
// 	asarray(gml.vcvs,1)
// 	asarray(gml.vcvs,2)
// // 	uhtred_p_init_ip(gml)
// }


 	//Done
	swap(*pGML,gml)
}

void uhtred_setup_core(`gml' gml, `SS' touse)
{
	gml.touse 	= touse
	gml.Nmodels	= strtoreal(st_local("EQ_n"))
	gml.fixedonly	= 0
	gml.iter 	= 0
	gml.survind	= 0
	gml.predict 	= st_local("predict")!=""
	if (gml.predict & st_local("npredict")!="") {
		gml.N = strtoreal(st_local("npredict"))	
	}
}	

void uhtred_setup_wrappers(`gml' gml)
{	
	if 	(st_local("bors")!="")	        gml.E.uhtred	= 4  //stuhtred
	else if (st_local("arthur")!="") 	gml.E.uhtred	= 3  //neuralnet
	else if (st_local("excalibur")!="")     gml.E.uhtred	= 2  //stmixed
	else if (st_local("sagramore")!="")     gml.E.uhtred	= 6  //jm
        else if (st_local("mordred")!="")       gml.E.uhtred    = 7  //stexcess
	else {
		if (gml.Nrelevels)		gml.E.uhtred	= 1  //uhtred
		else 				gml.E.uhtred	= 5  //uhtred
	}
	gml.morgana = st_global("c(prefix)")=="morgana"
}

void uhtred_setup_family_general(`gml' gml)
{
	gml.familys 	= tokens(st_local("familylist"))'
	gml.failures	= J(gml.Nmodels,1,"")
	gml.issurv	= J(gml.Nmodels,1,.)
	gml.ltruncated 	= J(gml.Nmodels,1,"")
	gml.hasltrunc	= J(gml.Nmodels,1,.)
	gml.linterval 	= J(gml.Nmodels,1,"")
	gml.haslint	= J(gml.Nmodels,1,.)
	gml.bhvarnames 	= J(gml.Nmodels,1,"")
	gml.hasbh	= J(gml.Nmodels,2,.)
	for (i=1;i<=gml.Nmodels;i++) {
		stri = strofreal(i)
		gml.failures[i] = st_local("failure"+stri)
		gml.issurv[i]	= gml.failures[i]!=""
		gml.ltruncated[i] = st_local("ltruncated"+stri)
		gml.hasltrunc[i] = gml.ltruncated[i]!=""
		gml.linterval[i] = st_local("linterval"+stri)
		gml.haslint[i] 	= gml.linterval[i]!=""
		gml.bhvarnames[i] = st_local("bhaz"+stri)
		gml.hasbh[i,1] 	= gml.bhvarnames[i]!=""	
		gml.hasbh[i,2] 	= st_local("bhfile"+stri)!=""
	}
	gml.hasanyltrunc = sum(gml.hasltrunc)
	gml.hasanylint = sum(gml.haslint)
	gml.ltflag	 	= 0	//gets updated when needed
	gml.distancb 		= asarray_create("real",2)
	uhtred_setup_daps(gml)
	uhtred_setup_ap(gml)
}

void uhtred_setup_ap(`gml' gml)
{
	gml.apxb = asarray_create("real",2)
	gml.Nap  = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		gml.Nap[i] = strtoreal(st_local("nap"+strofreal(i)))	
		apeqns	= J(1,0,"")
		for (a=1;a<=gml.Nap[i];a++) {
			apeqns = apeqns," " + 
				"(ap"+strofreal(i)+"_"+strofreal(a)+":)"
		}
		st_local("ap_eqns"+strofreal(i),invtokens(apeqns))
	}
}

`SR' uhtred_get_cmps(`RS' i)
{
	indepvars = J(1,0,"")
	st_local("rest",st_local("indepvars"+strofreal(i)))
	stata("gettoken lhs rest : rest, bind")						
	indepvars = indepvars,st_local("lhs")
	while (st_local("rest")!="") {
		stata("gettoken lhs rest : rest, bind")
		indepvars = indepvars,st_local("lhs")
	}
	return(indepvars)
}

void uhtred_get_ys(`gml' gml)
{
	gml.responses = J(0,1,"")
	gml.y = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
		gml.responses = gml.responses\st_local("response"+strofreal(i))
		if (gml.responses[i]!="") {
			asarray(gml.y,i,st_data(.,gml.responses[i],gml.modeltouses[i]))
		}
	}
}

void uhtred_get_timevars(`gml' gml)
{
	gml.tvarnames = J(gml.Nmodels,1,"")
	gml.timevar = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
		gml.tvarnames[i] = st_local("timevar"+strofreal(i))
		if (gml.issurv[i]) {
			gml.tvarnames[i] = tokens(gml.responses[i])[1]
		}
		if (gml.tvarnames[i]!="") {
			asarray(gml.timevar,i,
				st_data(.,gml.tvarnames[i],gml.modeltouses[i]))
		}
	}
}

void uhtred_init_vcvs(`gml' gml)
{
	gml.vcvs = asarray_create("real",1)
	for (i=1;i<=gml.Nrelevels;i++) {
		asarray(gml.vcvs,i,I(gml.Nres[i]))
	}
}

void uhtred_parse_covstructures(`gml' gml)
{
	`SR' covs
	
	covs = tokens(st_local("covariance2"))
	if (covs==J(1,0,"")) {
		//default independent
		gml.covariances = J(1,gml.Nlevels,(1\0\0))
	}
	else {
		//check not more than Nlevels
		ncovs = cols(covs)
		if (ncovs>1 & ncovs!=gml.Nlevels) {
			errprintf("Error in covariance()\n")		
			exit(198)
		}
		if (ncovs==1) {
			gml.covariances = J(1,gml.Nlevels,
				(covs:=="diagonal"\covs:=="exchangeable"\
				covs:=="unstructured"))
		}
		else {
			gml.covariances = covs:=="diagonal"\
			 covs:=="exchangeable"\covs:=="unstructured"
		}
	}
	stata("tempname vcvmat")
	st_matrix(st_local("vcvmat"),gml.covariances)
}

void uhtred_parse_vcv_eqns(`gml' gml, `RS' coef)
{
	vcveqns = J(1,0,"")
	gml.E.Nreparams = J(gml.Nlevels,1,0)
	gml.E.reeqns 	= J(gml.Nlevels,1,"")
	gml.E.reivscale = J(gml.Nlevels,1,"")
	gml.E.relabel 	= J(gml.Nlevels,1,"")

	for (lev=1;lev<=gml.Nlevels;lev++) {
		if (!gml.Nres[lev]) continue
	
		lats = asarray(gml.latlevs,lev)[,1]
		strlev = strofreal(lev)
	
		//ind or unstr
		if (gml.covariances[1,lev] | gml.covariances[3,lev]) {
			for (r=1;r<=gml.Nres[lev];r++) {
				todo = 1
				vars = "" //uhtred_get_vcv_vars(gml, lev)
				vcveqns = vcveqns,"(lns"+strlev+"_"+strofreal(r)+": "+vars+")"
				gml.E.Nreparams[lev] 	= gml.E.Nreparams[lev] + 1
				gml.E.reeqns[lev] 	= gml.E.reeqns[lev] + "lns"+strlev+"_"+strofreal(r)+" "
				gml.E.reivscale[lev]    = gml.E.reivscale[lev] + "exp "
				gml.E.relabel[lev] 	= gml.E.relabel[lev] + "sd("+lats[r]+") " 
			}
		}
		else {
			if (gml.Nres[lev]) {
				vars = "" //uhtred_get_vcv_vars(gml, lev)
				vcveqns = vcveqns,"(lns"+strlev+"_1: "+vars+")"
				gml.E.reeqns[lev] = gml.E.reeqns[lev] + "lns"+strlev+"_1 "
				gml.E.reivscale[lev] = gml.E.reivscale[lev] + "exp "
				gml.E.relabel[lev] 	= gml.E.relabel[lev] + "sd " 
				gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
			}
		}
		
		//exch and nres>1
		if (gml.covariances[2,lev] & gml.Nres[lev]>1) {
			vcveqns = vcveqns,"(art"+strlev+"_1_1: )"
			gml.E.reeqns[lev] = gml.E.reeqns[lev] + "art"+strlev+"_1_1 "
			gml.E.reivscale[lev] = gml.E.reivscale[lev] + "tanh "
			gml.E.relabel[lev] 	= gml.E.relabel[lev] + "corr " 
			gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
		}
		else if (gml.covariances[3,lev]) {
			ind = 1
			while (ind<gml.Nres[lev]) {
				for (r=ind+1;r<=gml.Nres[lev];r++) {
					vcveqns = vcveqns,"(art"+strlev+"_"+strofreal(ind)+"_"+strofreal(r)+": )"
					gml.E.reeqns[lev] = gml.E.reeqns[lev] + "art"+strlev+"_"+strofreal(ind)+"_"+strofreal(r)+" "
					gml.E.reivscale[lev] = gml.E.reivscale[lev] + "tanh "
					gml.E.relabel[lev] 	= gml.E.relabel[lev] + "corr("+lats[r]+","+lats[ind]+") " 
					gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
				}
				ind++
			}			
		}		
	}
	//post to Stata
	st_local("vcveqns",invtokens(vcveqns))
	
}

void uhtred_get_logl_p(`gml' gml)
{
	`RS' nm
	`PC' pmat
	
	familys 	= gml.familys
	nm 		= gml.Nmodels
	pmat 		= J(nm,1,NULL)
	pmathaz 	= pmat
	pmatloghaz 	= pmat
	pmatchaz 	= pmat
	pmatcdf 	= pmat
	pmatuserst	= J(nm,3,&nm)	//h,logh,ch
	
	gml.E.userfunctions = J(nm,4,"")
	
	for (i=1; i<=nm; i++) {
		if (familys[i]=="weibull") {
			pmat[i] = &_logl_weibull()
//			pmatloghaz[i] 	= &uhtred_weibull_logh()
// 			pmatchaz[i] 	= &uhtred_weibull_ch()
// 			pmatcdf[i] 	= &uhtred_weibull_cdf()
		}
		else if (familys[i]=="loghazard") {
			pmat[i] = &_logl_loghazard()
		}
		else if (familys[i]=="rp") {
			pmat[i] = &_logl_rp()
// 			pmatloghaz[i] 	= &uhtred_rp_logh()
// 			pmatchaz[i] 	= &uhtred_rp_ch()
// 			pmatcdf[i] 	= &uhtred_rp_cdf()
		}
		
	}
	
	gml.Plnl 	 = pmat
	gml.Ph 		= pmathaz
	gml.Plogh 	= pmatloghaz
	gml.Pch 	= pmatchaz
	gml.Pcdf	= pmatcdf
	gml.Puserst     = pmatuserst
}

void uhtred_setup_daps(`gml' gml)
{
	familys	 = gml.familys
	gml.Ndap = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		stri 	 = strofreal(i)
		if (familys[i]=="weibull") {
			gml.Ndap[i] = 1
			st_local("dap_eqns"+stri,"(dap"+stri+"_1:)")
		}
	}
}

void uhtred_get_noconstants(`gml' gml) 
{
	gml.hascons = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys!="cox") {
			gml.hascons[i] = 
				st_local("constant"+strofreal(i))!="noconstant"
		}
	}
}

void uhtred_setup_evaltype(`gml' gml)
{
	if (st_local("evaltype")=="") {
		gf = "gf0"
		if (!gml.Nrelevels) {
			check = J(gml.Nmodels,1,0)
			for (i=1;i<=gml.Nmodels;i++) {
				f = gml.familys[i]
				if (f=="weibull" | f=="rp" ) {
					check[i] = 2
				}
				if (gml.haslint[i]) {
					check[i] = 1
				}
				if (f=="weibull" & gml.hasbh[i,1]) {
					check[i] = 0
				}
				if (gml.hasbh[i,2]) {
					check[i] = 0
				}
			}
			gf = "gf"+strofreal(min(check))
		}
// 		else {
// 			if (gml.Nrelevels==1 & sum(gml.Nres)==1) gf = "gf2"
// 		}
		st_local("evaltype",gf)
	}
	gml.todo = strtoreal(substr(st_local("evaltype"),3,1))
	if (gml.Nrelevels & gml.todo==2) gml.Sb = asarray_create("real",1)
}

void uhtred_setup_mleqns(`gml' gml)
{
	space = " "
	for (i=1;i<=gml.Nmodels;i++) {
		mlspeci	= st_local("xb"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("xb"+strofreal(i)+"_time")
		mlspeci	= mlspeci+space+st_local("dap_eqns"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("ap_eqns"+strofreal(i))
		st_local("mlspec",st_local("mlspec")+space+mlspeci)
		for (j=1; j<=gml.Nrelevels; j++) {
			st_local("mlspec",st_local("mlspec") + space + 
				st_local("Z"+strofreal(i)+"_"+strofreal(j)))
		}
	}
	st_local("mlspec", st_local("mlspec")+space+st_local("vcveqns"))
}

void uhtred_cleanup(`SS' GML)
{
	rmexternal(GML)
}

end
