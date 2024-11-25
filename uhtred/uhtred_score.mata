*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct uhtred_struct scalar) scalar
local gml 	struct uhtred_struct scalar

version 17
mata:

/*
scores for multilevel model (without any ?EV[], ?XB[]

- one function for complex predictor
- one function for each distributional ancillary parameter
*/

void uhtred_score_panels(`TR' M, `RR' b, `RM' G, `RC' lnfi, `gml' gml)
{

	bind = 1
	for (model=1; model<=gml.Nmodels; model++) {
		
		gml.model = gml.modtoind = model
		hasxb 	= gml.hasxb[model]
		hastb 	= gml.hastb[model]
		haszb 	= gml.haszb[model]
		
		if (hasxb) {
			xeqn 	= gml.xeqn[model]
			sindex1 = uhtred_util_score_indices(M,xeqn)
			Nxbs 	= cols(sindex1)
			for (bi=1;bi<=Nxbs;bi++) {
				gml.survind 	= 0
				gindex 		= sindex1[bi]
				G[,gindex] 	= uhtred_panels(1,gml,
							&_score_rp_xb(),bi)
				bind++
			}
		}

		if (hastb) {
			teqn 	= gml.teqn[model]
			sindex1 = uhtred_util_score_indices(M,teqn)
			Nxbs 	= cols(sindex1)
			for (bi=1;bi<=Nxbs;bi++) {
				gml.survind 	= 0
				gindex 		= sindex1[bi]
				G[,gindex] 	= uhtred_panels(1,gml,
							&_score_rp_tb(),bi)
				bind++
			}
			
		}
		
		if (haszb) {
			//levels loop
			//!! needed for association structures later
			bind++
		}
		
		
	}
	
	//vcv
	Nres 	= gml.Nres
	covariances = gml.covariances
	
	for (i=1 ; i < gml.Nlevels ; i++) {

		if (Nres[i]) {
		
			if (covariances[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					gml.survind 	= 0
					G[,bind] = uhtred_panels(1,gml,
							&_score_rp_vcv(),bi)
					bind++
				}
			}
			
			
			
		}
		
	}
	
	
	//final scores
	G = G :/ exp(lnfi)
	
	
	
// 	}
// 	else {
//	
// 		for (mod=1;mod<=gml.Nmodels;mod++) {
//			
// 			gml.model 	= gml.modtoind = mod
// 			NHbs 		= asarray(gml.NHbs,mod)	
//			
// 			for (eqn=1;eqn<=gml.NHeqns[mod];eqn++) {
// 				if (eqn==1) {
// 					hasconstr = asarray(gml.hasconstraint,mod)
// 					for (el=1;el<=NHbs[eqn];el++) {
// 						if (!hasconstr[el]) G[,bindex] = uhtred_dscore_panels(gml,bindex)
// 						bindex++
// 					}
// 				}
// 				else {
// 					for (el=1;el<=NHbs[eqn];el++) {
// 						G[,bindex] = uhtred_dscore_panels(gml,bindex)
// 						bindex++
// 					}
// 				}
// 			}
// 		}
//		
// 	}
	
// 	//vcv - derivatives found numerically
// 	Nres 	= gml.Nres
// 	covariances = gml.covariances
//	
// 	for (i=1 ; i < gml.Nlevels ; i++) {
//
// 		if (Nres[i]) {
//
// 			if (covariances[1,i]) {
// 				for (j=1;j<=Nres[i];j++) {
// 					G[,eqnind] = uhtred_deriv(gml,bind)
// 					bind++
// 				}
// 			}
// 			else if (covariances[2,i]) {
// 				G[,eqnind] = uhtred_deriv(gml,bind)
// 					bind++
// 			}
// 			else if (covariances[3,i]) {
// 				for (j=1;j<=Nres[i];j++) {
// 					G[,eqnind] = uhtred_deriv(gml,bind)
// 					bind++
// 				}
// 				for (j=1;j<=Nres[i];j++) {
// 					ind = 1
// 					while (ind<j) {
// 						G[,eqnind] = uhtred_deriv(gml,bind)
// 						bind++
// 						ind++
// 					}
// 				}
// 			}
// 			else {
// 				G[,eqnind] = uhtred_deriv(gml,bind)
// 				bind++
// 			}			
// 		}
// 	}
	
}

`RC' uhtred_panels(`RS' index,		/// -level-
		   `gml' gml,		/// -uhtred object-
		   `PS' func,		/// -function to call-
		   | `RS' Xindex)	/// -design matrix element-
{
	`RS' index2
	`RM' res, panelindex
	
	if (args()==3) Xindex = 1
	
	index2 = index+1
	
	resq = J(gml.Nobs[index,1],gml.ndim[index],0)
	
	if (index<gml.Nrelevels) {
		panelindex = asarray(gml.panelindexes,(index,1))
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			resq[,q] = panelsum(uhtred_panels(index2,gml,func,Xindex),panelindex)
		}
	}
	else {
		for (j=1;j<=gml.Nmodels;j++) {
			gml.model = gml.modtoind = j
			resq2 = (*func)(gml,Xindex)
			resq = resq :+ panelsum(resq2,asarray(gml.panelindexes,(index,j)))
		}
	}

	resq = resq :* asarray(gml.Li_ip,gml.qind) 

	if (gml.usegh[index]) {			//GHQ
		return(resq * asarray(gml.baseGHweights,index))
	}
	else {					//MCI
		return(quadrowsum(resq):/gml.ndim[index])
	}
}

void uhtred_logl_panels_d(`RS' myb,`gml' gml, `RS' bind, `RC' lnf)
{
	gml.myb[bind] = myb
	uhtred_fillvcv(gml,gml.myb,bind)
	lnf = uhtred_logl_panels(1,gml)
}

`RM' uhtred_deriv(`gml' gml, `RS' bind)
{
	
	myb = gml.myb
	D   = deriv_init()
	deriv_init_evaluator(D, &uhtred_logl_panels_d())
	deriv_init_evaluatortype(D, "v")
// 	deriv_init_search(D,"off")
	deriv_init_argument(D, 1, gml)
	gml.myb = myb		
	deriv_init_params(D, myb[bind])
	deriv_init_argument(D, 2, bind)
	ddb = deriv(D, 1)
	gml.myb = myb
	return(deriv_result_scores(D))
}


`RC' uhtred_dscore_panels(`gml' gml, `RS' bind)
{
	b = gml.myb[bind]				
	hstep = c("epsdouble")^(1/3)
	if (abs(b)<=1 & abs(b)>0) {
		hstep = abs(b)*hstep
	}
	copymyb = gml.myb
	gml.myb[bind] = copymyb[bind] :+ hstep
	uhtred_xb(gml,gml.myb)
	lnf1 = uhtred_logl_panels(1,gml)
	gml.myb[bind] = copymyb[bind] :- hstep
	uhtred_xb(gml,gml.myb)
	lnf2 = uhtred_logl_panels(1,gml)
	gml.myb = copymyb
	return((lnf1:-lnf2):/(2:*hstep))
}

end
