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

void uhtred_hessian_panels(`TR' M, `RR' b, `RM' G, `RM' H, `RC' lnfi, `gml' gml)
{
	Nbs	= cols(b)
	NH	= 0
	for (i=1;i<=Nbs;i++) NH = NH + i
	H2 = J(gml.Nobs[1],NH,.)
	H2ind = 1
	bind 	= 1
	bind2 	= 1
	Sindex1 = Sindex2 = J(1,0,.)
	
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
				bi2 = 1
				while (bi2<=bi) {
					gml.survind 	= 0
					Hi1 		= sindex1[bi]
					Hi2		= sindex1[bi2]
					Sindex1 = Sindex1,Hi1
					Sindex2 = Sindex2,Hi2
					H2[,H2ind++] 	= uhtred_panels2(1,gml,
							&_hessian_rp_xb(),bi,bi2)
					bi2++
					bind2++
				}
				bind++
			}
		}

		if (hastb) {
			teqn 	= gml.teqn[model]
			sindex2 = uhtred_util_score_indices(M,teqn)
			Ntbs 	= cols(sindex2)
			
			if (hasxb) {
				for (bi=1;bi<=Nxbs;bi++) {
					for (bi2=1;bi2<=Ntbs;bi2++) {
						gml.survind 	= 0
						Hi1 		= sindex1[bi]
						Hi2		= sindex2[bi2]
						Sindex1 = Sindex1,Hi2
						Sindex2 = Sindex2,Hi1
						H2[,H2ind++] 	= uhtred_panels2(1,gml,
									&_hessian_rp_xb_tb(),bi,bi2)
					}
				}
			}
			
			for (bi=1;bi<=Ntbs;bi++) {
				bi2 = 1
				while (bi2<=bi) {
					gml.survind 	= 0
					Hi1 		= sindex2[bi]
					Hi2		= sindex2[bi2]
					Sindex1 = Sindex1,Hi1
					Sindex2 = Sindex2,Hi2
					H2[,H2ind++] 	= uhtred_panels2(1,gml,
								&_hessian_rp_tb(),bi,bi2)
					bi2++
					bind2++
				}
				bind++
			}

			
		}
		
		if (haszb) {
			zeqn 	= gml.zeqn[model]
			sindex3 = uhtred_util_score_indices(M,zeqn)
			Nzbs 	= cols(sindex3)
			
			if (hasxb) {
				for (bi=1;bi<=Nxbs;bi++) {
					for (bi2=1;bi2<=Nzbs;bi2++) {
						Sindex1 = Sindex1,sindex3[bi2]
						Sindex2 = Sindex2,sindex1[bi]
						H2ind++
					}
				}
			}
			
			if (hastb) {
				for (bi=1;bi<=Ntbs;bi++) {
					for (bi2=1;bi2<=Nzbs;bi2++) {
						Sindex1 = Sindex1,sindex3[bi2]
						Sindex2 = Sindex2,sindex2[bi]
						H2ind++
					}
				}
			}
			
			
			//levels loop
			//!! needed for association structures later
			
			for (bi2=1;bi2<=Nzbs;bi2++) {
				Sindex1 = Sindex1,sindex3[bi2]
				Sindex2 = Sindex2,sindex3[bi2]
				H2ind++
				bind++
				bind2++
			}
			
		}
		
	}
	
	//vcv
	Nres 	= gml.Nres
	covariances = gml.covariances

	for (i=1 ; i < gml.Nlevels ; i++) {

		if (Nres[i]) {
		
			if (covariances[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					
					if (hasxb) {
						for (bi=1;bi<=Nxbs;bi++) {
							Sindex1 = Sindex1,bind
							Sindex2 = Sindex2,sindex1[bi]
							H2ind++
						}
					}
					if (hastb) {
						for (bi=1;bi<=Ntbs;bi++) {
							Sindex1 = Sindex1,bind
							Sindex2 = Sindex2,sindex2[bi]
							H2ind++
						}
					}
					if (haszb) {
						for (bi=1;bi<=Nzbs;bi++) {
							Sindex1 = Sindex1,bind
							Sindex2 = Sindex2,sindex3[bi]
							H2ind++
						}
					}
					
					Sindex1 = Sindex1,bind
					Sindex2 = Sindex2,bind
					bind++
					bind2++
					
					gml.survind 	= 0
					
					H2[,H2ind++] = uhtred_panels2(1,gml,
							&_hessian_rp_vcv(),bi,bi2)
				}
			}
			
		}
		
	}
	
	
	//final hessian
	
// 	dlogl /dbeta = (d logl/ d beta) / L
//	
// 	1/L (dd2 - d1 * d2)
//	
// 	1/L * dd2 - d1 * d2 / (L)
//	
// 	G = G :/ exp(lnfi)

	HH = quadcolsum((H2 :/ exp(lnfi) :- G[,Sindex1] :* G[,Sindex2]))
	
	bi = 1
	for (i=1;i<=Nbs;i++) {
		j=1
		while (j<=i) {
			H[i,j] = H[j,i] = HH[bi++]
			j++
		}
	}
	
	H
	
}

`RC' uhtred_panels2(`RS' index,		/// -level-
		   `gml' gml,		/// -uhtred object-
		   `PS' func,		/// -function to call-
		   `RS' Xindex, 	///
		   `RS' Xindex2)	/// -design matrix element-
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
			resq[,q] = panelsum(uhtred_panels(index2,gml,func,Xindex,Xindex2),panelindex)
		}
	}
	else {
		for (j=1;j<=gml.Nmodels;j++) {
			gml.model = gml.modtoind = j
			resq2 = (*func)(gml,Xindex,Xindex2)
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


end
