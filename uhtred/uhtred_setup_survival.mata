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

void uhtred_setup_survival(`gml' gml)
{	
	gml.survind	= 0
	gml.Nsurv 	= J(gml.Nmodels,7,0) 
	gml.surv_index 	= asarray_create("real",2)
	familys 	= gml.familys
	gml.NI 		= J(gml.Nmodels,1,0)
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		gml.model = gml.modtoind = mod
		fam = familys[mod]

		if (fam=="exponential" | fam=="weibull" | 
			fam=="gompertz" | (fam=="user" & 
			st_local("hazfunction"+strofreal(mod))!="") | 
			(fam=="user" & 
			st_local("loghazfunction"+strofreal(mod))!="") | 
			(fam=="user" & 
			st_local("failure"+strofreal(mod))!="")) {
			if (gml.hastb[mod,1]) gml.NI[mod] = 1
		}
		if (fam=="loghazard" | fam=="addhazard") {
			gml.NI[mod] = 1
		}
		if (gml.hasbh[mod,2]) gml.NI[mod] = 1

		if (sum(gml.hasbh[,2])) {
			printf("Setting up expected rates...")
			displayflush()
			gml.bhazard = asarray_create("real",2)
		}
		
		if (st_local("pchintpoints")!="") {
			gml.chip = strtoreal(st_local("pchintpoints"))
		}
		else {
			chipstr = st_local("chintpoints")
			if (chipstr=="") gml.chip = 30
			else 		 gml.chip = strtoreal(chipstr)
		}
			
		fam = familys[mod]
		
		if (gml.failures[mod]!="") {
			
			y = uhtred_util_depvar(gml)

			if (gml.hasbh[mod,2]) {
				uhtred_setup_check_expratevars(gml)
				uhtred_setup_get_expratefile(gml)
			}
			
			//exactly observed events
			//-> hazard function	
			index 	= selectindex(y[,2]:==1)
			Ns 	= rows_cols(index)
			if (Ns[1] & Ns[2]) {
				asarray(gml.surv_index,(mod,1),index)
				gml.Nsurv[mod,1] = Ns[1]
				
				//expected hazard at event times
				if (gml.hasbh[mod,1] | gml.hasbh[mod,2]) {
					uhtred_setup_exph1(gml)
				}
			}

			//exactly observed events and/or right censoring
			//-> survival function
			index 	= selectindex(y[,2]:<2)
			Ns 	= rows_cols(index)
			if (Ns[1] & Ns[2]) {
				asarray(gml.surv_index,(mod,2),index)
				gml.Nsurv[mod,2] = Ns[1]
				if (gml.hasbh[mod,2]) uhtred_setup_exph2(gml)
			}
			
			//interval censoring
			//-> cdf function
			//-> left interval handled separately in case of 0s
			if (gml.haslint[mod]) {
				index 	= selectindex(y[,2]:==2)
				Ns 	= rows_cols(index)
				if (Ns[1] & Ns[2]) {
					asarray(gml.surv_index,(mod,3),index)
					gml.Nsurv[mod,3] = Ns[1]
				}
				ind2 = 3
				if (gml.hasltrunc[mod]) ind2 = 4
				index 	= selectindex((y[,2]:==2) :* (y[,ind2]:>0))
				Ns 	= rows_cols(index)
				if (Ns[1] & Ns[2]) {
					asarray(gml.surv_index,(mod,5),index)
					gml.Nsurv[mod,5] = Ns[1]
				}
			}

			//left truncation
			//-> survival function
			if (gml.hasltrunc[mod]) {
				index 	= selectindex(y[,3]:>0)
				Ns 	= rows_cols(index)
				if (Ns[1] & Ns[2]) {
					asarray(gml.surv_index,(mod,4),index)
					gml.Nsurv[mod,4] = Ns[1]
					
				}
				if (!gml.Nsurv[mod,4]) {
					gml.hasltrunc[mod] = 0
				}
			}
			
			//right censored - needed for generalised gamma/log normal/ log logistic
			if (gml.familys[mod]=="ggamma" | 
				gml.familys[mod]=="lognormal" | 
				gml.familys[mod]=="loglogistic") {
				
				index = selectindex(y[,2]:==0)
				Ns = rows_cols(index)
				if (Ns[1] & Ns[2]) {
					asarray(gml.surv_index,(mod,6),index)	
					gml.Nsurv[mod,6] = Ns[1]
				}
			}
		
			// CH quadrature nodes
			if (gml.NI[mod]) {
				gml.survind 	= 2
				Nobs2 		= uhtred_util_nobs(gml)
				if (Nobs2) {
					index2 	= uhtred_get_surv_index(gml)
					Ngq 	= gml.chip
					chq2 	= J(Nobs2,Ngq,0)
					gq 	= uhtred_glegendre(Ngq)
					y2	= y[index2,1] :/ 2
					if (gml.hasltrunc) {
						y02   = y[index2,3] :/ 2
						gml.y2 = (y2 :- y02)
						gml.chqp    = gml.y2 :* J(Nobs2,1,gq[,1]') :+ y2 :+ y02
						gml.chqw    = gq[,2]
					}
					else {
						gml.y2 = y2
						gml.chqp = y2 :* J(Nobs2,1,gq[,1]') :+ y2
						gml.chqw = gq[,2]
					}
				}
			}
		
		} //endif

	}	
		
	if (sum(gml.hasbh[,2])) {
		printf("done\n")
		displayflush()
	}
	
}

end

