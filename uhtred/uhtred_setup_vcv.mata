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

void uhtred_setup_vcv(`gml' gml, `RS' coef)
{
	if (gml.Nrelevels) {
		uhtred_init_vcvs(gml)
		uhtred_parse_covstructures(gml)
		uhtred_parse_vcv_eqns(gml,coef)
		uhtred_init_integration(gml)
		prolog = "derivprolog(uhtred_prolog())"
		st_local("mlprolog",prolog)
	}
}

void uhtred_init_integration(`gml' gml)
{	
	//setup level specific integration techniques
	
	if (st_local("intmethod")=="") {
		gml.usegh = J(gml.Nrelevels,1,1)
		gml.adapt = J(gml.Nrelevels,1,1)
	}
	else {
		intmethod 	= tokens(st_local("intmethod"))'
		if (rows(intmethod)!=gml.Nrelevels) intmethod = J(gml.Nrelevels,1,intmethod)
		strlen 		= strlen(intmethod)
		gml.usegh	= gml.adapt = J(gml.Nrelevels,1,.)
		for (i=1;i<=gml.Nrelevels;i++) {
			if (intmethod[i]==substr("ghermite",1,max((strlen[i],2)))) {
				gml.usegh[i] = 1
				gml.adapt[i] = 0
			}
			else if (intmethod[i]==substr("mvaghermite",1,max((strlen[i],5)))) {
				gml.usegh[i] = 1
				gml.adapt[i] = 1
			} 
			else if (intmethod[i]==substr("mcarlo",1,max((strlen[i],2)))) {
				gml.usegh[i] = 0
				gml.adapt[i] = 0
			}
			else {
				errprintf("Invalid intmethod()\n")
				exit(198)
			}
		}
		if (rows(intmethod)>1 & rows(intmethod)!=gml.Nrelevels) {
			errprintf("Invalid intmethod()\n")
			exit(198)
		}
	}
	
	//setup integration points
	
	gml.ip = J(gml.Nrelevels,1,0)
	
	if (st_local("intpoints")=="") {							//defaults
		
		for (i=1;i<=gml.Nrelevels;i++) {
			if (gml.usegh[i]) 	gml.ip[i] = 7
			else 			gml.ip[i] = gml.Nres[i] * 50
		}
		
	} 
	else {

		ip = strtoreal(tokens(st_local("intpoints")))'
		if (rows(ip)>1 & rows(ip)!=gml.Nrelevels) {
			errprintf("Invalid intpoints()\n")
			exit(198)
		}
		if (rows(ip)==1) {
			for (i=1;i<=gml.Nrelevels;i++) {
				gml.ip[i] = ip
			}
		}
		else gml.ip = ip

		for (i=1;i<=gml.Nrelevels;i++) {
			if (gml.usegh[i]) {
				if (gml.ip[i]<3& gml.ip[i]>150) {
					errprintf("intpoints() must be between 3 and 150 with Gauss-Hermite quadrature\n")
					exit(198)
				}
			}
			else {
				if (gml.ip[i]<10) {
					errprintf("intpoints() must be >=10 with Monte Carlo integration\n")
					exit(198)
				}
			}
		}
	}
	
	//core index that gets updated through the recursive function
	
	gml.qind = 1,J(1,gml.Nrelevels,0)
	
	if (sum(gml.adapt) | gml.todo | gml.predict) {
		gml.Li_ip = asarray_create("real",cols(gml.qind))
		gml.c 	  = asarray_create("real",cols(gml.qind))
	}
	
	if (max(gml.usegh)==1) {
		if (max(gml.adapt :* gml.usegh)==1) {
			gml.atol = 1e-08
			gml.Pupdateip = &uhtred_gh_update_ip_alllevs()
			gml.showadapt = st_local("showadapt")!=""
			if (st_local("adaptiterations")=="") gml.adaptit = 1001
			else gml.adaptit = strtoreal(st_local("adaptiterations"))
		}
	}
	if (min(gml.usegh)==0) {	
		gml.seed = rseed()
		gml.bdraws = asarray_create("real",1)
		gml.atol = 1e-08
		if (max((gml.usegh:==0) :* gml.adapt)==1) {
			gml.Pupdateip = &uhtred_mc_update_ip()
			gml.showadapt = st_local("showadapt")!=""
			if (st_local("adaptiterations")=="") gml.adaptit = 1001
			else gml.adaptit = strtoreal(st_local("adaptiterations"))
		}
	}
	
	//setup NI
	uhtred_init_ip(gml)
}

void uhtred_init_ip(`gml' gml)
{

	//initialise stuff needed if at least one level's integ is adaptive
	if (sum(gml.adapt)) {
		gml.aghip 			= asarray_create("real",2)
		gml.aghip2 			= asarray_create("real",2)
		gml.aghlogl			= asarray_create("real",1)	
		gml.stackednodes 		= asarray_create("real",1)
		gml.adpanelindexes 		= asarray_create("real",2)
		for (j=1;j<=1;j++) {
			ids = st_data(.,gml.levelvars[1,1..gml.Nrelevels],gml.touse)
			for (i=1;i<gml.Nlevels;i++) {
				if (gml.adapt[i]) {
					asarray(gml.adpanelindexes,(i,j),ids[,i])
				}
			}
		}
	}
	else {
		gml.b = asarray_create("real",1)
	}

	//initialise stuff for GHQ
	if (sum(gml.usegh)) {
		gml.baseGHnodes = gml.baseGHweights = asarray_create("real",1)
	}
	
	gml.ndim = J(gml.Nrelevels,1,0)
	
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
