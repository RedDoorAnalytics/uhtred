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

`TR' uhtred_setup_el_rcs(`gml' gml, `SS' elsyn, `SS' stub, 
			 `RS' model, `RS' cmp, `RS' el)
{
	pos1 = strpos(elsyn,"(")
	pos2 = strlen(elsyn)
	synt = substr(elsyn,pos1+1,pos2-pos1-1)
	resps = tokens(gml.responses[model])
	touse = gml.modeltouses[model]
	
	stata("capture drop "+stub+"_*")
	stcall = "uhtred_rcs " + synt
	stcall = stcall + " touse(" + touse + ")"
	stcall = stcall + " stub(" + stub + ")"
	if (gml.predict) {
		stcall = stcall + " predict"
		stcall = stcall + " model(" + strofreal(model) + ")"
		stcall = stcall + " cmp(" + strofreal(cmp) + ")"
		stcall = stcall + " el(" + strofreal(el) + ")"
	}

	if (gml.issurv[model]) {
		stcall = stcall + " eventvar(" + resps[2] + ")"
	}
	stata(stcall)

	//store elinfo
	yvar 	= st_global("r(yvar)")
	logopt 	= st_global("r(log)")
	knots 	= st_global("r(knots)")
	rmat 	= st_matrix("r(rmatrix)")
	hasrmat = rmat!=J(0,0,.)

	offset		= st_global("r(offset)")
	hasoffset 	= offset!=""
	offsetist0var = 0
	if (gml.hasltrunc[model]) {
		offsetist0var = offset==resps[3]
	}
	moffset		= st_global("r(moffset)")
	hasmoffset 	= moffset!=""
	moffsetist0var = 0
	if (gml.hasltrunc[model]) {
		moffsetist0var = moffset==resps[3]
	}

	info = asarray_create("real",1)
	asarray(info,3,knots)
	asarray(info,4,logopt!="")
	asarray(info,5,hasrmat)
	if (hasrmat) asarray(info,6,rmat)
	asarray(info,7,(hasoffset\offsetist0var))
	if (hasoffset & !offsetist0var) {
		toff = st_data(.,offset,touse)
		asarray(info,8,toff)
		asarray(info,11,offset)
	}
	asarray(info,9,(hasmoffset\moffsetist0var))
	if (hasmoffset & !moffsetist0var) {
		toff = st_data(.,moffset,touse)
		asarray(info,10,toff)
		asarray(info,12,moffset)
	}

	asarray(gml.elinfo,(model,cmp,el),info)
}


`SS' uhtred_stata_setup_rcs(`SS' touse, `SS' yvar, `SS' options, 
			`SR' knots, `RM' rmat, `SS' stub)
{
	stata("cap drop " + stub + "_?")
	synt = "," + options
	stcall = "uhtred_rcs "+yvar + synt
	stcall = stcall + " touse("+touse+")"
	stata(stcall + " stub(" + stub +")")
	knots = st_global("r(knots)")
	rmat = st_matrix("r(rmatrix)")
	hasrmat = rmat!=J(0,0,.)
	return(st_global("r(varnames)"))
}

`SS' uhtred_stata_rcs(`SS' touse, `SS' yvar, `SR' knots, `RM' rmat, 
			`SS' opts, `SS' stub)
{
	hasrmat = rmat!=J(0,0,.)
	stata("cap drop " + stub + "_?")
	stcall = "uhtred_rcs " + yvar + ", " + opts
	stcall = stcall + " touse("+touse+")"
	stcall = stcall + " knots("+knots+")"
	if (hasrmat) {
		st_matrix("_rmatrix",rmat)
		stcall = stcall + " rmatrix(_rmatrix)"
	}		
	stata(stcall + " stub("+ stub +")")
	return(st_global("r(varnames)"))
}

`RM' uhtred_rcs(real colvector variable,	///
                real rowvector knots,| 		///
                real scalar deriv,		///
                real matrix rmatrix		///
        )
{
	if (args()==2) {
		deriv = 0
	}
	hasrmat = args()==4

	// extract knot locations //

	Nobs 	= rows(variable)
	Nknots 	= cols(knots)
	kmin 	= knots[1,1]
	kmax 	= knots[1,Nknots]

	if (Nknots==2) 	interior = 0
	else 		interior = Nknots - 2
	Nparams = interior + 1
		
	// calculate splines //

	splines = J(Nobs,Nparams,.)
	ind1 = 1..Nparams
	ind2 = 2..Nparams
	
	if (Nparams>1) {
		lambda = J(Nobs,1,(kmax:-knots[,ind2]):/(kmax:-kmin))
		knots2 = J(Nobs,1,knots[,ind2])
	}

	if (deriv==0) {
		splines[,1] = variable
		if (Nparams>1) {
			splines[,ind2] = (variable:-knots2):^3 :* 
				(variable:>knots2) :- 
				lambda :* 
				((variable:-kmin):^3):*(variable:>kmin) :- 
				(1:-lambda) :* 
				((variable:-kmax):^3) :* (variable:>kmax) 
		}
	}
	else if (deriv==1) {
		splines[,1] = J(Nobs,1,1)
		if (Nparams>1) {
			splines[,ind2] = 3 :* (variable:-knots2):^2 :* 
				(variable:>knots2) :- 
				lambda:*(3:*(variable:-kmin):^2) :* 
				(variable:>kmin) :- 
				(1:-lambda) :* (3:*(variable:-kmax):^2) :* 
				(variable:>kmax) 	
		}
	}
	else if (deriv==2) {
		splines[,1] = J(Nobs,1,0)
		if (Nparams>1) {
			splines[,ind2] = 6:*(variable:-knots2) :* 
				(variable:>knots2) :- 
				lambda:*(6:*(variable:-kmin)) :* 
				(variable:>kmin) :- 
				(1:-lambda):*(6:*(variable:-kmax)) :* 
				(variable:>kmax) 	
		}
	}
	else if (deriv==3) {
		splines[,1] = J(Nobs,1,0)
		if (Nparams>1) {
			splines[,ind2] = 6 :* (variable:>knots2) :- 
				lambda :* 6 :* (variable:>kmin) :- 
				(1:-lambda) :* 6 :* (variable:>kmax)
		}
	}
	
	//orthog
	if (hasrmat) {
		rmat = luinv(rmatrix)
		if (deriv==0) {
			splines = (splines,J(Nobs,1,1)) * rmat[,ind1]
		}
		else {
			splines = splines * rmat[ind1,ind1]
		}
	}
	return(splines)
}

`RM' uhtred_orthog(`RM' x)
{
	meanx = mean(x)
	colsx = cols(x)
	v = x :- meanx ,J(rows(x),1,1) 
	rowsv = rows(v)
	colsv = cols(v)
	q = J(rowsv,0,.)
	R = J(colsv,colsv,0)
	R[colsv,] = (meanx,1)
	for (i=1;i<=colsx;i++){
                r = norm(v[,i])/sqrt(rowsv)
                q = q, (v[,i]:/ r)
                R[i,i] = r
                for (j = i + 1; j<=colsx; j++){
                        r = (q[,i]' * v[,j])/rowsv
                        v[,j] = v[,j] - r*q[,i]
                        R[i,j] = r 
                }
	}
	return(R)
}

end
