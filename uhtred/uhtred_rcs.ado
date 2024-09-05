*! version 1.0.0  ?????2024

program uhtred_rcs, rclass
	version 17
	
	
	syntax varname, [				///
				DF(string) 		///
				KNOTS(numlist asc) 	///
				LOG 			///
				ORTHog 			///
				NOORTHog 		///
				RMATRIX(string)		///
				EVent 			///
				EVENTVAR(varname)	///
				OFFset(varname) 	///
				MOFFset(varname)	///
							///
				TOUSE(varname)		///
				MODEL(string)		///
				STUB(string)		///
				DERIV(real 0)		///
				TIMEST			///
			]
	
	if "`df'"=="" & "`knots'"=="" {
		di as error "option {bf:df()} or option " ///
			"{bf:knots()} required"
		exit 198
	}
	if "`df'"!="" & "`knots'"!="" {
		di as error "only one of options {bf:df()} and " ///
			"{bf:knots()} may be specified"
		exit 198
	}
	
	mata: _uhtred_rcs()
	
	return local yvar `varlist'
	return local varnames `varnames'
	return local varnames_stub `varnames_stub'
	return local knots `knots'
	return local log `log'
	return local event event 
	return local eventvar `eventvar'
	return local offset `offset'
	return local moffset `moffset'
	cap confirm matrix _rmatrix
	if !c(rc) {
		return matrix rmatrix _rmatrix
	}
end

version 17

local RS 	real scalar
local SS 	string scalar
local RM 	real matrix
local SM	string matrix
local PC 	pointer colvector
local SC 	string colvector
local SR	string rowvector
local TR	transmorphic
local RC	real colvector
local RR	real rowvector

mata:
`RM' _uhtred_rcs()
{
	varname = st_local("varlist")
	touse   = st_local("touse")
	st_view(rcsvar=.,.,varname,touse)
	
	deriv 		= strtoreal(st_local("deriv"))
	timest 		= st_local("timest")!=""
	islog 		= st_local("log")!=""
	orth 		= st_local("orthog")!=""
	offsetvar 	= st_local("offset")
	hasoffset 	= offsetvar!=""
	moffsetvar	= st_local("moffset")
	hasmoffset	= moffsetvar!=""
	offsetist0var 	= .
	moffsetist0var 	= .

	if (hasoffset)  offset  =  st_data(.,offsetvar,touse)
	if (hasmoffset) moffset = -st_data(.,moffsetvar,touse)

	//get element info
	
// 	if (gml.predict) { //called from predict
// 		index 	= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(Nels)
// 		knots 	= strtoreal(tokens(st_global("e(knots_"+index+")")))
// 		if (orth) rmat 	= st_matrix("e(rmat_"+index+")")
// 	}
// 	else {
		
	
		if (st_local("df")!="") {
			df = strtoreal(st_local("df"))	
			tv = rcsvar
			if (hasoffset) 		tv = tv :+ offset
			if (hasmoffset) 	tv = tv :+ moffset
			if (st_local("event")!="")  {
				ev = st_data(.,st_local("eventvar"),touse)
				tv = select(tv,ev)
			}
			if (islog) tv = log(tv)

			tv	= sort(tv,1)
			nrows 	= rows(tv)
			if (df==1) 			index = 1\nrows
			else {
				if (df==2) 		index = 50
				else if (df==3) index = 33.3333333333\66.66666666
				else if (df==4) index = 25\50\75
				else if (df==5) index = 20\40\60\80
				else if (df==6) index = 17\33\50\67\83
				else if (df==7) index = 14\29\43\57\71\86
				else if (df==8) index = 12.5\25\37.5\50\62.5\75\87.5
				else if (df==9) index = 11.1\22.2\33.3\44.4\55.6\66.7\77.8\88.9
				else if (df==10) index = 10\20\30\40\50\60\70\80\90 
				index = 1\round(index :/100 :* nrows)\nrows
			}
			knots = tv[index]'
			if (orth) {
				tv = rcsvar
				if (hasoffset)  tv = tv :+ offset
				if (hasmoffset) tv = tv :+ moffset
				if (islog) 	tv = log(tv)
				rmat = uhtred_orthog(uhtred_rcs(tv,knots)) 
			}
		}
		else {
			knots = strtoreal(tokens(st_local("knots")))
			if (st_local("rmatrix")!="") {
				rmat = st_matrix(st_local("rmatrix"))
			}
			else if (orth) {
				tv = rcsvar
				if (hasoffset)  tv = tv :+ offset
				if (hasmoffset) tv = tv :+ moffset
				if (islog)	rmat = uhtred_orthog(uhtred_rcs(log(tv),knots)) 
				else 		rmat = uhtred_orthog(uhtred_rcs(tv,knots)) 
			}
		}
	
// 	}

        //check knots
        if (cols(knots)!=rows(uniqrows(knots'))) {
                errprintf("Knot locations not unique\n")
                exit(198)
        }
        
	if (timest & !islog)  rcsvar2 = rcsvar
	
	if (hasoffset) 	rcsvar = rcsvar :+ offset
	if (hasmoffset) rcsvar = rcsvar :+ moffset
	if (orth | st_local("rmatrix")!="") {
		if (islog)      rcsvars = uhtred_rcs(log(rcsvar),knots,deriv,rmat)
		else 		rcsvars = uhtred_rcs(rcsvar,knots,deriv,rmat)
	}
	else {
		if (islog) 	rcsvars = uhtred_rcs(log(rcsvar),knots,deriv)
		else 		rcsvars = uhtred_rcs(rcsvar,knots,deriv)
	}
	
	if (timest & !islog)  rcsvars = rcsvars :* rcsvar2
	
	//build stata vars
	Nvars = cols(rcsvars)
	names = J(1,0,"")
	stub  = st_local("stub")
	
	for (k=1;k<=Nvars;k++) {
		names = names,(stub+"_"+strofreal(k))
	}

	id = st_addvar("double",names)
	st_store(.,id,touse,rcsvars)

	st_local("varnames",invtokens(names))
	st_local("varnames_stub",stub+"_*")
	st_local("knots",invtokens(strofreal(knots)))
	st_matrix("_rmatrix",rmat)
}

end
