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

void uhtred_build_clp(`gml' gml, `RS' coef)
{
	eqn = coef = 1
	gml.hasxb   = gml.hastb = gml.haszb = J(gml.Nmodels,1,0)
	gml.xeqn    = gml.teqn = gml.hasxb
	gml.zeqn    = J(1,gml.Nrelevels,gml.hasxb)
	gml.Ncmps   = 0
	gml.Nels    = asarray_create("real",1)
	gml.elindex = asarray_create("real",2)
	gml.eltvar  = asarray_create("real",2)
	gml.bindices = asarray_create("real",2)
	gml.Hindices = asarray_create("real",3)
	gml.tvarlist = J(gml.Nmodels,1,"")

// 	gml.X = gml.XT = asarray_create("real",1)
// 	if (gml.hasanyltrunc) gml.XT0 = gml.X
// 	if (gml.hasanylint) gml.XTL = gml.X
	if (gml.Nrelevels) {
		gml.Z = asarray_create("real",2)
		gml.Zbindex = asarray_create("real",2)
	}
// 	if (sum(gml.familys:=="rp")) {
// 		gml.dXT = gml.X
// 	}
	
	for (i=1;i<=gml.Nmodels;i++) {
		uhtred_build_xz(gml,i,eqn,coef)
	}
	
}

void uhtred_build_xz(`gml' gml, `RS' model, `RS' eqn, `RS' coef)
{
	xbsyn     = uhtred_get_cmps(model)
	if (xbsyn=="") {
		gml.Ncmps = Ncmps = 0
	}
	else gml.Ncmps = Ncmps = cols(xbsyn)

	hascons   = gml.hascons[model]
	haslt 	  = gml.hasltrunc[model]
	hasic 	  = gml.haslint[model]
	yvars     = tokens(gml.responses[model])
	yvar	  = yvars[1]
	strmod	  = strofreal(model)
	touse 	  = gml.modeltouses[model]
	Nobs	  = gml.Nobs[gml.Nlevels,model]
	
	if (gml.issurv[model]) {
		dvar = yvars[2]
		if (haslt) {
			y0var = yvars[3]
			if (hasic) ylvar = yvars[4]
		}
		else {
			if (hasic) ylvar = yvars[3]
		}
	}
	
	/*
	xbsyntax is a SR with each element a component
	- loop over and check for rcs(), fp() etc.
	  -> for each create basis functions
	  -> subinstr rcs() etc. for new basis functions
	- now rhs contains a rhs syntax which varlist can handle directly
	- create equations for ml model
	
	timexbsyn is an extra linear predictor containing only 
	time-dependent functions which get extracted out and then 
	processed in the same way as xbsyntax
	*/
	
	Nelsmat  = J(0,1,.)
	newxbsyn = timexbsyn = J(1,0,"")
	getdt    = gml.familys[model]=="rp" | 
			gml.familys[model]=="logchazard" 
	
	if (getdt) dtimexbsyn = newxbsyn
	if (haslt) time0xbsyn = newxbsyn
	if (hasic) timelxbsyn = newxbsyn
	if (gml.Nrelevels) {
		rexbsyn     = J(gml.Nrelevels,1,"")
		renocons    = J(gml.Nrelevels,1,", noconstant")
	}
	
	for (i=1;i<=Ncmps;i++) {

		//check for at
		//-> extract constraint value if there
		//-> strip off at
		atpos	= strpos(xbsyn[i],"@")
		if (atpos) {
			at = substr(xbsyn[i],atpos+1,.)
			work = substr(xbsyn[i],1,atpos-1)
		}
		else {
			at = ""
			work = xbsyn[i]
		}
		
		mainwork = mainworkt = mainworkt0 = mainworkl = work
		eltype = J(0,1,.)
		eltvar = J(0,2,.)
		elsyn  = J(0,1,"")
		Nels   = 0

		//if needed, split up based on # or ##
		pos 	= 1	
		while (pos) {
			pos1 = strpos(work,"##")
			pos2 = strpos(work,"#")
			pos  = max((pos1,pos2))
			if (pos) {
				step = 1
				if (pos1==pos) step = 2
				work2 = substr(work,1,pos-1)
				work  = substr(work,pos+step,.) 
			}
			else work2 = work

			if (work2!="") {
				Nels++
				eltype 	= eltype\uhtred_get_element_codes(gml,1,work2)
				eltvar	= eltvar\uhtred_get_element_istimevar(gml,1,work2)
				elsyn 	= elsyn\work2
			}
		}

		istimedep = max(eltvar[,1])
		hasre 	  = any(eltype:==2)
		Nelsmat   = Nelsmat\Nels
		
		asarray(gml.elindex,(model,i),eltype)
		asarray(gml.eltvar,(model,i),eltvar)

		//now build element vars if needed
		for (j=1;j<=Nels;j++) {
			if (eltype[j]==8) {
				stub  = "_rcs" + strofreal(model) +
					"_" + strofreal(i) + 
					"_"+strofreal(j)
					
				uhtred_setup_el_rcs(gml,elsyn[j],stub,
							model,i,j)
				//now strsub in new stub 
				mainwork = subinstr(mainwork,elsyn[j],
						"c."+stub+"_*",1)

				if (getdt & eltvar[j,1]) {
					elinfo = asarray(gml.elinfo,(model,i,j))
					knots  = asarray(elinfo,3)
					rmat   = asarray(elinfo,6)
					opts   = ""
					if (asarray(elinfo,4)) opts = "log"
					offset = asarray(elinfo,7)
					if (offset[1]) {
						opts = opts + " offset(" + 
							asarray(elinfo,11) + ")"
					}
					moffset = asarray(elinfo,9)
					if (moffset[1]) {
						opts = opts + " moffset(" + 
							asarray(elinfo,12) + ")"
					}

					//get deriv splines
					stub  = "_d_rcs"+ strofreal(model) +
					"_"+strofreal(i)+"_"+strofreal(j)
					svars = uhtred_stata_rcs(touse,yvar,knots,
						rmat,opts+" deriv(1) timest",
						stub)

					//now strsub in new stub 
					mainworkt = subinstr(mainworkt,
							elsyn[j],
							"c."+stub+"_*",1)
					//linterval
					if (hasic) {
						stub  = "_l_rcs"+ strofreal(model) +
					"_"+strofreal(i)+"_"+strofreal(j)
						svars = uhtred_stata_rcs(touse,ylvar,knots,
							rmat,opts,stub)
						mainworkl = subinstr(mainworkl,
								elsyn[j],
								"c."+stub+"_*",
								1)
					}
					//ltruncation
					if (haslt) {
						stub  = "_s0_rcs"+ strofreal(model) +
					"_"+strofreal(i)+"_"+strofreal(j)
						svars = uhtred_stata_rcs(touse,y0var,knots,
							rmat,opts,stub)
						mainworkt0 = subinstr(
								mainworkt0,
								elsyn[j],
								"c."+stub+"_*",
								1)
					}
				}
			}
			else if (eltype[j]==2) {
				elinf = asarray(gml.elinfo,(model,i,j))
				relev = asarray(elinf,1)
				reel = j
			}
			else if (eltype[j]==1) {
				if (hasre & !strpos(elsyn[j],"c.")) {
					mainwork = subinstr(mainwork,elsyn[j],"c."+elsyn[j],1)
				}
			}

		}

		if (hasre) {
			
			//now process re
			//note: they must be done after any basis variables 
			//      are created in above el loop
			allres = asarray(gml.latlevs,relev)
			
			newrework = uhtred_remove_re(mainwork,allres,lat="")
			//!!fvexpand, get columns and match index
			reindex = selectindex(lat:==allres[,1])
			
			oldreindex = asarray(gml.Zbindex,(model,relev))
			if (oldreindex==J(0,0,.)) {
				asarray(gml.Zbindex,(model,relev),reindex)
			}
			else {
				asarray(gml.Zbindex,(model,relev),
					asarray(gml.Zbindex,(model,relev))\reindex)
			}

			if (istimedep) {
				//!!
			}
			else {
				reconsname = ""
				if (newrework=="") {	//random intercept
					reconsname = "_re_"+lat
					stata("cap drop "+reconsname)
					stata("qui gen byte "+reconsname+" = 1 if "+touse)
					rexbsyn[relev] =  rexbsyn[relev] + " " + reconsname
					renocons[relev] = ""
				}		
				else {	//random coefficient
					rexbsyn[relev] =  rexbsyn[relev] + " " + invtokens(newrework)
				}
			}
			
		}
		else {
			if (istimedep) {
				timexbsyn = timexbsyn,mainwork
				if (getdt) {
					dtimexbsyn = dtimexbsyn,mainworkt
					if (haslt) time0xbsyn = 
							time0xbsyn,mainworkt0
					if (hasic) timelxbsyn = 
							timelxbsyn,mainworkl
				}
			}		
			else 	newxbsyn  = newxbsyn,mainwork
		}
		
		//handle at
		if (at!="") {
			if (hasre) {
				_uhtred_build_at_on_re(at,model,relev,
							newrework,reconsname)
			}
			else {
				if (istimedep) {
					ateqn = "[tb" + strofreal(model) + "]"
				}
				else 	ateqn = "[xb" + strofreal(model) + "]"
				_uhtred_build_at_on_re(at,model,ateqn,mainwork)
			}
		}
		
	}

	asarray(gml.Nels,model,Nelsmat)
	
	//baseline time splines if needed
	if (gml.familys[model]=="rp") {
		rcsopts = st_local("rcsopts"+strmod)
		opts = st_local("rcsopts"+strmod) + 
			" log event" +
			" eventvar("+dvar+")"
			
		if (gml.predict) {
			opts = opts + " predict baseline model("+strmod+")"
		}
		
		if (!strpos(rcsopts,"noorthog")) opts = opts + " orthog"
		timexbsyn = timexbsyn,
				uhtred_stata_setup_rcs(touse,yvar,opts,
					knots="",rmat=J(0,0,.,),
					"_brcs" + strmod)
		//store for eret list
		asarray(gml.distancb,(model,3),strtoreal(tokens(knots)))
		if (rmat==J(0,0,.)) {
			asarray(gml.distancb,(model,4),0)
		}
		else {
			asarray(gml.distancb,(model,4),1)
			asarray(gml.distancb,(model,5),rmat)
		}
		
		//deriv splines
		dtimexbsyn = dtimexbsyn,
				uhtred_stata_rcs(touse,yvar,knots,
						rmat,"log deriv(1)",
						"_d_brcs" + strmod)
		//ltruncation
		if (haslt) {
			time0xbsyn = time0xbsyn,
					uhtred_stata_rcs(touse,y0var,knots,
							rmat,"log",
							"_s0_brcs" + strmod)
		}
		//linterval
		if (hasic) {
			timelxbsyn = timelxbsyn,
					uhtred_stata_rcs(touse,ylvar,knots,
							rmat,"log",
							"_l_brcs" + strmod)
		}
	}	

	//get equations and design matrices //
	
	hasxb 	  = newxbsyn!=J(1,0,"")
	hastimexb = timexbsyn!=J(1,0,"")
	if (hascons) {
		hasxb = 1
		_cons = st_local("_cons")
	}

	//xb equation:
	
	if (hasxb) {
		xb = invtokens(newxbsyn)
		xbeqnlocal = "xb"+strmod
		if (!hascons) xb = xb + ", noconstant" 
		st_local(xbeqnlocal,"("+xbeqnlocal+": "+xb+")")
		gml.hasxb[model] = 1
		gml.xeqn[model] = eqn
		
		//store design matrix for xb
		if (hascons) newxbsyn = newxbsyn,_cons
		gml.X = st_data(.,newxbsyn,touse)
// 		asarray(gml.X,model,X)
		
		//get coefficient indices
		coefx1 = coef
		coefx2 = coef + cols(gml.X) - 1
		asarray(gml.bindices,(model,eqn),(1,coefx1\1,coefx2))

		if (gml.todo==2) {
			asarray(gml.Hindices,(model,eqn,eqn),
				(coefx1,coefx1\coefx2,coefx2))
		}
		
		coef = coefx2 + 1
		eqn++
	}

	// time xb:
	
	if (hastimexb) { 
		xb_time = invtokens(timexbsyn)
		gml.tvarlist[model] = xb_time
		xbtimeeqnlocal = "xb"+strmod+"_time"
		xb_time = xb_time + ", noconstant"
		st_local(xbtimeeqnlocal,"(tb" + strmod + ": " + xb_time+")")
		gml.hastb[model] = 1
		gml.teqn[model] = eqn
		
		//store design matrix for time xb 
		gml.XT = st_data(.,timexbsyn,touse)
// 		asarray(gml.XT,model,XT)
		
		//get coefficient indices
		coeft1 = coef
		coeft2 = coef + cols(gml.XT) - 1
		asarray(gml.bindices,(model,eqn),(1,coeft1\1,coeft2))
		
		if (gml.todo==2) {
			asarray(gml.Hindices,(model,gml.xeqn[model],eqn),
				(coefx1,coeft1\coefx2,coeft2))
			asarray(gml.Hindices,(model,eqn,gml.xeqn[model]),
				(coeft1,coefx1\coeft2,coefx2))	
			asarray(gml.Hindices,(model,eqn,eqn),
				(coeft1,coeft1\coeft2,coeft2))
		}
		
		coef = coeft2 + 1
		eqn++
	}

	if (getdt) {
		//store deriv of XT
		gml.dXT = st_data(.,dtimexbsyn,touse)
// 		asarray(gml.dXT,model,dXT)
		
		if (haslt) {
			//asarray(gml.XT0,model,st_data(.,time0xbsyn,touse))
			gml.XT0 = st_data(.,time0xbsyn,touse)
		}
		if (hasic) {
// 			asarray(gml.XTL,model,st_data(.,timelxbsyn,touse))
			gml.XTL = st_data(.,timelxbsyn,touse)
		}
	}
	
	// re xb :
	
	if (gml.Nrelevels) {
		for (lev=1;lev<=gml.Nrelevels;lev++) {
			Z = st_data(.,rexbsyn[lev],touse)
			asarray(gml.Z,(model,lev),Z)
			st_local("Z"+strmod+"_"+strofreal(lev),"(zb"+strmod+"_" + 
				strofreal(lev) + ":" + rexbsyn[lev]+",nocons)")
			gml.zeqn[model,lev] = eqn
			
			//get coefficient indices
			coef2 = coef + cols(Z) - 1
			asarray(gml.bindices,(model,eqn),(1,coef\1,coef2))
			
			coef = coef2 + 1
			eqn++
		}
		gml.haszb[model] = 1
	}
	
}

/*
uhtred_build_els()
-> get return code for type of element:
	1 - variable
	2 - random effect
	3 - time/variable function
	8 - rcs()
	9 - fp()
	14 - bs()
	15 - ps()
	16 - gp()
	17 - pc()
*/

`RS' uhtred_get_element_codes(`gml' gml, `RS' mod, `SS' dv2)
{
	el 			= dv2
	hassquareb	= strpos(el,"[")
	hasroundb	= strpos(el,"(")
	
	if (hassquareb & !hasroundb) 		return(2)	//random effect
	else {	
		if 	(strpos(el,"mf("))  	return(3)	//user-defined
		else if (strpos(el,"rcs(")) 	return(8)	//rcs function
		else if (strpos(el,"fp(")) 	return(9)	//fp function
		else if (strpos(el,"bs(")) 	return(14)	//bs function
		else if (strpos(el,"pc(")) 	return(17)	//pc function
		else 				return(1)	//varname
	}
}

//flag whether varname, or function input var is timevar()
`RR' uhtred_get_element_istimevar(`gml' gml, `RS' mod, `SS' dv2)
{
	el 		= dv2
	hassquareb	= strpos(el,"[")
	hasroundb	= strpos(el,"(")
	
	if (hassquareb & !hasroundb) {
		varname = substr(el,hassquareb+1,strpos(el,"]")-hassquareb+1)
	}
	else if (hasroundb) {
		varname = substr(el,hasroundb+1,strpos(el,",")-hasroundb-1)
	}
	else 	varname = el
	varname = strtrim(varname)
	
	if (strpos(el,"mf(")) {
		//override for mf()
		varname = gml.tvarnames[mod] 
	}

	return(varname==gml.tvarnames[mod],varname==gml.ltruncated[mod])
}

void _uhtred_build_at_on_re(`SS' at, `RS' model, `RS' relev, `SM' lats, `SS' reconsname)
{
	cneqn 	= "[zb"+strofreal(model)+"_"+strofreal(relev)+"]"
	if (lats=="") {
		stata("constraint free")	
		cn = st_global("r(free)")
		stata("constraint "+cn+cneqn+" : "+reconsname+" = "+at)
		stata("local constraints "+st_local("constraints")+" "+cn)
	}
	else {
		for (r=1;r<=cols(lats);r++) {
			if (!strpos(lats[1,r],"0.")) {
				stata("constraint free")	
				cn = st_global("r(free)")
				stata("constraint "+cn+cneqn+" : "+lats[1,r]+" = "+at)
				stata("local constraints "+st_local("constraints")+" "+cn)
			}
		}
	}
}

void _uhtred_build_at(`SS' at, `SS' ateqn, `SR' mainwork)
{
	mw = tokens(mainwork)
	for (a=1;a<=cols(mw);a++) {
		stata("constraint free")	
		cn = st_global("r(free)")
		stata("constraint "+cn+ateqn+" : "+mw[a]+" = "+at)
		stata("local constraints "+st_local("constraints")+" "+cn)
	}
}


`SS' uhtred_remove_re(`SS' mainwork,`SM' res, `SS' lat)
{
	newwork = mainwork
	done = 0
	for (r=1;r<=rows(res);r++) {
		newwork = subinstr(newwork,"##"+res[r,2],"")
		newwork = subinstr(newwork,res[r,2]+"##","")
		newwork = subinstr(newwork,"#"+res[r,2],"")
		newwork = subinstr(newwork,res[r,2]+"#","")
		newwork = subinstr(newwork,res[r,2],"")
		if (strlen(newwork)!=strlen(mainwork) & !done) {
			lat = res[r,1]
			done = 1
		}
	}
	return(newwork)
}

end
