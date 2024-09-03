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

void uhtred_setup_touses(`gml' gml)
{
	//get any variable names for touses
	allvars	= J(1,0,"")

	for (i=1;i<=gml.Nmodels;i++) {

		modelvars = J(1,0,"")
		cmps 	  = uhtred_get_cmps(i)
		Ncmps 	  = cols(cmps)

		//start with any level vars if there are some
		if (gml.Nrelevels) modelvars = modelvars,gml.levelvars[1,]
		
		//extract all variable names from complex predictor
		for (j=1;j<=Ncmps;j++) {
			
			//prep
			dv = strtrim(cmps[j])
			_uhtred_remove_at(dv)
			
			//one element, and its a variable
			pos = 1	
			if (!strpos(dv,"[") & !strpos(dv,"(") & 
				!strpos(dv,"#")) {
				pos = 0
				if (!(gml.predict & 
					dv==st_global("e(timevar" + 
						strofreal(i)+")"))) {
					modelvars = modelvars,dv
				}
			}

			//hash it up
			while (pos) {
				pos1 = strpos(dv,"##")
				pos2 = strpos(dv,"#")
				pos  = max((pos1,pos2))
				if (pos) {
					step = 1
					if (pos1==pos) step = 2
					dv2 = substr(dv,1,pos-1)
					dv  = substr(dv,pos+step,.) 
				}
				else dv2 = dv
				if (!strpos(dv2,"[") & !strpos(dv2,"(")) { 
					//variable name
					if (!(gml.predict & 
						dv2==st_global("e(timevar" + 
							strofreal(i)+")"))) {
						modelvars = modelvars,dv2
					}
				}
				//rcs()/fp()/bs()/ps()/gp() varname
				if (strpos(dv2,"pc(") | strpos(dv2,"fp(") | 
					strpos(dv2,"rcs(") | strpos(dv2,"bs(") | 
					strpos(dv2,"ps(") | strpos(dv2,"gp(")) {
					p1 = strpos(dv2,"(")+1
					len = strpos(dv2,",") - p1
					tpvar = strtrim(substr(dv2,p1,len))
					if (!(gml.predict & 
						tpvar==st_global("e(timevar" + 
							strofreal(i)+")"))) {
						//if predicting and the var is 
						//equal to timevar then skip
						modelvars = modelvars,tpvar
					}
					//get offset/moffset
					len = strlen(dv2) - p1
					synt = "local 0 "+substr(dv2,p1,len)
					stata(synt)
					rc = _stata("syntax varname, [offset(varname) moffset(varname) *]")
					if (rc) exit(rc)
					
					if (st_local("offset")!="") {
						modelvars = modelvars,
							st_local("offset")
					}
					if (st_local("moffset")!="") {
						modelvars = modelvars,
							st_local("moffset")
					}
				}
			}
		}

		//timevar
		//-> if not predicting as it can overide predicts timevar
		if (!gml.predict) {
			modelvars = modelvars,
				st_local("timevar"+strofreal(i))
		}
		//markout on predicts timevar, if there
		//unless obtaining a standardised prediction
		if (st_local("ptvar")!="" & st_local("standardise")=="") {
			modelvars = modelvars,st_local("ptvar")
		}

		//response vars
		//-> if survival model, can override predicts timevar
		//-> but needed if reffects specified
		predandsurv = gml.predict & st_local("failure" + 
			strofreal(i))!="" & st_local("reffects")==""
		if (!gml.predict) {
			modelvars = modelvars,tokens(st_local("response"+strofreal(1)))
			//remove linterval if there, as missings are allowed
			if (gml.haslint[i]) {
				allvars = allvars,modelvars[1,cols(modelvars)]
				modelvars[1,cols(modelvars)] = ""
			}
		}

		//collect all variables
		allvars = allvars,modelvars

		//post complete case model specific touse var
		modeltousei = gml.touse+strofreal(i)
		st_local(modeltousei,modeltousei)
		
		//avoid model specific ifs/ins if being called from predict
		if (gml.predict) {
			stata("mark "+st_local(modeltousei))
		}
		else {
			stata("mark "+st_local(modeltousei) + " " + 
				st_local("if"+strofreal(i)) + " " + 
				st_local("in"+strofreal(i)))
		}

		//markout vars, and update based on main touse 
		//(markout only marks out missing obs of vars)
		stata("markout "+st_local(modeltousei)+" "+invtokens(modelvars))
		stata("qui replace "+st_local(modeltousei)+"=0 if !"+gml.touse)
		gml.modeltouses = gml.modeltouses,st_local(modeltousei)

		//to post for predict(ms)
		gml.allvars = invtokens(uniqrows(allvars')')
		
		//if being called from predict(ms), update on N
		if (gml.predict) {
			if (st_local("npredict")!="") {
				stata("qui replace "+st_local(modeltousei)+
					" = _n <= "+st_local("npredict"))
			}
		}

		//update global touse based on modeltouses, as models can have 
		//different # of obs, but everything must get read in
		//for possible indexing from other models
		if (i==gml.Nmodels) st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 ")
		else st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 | ")
	
	}
	
	stata("qui replace "+gml.touse+" = ("+st_local("touseup")+")")
	stata("qui count if "+gml.touse+"==1")
	gml.N = st_numscalar("r(N)")		
		
}

end
