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

void uhtred_setup_check_clp(`gml' gml)
{
	for (i=1;i<=gml.Nmodels;i++) {
		st_local("rest",st_local("indepvars"+strofreal(i)))
		stata("gettoken lhs rest : rest, bind")
		uhtred_cmp_errorchecks(gml,strtrim(st_local("lhs")))
		while (st_local("rest")!="") {
			stata("gettoken lhs rest : rest, bind")
			uhtred_cmp_errorchecks(gml,strtrim(st_local("lhs")))
		}
	}
}

void uhtred_cmp_errorchecks(`gml' gml,`SS' dv)
{
	//check each component
	//-> note; not effected by @'s
	pos 	= 1	
	varind 	= 0
	reind 	= 0
	fvi 	= strpos(dv,"i.") | strpos(dv,"i(")
	h2	= strpos(dv,"##")
	hasrcs 	= strpos(dv,"rcs(")

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
		uhtred_dvcheck(gml,dv2,varind,reind)
	}

	if (reind) {
		if (fvi) {
			errprintf("i. notation cannot be combined with a random effect\n")
			exit(1986)
		}
		if (h2) {
			errprintf("## notation cannot be combined with a random effect\n")
			exit(1986)
		}
		if (hasrcs) {
			errprintf("rcs() notation cannot be combined with a random effect\n")
			exit(1986)
		}
	}
	
	if (reind>1) {
		txt = "Only one random effect element is "
		txt = txt + "allowed in each component\n"
		errprintf(txt)
		exit(1986)
	}
}

void uhtred_dvcheck(`gml' gml,`SS' dv2, `RS' varind, reind)
{
	//check parentheses
	pos1 = strpos(dv2,"[")
	pos2 = strpos(dv2,"]")	
	if ((!pos1 & pos2) | (pos1 & !pos2)) {
		errprintf("Missing [ or ]\n")
		exit(1986)
	}

	pos5 = strpos(dv2,"(")
	pos6 = strpos(dv2,")")	
	if ((!pos5 & pos6) | (pos5 & !pos6)) {
		errprintf("Missing ( or )\n")
		exit(1986)
	}
	
	//check at
	if (strpos(dv2,"@")) {			
		at = substr(dv2,strpos(dv2,"@")+1,.)
		if (strtoreal(at)==.) {
			uhtred_error("@"+at+" must be a real number")
		}
		//strip of @ to continue
		dv2 = uhtred_strip_at(dv2)
	}
	
	hasrcs	= substr(dv2,1,4)=="rcs(" 

	//random effect or ?EV[outcome]
	if (pos1 & pos2) {
		//random effect
		uhtred_check_name(dv2)
		uhtred_check_pipes(dv2,pos1,pos2)		
		uhtred_confirm_idvars(dv2,pos1,pos2)	
		reind++
	}
	//rcs()/fp()/mf() syntax checks
	else if (pos5 & pos6 & hasrcs) {
		hasfp	= substr(dv2,1,3)=="fp(" 
		hasmf	= substr(dv2,1,3)=="mf(" 
		hasbs	= substr(dv2,1,3)=="bs(" 
		hasps	= substr(dv2,1,3)=="ps(" 
		haspc	= substr(dv2,1,3)=="pc(" 
		hasfn 	= hasrcs | hasfp | hasmf | hasbs | hasps | haspc
		if (!hasfn) uhtred_error("Invalid function "+dv2)
		pos1 	= strpos(dv2,"(") + 1
		len 	= strlen(dv2) - pos1 
		//-> note: can't match on ")" as it's in rcs(...df()... )
		
		if (hasrcs) {
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) stata("local 0 "+synt)
			else stata("local 0 , "+synt)
			cmdsynt = "syntax anything , [DF(string) "
			cmdsynt = cmdsynt + "KNOTS(numlist asc) LOG ORTHog "
			cmdsynt = cmdsynt + "EVent OFFset(varname) "
			cmdsynt = cmdsynt + "MOFFset(varname)]"
			stata(cmdsynt)
			anything = st_local("anything")
			hasrcsEV = strpos(anything,"EV[")
			
			if (st_local("df")!="") {
				if (st_local("knots")!="") {
					errtxt = "Cannot specify both df() "
					errtxt = errtxt + "and knots() within "
					errtxt = errtxt + "rcs()"
					uhtred_error(errtxt)
				}
				df = strtoreal(st_local("df"))
				if (df<1 | df>10) {
					errtxt = "df() must be between "
					errtxt = errtxt + "1 and 10"
					uhtred_error(errtxt)
				}
				if (df!=round(df)) {
					uhtred_error("df() must be an integer")
				}
			}
			else {
				if (st_local("knots")=="") {
					errtxt = "One of df() or knots() "
					errtxt = errtxt + "needed in rcs()"
					uhtred_error(errtxt)
				}
			}
		}
	}
	else {  //variable
		uhtred_confirmvars(dv2)				
		varind++
	}
}

void uhtred_check_varname(`SS' dv)
{
	pos = strpos(dv,"#")					
	if (pos) {
		uhtred_confirmvars(substr(dv,1,pos-1))	
		dv = substr(dv,pos+1,.)
	}	
}

void uhtred_check_name(`SS' dv)
{	
	`SS' rest
	`RS' hasm
	
	hasm 	= substr(dv,1,1)=="M"
	if (!hasm) {
		errprintf("invalid latent variable specification;\n")
		errprintf("random effect names must start with M\n")
		exit(1986)	
	}
	num = substr(dv,2,1)
	if (strtoreal(num)==.) {
		errprintf("invalid latent variable specification;\n")
		errprintf("M must be followed by an integer\n")
		exit(1986)	
	}
}

void uhtred_check_pipes(`SS' dv, `RS' pos1, `RS' pos2)
{
	`SS' newvar
	newvar = substr(dv,pos1+1,pos2-pos1-1)
	if (strpos(newvar,"<") & strpos(newvar,">")) {
		errprintf("invalid latent variable specification;\n")
		errprintf("both '<' and '>' detected within the ")
		errprintf("same specification\n")
		exit(1986)
	}
}

void uhtred_confirm_idvars(`SS' dv, `RS' pos1, `RS' pos2)
{
	`SS' newvar,copyvar, tempidvar
	`RS' pipepos
	
	ind1 = ind2 = 0
	newvar = substr(dv,pos1+1,pos2-pos1-1)
	copyvar = newvar
	pipepos = strpos(copyvar,"<")
	while (pipepos) {
		ind1 = 1
		tempidvar = substr(copyvar,1,pipepos-1)
		uhtred_confirmvars(tempidvar)
		copyvar = substr(copyvar,pipepos+1,.)
		pipepos = strpos(copyvar,"<")
	}
	copyvar = newvar
	pipepos = strpos(copyvar,">")
	while (pipepos) {
		ind2 = 1
		tempidvar = substr(copyvar,1,pipepos-1)
		uhtred_confirmvars(tempidvar)
		copyvar = substr(copyvar,pipepos+1,.)
		pipepos = strpos(copyvar,">")
	}	
	if (ind1 & ind2) {
		errprintf("Can't have both > and < in level specifier\n")
		exit(1986)
	}
}

void uhtred_confirmvars(`SR' var1)
{
	rc = _stata("fvexpand "+var1)
	if (rc) exit(rc)
	true = st_global("r(fvops)")
	if (true!="true") {	
		rc = _stata("confirm numeric variable " + 
			st_global("r(varlist)"))
		if (rc) exit(rc)
	}
}

void uhtred_confirm_expval(`gml' gml, `RS' k)
{
	f = gml.familys[k]
	if (f=="gompertz" | f=="lquantile" | f=="loghazard" | 
		f=="user" | f=="ordinal") {
		errtxt = "EV[]/XB[] of a family("
		errtxt = errtxt +f+") not currently supported"
		uhtred_error(errtxt)
	}
}

`SS' uhtred_strip_at(`ss' txt)
{
	if (strpos(txt,"@")) {
		return(substr(txt,1,strpos(txt,"@")-1))
	}
	return(txt)
}

end
