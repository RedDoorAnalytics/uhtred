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

void uhtred_setup_check_expratevars(`gml' gml)
{
	matchvars = tokens(st_local("matchby"+strofreal(1)))
	if (matchvars=="") uhtred_error("matchby() required with bhfile()")
	
	for (i=1;i<=cols(matchvars);i++) {
		stata("confirm var "+matchvars[i],1)
	}
}

void uhtred_setup_get_expratefile(`gml' gml)
{
	matchvars = tokens(st_local("matchby"+strofreal(1)))
	
	stata("preserve")
	
	//read in expected rate file
	filename = st_local("bhfile"+strofreal(1))
	stata("confirm file "+filename+".dta")
	stata("qui use "+filename+",clear")
	
	//check core variables exist
	for (i=1;i<=cols(matchvars);i++) {
		stata("confirm var "+matchvars[i])
	}
	stata("confirm var _rate")

	asarray(gml.bhazard,(1,0),st_data(.,("_rate",matchvars)))
	
	stata("restore")
}

void uhtred_setup_exph1(`gml' gml)
{
	gml.survind = 1
	mod			= 1
	if (gml.hasbh[mod,1]) {
		gml.bhazard = st_data(.,st_local("bhaz"+strofreal(mod)),gml.touse[mod])
	}
	else {
		t 		= uhtred_util_depvar(gml)[,1]
		N 		= uhtred_get_nobs(gml)
		matchvars 	= tokens(st_local("matchby"+strofreal(1)))
		Nmatchvars	= cols(matchvars)
		exprate 	= asarray(gml.bhazard,(mod,0))
		x		= st_data(.,matchvars,gml.touse[mod])[uhtred_get_index(gml),]
		//need to add t to -> age and year
		//then floor it to get the expected rate at the age, year previous
		ageyearindex = selectindex(matchvars:==("_age")),selectindex(matchvars:==("_year"))
		x[,ageyearindex] = floor(x[,ageyearindex] :+ t)
		//now get matched rate

		matchedrate = J(N,1,.)
		for (i=1;i<=N;i++) {
			index 		= match_row_index(exprate[|1,2\.,.|],x[i,])
			matchedrate[i] 	= exprate[index,1]
		}
		
		asarray(gml.bhazard,(mod,1),matchedrate)
		
		if (gml.familys[mod]=="rp" | gml.familys[mod]=="logchazard") {
			asarray(gml.bhazard,(mod,3),uhtred_setup_expH(gml,exprate,x))
		}
	}
}

`RS' match_row_index(`RM' xmat, `RR' row)
{
	index = selectindex(rowsum(xmat:==row):==cols(row))
	if (index==J(0,1,.)) {
		printf("\n")
		uhtred_error("No match found for some observations")
	}
	return(index)
}

`RC' uhtred_setup_expH(`gml' gml, `RM' exprate, `RM' x)
{
	mod				= 1
	t 				= uhtred_util_depvar(gml)[,1]
	N 				= uhtred_get_nobs(gml)
	
	matchedH		= J(N,1,0)
	matchvars 		= tokens(st_local("matchby"+strofreal(mod)))
	ageyearindex 	= selectindex(matchvars:==("_age")),selectindex(matchvars:==("_year"))
	
	for (i=1;i<=N;i++) {

 		//build rate at each year from entry to exit
		index = match_row_index(exprate[|1,2\.,.|],x[i,])
 		//how many years
		ft 		= floor(t[i])
		extrat	= t[i]:-ft

		if (ft) {
			xcopy = x[i,]
			for (f=1;f<=ft;f++) {
				matchedH[i] = matchedH[i] + exprate[index,1]
				xcopy[,ageyearindex] = xcopy[,ageyearindex] :+ 1	//increase age and year
				index = match_row_index(exprate[|1,2\.,.|],xcopy)
			} 
			if (extrat>0) matchedH[i] = matchedH[i] + extrat * exprate[index,1]
			
		}
		else matchedH[i] = exprate[index,1] * extrat
	}
	
	return(matchedH)

}

void uhtred_setup_exph2(`gml' gml)
{
	gml.survind 	= 2
	mod				= 1
	exprate 		= asarray(gml.bhazard,(mod,0))
	matchvars 		= tokens(st_local("matchby"+strofreal(1)))
	ageyearindex 	= selectindex(matchvars:==("_age")),selectindex(matchvars:==("_year"))
	x				= st_data(.,matchvars,gml.touse[mod])[uhtred_get_index(gml),]
	
	if (gml.familys[mod]=="rp" | gml.familys[mod]=="logchazard") {
		asarray(gml.bhazard,(mod,4),uhtred_setup_expH(gml,exprate,x))
	}
	else {
	
		t 			= uhtred_util_depvar(gml)
		N 			= uhtred_get_nobs(gml)
		
		//need to add t to x[,1] and x[,2] -> age and year
		//then floor it to get the expected rate at the age, year previous

		Ngq 	= gml.chip
		chq2 	= J(N,Ngq,0)
		gq 		= uhtred_gq(Ngq,"legendre")
		if (gml.hasltrunc[1]) 	qp2 = (t[,1]:-t[,3]) :/ 2 :* J(N,1,gq[,1]') :+ (t[,1]:+t[,3]):/2
		else							qp2 = t[,1] :/ 2 :* J(N,1,gq[,1]') :+ t[,1]:/2
		//now get matched rate

		matchedrate = J(N,Ngq,.)
		for (q=1;q<=Ngq;q++) {
			xq 					= x
			xq[,ageyearindex] 	= floor(xq[,ageyearindex] :+ qp2[,q])
			for (i=1;i<=N;i++) {
				matchedrate[i,q] = exprate[match_row_index(exprate[|1,2\.,.|],xq),1]
			}
		}

		asarray(gml.bhazard,(mod,2),matchedrate)
	}
	
}

`RC' survsim_get_exprate(`RM' exprate, `RM' x, `RR' t)
{
	N 		= rows(x)
	x[,(1,2)] 	= floor(x[,(1,2)] :+ t)
	//now get matched rate
	matchedrate = J(N,1,.)
	for (i=1;i<=N;i++) {
		matchedrate[i] = exprate[selectindex(rowsum(x[i,]:==exprate[,(2,3)]):==2),1]
	}
	return(matchedrate)
}

`RR' survsim_get_exprate_row(`RR' t, `RS' i)
{
	Nq = cols(t)
	//read in expected rate file
	stata("preserve")
	stata("qui use exprate,clear")
	exprate = st_data(.,("_rate","_age","_year"))
	stata("restore")
	x = st_data(i,("_age","_year"))
	//now get matched rate
	matchedrate = J(1,Nq,.)
	for (q=1;q<=Nq;q++) {
		xq = x
		xq[,(1,2)] 	= floor(xq[,(1,2)] :+ t[1,q])
		matchedrate[1,q] = exprate[selectindex(rowsum(xq:==exprate[,(2,3)]):==2),1]
	}
	return(matchedrate)
}


end

