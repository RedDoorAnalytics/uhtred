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

void uhtred_setup_levelvars(`gml' gml)
{
	gml.Nobs	= J(gml.Nlevels,gml.Nmodels,.)
	
	//build new idvars
	if (gml.Nlevels>1) {
		if (!gml.predict) stata("qui sort "+invtokens(gml.levelvars[1,])+",stable")
		//now create temps of them to account for gaps in id indexes
		newidvars = J(1,0,"")
		for (i=1;i<gml.Nlevels;i++) {
			tempid = "_tempidindex"+strofreal(i)
			stata("tempvar "+tempid)
			newidvars = newidvars,st_local(tempid)
			stata("qui egen "+st_local(tempid)+
                                "= group("+invtokens(gml.levelvars[1,1..i])+
                                ") if "+gml.touse)
		}
		gml.levelvars[1,] = newidvars
	}
	
	//ob level, now full touse is built
	stata("tempvar coreindex")
	coreindex = st_local("coreindex")
	stata("qui egen "+coreindex+" = seq() if "+gml.touse)
		
	//need overall ids built with | model touses, which are then applied to the ids and stored
	gml.xbindex = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
		asarray(gml.xbindex,i,
                        st_data(.,coreindex,gml.modeltouses[1,i]))
		gml.Nobs[gml.Nlevels,i] = rows(asarray(gml.xbindex,i))	
	}

	gml.levelvars = gml.levelvars,(st_local("coreindex")\st_local("coreindex"))
	
	//id vars setup
	gml.panelindexes = asarray_create("real",2)
	if (gml.Nlevels>1) {
		for (j=1;j<=gml.Nmodels;j++) {
			ids = st_data(.,gml.levelvars[1,],gml.modeltouses[1,j])
			for (i=1;i<gml.Nlevels;i++) {
				psetup 		= panelsetup(uniqrows(ids[,1..(i+1)]),i)
				gml.Nobs[i,j]	= rows(psetup)
				asarray(gml.panelindexes,(i,j),psetup)
			}
		}
	}
}

end
