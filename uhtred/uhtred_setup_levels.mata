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

void uhtred_setup_levels(`gml' gml)
{
	uhtred_get_cluster_varnames(gml)
	gml.E.levelvars = gml.levelvars[1,]
	gml.Nrelevels	= cols(gml.levelvars)
	gml.Nlevels 	= gml.Nrelevels + 1
	gml.Nobs	= J(gml.Nlevels,1,.)
}

void uhtred_get_cluster_varnames(`gml' gml)
{
	//get unique cluster spec's within []
	//-> first row contains level name
	//-> second row contains [] as well
	gml.levelvars = J(2,0,"")		
	//get any level specs within [] - includes pipes
	idspec = J(0,1,"")
	i = 1
	vars 	= uhtred_get_cmps(i)
	Ndv 	= cols(vars)
	for (j=1;j<=Ndv;j++) {
		for (k=1;k<=Ndv;k++) {
			dv = strtrim(vars[1,k])
			pos = 1	
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
				pos1 = strpos(dv2,"[") 
				pos2 = strpos(dv2,"]")
				if (pos1) {
					idspec = idspec\substr(dv2,pos1+1,pos2-pos1-1)
				}	
			}	
		}		
	}
	
	if (idspec!=J(0,1,"")) {
		//find unique ones to get # levels
		idspec = uniqrows(idspec)
		Nlevs  = rows(idspec)
		//get id vars, stripping out pipes
		idspecind = J(Nlevs,1,1)
		for (i=1;i<=Nlevs;i++) {
			ind 	= 1
			tempid 	= idspec[i]
			pos 	= strpos(tempid,">")
			while (pos) {
				idspecind[i] 	= idspecind[i] + 1
				tempid 		= substr(tempid,pos+1,.)
				pos 		= strpos(tempid,">")
			}
			pos 	= strpos(tempid,"<")
			while (pos) {
				idspecind[i] 	= idspecind[i] + 1
				tempid 		= substr(tempid,pos+1,.)
				pos 		= strpos(tempid,"<")
			}
			gml.levelvars = gml.levelvars,(tempid\("["+idspec[i]+"]"))
		}
		//make sure in order
		gml.levelvars = gml.levelvars[,idspecind']
	}
	
}
end
