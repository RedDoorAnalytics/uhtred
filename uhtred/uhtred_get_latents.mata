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

void uhtred_get_latents(`gml' gml)
{
	gml.elinfo = asarray_create("real",3)

	//notes; not effected by @'s

	//unique latents at each level
	gml.latlevs 	= asarray_create("real",1)
	gml.Nres 	= J(gml.Nlevels,1,0)
	cmps		= uhtred_get_cmps(1)
	Ncmps		= cols(cmps)
	
	for (i=1;i<=gml.Nlevels;i++) {
		lats = J(0,1,"")
		//go through everything for each level
		j = 1
		for (k=1;k<=Ncmps;k++) {
			dv  = strtrim(cmps[k])
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
				posid = strpos(dv2,gml.levelvars[2,i])
				if (posid) lats = lats\substr(dv2,1,posid-1)
			}	
		}
		
		if (rows(lats)) {
			//they could appear in multiple models and multiple times in same model
			lats 	= uniqrows(lats)			
			//now check that if a, b, c is specified, then its core random effect is also specified
			Nlats 	= rows(lats)
			len 	= strlen(lats)
			last 	= substr(lats,len,.)
			for (p=1;p<=Nlats;p++) {
				if (strtoreal(last[p])==.) {
					//strip off a, b to match
					core = substr(lats[p],1,len[p]-1)	
					if (!sum(core:==lats)) uhtred_error("random effect "+lats[p]+" requires "+core)
				}
			}
			
			//now go through elements, expanding any for design matrices
			latsup = J(0,2,"")
			depvars = uhtred_get_cmps(j)
			Ncmp 	= cols(depvars)
			for (c=1;c<=Ncmps;c++) {
				dv 	= strtrim(cmps[1,c])
				el  = 1
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
					posid = strpos(dv2,gml.levelvars[2,i])
					if (posid) {		
						posrb = strpos(dv2,"]")
						latsup = latsup\(substr(dv2,1,posid-1),substr(dv2,1,posrb))
						//store level
						//store re index(es)
						asarray(gml.elinfo,(j,c,el),i)
					}
					el++
				}	
			}
			
			uniqlats 	= uniqrows(latsup)
			gml.Nres[i]	= rows(uniqlats)
			
			//now go through elements, storing level and random effect index 
			for (c=1;c<=Ncmps;c++) {
				dv 	= strtrim(cmps[1,c])
				el  = 1
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
					posid = strpos(dv2,gml.levelvars[2,i])
					if (posid) {				
						ellatinfo = asarray_create("real",1)
						asarray(ellatinfo,1,i)		//store level
						reindex = J(0,1,.)
						lat = substr(dv2,1,posid-1)
						reindex = reindex\selectindex(lat:==uniqlats)
						asarray(ellatinfo,2,reindex)	//store re index(es)
						asarray(gml.elinfo,(j,c,el),ellatinfo)
					}
					el++
				}	
			}
			asarray(gml.latlevs,i,uniqlats)
		}
	}
}

end
