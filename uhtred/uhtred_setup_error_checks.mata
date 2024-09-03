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

void uhtred_setup_error_checks(`gml' gml)
{
	if (st_local("timevar"+strofreal(1))!="") {
		if (gml.familys=="lognormal") {
			uhtred_error("timevar() not currently supported with family(lognormal)")
		}
		if (gml.familys=="loglogistic") {
			uhtred_error("timevar() not currently supported with family(loglogistic)")
		}
	}
	
	if (gml.familys=="addhazard" & gml.Nlevels>1) {
		uhtred_error("family(addhazard) doesn't currently allow random effects")
	}
}

end
