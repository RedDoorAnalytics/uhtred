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


/*
-> functions for family(rp, ...)
*/

mata:

`RM' uhtred_p_rp_logch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(uhtred_util_xtzb(gml,gml.myb,t))
}

`RM' uhtred_p_rp_ch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(exp(uhtred_util_xtzb(gml,gml.myb,t)))
}

`RM' uhtred_p_rp_s(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(exp(-exp(uhtred_util_xtzb(gml,gml.myb,t))))
}

`RM' uhtred_p_rp_logh(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)

	dXT		= asarray(gml.dXT,model)[index1,]
	dxb 		= dXT * brcs
	xtzb[index1,] :+ log(dxb) :- log(t)
	
	return(logh)
}

`RM' uhtred_p_rp_h(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(exp(uhtred_p_rp_logh(gml,t)))
}


`RM' uhtred_p_rp_f(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(1 :- uhtred_p_rp_s(gml,t))
}

`RM' uhtred_p_rp_dens(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = uhtred_util_timevar(gml)
	return(uhtred_p_rp_h(gml,t):*uhtred_p_rp_s(gml,t))
}

end
