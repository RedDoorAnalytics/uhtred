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

void uhtred_error(`SS' text)
{
	errprintf(text+"\n")
	exit(1986)
} 

/*
- removes @ and subsequent text from a string
- if not there it's left unchanged
*/
void _uhtred_remove_at(`SS' dv)
{
	if (strpos(dv,"@")) dv = substr(dv,1,strpos(dv,"@")-1)
}

/*
Multivariate normal probability density function
*/

`RM' uhtred_dmvnorm(`RM' X, `RM' V)
{
	`RM' invsigma, Vvecs, Vvals, q2, denom
	
	symeigensystem(V,Vvecs,Vvals)
	invsigma = Vvecs * (Vvecs' :/ Vvals')
	q2 = 0.5 :* rowsum((X*invsigma):*X)
	d1 = -0.5 :* (rows(V)*log(2*pi()) :+ sum(log(Vvals)) )
	return(exp(d1:-q2))
}

`RM' uhtred_lndmvnorm(`RM' X, `RM' V)
{
	`RM' invsigma, Vvecs, Vvals, q2, denom
	
	symeigensystem(V,Vvecs,Vvals)
	invsigma = Vvecs * (Vvecs' :/ Vvals')
	q2 = 0.5 :* quadrowsum((X*invsigma):*X)
	d1 = -0.5 :* (rows(V)*log(2*pi()) :+ quadsum(log(Vvals)) )
	return(d1:-q2)
}


`RM' uhtred_glegendre(real scalar n)
{
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	alpha = 0
			
	muzero = 2
	a = J(1,n,0)
	b = i1:/sqrt(4 :* i1:^2 :- 1)

	A= diag(a)
	for(j=1;j<=n-1;j++){
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	return(nodes,weights)
}

`RM' _uhtred_get_hessian_indices(`RS' Nb, | `RS' Nb2)
{
	xindex = J(2,0,.)
	
	if (args()==1) {
		for (i=1; i<=Nb; i++) {
			refind = 1
			while (refind<=i) {
				xindex = xindex,(i\refind)
				refind++
			}
		}
	}
	else {
		for (i=1; i<=Nb; i++) {
			for (j=1;j<=Nb2;j++) {
				xindex = xindex,(i\j)
			}
		}
	}
	return(xindex)
}


`RM' uhtred_xz_rcs(`gml' gml, `RS' c, `RS' el, `RS' deriv, | `RC' t, `RC' t0)
{
	hast		= args()==5
	hast0		= args()==6
	rcsinfo 	= asarray(gml.elinfo,(1,c,el))

	istvars 	= uhtred_util_istimevar(gml,c,el) 	//asarray(rcsinfo,1)
	index		= uhtred_util_index(gml)

	if (hast) {
		if 	(istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = uhtred_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else if (hast0) {
		if 	(istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = t0
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else {
		if 	(istvars[1]) 	rcsvar = uhtred_util_timevar(gml)
		else if (istvars[2]) 	rcsvar = uhtred_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}

	hasoffset = asarray(rcsinfo,7) //offset

	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :+ t0
			else 		rcsvar = rcsvar :+ uhtred_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,8)[index]
	}
	hasoffset 	= asarray(rcsinfo,9) //moffset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :- t0
			else 		rcsvar = rcsvar :- uhtred_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,10)[index]
	}

	islog 	= asarray(rcsinfo,4)
	knots 	= strtoreal(tokens(asarray(rcsinfo,3)))
	isorth 	= asarray(rcsinfo,5)

	if (isorth) {
		rmat = asarray(rcsinfo,6)
		if (islog)  result = uhtred_rcs(log(rcsvar),knots,deriv,rmat)
		else 	    result = uhtred_rcs(rcsvar,knots,deriv,rmat)		
	}
	else {
		if (islog)  result = uhtred_rcs(log(rcsvar),knots,deriv)
		else 	    result = uhtred_rcs(rcsvar,knots,deriv)
	}

	if (deriv & islog)  return(result:/rcsvar)
	else 		    return(result)
	
}

/*
Expand a matrix into all permutations
*/

`RM' uhtred_expand_matrix(`RM' x,| `RS' weights )
{
	`RS' nrows, ncols, nexp, index
	`RM' newx
	
	if (rows(x)==1) return(x)
	else {
		nrows = rows(x)
		ncols = cols(x)	
		nexp  = ncols^(nrows-1)
		newx = J(1,nexp,x)
		index = 1
		for (i=1; i<=nrows-1; i++) {
			reps = ncols^(nrows-i)
			xrep = J(reps,1,x[i,])'
			xrep = rowshape(xrep,1)
			newx[i,] = J(1,index,xrep)
			index = index * ncols
		}
		if (args()>1) {
			for (i=2; i<=nrows; i++) newx[1,] = newx[1,] :* newx[i,]
			return(newx[1,])
		}
		else return(newx)
	}	
}

/*
Outer product by column
*/

`RM' uhtred_outerprod_by_col(`RM' x)
{
	`RS' ncols, nrows
	`RM' newx
	nrows = rows(x)
	ncols = cols(x)
	newx = J(nrows,ncols^2,.)
	for (i=1;i<=nrows;i++) newx[i,] = x[i,] # x[i,] //rowshape(cross(x[i,],x[i,]),1)
	return(newx)
}

/*
error check for se(blups) when 0
*/

void uhtred_seblup_check(`RS' x)
{
	if (x==0 | x==.) {
		errprintf("BLUP calculation failed in adaptive quadrature algorithm\n")
		errprintf("Try increasing gh()\n")
		exit(1986)
	}
}

`RC' uhtred_get_adpanelindex(`gml' gml, `RS' lev)
{
	index = asarray(gml.adpanelindexes,(lev,1))
	if (gml.survind==0) return(index)
	else 		    return(index[merlin_get_surv_index(gml)])
}

end
