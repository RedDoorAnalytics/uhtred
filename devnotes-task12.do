/*

//============================================================================//
// notes



//============================================================================//
// to do

- tidy up touses
- random effects on time-dep effects
- help file
- multiple time interactions for rp -> chain rule
- predictions

- swap to views if possible 
  -> currently st_data() for gf1/2 as variable order gets changed (don't know 
     why)
  -> test doing views each time within the evaluator
- rp
  -> gf0 ok
  -> gf1 ok
  -> gf2 ok
  - ltruncation ok
  - interval censoring
	-> gf1 done, removing gf2 for now
  - bhazard() ok
      -> added gf2

- weibull
- loghazard

- stmixed
- morgana

//============================================================================//
// done

- sync CH quad time design matrices being stored
- new subroutine for extracting design matrix on random effects
- fix spline names with model in index
- block any rcs() with random effect, only allow varname or c.varname
- error check for ## with a random effect
- error check for i. interacted with random effect
- strip out left truncation with Nlevels>1 
- strip out predictms
- store CH quad points
- store eqn indices for each xb, tb, zb, vcvs
- constraint parsing with @
- gridsearch synced
- adaptive quad
- remove random option
- fix display subroutine
- parse rcs() and replace with indexes vars
- allow factor variable notation
- sync starting values
- sync time-dependent matching
  -> all moved to a second "time" equation
- check_clp improved and tidied
  -> now allows factor variables for anything
- move constant to xb equation


//============================================================================//
// next steps

- stmerlin
- predictms
- cox
- sync fp()
- sync bs()

*/
