//============================================================================//
// cert. script for uhtred

//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

cscript uhtred

//============================================================================//

//survival

do ./cert/cert_rp_1level

//gf2

//multilevel survival

do ./cert/cert_rp_2level

//cr 


//============================================================================//

di "certs complete - no errors"

//============================================================================//
