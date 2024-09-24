//============================================================================//
// cert. script for uhtred

//source paths
local drive C:\Users\Hannah Bower\Documents\GitHub
cd "`drive'\uhtred"
adopath ++ "`drive'\uhtred"
adopath ++ "`drive'\uhtred/uhtred"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

cscript uhtred

//============================================================================//

//survival- 1 level - predictions
do ./cert/cert_rp_1level
do ./cert/cert_rp_1level_lt 

do ./cert/cert_er_rp_1level

*do ./cert/cert_rp_1level_mts
*do ./cert/cert_er_rp_1level_mts

*do ./cert/cert_rp_cr

//survival- 1 level - predictions
do ./cert/cert_rp_1level_pred
do ./cert/cert_er_rp_1level_pred


//survival - 2 level 
do ./cert/cert_rp_2level_int_slope
do ./cert/cert_er_rp_2level


//MC old code
*do ./cert/cert_rp_2level_int_slope




//weib 3 level
do ./cert/cert_weib_3level_int

//gf2

//multilevel survival

*do ./cert/cert_rp_2level





//============================================================================//

di "certs complete - no errors"

//============================================================================//
