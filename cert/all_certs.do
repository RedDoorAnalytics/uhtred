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

//survival- 1 level 
do ./cert/cert_rp_1level
do ./cert/cert_rp_1level_lt 

do ./cert/cert_er_rp_1level


//survival- 1 level - predictions
do ./cert/cert_rp_1level_pred
do ./cert/cert_er_rp_1level_pred


//survival - 1 level multiple timescales
*do ./cert/cert_rp_1level_mts
*do ./cert/cert_er_rp_1level_mts


//survival - 1 level competing risks 
*do ./cert/cert_rp_cr


//survival - 2 level weibull
do ./cert/cert_weibull_2level

do ./cert/cert_er_weibull_2level


//survival - 3 level weibull
do ./cert/cert_weibull_3level







//============================================================================//

di "certs complete - no errors"

//============================================================================//
