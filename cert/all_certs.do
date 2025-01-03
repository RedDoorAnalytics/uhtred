//============================================================================//
// cert. script for uhtred

//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/uhtred"
adopath ++ "`drive'/uhtred"
adopath ++ "`drive'/uhtred/uhtred"

//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

cscript uhtred

//============================================================================//
// ob level models

do ./cert/cert_rp_1level_errors
do ./cert/cert_rp_1level
do ./cert/cert_rp_1level_ltrunc 
// do ./cert/cert_rp_1level_mts
// do ./cert/cert_er_rp_1level_mts

// competing risks
// do ./cert/cert_rp_cr

//============================================================================//
// 2-level models

do ./cert/cert_rp_2level
do ./cert/cert_rp_2level_gf12


//============================================================================//
// 3-level models

do ./cert/cert_rp_3level


//============================================================================//

di "certs complete - no errors"

//============================================================================//
