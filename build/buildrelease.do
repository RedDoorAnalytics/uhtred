// build new version of uhtred

local drive /Users/Michael
cd "`drive'/My Drive/software/uhtred/"

local includemata       = 0

//=======================================================================================================================//

//build new release -> current version up is 1.1.0
local newversion 1_1_0
if `includemata' {
        local newversion `newversion'_mata
}
cap mkdir ./release/version_`newversion'
local fdir ./release/version_`newversion'/

//=======================================================================================================================//

//pkg files
// copy ./build/uhtred_details.txt `fdir', replace
	
//=======================================================================================================================//

//uhtred

	copy ./uhtred/uhtred.ado `fdir', replace
	copy ./uhtred/uhtred_parse.ado `fdir', replace
	copy ./uhtred/uhtred_p.ado `fdir', replace
	copy ./uhtred/uhtred_rcs.ado `fdir', replace
	copy ./uhtred/stuhtred.ado `fdir', replace	

	//help files
	copy ./uhtred/uhtred.sthlp `fdir', replace
	copy ./uhtred/uhtred_estimation.sthlp `fdir', replace
	copy ./uhtred/uhtred_model_options.sthlp `fdir', replace
	copy ./uhtred/uhtred_models.sthlp `fdir', replace
	copy ./uhtred/uhtred_reporting.sthlp `fdir', replace
	copy ./uhtred/uhtred_postestimation.sthlp `fdir', replace
	
//mlib
cap erase `fdir'luhtred.mlib
copy ./luhtred.mlib `fdir', replace
