version 17

capture erase luhtred.mlib
mata: mata set matastrict off
qui {
	do "./uhtred/uhtred_setup.mata"
	do "./uhtred/uhtred_setup_check_clp.mata"
	do "./uhtred/uhtred_setup_levels.mata"
	do "./uhtred/uhtred_setup_error_checks.mata"
	do "./uhtred/uhtred_setup_touses.mata"
	do "./uhtred/uhtred_setup_levelvars.mata"
	do "./uhtred/uhtred_get_latents.mata"
	do "./uhtred/uhtred_rcs.mata"
	do "./uhtred/uhtred_setup_build_clp.mata"
	do "./uhtred/uhtred_setup_survival.mata"
	do "./uhtred/uhtred_setup_exprate.mata"
	do "./uhtred/uhtred_setup_vcv.mata"
	do "./uhtred/uhtred_starting_values.mata"
	do "./uhtred/uhtred_setup_build_clpq.mata"
	
	//functions
	do "./uhtred/uhtred_xb.mata"
	do "./uhtred/uhtred_gf.mata"
	do "./uhtred/uhtred_prolog.mata"
	do "./uhtred/uhtred_utils.mata"
	do "./uhtred/uhtred_functions.mata"
	
	do "./uhtred/uhtred_ereturn.mata"
	
	//distributions
	do "./uhtred/uhtred_logl_weibull.mata"
	do "./uhtred/uhtred_logl_loghazard.mata"
	do "./uhtred/uhtred_logl_rp.mata"
	
	//predictions
	do "./uhtred/uhtred_predict.mata"
	do "./uhtred/uhtred_p_rp.mata"
	
	mata: mata mlib create luhtred, dir(.)
	mata: mata mlib add    luhtred *(), dir(.)
	mata: mata d *()  
	mata mata clear
	mata mata mlib index
}
