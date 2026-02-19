*! version 1.6.0  12feb2026

/*
History
v1.6.0: 12feb2026
- bug fix in predictions, timevar() was ignored when rebuilding splines
v1.5.1: 23may2025 
- bug fix in eta predictions
v1.5.0: 28mar2025
v1.4.0: 25mar2025
- fully adaptive quadrature
v1.0.0: 02sep2024 
- initial release
*/

program uhtred, eclass 
        version 17

        if replay() {
                if "`e(cmd)'" != "uhtred" {
                        error 301
                }
                Display `0'
                exit
        }
		
        if substr("`0'",1,1)!="(" {
                di as error "missing bracket"
                exit 198
        }
        
        tempname GML
        capture noisily Estimate `GML' `0'
        local rc = c(rc)
        if "`c(prefix)'"!="morgana" {
		di "here"
                capture n mata: uhtred_cleanup("`GML'")
                capture drop `GML'*
                capture mata: mata drop chazf hazf loglf
        }
        else ereturn local object "`GML'"
		
        if (`rc') {
		capture drop _uhtred_*
                capture mata: uhtred_cleanup("`GML'")
                capture drop `GML'*
                capture mata: mata drop chazf hazf loglf
                exit `rc'
        }
	capture drop _uhtred_*
        ereturn local cmdline `"uhtred `0'"'
end

program Estimate, eclass
        version 17
        gettoken GML : 0

        Fit `0'		//!! should leave behind diopts

	mata: uhtred_ereturn("`GML'")
		
        if "`c(prefix)'"!="morgana" {
                Display, `diopts'
        }
end

program Fit, eclass sortpreserve
        version 17
        gettoken GML 0 : 0

        tempname touse b
        uhtred_parse `GML', touse(`touse') : `0'
        if "`r(predict)'"!="" | "`c(prefix)'"=="morgana" | "`c(prefix)'"=="crossval" {
                exit
        }
        local hasopts   `"`r(hasopts)'"'
        local mltype    `"`r(mltype)'"'
        local mleval    `"`r(mleval)'"'
        local mlspec    `"`r(mlspec)'"'
        local mlopts    `"`r(mlopts)'"'
        local mlvce     `"`r(mlvce)'"'
        local mlwgt     `"`r(mlwgt)'"'
        local mlcns	`"`r(constr)'"'
        local mlinitcns	`"`r(initconstr)'"'
        local nolog     `"`r(nolog)'"'
        local mlprolog  `"`r(mlprolog)'"'
	local mftodrop  `r(mftodrop)'
        c_local diopts  `"`r(diopts)'"'
        local mlfrom    `"`r(mlfrom)'"'
	local initextra `"`r(initextra)'"'
        local mlzeros	`"`r(mlzeros)'"'
        local modellabels `"`r(modellabels)'"'
	local firth 	"`r(firth)'"
        if "`mlcns'" != "" {
                local cnsopt constraint(`mlcns')
        }
		
        if "`mlfrom'"!="" {
                matrix `b' = r(b)
                local mlinit init(`b')
        }
	if `"`initextra'"'!="" {
		local initextra init(`initextra')
	}
        
        di
        di as txt "Fitting full model:"

        cap n ml model `mltype' `mleval'              		///
                                `mlspec'                        ///
                                `mlwgt'                         ///
                                if `touse',                     ///
                                `mlopts'                        ///
                                `mlvce'                         ///
                                `mlprolog'                      ///
                                `cnsopt'                        ///
                                `mlinit'			///
				`initextra' 			///
                                collinear                       ///
                                maximize     	                ///
                                missing    		        ///
                                nopreserve   	                ///
                                search(off)                     ///
                                userinfo(`GML')	                ///
                                wald(0)
                                                                                                         
        if _rc>0 {
        
                if _rc==1400 & "`mlzeros'"=="" {
                
                        di ""
                        di as text "-> Starting values failed - trying zero vector"
                
                        ml model `mltype' 	`mleval'              		///
                                                `mlspec'                        ///
                                                `mlwgt'                         ///
                                                if `touse',                     ///
                                                `mlopts'                        ///
                                                `mlvce'                         ///
                                                `mlprolog'                      ///
                                                `cnsopt'                        ///
                                                collinear                       ///
                                                maximize     	                ///
                                                missing    		        ///
                                                nopreserve   	                ///
                                                search(off)                     ///
                                                userinfo(`GML')	                ///
                                                wald(0)   	                    
                }
                else {
                        exit _rc
                }
                
        }
        
	if "`firth'"!="" {
		di ""
		di as text "Refining variance-covariance matrix"
		ereturn repost V = firthV
		capture mat drop firthV
	}
	
        ereturn local predict   uhtred_p
        ereturn local from = "`mlfrom'"!=""
        ereturn local hasopts = `hasopts'
        ereturn local cmd uhtred
        ereturn local modellabels "`modellabels'"
        cap mata: mata drop `mftodrop'
        if "`mlcns'" != "" {
                cap constraint drop `mlcns'
        }
end

program Display
        syntax [,       noHeader        ///
                        noDVHeader      ///
                        noLegend        ///
                        notable         ///
                        EFORM           ///
                        *               ///
        ]

        _get_diopts diopts, `options'
        if e(estimates) == 0 {
                local coefl coeflegend selegend
                local coefl : list diopts & coefl
                if `"`coefl'"' == "" {
                        local diopts `diopts' coeflegend
                }
        }
		
        local Nrelevels = `e(Nlevels)'-1
        
        if "`eform'"!="" {
                local exp exp
        }
        
        local plus
        if "`e(Nres1)'"!="" {
                local plus plus
                local neq = `e(k_eq)'
                forval l=1/`Nrelevels' {
                        local neq = `neq' - `e(Nreparams`l')'
                }
                local neq neq(`neq')
        }

        _coef_table_header
        if "`exp'"!="" {
                local coeftitle exp(b)
        }
        else    local coeftitle 
        
	//get neq
	local neq = 0
	forvalues i=1/`e(Nmodels)' {
		if `: word `i' of `e(hasxb)'' {
			local `neq++'
		}
		if `: word `i' of `e(hastb)'' {
			local `neq++'
		}
		if `: word `i' of `e(haszb)'' {
			local neq = `neq' + `Nrelevels'
		}
		if `e(ndistap`i')'>0 {
			local neq = `neq' + `e(ndistap`i')'
		}
		if `e(nap`i')'>0 {
			local neq = `neq' + `e(nap`i')'
		}
	}

	if `Nrelevels'>0 {
		local plus plus
	}
	
	_coef_table, neq(`neq') `plus' nocnsreport coeftitle("`coeftitle'")
	
	
	
        //VCV display
        if "`e(Nres1)'"!="" {
                
                forval i=1/`Nrelevels' {
                        local lev : word `i' of `e(levelvars)'
                        _diparm __lab__ , label("`lev':") eqlabel
                        
                        forval j=1/`e(Nreparams`i')' {
                                local param : word `j' of `e(re_eqns`i')'
                                local scale : word `j' of `e(re_ivscale`i')'
                                local label : word `j' of `e(re_label`i')'
                                _diparm `param', `scale' label(`label')
                        }
                        if `i'<`Nrelevels' {
                                _diparm __sep__
                        }
                }
                
                _diparm __bot__
        }
        
        if "`e(penalty)'"!="" {
                di as text " Estimation: Maximum Penalised Likelihood"
                di as text "    Penalty: `e(penalty)' with lambda = `e(lambda)'"
        }
        
end

exit
