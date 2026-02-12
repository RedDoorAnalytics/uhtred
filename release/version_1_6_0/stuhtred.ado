*! version 1.0.0  06sep2024 MJC

/*
History
06sep2024: version 1.0.0
*/

program stuhtred, eclass sortpreserve properties(st)
	version 15.1
	if replay() {
		if (`"`e(cmd2)'"'!="stmixed") error 301
		merlin `0'
	}
	else Estimate `0'
end

program Estimate, eclass
	st_is 2 analysis
	
	cap which uhtred
	if _rc {
		di as error "stuhtred requires uhtred"
		exit 198
	}
	
	// model opts
		local modelopts `"Distribution(string)"'
		local modelopts `"`modelopts' BHazard(passthru)"'
		local modelopts `"`modelopts' COVariance(passthru)"'
		local modelopts `"`modelopts' DF(passthru)"'
		local modelopts `"`modelopts' KNOTS(passthru)"'
		local modelopts `"`modelopts' TVC(varlist)"'
		local modelopts `"`modelopts' DFTvc(string)"'
		local modelopts `"`modelopts' KNOTSTvc(string)"'
		local modelopts `"`modelopts' NOORTHog"'
		
		
		local modelopts `"`modelopts' LLFunction(passthru)"'
		local modelopts `"`modelopts' LOGHFunction(passthru)"'
		local modelopts `"`modelopts' HFunction(passthru)"'
		local modelopts `"`modelopts' CHFunction(passthru)"'
		local modelopts `"`modelopts' NAP(passthru)"'
		
	//Max opts
		local modelopts `"`modelopts' INTPoints(passthru)"'
		local modelopts `"`modelopts' INTMethod(passthru)"'	
		local modelopts `"`modelopts' FROM(passthru)"'				
		local modelopts `"`modelopts' RESTARTVALues(passthru)"'
		local modelopts `"`modelopts' APSTARTValues(passthru)"'
		local modelopts `"`modelopts' ZEROS"'
		local modelopts `"`modelopts' ADAPTopts(passthru)"'
		
	// mlopts
        local mlopts `"NONRTOLerance NRTOLerance(string)"'
        local mlopts `"`mlopts' TRace GRADient HESSian showstep"'
        local mlopts `"`mlopts' TECHnique(string) SHOWNRtolerance"'
        local mlopts `"`mlopts' ITERate(string) TOLerance(string)"'
        local mlopts `"`mlopts' LTOLerance(string) GTOLerance(string)"'
        local mlopts `"`mlopts' DIFficult dots depsilon(passthru) showh"'
	
	//Display opts
		local modelopts `"`modelopts' Level(passthru)"'
		local modelopts `"`modelopts' NOLOG"'
		local modelopts `"`modelopts' SHOWMERLIN"'
			
	//Undocumented for use by predict and parsing check
		local modelopts `"`modelopts' debug"'
		
		local globallow `"`modelopts' `mlopts'"'
		local cmdline `0'
		
	//parse
	
        _parse expand cmd glob : 0 , common(`globallow')
	
        if `cmd_n'==1 {
                if strpos(`"`cmd_1'"',":") {
                        local cmd_2 `cmd_1'
                        local cmd_1 
                        local cmd_n = 2
                }
                else {
                        di as error "Missing level variable"
                        exit 198
                }
        }
        
        forvalues k = 1/`cmd_n' {              
			local cmds `"`cmds' `"`cmd_`k''"'"'
        }
        _mixed_parseifin stuhtred `=`cmd_n'+1' `cmds' `"`glob_if' `glob_in'"'

        local ifin if _st==1 
        if "`glob_if'"!="" local ifin `ifin' & `glob_if'
        if "`glob_in'"!="" local ifin `ifin' & `glob_in'
		
	//global options are in glob_op 
	//parse model options

        local 0 `", `glob_op'"'
        syntax [ , `modelopts' *]

        //family
        
        local family "`distribution'"
        if "rp"=="`family'" {
                local family rp
        }
        else {
                di as error "distribution(`distribution') not supported"
                exit 198
        }
                
        //delayed entry
        qui su _t0 `ifin', meanonly
        if (`r(min)'==0 & `r(max)'>0) {
                di as error "{p}delayed entry detected for a subset of " ///
                        "observations; must be all or none{p_end}"
                exit 198
        }
        
        if `r(max)'>0 & `r(min)'>0 {
                local ltruncated ltruncated(_t0)
                di as text "note; a delayed entry model is being fitted"
        }
        
        //parse mlopts
        mlopts mlopts , `options'

        //final family
        local family family(`family', failure(_d) ///
              `userfunc' `bhazard' `ltruncated' `df' `knots' `noorthog')		
        
        //parse mlopts
        mlopts mlopts , `options'		
							
	//===================================================================//
	// build complex predictor
		
        local Mind = 1
        
        local merlincp `cmd_1)'
        
        //tvcs
        if "`tvc'"!="" {
                
                if "`dftvc'"=="" & "`knotstvc'"=="" {
                        di as error "dftvc() or knotstvc() required"
                        exit 198
                }
                
                if "`dftvc'"!="" local tvck df(`dftvc')
                else 		local tvck knots(`knotstvc')
                
                foreach var in `tvc' {
                        local merlincp `merlincp' c.`var'#rcs(_t, `tvck' log orthog)
                }
        
        }
                
        //random effects
        forval lev=2/`cmd_n' {
        
                //extract level var
                gettoken levvar rest : cmd_`lev', parse(":")
                gettoken colon rest : rest 		, parse(":")
                if `"`colon'"' == `":"' {
                        local lev = trim("`levvar'")
                        if "`levvar'" != "_all" {
                                confirm variable `levvar'
                        }
                        else {
                                di as error "_all not supported"
                                exit 198
                        }
                }
                else {
                        di as err "{p 0 6 2}invalid random-effects "
                        di as err "specification; perhaps you omitted the "
                        di as err "colon after the level "
                        di as err "variable{p_end}"
                        exit 198
                }
        
                local levelvars `levelvars' `levvar'
                local lev`lev' = subinstr("`levelvars'"," ",">",.)
                
                //random effects
                local 0 `rest'
                syntax [varlist(numeric default=none)] [, NOConstant]

                if ("`noconstant'"!="" & "`varlist'"=="") {
                        di as error "No random effects have been specified at level `levvar':"
                        exit 198
                }
                if "`noconstant'"=="" {
                        local merlincp `merlincp' M`Mind'[`lev`lev'']@1
                        di as text "Random effect M`Mind': Intercept at level `levvar'"
                        local Mind = `Mind' + 1
                }
                
                local Nvars : word count of `varlist'
                foreach var in `varlist' {
                        local merlincp `merlincp' `var'#M`Mind'[`lev`lev'']@1
                        di as text "Random effect M`Mind': `var' at level `levvar'"
                        local Mind = `Mind' + 1
                }
                
        }

        //===================================================================//
		// merlin
			
        if "`debug'"!="" {
                di "`merlincp'"
                di "`family'"
        }
        
        if "`showmerlin'"!="" {
                di "uhtred (_t `merlincp', `family' `timevar') `ifin' , `mlopts' `covariance' `intmethod' `intpoints' `from' `debug' `restartvalues' `nolog' `adaptopts' `level'"
        }

        uhtred 	(_t	`merlincp', `family' `timevar')		///
                         `ifin' ,	 			///
                         excalibur				///
                        `mlopts'				///
                        `covariance' `intmethod' 		///
                        `intpoints' `from' `debug'		///
                        `restartvalues'	`nolog'	`adaptopts'	///
                        `level' `zeros' `restartvalues' 	///
                        `apstartvalues'			

        ereturn local cmd2 stuhtred
        ereturn local cmdline2 stuhtred `cmdline'

end
