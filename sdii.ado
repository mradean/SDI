
program sdii
	version 14.2

	syntax	[varlist(default=none min=0 max=2 numeric)] [if] [in] [fw aw pw iw] ///
			[, by(str) CORRelation(numlist max=1 >=-1 <=1) DIFFerence DISTribution Level(cilevel) median Mvalue(real 0) ///
			   NOLegend UNEqual UNPaired percentile PERCENTILE_s(str) pctile PCTILE_s(str) PRECision(integer 1) power POWER_n(numlist max=1 >0 <1) ///
			   REVerse sample1(numlist min=2 max=3) s1(numlist min=2 max=3) sample2(numlist min=2 max=3) s2(numlist min=2 max=3) var2(str) Welch ]
	
	cap noi {
		qui {
			local cmdl `"sdii `*'"'

			cap mata: st_local("sname",strofreal(max(substr(direxternal("*"),1,4):=="__s_"))) 
			if (!_rc & "`sname'"!="0") {
				di as err "There are Mata matrixes that start with '__s_' that you must either drop or rename before proceeding." ///
				_n "Use either {it:mata drop} or {it:mata rename}; see {ul:[M-3] Commands for controlling Mata}." ///
				_n "(To drop all, you can write {cmd:mata: mata drop __s_*} in Stata's command window.)"
				exit 198
			}

			local unpair = "`unpaired'"!=""
			if (`unpair' & "`correlation'"!="") {
				di as err "option {bf:correlation()} not allowed with unpaired data"
				exit 198
			}
			
			if (!`unpair' & ("`by'"!="" | "`welch'"!="")) local unpair 1
			local options `unequal'
			if ("`welch'"!="") {
				if ("`varlist'"=="") local options "unequal welch"
				else  local options "unequal unpaired welch"
				if ("`unequal'"=="") local unequal "unequal"
			}
			else if ("`unequal'"=="") local unequal "equal"
					
			if ("`power'"!="") local power .8
			else if ("`power_n'"!="") local power `power_n'

			tempname rho
			if ("`correlation'"=="") scalar `rho' = 0
			else {
				cap confirm number `correlation'
				if (_rc) {
					di as err "invalid option {bf:correlation()};" _n " {it:numlist} not allowed"
					exit(198)
				}
				else scalar `rho' = `correlation'
			}

			local inf `in' `if'

			if ("`weight'"!="") {
				if ("`varlist'"=="") {
					di as err "weights require {it:var1 var2}"
					exit 198
				}
				local wght "[`weight' `exp']"
			}

			
			if ("`varlist'"=="") {
				if ("`distribution'"!="") {
					di as err "option {bf:distribution} only allowed with {it:var1 var2}"
					exit 198
				}
				if ("`sample1'"!="" & "`s1'"!="") { {
					di as err "option {bf:s1()} not allowed"
					exit 198
				}
				else if ("`sample1'"=="" & "`s1'"!="") local sample1 `s1'
				if ("`sample2'"!="" & "`s2'"!="") { {
					di as err "option {bf:s2()} not allowed"
					exit 198
				}
				else if ("`sample2'"=="" & "`s2'"!="") local sample2 `s2'
				if ("`sample1'"=="" & "`sample2'"=="") {
					di as err "one of ({it:var1 var2}) or ({bf:sample1()} {bf:sample2()}) is required"
					exit 198
				}
				else if ("`sample1'"=="") {
					di as err "{bf:sample1()} is missing"
					exit 198
				}
				else if ("`sample2'"=="") {
					di as err "{bf:sample2()} is missing"
					exit 198
				}
				
				if ("`reverse'"!="") {
					local switch `sample1'
					local sample1 `sample2'
					local sample2 `switch'
					local mvalue `=-`mvalue''
					local switch
				}				
			}
			else {
				if ("`sample1'"!="" | "`sample2'"!="") {
					di as err "only one of {it:var1 var2} or {bf:sample1() sample2()} is allowed"
					exit 198
				}
				if ("`correlation'"!="") {
					di as err "option {bf:correlation()} not allowed with {it:var1 var2}"
					exit 198
				}
			}
			
			if (`precision'<0 | `precision'>6) {
				di as err "{bf:precision()} invalid -- invalid number, outside of allowed range"
				exit 198
			}
			else local form = "%`=`precision'+4'.`precision'f"

			if ("`percentile_s'"!="" | "`pctile_s'"!="") {
				if ("`varlist'"=="") {
					di as err "option {bf:percentile} only allowed with {it:var1 var2}"
					exit 198
				}
				if (`unpair') {
					di as err "option {bf:percentile} not allowed with unpaired data"
					exit 198
				}

				local ps `percentile_s' `pctile_s'
				local ps : list uniq ps
				if (!inlist("`ps'","alt","altd","altde","altdef")) {
					di as err "invalid {bf:percentile()} suboption"
					exit 198
				}
				else {
					if ("`percentile'"!="" | "`pctile'"!="") {
						di as err "only one of {bf:percentile} or {bf:percentile(altdef)} is allowed"
						exit 198
					}
					if ("`wght'"!="") {
						di as err "option {bf:altdef} not allowed with weights"
						exit 198
					}
					
					local alt "altdef"
				}
				local pct 1
			}
			else if ("`percentile'"!="" | "`pctile'"!="") {
				if ("`varlist'"=="") {
					di as err "option {bf:percentile} only allowed with {it:var1 var2}"
					exit 198
				}
				if (`unpair') {
					di as err "option {bf:percentile} not allowed with unpaired data"
					exit 198
				}
				local pct 2
			}
			else if ("`median'"!="") {
				di as err "option {bf:median} only allowed with {bf:percentile}"
				exit 198
			}

			if (_N<7) {
				local N1 `=_N+1'
				set obs 7
			}
			tempvar diff quant statis std_err zt_d lo_sdi hi_sdi levtype df
			gen `df' = .
			tempname lbl
			label def `lbl' 1 "(1)" 2 "(2)"
			local nlab 2

			if ("`varlist'"!="") {
				tokenize "`varlist'"
				if ("`2'"=="" & "`var2'"=="" & "`by'"=="") {
					di as err "'var2' is required"
					exit 198
				}
				else if (("`2'"!="" | "`var2'"!="") & "`by'"!="") {
					di as err "option {bf:by()} not allowed with 'var2'"
					exit 198
				}
				if ("`2'"!="" & "`var2'"!="") {
					di as err "option {bf:var2} not allowed;" _n "'var2' already defined (ie, {it:`2'})"
					exit 198
				}
				else if ("`2'"=="" & "`by'"=="") {
					tempvar 2
					__s_var2 `2' = `var2', ifn(`inf')
				}
				
				if ("`reverse'"!="" & "`by'"=="") {
					local switch `1'
					local 1 `2'
					local 2 `switch'
					local mvalue `=-`mvalue''
					local switch
				}				

				if ("`pct'"=="") {
					if ("`distribution'"=="") {
						if ("`by'"=="") {
							if (!`unpair') {
								corr `1' `2' `inf'
								scalar `rho' = r(rho)
							}
							ttest `1'==`2' `inf', `unpaired' `options'
							mata: __s_sd = st_numscalar("r(se)")*sqrt(st_numscalar("r(N_1)"))
						}
						else {
							if ("`reverse'"=="") {
								cap sum `by', meanonly
								local bymi `r(min)'
								local bymx `r(max)'
								ttest `1' `inf', by(`by') `options'
								if (r(N_1)==r(N_2)) local equnp 1
							}
							else {
								local mvalue `=-`mvalue''
								tempvar byrev
								gen `byrev' = `by'
								cap sum `byrev', meanonly
								local bymi `r(max)'
								local bymx `r(min)'
								recode `byrev' (`r(min)'=`r(max)') (`r(max)'=`r(min)')
								ttest `1' `inf', by(`byrev') `options'
								if (r(N_1)==r(N_2)) local equnp 1
							}
						}
					}
					else {
						if (`unpair') {
							di as err "option {bf:distribution} only allowed with paired tests"								
							exit 198
						}
						if ("`by'"!="") {
							di as err "option {bf:by()}, which assumes unpaired data, not allowed with {bf:distribution}"
							exit 198
						}

						sum `1' `inf'
						local sample1 `r(mean)' `r(sd)'
						sum `2' `inf'
						local sample2 `r(mean)' `r(sd)'
						corr `1' `2' `inf'
						scalar `rho' = r(rho)
						
						local varlist
					}					
				}
				else {
					if ("`if'"!="") local inf `inf' & `1'!=. & `2'!=.
					else local inf `inf' if `1'!=. & `2'!=.

					sum `1' `inf' `wght'
					replace `df' = r(N)-1 in 1/2
					mata: __s_s1 = st_numscalar("r(sd)")
					gen double `statis' = r(mean) in 1
					gen double `std_err' = r(sd)/sqrt(r(N)) in 1
					sum `2' `inf' `wght'
					mata: __s_s2 = st_numscalar("r(sd)")
					replace `statis' = r(mean) in 2
					replace `std_err' = r(sd)/sqrt(r(N)) in 2
					if (`statis'[1]<`statis'[2]) {
						local switch `1'
						local 1 `2'
						local 2 `switch'
						local mvalue `=-`mvalue''					
					}
				
					gen `diff' = `1'-`2'
					local dlb = .5*(100-`level')
					_pctile `diff' `inf' `wght', p(`dlb' `=100-`dlb'') `alt'
					gen double `lo_sdi' = r(r1) in 3
					gen double `hi_sdi' = r(r2) in 3
					mata: __s_rhs = st_numscalar("r(r1)")-strtoreal(st_local("mvalue"))

					numlist ".5(.5)50"
					_pctile `1' `inf' `wght', p(`r(numlist)') `alt'
					mata: __s_pct = J(100,2,.)
					mata: for (i=1; i<=100; i++) __s_pct[i,1] = st_numscalar("r(r"+strofreal(i)+")")		
					numlist "50(.5)99.5"
					_pctile `2' `inf' `wght', p(`r(numlist)') `alt'
					mata: for (i=1; i<=100; i++) __s_pct[i,2] = st_numscalar("r(r"+strofreal(101-i)+")")
					sum `1' `inf' `wght', meanonly
					mata: __s_pct = (st_numscalar("r(min)"),.) \ __s_pct
					sum `2' `inf' `wght', meanonly
					mata: __s_pct[1,2] = st_numscalar("r(max)")

					gen str14 `levtype' = ""
					tempname kp
					mata: __s_k=__s_lev=__s_lhs = .
					mata: __s_centile(__s_k,__s_lev,__s_lhs,__s_pct,__s_rhs,"`1'","`2'","`kp'","`levtype'")
					if ("`switch'"!="") {
						local switch `1'
						local 1 `2'
						local 2 `switch'
						local mvalue `=-`mvalue''					
					}

					if (inlist(`kp',0,50)) {
						if ("`median'"!="") {
							replace `std_err' = `std_err'*1.2533 in 1/2
							_pctile `1' `inf' `wght', p(50) `alt'
							replace `statis' = r(r1) in 1
							_pctile `2' `inf' `wght', p(50) `alt'
							replace `statis' = r(r1) in 2
						}
						
						if (`kp'==0) {
							sum `1' `inf' `wght', meanonly
							replace `lo_sdi' = r(min) in 1
							replace `hi_sdi' = r(max) in 1
							sum `2' `inf' `wght', meanonly
							replace `lo_sdi' = r(min) in 2
							replace `hi_sdi' = r(max) in 2
						}
						else {
							replace `lo_sdi' = `statis' in 1/2
							replace `hi_sdi' = `statis' in 1/2 
						}						
					}
					else {
						if ("`median'"=="") {
							_pctile `1' `inf' `wght', p(`=`kp'' `=100-`kp'') `alt'
							replace `lo_sdi' = r(r1) in 1
							replace `hi_sdi' = r(r2) in 1
							_pctile `2' `inf' `wght', p(`=`kp'' `=100-`kp'') `alt'
							replace `lo_sdi' = r(r1) in 2
							replace `hi_sdi' = r(r2) in 2
						}
						else {
							replace `std_err' = `std_err'*1.2533 in 1/2
							_pctile `1' `inf' `wght', p(`=`kp'' 50 `=100-`kp'') `alt'
							replace `statis' = r(r2) in 1
							replace `lo_sdi' = r(r1) in 1
							replace `hi_sdi' = r(r3) in 1
							_pctile `2' `inf' `wght', p(`=`kp'' 50 `=100-`kp'') `alt'
							replace `statis' = r(r2) in 2
							replace `lo_sdi' = r(r1) in 2
							replace `hi_sdi' = r(r3) in 2
						}
					}

					if ("`difference'"!="") {
						sum `diff' `inf' `wght'
						mata: __s_sd = st_numscalar("r(sd)")
						replace `std_err' = r(sd)/sqrt(r(N)) in 3
						if ("`switch'"=="") {
							if ("`median'"=="") replace `statis' = r(mean) in 3
							else {
								replace `std_err' = `std_err'*1.2533 in 3
								_pctile `diff' `inf' `wght', p(50) `alt'
								replace `statis' = r(r1) in 3
							}
						}
						else {
							if ("`median'"=="") replace `statis' = -r(mean) in 3
							else {
								replace `std_err' = `std_err'*1.2533 in 3
								_pctile `diff' `inf' `wght', p(50) `alt'
								replace `statis' = -r(r1) in 3
							}
							tempname lo
							scalar `lo' = -`lo_sdi'[3]
							replace `lo_sdi' = -`hi_sdi' in 3
							replace `hi_sdi' = `lo' in 3
						}
						
						if (`mvalue' & `statis'[3] & (sign(`mvalue')!=sign(`statis'[3]))) local mtext = sign(`mvalue')
						replace `levtype' = "`level'  CI" in 3
						mata: __s_lev = __s_lev\strtoreal(st_local("level"))
						label def `lbl' 3 "(1-2)", add
						local nlab 3
						local dfci " & CI"
					}
					else if (`mvalue') {
						if ("`median'"=="") {
							sum `diff' `inf' `wght', meanonly
							local rval `r(mean)'
						}
						else {
							_pctile `diff' `inf' `wght', p(50) `alt'
							local rval `r(r1)'
						}
						if (`rval' & sign(`mvalue')!=sign(`rval')) local mtext = sign(`mvalue')
					}
					
					putmata __s_matn=(`statis' `std_err' . . . . `lo_sdi' `hi_sdi') in 1/`nlab'
				}
			}
	
			if ("`varlist'"=="") {
				if ("`inf'"!="") {
					di as err "{bf:in} and {bf:if} qualifiers require {it:varlist}"
					exit 198
				}				
				if (`: word count `sample1''!=`: word count `sample2'') {
					di as err "{bf:sample1()} and {bf:sample2()} must have the same number of arguments"
					exit 198
				}
				else {
					if (`: word count `sample1''==2) {
						if (`unpair') {
							di as err "option {bf:unpaired} not allowed with distributions"
							exit 198			
						}
						local nosamp 1
					}
					else if (`:word 1 of `sample1''==`:word 1 of `sample2'') {
						if (`unpair') local equnp 1 
					}
					else if (!`unpair') local unpair 1

					if ("`nosamp'"!="" & "`power'"!="") {
						di as err "option {bf:power()} not allowed with distributions"
						exit 198			
					}
				}
			}

			if ("`pct'"=="") {
				tempname int prec sd zts ztd
				scalar `prec' = 1/10^`precision'

				if ("`nosamp'"!="") {
					ttesti 100 `sample1' 100 `sample2', `options' 
					scalar `zts' = invnormal(1-.5*(1-(`level'/100)))
					scalar `sd' = sqrt(r(sd_1)^2+r(sd_2)^2-2*`rho'*r(sd_1)*r(sd_2))
					scalar `ztd' = (`zts'*`sd'+`mvalue')/(r(sd_1)+r(sd_2))		
					if (`ztd'<0) scalar `ztd' = 0
					gen double `zt_d' = `ztd' in 1/2

					gen	double `std_err' = r(sd_1) in 1
					replace `std_err' = r(sd_2) in 2
					mata: __s_lev = J(2,1,(1-2*(1-(normal(st_numscalar("`ztd'")))))*100)
					scalar `int' = (1-2*(1-(normal(`ztd'))))*100/`prec'
					if (`int'==. & `ztd'!=.) {
						mata: _editmissing(__s_lev,100)
						scalar `int' = 100
					}
					if (int(`int')==100) gen `levtype' = "100 SDI" in 1/2
					else if (`ztd'==0) gen `levtype' = "0 SDI" in 1/2
					else if (int(`int')!=`int') gen `levtype' = strofreal(`prec'*floor(`int')+`prec',"`form'")+" SDI" in 1/2
					else gen `levtype' = strofreal(`prec'*floor(`int'),"`form'")+" SDI" in 1/2
				}
				else if (!`unpair' | "`equnp'"!="") {
					if ("`varlist'"=="") {
						ttesti `sample1' `sample2', `options'
						scalar `sd' = sqrt(r(sd_1)^2+r(sd_2)^2-2*`rho'*r(sd_1)*r(sd_2))/sqrt(r(N_1))
					}
					else scalar `sd' = r(se)				

					gen	double `std_err' = r(sd_1)/sqrt(r(N_1)) in 1
					replace `std_err' = r(sd_2)/sqrt(r(N_2)) in 2
					if ("`equnp'"=="") scalar `zts' = invttail(`=r(N_1)-1',.5*(1-(`level'/100)))
					else scalar `zts' = invttail(r(df_t),.5*(1-(`level'/100)))
					scalar `ztd' = (`zts'*`sd'+`mvalue')/(`std_err'[1]+`std_err'[2])
					if (`ztd'<0) scalar `ztd' = 0
					gen double `zt_d' = `ztd' in 1/2
				
					if ("`varlist'"=="" & "`equnp'"=="") {
						mata: __s_lev = J(2,1,(1-2*ttail(st_numscalar("r(N_1)")-1,st_numscalar("`ztd'")))*100)
						scalar `int' = (1-2*ttail(r(N_1)-1,`ztd'))*100/`prec'
						if (`int'==.) {
							mata: _editmissing(__s_lev,100)
							scalar `int' = 100
						}
					}
					else {
						mata: __s_lev = J(2,1,(1-2*ttail(st_numscalar("r(df_t)"),st_numscalar("`ztd'")))*100)
						scalar `int' = (1-2*ttail(r(df_t),`ztd'))*100/`prec'
						if (`int'==. & `ztd'!=.) {
							mata: _editmissing(__s_lev,100)
							scalar `int' = 100
						}
					}
					if (int(`int')==100) gen `levtype' = "100 SDI" in 1/2
					else if (`ztd'==0) gen `levtype' = "0 SDI" in 1/2
					else if (int(`int')!=`int') gen `levtype' = strofreal(`prec'*floor(`int')+`prec',"`form'")+" SDI" in 1/2
					else gen `levtype' = strofreal(`prec'*floor(`int'),"`form'")+" SDI" in 1/2
				}
				else {
					if ("`varlist'"=="") ttesti `sample1' `sample2', `options'
					
					scalar `zts' = invttail(r(df_t),.5*(1-`level'/100))
					gen	double `std_err' = .
					gen	double `zt_d' = .
					gen	str14 `levtype' = ""
					mata: __s_lev = .
					mata: __s_tlevel(__s_lev,"`zts'","`zt_d'","`std_err'","`levtype'")
				}
				
				if (inlist(r(sd_1),0,.) | inlist(r(sd_2),0,.)) {
					di as err "invalid argument; one of the samples is a constant"
					exit 198		
				}

				gen	double `statis' = r(mu_1) in 1
				replace `statis' = r(mu_2) in 2
				gen double `lo_sdi' = `statis'-`zt_d'*`std_err' in 1/2
				gen double `hi_sdi' = `statis'+`zt_d'*`std_err' in 1/2

				if ("`difference'"!="") {
					replace `statis' = `statis'[1]-`statis'[2] in 3
					if (`mvalue' & `statis'[3] & (sign(`mvalue')!=sign(`statis'[3]))) local mtext = sign(`mvalue')
					if (`unpair') replace `std_err' = r(se) in 3
					else replace `std_err' = `sd' in 3
					replace `zt_d' = `zts' in 3
					replace `lo_sdi' = `statis'-`zt_d'*`std_err' in 3
					replace `hi_sdi' = `statis'+`zt_d'*`std_err' in 3
					replace `levtype' = "`level'  CI" in 3
					mata: __s_lev = __s_lev\strtoreal(st_local("level"))
					label def `lbl' 3 "(1-2)", add
					local nlab 3
					local dfci " & CI"

					if ("`nosamp'"=="") {
						if (`unpair') replace `df' = r(df_t) in 3
						else replace `df' = r(N_1)-1
					}
					if ("`r(sd)'"!="") mata: __s_sd = st_numscalar("r(sd)")
				}
				else if (`mvalue' & (`statis'[1]-`statis'[2]) & (sign(`mvalue')!=sign(`statis'[1]-`statis'[2]))) local mtext = sign(`mvalue')

				if ("`nosamp'"=="") {
					replace `df' = r(N_1)-1 in 1
					replace `df' = r(N_2)-1 in 2
				}
				mata: __s_s1 = st_numscalar("r(sd_1)")
				mata: __s_s2 = st_numscalar("r(sd_2)")

				if ("`power'"!="") {
					if ("`varlist'"!="") {
						local sample1 "`r(N_1)' `r(mu_1)' `r(sd_1)'"
						local sample2 "`r(N_2)' `r(mu_2)' `r(sd_2)'"
					}
					if (!`unpair') {
						if (abs(`rho')!=1) cap power pairedmeans `: word 2 of `sample1'', n(`=2* `: word 1 of `sample1''') sd1(`: word 3 of `sample1'') sd2(`: word 3 of `sample2'') corr(`=`rho'') alpha(`=1-`level'/100') power(`power')
						else cap power pairedmeans `: word 2 of `sample1'', n(`=2* `: word 1 of `sample1''') sd1(`: word 3 of `sample1'') sd2(`: word 3 of `sample2'') corr(`=`rho'-sign(`rho')*.0000000000001') alpha(`=1-`level'/100') power(`power')
						if (!_rc & r(converged)) {
							if (!strpos("`=`power'*100'",".")) local pow = string(`=`power'*100')
							else if (strlen(substr("`=`power'*100'",strpos("`=`power'*100'","."),.))==2) local pow = string(`=`power'*100',"%4.1f")
							else local pow = string(round(`=`power'*100',0.01),"%5.2f")
							local del = string(round(`r(da)',0.0001),"%`=strlen("`=int(`r(da)')'")+5'.4f")
							local dtext "The minimum detectable effect size with a power of `pow'% is 𝛿= `del'."
						}
						else local dtext "Minimum detectable value of the effect size could not be estimated."
						mata __s_del = st_numscalar("r(da)")
						if (substr("`dtext'",1,1)=="T") {
							mata __s_alp = st_numscalar("r(alpha)")
							mata __s_pow = st_numscalar("r(power)")
						}
					}
					else {
						cap power twomeans `: word 2 of `sample1'', n1(`: word 1 of `sample1'') n2(`: word 1 of `sample2'') sd1(`: word 3 of `sample1'') sd2(`: word 3 of `sample2'') alpha(`=1-`level'/100') power(`power')
						if (!_rc & r(converged)) {
							if (!strpos("`=`power'*100'",".")) local pow = string(`=`power'*100')
							else if (strlen(substr("`=`power'*100'",strpos("`=`power'*100'","."),.))==2) local pow = string(`=`power'*100',"%4.1f")
							else local pow = string(round(`=`power'*100',0.01),"%5.2f")
							local del = string(round(`r(delta)',0.0001),"%`=strlen("`=int(`r(delta)')'")+5'.4f")
							local dtext "The minimum detectable effect size with a power of `pow'% is 𝛿= `del'."
						}
						else local dtext "Minimum detectable value of the effect size could not be estimated."
						mata __s_del = st_numscalar("r(delta)")
						if (substr("`dtext'",1,1)=="T") {
							mata __s_alp = st_numscalar("r(alpha)")
							mata __s_pow = st_numscalar("r(power)")
						}
					}
				}

				putmata __s_matn=(`statis' `std_err' . . `df' `zt_d' `lo_sdi' `hi_sdi') in 1/`nlab'		
				mata: __s_matn[,3] = __s_matn[,1]:/__s_matn[,2]				
				if ("`nosamp'"=="") mata: __s_matn[,4] = 2*ttail(__s_matn[,5],abs(__s_matn[,2]))
				else mata: __s_matn[,4] = 2*(normal(-abs(__s_matn[,2])))
			}
			else if ("`power'"!="") {
				di as err "option {bf:power()} not allowed with {bf:percentile}"
				exit 198			
			}
	
			
			noi di as text _newline "{bf:SDI} Results"
			if ("`nolegend'"=="") {
				tempvar ra rb rc
				gen `ra' = "Comparison" in 1
				gen `rb' = ""
				gen `rc' = ""
				format `ra' %-19s
				format `rb' %-1s
				format `rc' %-45s

				if (`unpair') {		
					replace `rc' = "Two unpaired samples with `unequal' population variances" in 1			
					replace `ra' = "(1)" in 2
					if ("`varlist'"!="") {
						if ("`by'"=="") {
							replace `rc' = "`1' (variable)" in 2
							replace `rc' = "`2' (variable)" in 4
						}
						else {
							replace `rc' = "`1' if `by'==`bymi' (variable)" in 2
							replace `rc' = "`1' if `by'==`bymx' (variable)" in 4
						}
					}
					else if ("`reverse'"=="") {
						replace `rc' = "sample1" in 2
						replace `rc' = "sample2" in 4
					}
					else {
						replace `rc' = "sample2" in 2
						replace `rc' = "sample1" in 4
					}
					replace `ra' = "Number of obs (1)" in 3
					local len = strlen(string(`df'[1]+1))
					if (`len'>3) replace `rc' = string(`df'[1]+1,"%-"+"`=`len'+int(`len'/3)+2'"+".0gc") in 3
					else replace `rc' = string(`df'[1]+1) in 3
					replace `ra' = "(2)" in 4
					replace `ra' = "Number of obs (2)" in 5
					local len = strlen(string(`df'[2]+1))
					if (`len'>3) replace `rc' = string(`df'[2]+1,"%-"+"`=`len'+int(`len'/3)+2'"+".0gc") in 5
					else replace `rc' = string(`df'[2]+1) in 5
					replace `ra' = "Statistic" in 6
					if ("`median'"=="" & "`difference'"!="") replace `rc' = "mean; (1-2) = mean(1) - mean(2)" in 6
					else if ("`median'"=="") replace `rc' = "mean" in 6
					else replace `rc' = "median" in 6
					replace `ra' = "Interval type" in 7
					replace `rc' = "Standard error-based SDIs`dfci'" in 7
					replace `rb' = ":" in 1/7	
					replace `rb' = "=" if inlist(_n,3,5)
					local nra 7
				}
				else {
					replace `ra' = "Statistic" in 4
					if ("`median'"!="") replace `rc' = "median" in 4
					else if ("`difference'"!="") {
						if ("`varlist'"!="") replace `rc' = "mean; (1-2) = mean(1 - 2)" in 4
						else replace `rc' = "mean; (1-2) = mean(1) - mean(2)" in 4
					}
					else replace `rc' = "mean" in 4
					replace `rb' = ":" in 1/6

					if ("`nosamp'"!="") {
						if ("`distribution'"=="") {
							if (int(`rho')!=`rho') replace `rc' = "Two distributions with correlation ρ = "+string(round(`rho',0.01),"%4.2f") in 1
							else replace `rc' = "Two distributions with correlation ρ = `=`rho''" in 1
						}
						else {
							if (int(`rho')!=`rho') replace `rc' = "Two sampling distributions with correlation ρ = "+string(round(`rho',0.01),"%4.2f") in 1
							else replace `rc' = "Two sampling distributions with correlation ρ = `=`rho''" in 1
						}
						replace `ra' = "(1)" in 2
						replace `ra' = "(2)" in 3
						if ("`distribution'"!="") {
							replace `rc' = "`1' (variable)" in 2
							replace `rc' = "`2' (variable)" in 3
						}
						else if ("`reverse'"=="") {
							replace `rc' = "sample1" in 2
							replace `rc' = "sample2" in 3
						}
						else {
							replace `rc' = "sample2" in 2
							replace `rc' = "sample1" in 3
						}
						replace `ra' = "Interval type" in 5
						replace `rc' = "Standard error-based SDIs`dfci'" in 5
						local nra 5
					}
					else {
						if ("`varlist'"=="") {
							if (int(`rho')!=`rho') replace `rc' = "Two paired samples with correlation ρ = "+string(round(`rho',0.01),"%4.2f") in 1
							else replace `rc' = "Two paired samples with correlation ρ = `=`rho''" in 1
							if ("`reverse'"=="") {
								replace `rc' = "sample1" in 2
								replace `rc' = "sample2" in 3
							}
							else {
								replace `rc' = "sample2" in 2
								replace `rc' = "sample1" in 3
							}
						}
						else {
							if ("`distribution'"=="") replace `rc' = "Two paired samples" in 1
							else replace `rc' = "Two sampling distributions" in 1
							replace `rc' = "`1' (variable)" in 2
							replace `rc' = "`2' (variable)" in 3
						}
						replace `ra' = "(1)" in 2
						replace `ra' = "(2)" in 3
						replace `ra' = `ra'[4] in 5
						replace `rb' = `rb'[4] in 5
						replace `rc' = `rc'[4] in 5
						replace `ra' = "Interval type" in 6
						if ("`pct'"=="") replace `rc' = "Standard error-based SDIs`dfci'" in 6
						else replace `rc' = "Percentile-based SDIs`dfci'" in 6

						replace `ra' = "Number of obs" in 4
						local len = strlen(string(`df'[1]+1))
						if (`len'>3) replace `rc' = string(`df'[1]+1,"%-"+"`=`len'+int(`len'/3)+2'"+".0gc") in 4
						else replace `rc' = string(`df'[1]+1) in 4
						replace `rb' = "=" in 4
						local nra 6
					}
				}

				noi di as text ""
				noi flist `ra' `rb' `rc' in 1/`nra', noobs noh clean
			}

			gen `quant' = _n in 1/`nlab'
			label var `quant' "`=char(13)'"
			label values `quant' `lbl'
			label var `statis' "Statistic"
			label var `std_err' "Std. Err."
			label var `lo_sdi' "[Interval"
			label var `hi_sdi' "Bounds]"
			label var `levtype' "(%) Type"
			noi tabdis `quant' in 1/`nlab', c(`statis' `std_err' `lo_sdi' `hi_sdi' `levtype') left cellwidth(13)
			if ("`nolegend'"=="") {
				if ("`mtext'"=="" & "`power'"=="") noi di as text "Note: SDIs indicate significance of difference from `mvalue'."
				else {
					noi di as text "Notes: SDIs indicate significance of difference from `mvalue'."
					if ("`mtext'"=="-1") noi di as text _skip(7) "Meaningful value {it:m} is negative whereas the difference in estimates is positive."
					else if ("`mtext'"=="1") noi di as text _skip(7) "Meaningful value {it:m} is positive whereas the difference in estimates is negative."
					if ("`power'"!="") {
						if (substr("`dtext'",1,1)=="T" & (`statis'[1]<`statis'[2])) local dtext = subinstr("`dtext'","= ","= -",1)
						noi di as text _skip(7) "`dtext'"
					}
				}
			}

			mata: st_rclear()
			mata: st_global("r(cmd)","sdii")
			mata: st_global("r(cmdline)",st_local("cmdl"))
			mata: st_numscalar("r(level)",strtoreal(st_local("level")))
			mata: st_numscalar("r(mvalue)",strtoreal(st_local("mvalue")))
			mata: st_numscalar("r(sd_1)",__s_s1)
			mata: st_numscalar("r(sd_2)",__s_s2)
			
			if ("`nosamp'"=="") {
				mata: st_numscalar("r(N_1)",st_data(1,"`df'")+1)
				mata: st_numscalar("r(N_2)",st_data(2,"`df'")+1)
				if ("`difference'"!="") {
					if (!`unpair') mata: st_numscalar("r(N_d)",st_numscalar("r(N_1)"))
					else mata: st_numscalar("r(N_d)",st_numscalar("r(N_1)")+st_numscalar("r(N_2)"))
				}
				
				if ("`power'"!="") {
					mata: st_numscalar("r(delta)", __s_del)
					if (substr("`dtext'",1,1)=="T") {
						mata: st_numscalar("r(alpha)", __s_alp)
						mata: st_numscalar("r(power)", __s_pow)
					}
				}
			}

			mata: st_matrix("r(sdii)",(__s_matn,__s_lev))
			if ("`pct'"=="") {
				if ("`nosamp'"=="") mata: st_matrixcolstripe("r(sdii)",(J(9,1,""),("b"\"se"\"t"\"pvalue"\"df"\"crit"\"ll"\"ul"\"level")))
				else mata: st_matrixcolstripe("r(sdii)",(J(9,1,""),("b"\"se"\"z"\"pvalue"\"df"\"crit"\"ll"\"ul"\"level")))
				if ("`varlist'"=="" & !`unpair') mata: st_numscalar("r(rho)",st_numscalar("`rho'"))
			}
			else mata: st_matrixcolstripe("r(sdii)",(J(9,1,""),("b"\"se"\"t"\"pvalue"\"df"\"crit"\"ll"\"ul"\"level")))

			if ("`difference'"!="") {
				mata: st_numscalar("r(sd_d)",__s_sd)
				mata: st_global("r(type)","SDI,SDI,CI")
				mata: st_matrixrowstripe("r(sdii)",(J(strtoreal(st_local("nlab")),1,""),("(1)"\"(2)"\"(1-2)")))
			}
			else {
				mata: st_global("r(type)","SDI,SDI")
				mata: st_matrixrowstripe("r(sdii)",(J(strtoreal(st_local("nlab")),1,""),("(1)"\"(2)")))
			}
		}
	}

	if (_rc) {
		local rc = _rc
		cap if ("`N1'"!="") drop in `N1'/`=_N'
		cap mata: mata drop __s_*
		exit(`rc')
	}
	else {
		cap if ("`N1'"!="") drop in `N1'/`=_N'
		cap mata: mata drop __s_*
	}

end


program __s_var2
	version 14.2

	syntax newvarname(max=1 numeric generate)=exp [, ifn(str) ]
	replace `varlist' `exp' `ifn'
end


version 14.2
mata:

	void __s_centile(__s_k,__s_lev,__s_lhs,__s_pct,__s_rhs,string scalar v1,string scalar v2,string scalar pk,string scalar typ) {
		__s_lhs = __s_pct[,1]:-__s_pct[,2]
		perc = range(0,50,.5) , range(100,50,.5)		
		prec = 1/10^strtoreal(st_local("precision"))

		plus = 0
		for (i=101; i>0; i--) {
			if (__s_lhs[i]>__s_rhs & i!=1) continue
			else if (__s_lhs[i]<__s_rhs & i!=101) {
				__s_k = perc[i,1] , .05
				__s_lefths(__s_k,__s_lhs,v1,v2)

				for (j=10; j>0; j--) {
					if (__s_lhs[j]>__s_rhs) continue
					else if (__s_lhs[j]<__s_rhs) {
						__s_k = __s_k[1]+(j-1)*.05 , .005
						__s_lefths(__s_k,__s_lhs,v1,v2)
					
						for (k=10; k>0; k--) {
							if (__s_lhs[k]>__s_rhs) continue
							else if (__s_lhs[k]<__s_rhs) {
								__s_k = __s_k[1]+(k-1)*.005 , .0005
								__s_lefths(__s_k,__s_lhs,v1,v2)

								for (l=10; l>0; l--) {
									if (__s_lhs[l]>__s_rhs) continue
									else if (__s_lhs[l]<__s_rhs) {
										__s_k = __s_k[1]+(l-1)*.0005 , .00005
										__s_lefths(__s_k,__s_lhs,v1,v2)
										
										for (m=10; m>0; m--) {
											if (__s_lhs[m]>__s_rhs) continue
											else if (__s_lhs[m]<__s_rhs) {
												__s_k = __s_k[1]+(m-1)*.00005 , .000005
												__s_lefths(__s_k,__s_lhs,v1,v2)

												for (n=10; n>0; n--) {
													if (__s_lhs[n]>__s_rhs) continue
													else if (__s_lhs[n]<__s_rhs) {
														__s_k = __s_k[1]+(n-1)*.000005 , .0000005
														__s_lefths(__s_k,__s_lhs,v1,v2)

														for (o=10; o>0; o--) {
															if (__s_lhs[o]>__s_rhs) continue
															else if (sign(__s_rhs)==1) {
																if (sign(__s_lhs[o])==1) kp = __s_k[1]+(o-1)*.0000005
																else kp = __s_k[1]+o*.0000005
															}
															else kp = __s_k[1]+(o-1)*.0000005
															if (__s_lhs[o]<__s_rhs) plus = 1 ;;
															break
														}
													}
													else kp = __s_k[1]+(n-1)*.000005
													break
												}
											}
											else kp = __s_k[1]+(m-1)*.00005
											break
										}
									}
									else kp = __s_k[1]+(l-1)*.0005
									break
								}
							}
							else kp = __s_k[1]+(k-1)*.005
							break
						}
					}
					else kp = __s_k[1]+(j-1)*.05
					break
				}
			}
			else kp = perc[i,1]																	
			break
		}
		
		st_numscalar(pk,kp)
		__s_lev = J(2,1,100-2*kp)
		if (kp==0) sdi = J(2,1,"100 SDI")
		if (kp==50) sdi = J(2,1,"0 SDI")
		else if (plus) sdi = J(2,1,strofreal(prec*floor((100-2*kp)/prec)+prec,st_local("form"))):+" SDI"
		else sdi = J(2,1,strofreal(prec*floor((100-2*kp)/prec),st_local("form"))):+" SDI"
		st_sstore((1,2),typ,sdi)
	}

	void __s_lefths(__s_k,__s_lhs,v1,v2) {

			perc = range(__s_k[1],__s_k[1]+(10*__s_k[2]-__s_k[2]),__s_k[2]) , range(100-__s_k[1],100-(__s_k[1]+(10*__s_k[2]-__s_k[2])),__s_k[2])
			pct = J(10,2,.)
			st_local("numl",invtokens(strofreal(perc[,1]',"%10.7f")))
			stata("_pctile "+v1+st_local(" inf")+st_local("wght")+", p("+st_local("numl")+")"+st_local("alt"))
			for (i=1; i<=10; i++) pct[i,1] = st_numscalar("r(r"+strofreal(i)+")")
			st_local("numl",invtokens(strofreal(sort(perc[,2],1)',"%10.7f")))
			stata("_pctile "+v2+st_local(" inf")+st_local("wght")+", p("+st_local("numl")+")"+st_local("alt"))
			for (i=1; i<=10; i++) pct[i,2] = st_numscalar("r(r"+strofreal(11-i)+")")
			
			__s_lhs = pct[,1]:-pct[,2]
	}

	void __s_tlevel(__s_lev,string scalar ts,string scalar td,string scalar se,string scalar typ) {

		prec = 1/10^strtoreal(st_local("precision"))
		ts = st_numscalar(ts)
		rhs = ts*st_numscalar("r(se)")+strtoreal(st_local("mvalue"))

		pmax = 99.99999999999999 // maximum percentile
		perc = (0::99)\pmax
		sub = range(.9,0,.1)
		dif = st_numscalar("r(mu_1)")-st_numscalar("r(mu_2)")
		df1 = st_numscalar("r(N_1)")-1
		se1 = st_numscalar("r(sd_1)")/sqrt(st_numscalar("r(N_1)"))
		df2 = st_numscalar("r(N_2)")-1
		se2 = st_numscalar("r(sd_2)")/sqrt(st_numscalar("r(N_2)"))
		lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2
	
		plus = 0
		for (i=1; i<=101; i++) {

			if (lhs[i]<rhs & i!=101) continue
			else if (lhs[i]>rhs & i!=1) {
				if (i!=101) perc = perc[i]:-sub
				else perc = 100:-sub
				lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2
	
				for (j=1; j<=10; j++) {
					if (lhs[j]<rhs) continue
					else if (lhs[j]>rhs) {
						perc = perc[j]:-(sub/10)
						lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2
					
						for (k=1; k<=10; k++) {
							if (lhs[k]<rhs) continue
							else if (lhs[k]>rhs) {
								perc = perc[k]:-(sub/100)
								lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2

								for (l=1; l<=10; l++) {
									if (lhs[l]<rhs) continue
									else if (lhs[l]>rhs) {
										perc = perc[l]:-(sub/1000)
										lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2
										
										for (m=1; m<=10; m++) {
											if (lhs[m]<rhs) continue
											else if (lhs[m]>rhs) {
												perc = perc[m]:-(sub/10000)
												lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2

												for (n=1; n<=10; n++) {
													if (lhs[n]<rhs) continue
													else if (lhs[n]>rhs) {
														perc = perc[n]:-(sub/100000)
														lhs = invttail(df1,.5*(1:-perc/100))*se1 :+ invttail(df2,.5*(1:-perc/100))*se2
														
														for (o=1; o<=10; o++) {
															if (lhs[o]<rhs & o!=10) continue
															else if (sign(dif-rhs)==1) {
																if (sign(dif-lhs[o])==1) {
																	if (o==10 & i==101) perc = pmax
																	else perc = perc[o]
																}
																else perc = perc[o]-.000001
															}
															else {
																if (o==10 & i==101) perc = pmax
																else perc = perc[o]
															}
															if (lhs[o]>rhs & i!=101) plus = 1 ;;
															break
														}
													}
													else perc = perc[n]
													break
												}
											}
											else perc = perc[m]
											break
										}
									}
									else perc = perc[l]
									break
								}
							}
							else perc = perc[k]
							break
						}
					}
					else perc = perc[j]
					break
				}
			}
			else perc = perc[i]
			break
		}

		t1_2 = invttail((df1\df2),(.5*(1-perc/100)))
		st_store((1,2),td,t1_2)
		if (perc>99.999999) perc = 100 ;; 
		__s_lev = J(2,1,perc)
		if (perc==0) sdi = J(2,1,"0 SDI")
		else if (perc==100) sdi = J(2,1,"100 SDI")
		else if (plus) sdi = J(2,1,strofreal(prec*floor(perc/prec)+prec,st_local("form"))):+" SDI"
		else sdi = J(2,1,strofreal(prec*floor(perc/prec),st_local("form"))):+" SDI"
		st_store((1,2),se,(se1\se2))
		st_sstore((1,2),typ,sdi)
	}

end
