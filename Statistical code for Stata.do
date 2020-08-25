version 16


****** Descriptive analysis and definition of variables *************

tab ensidementia     // byte    %8.0g                any dementia diagnosis (yes/no)
tab dementiatyyppi	//								 type of dementia
sum anydementiapvm  // int     %9.0g                 date of first dementia diagnosis

tab anyanyinf       // byte    %8.0g      yesno      Hospitalisation for any infectious disease

tab anybactinf      // byte    %8.0g      yesno      Hospitalisation for any bacterial infection

tab anysysbact      // byte    %8.0g      yesno      Hospitalisation for any potentially invasive bacterial infections
tab anylocalbact    // byte    %8.0g      yesno      Hospitalisation for any localised bacterial infection

tab anygramplus     // byte    %8.0g      yesno      Hospitalisation for any Gram-positive bacterial infection

tab anygramminus    // byte    %8.0g      yesno      Hospitalisation for any Gram-negative bacterial infection

tab anymycobact     // byte    %8.0g      yesno      Hospitalisation for any mycobacterial infection

tab anymycoplasma   // byte    %8.0g      yesno      Hospitalisation for any mycoplasma infection

tab anyextrac       // byte    %8.0g      yesno      Hospitalisation for any extracellular bacterial infection
tab anyintrac       // byte    %8.0g      yesno      Hospitalisation for any intracellular bacterial infection

tab anyviralinf     // byte    %8.0g      yesno      Hospitalisation for any viral infection

tab anyacuteviral   // byte    %8.0g      yesno      Hospitalisation for any acute viral infection
tab anyherpes       // byte    %8.0g      yesno      Hospitalisation for any herpesvirus infection
tab anypersviral    // byte    %8.0g      yesno      Hospitalisation for any potentially persistent viral infection (excluding herpesvirus infections)

tab anyparasite     // byte    %8.0g      yesno      Hospitalisation for any parasitic infection

tab anymycose       // byte    %8.0g      yesno      Hospitalisation for any fungal infection (mycosis)

tab anybactsepsis		// byte    %8.0g      yesno      Hospitalisation for bacterial infection with sepsis
tab anybact_nonsepsis	// byte    %8.0g      yesno      Hospitalisation for bacterial infection without sepsis
tab anyanycns       	// byte    %8.0g      yesno      Hospitalisation for central nervous system infection
tab anynoncsn       	// byte    %8.0g      yesno      Hospitalisation for extra-central nervous system infection
tab anycnspredisp       // byte    %8.0g      yesno      Hospitalisation for an infection predisposed to affect the central nervous system
tab anynoncnspredisp	// byte    %8.0g      yesno      Hospitalisation for an infection not predisposed to affect the central nervous system
tab anynonpers       	// byte    %8.0g      yesno      Hospitalisation for acute infection
tab anypers     		// byte    %8.0g      yesno      Hospitalisation for chronic infection

*For all categories of infections listed above, there exists a variable for the date of first hospitalisation for an infection in that category.
*The variable is named "ensi(name_of_disease)pvm" where "name_of_disease" is that of the variable above when the prefix "any" is removed. 

tab anyhiv			//								Known HIV (yes/no)

sum entrypvm			//							date of entry to cohort
sum exitpvm				//							date of end of follow-up
sum syntpvm             //							date of birth

tab ses				// 								socioeconomic status (low, intermediate, high)
tab supu			//								sex (man/woman)
tab alcocl			//								alcohol drinking (no, light, moderate, heavy)
tab tupakka			//								smoking (never, ex, current)
tab bmi_who				//							body mass index (normal weight, overweight, obese)	
sum hypertensiondate	//							date of hypertension diagnosis
sum diabetesdate		//							date of diabetes diagnosis



*** PRIMARY ANALYSES USING THREE FINNISH COHORTS ***


*FIGURE 2.
capture log close
log using "FIGURE 2.log", replace
local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral
local row 7 10 13 14 16 17 20 21 24 25 28 31 32 33
local n_dg: word count `dg'
local dg_name `""Any infectious disease vs not" "Any bacterial infection vs not" "Potentially invasive bacterial infection vs not"  "Localised bacterial infection vs not" "Bacterial infection with sepsis vs not" "Bacterial infection without sepsis vs not" "Extracellular bacterial infection vs not" "Intracellular bacterial infection vs not" "Gram-positive bacterial infection vs not" "Gram-negative bacterial infection vs not" "Any viral infection vs not" "Herpesvirus (persistent) infection vs not" "Other persistent viral infection vs not" "Acute viral infection vs not""'
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local irow : word `i' of `row'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local irow = `irow' + 29 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		}
log close

*FIGURE 3.
capture log close
log using "FIGURE 3.log", replace
local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral
local row 7 11 15 18 21 24 28 31 35 38 42 46 49 52
local n_dg: word count `dg'
local dg_name `""Any infectious disease vs not" "Any bacterial infection vs not" "Potentially invasive bacterial infection vs not"  "Localised bacterial infection vs not" "Bacterial infection with sepsis vs not" "Bacterial infection without sepsis vs not" "Extracellular bacterial infection vs not" "Intracellular bacterial infection vs not" "Gram-positive bacterial infection vs not" "Gram-negative bacterial infection vs not" "Any viral infection vs not" "Herpesvirus (persistent) infection vs not" "Other persistent viral infection vs not" "Acute viral infection vs not""'
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local irow : word `i' of `row'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		*Alzheimer's disease dementia
		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		quietly: stset exitpvm, id(id) failure(dementiatyyppi == 1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'_AD")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		
		*Non-Alzheimer's disease dementia
		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		quietly: stset exitpvm, id(id) failure(dementiatyyppi == 2 3 4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'_nonAD")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		}
log close

	
*FIGURE sensitivity, acute and chronic infections
capture log close
log using "Sensitivity, acute, chronic, coinfections.log", replace
local dg nonpers anypers bactviralcoinf
local n_dg: word count `dg'
local dg_name `""Any acute infectious disease" "Any chronic infectious disease""'

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "Coinfections_trimmattu"
		drop _merge
		drop if ses==.
		
quietly: local irow = 7
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
log close

		
***** SINGLE VS MULTIPLE INFECTIONS *******		

capture log close
log using "Single vs multiple infections.log", replace

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "Coinfections_trimmattu"
		drop _merge
		drop if ses==.
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id anyanyinf supu cohort exitpvm ensidementia syntpvm entrypvm ensi*anyinfpvm *toinen* *kolmas*  _st _d _origin _t _t0 ses
		sum ensi*anyinfpvm
		
		stsplit first_split, at(0) after(time=ensianyinfpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anyanyinf==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=ensitoinenanyinfpvm)
		replace multiinf_timedep = 2 if toinenanyinf==1 & second_split==0
		tab multiinf_timedep
		stsplit third_split, at(0) after(time=ensikolmasanyinfpvm)
		replace multiinf_timedep = 3 if kolmasanyinf==1 &third_split==0
		tab multiinf_timedep
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/3 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 & third_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 & third_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		
		
		*Dementia occurring from year 10 onwards
		local i 1
		foreach idg in anyinf toinenanyinf kolmasanyinf {
		display `i'
		replace _t0 = _t0 + 10 if multiinf_timedep==0 & `i'==1
		drop if _t0>=_t
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=. & multiinf_timedep==`i'
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==. & multiinf_timedep==`i'
		drop if _t0>=_t
		local ++i
		}
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		local irow 17
		forvalues i =1/3 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Dementia occurring more than 10 years after the hospitalisation, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 & third_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 & third_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow 
		local ++irow 
		}
		
		*P for trend
		local ++prow
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)

log close

***** One or several different pathogens *******		

capture log close
log using "One vs at least two different pathogens.log", replace

local dg_name `""Multiple organisms""'
local irow 82
local prow = `irow'

		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		egen ensibactviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
		gen toinenbactviralpvm = ensibactviralpvm if ensibactinfpvm==ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensibactinfpvm if ensibactinfpvm>ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensiviralinfpvm if ensiviralinfpvm>ensibactinfpvm & ensiviralinfpvm!=.
		gen anybactviral = 1 if ensibactviralpvm!=.
		gen toinenbactviral = 1 if toinenbactviralpvm!=.
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral
		tab toinenbactviral
		
		stsplit first_split, at(0) after(time=ensibactviralpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anybactviral==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=toinenbactviralpvm)
		replace multiinf_timedep = 2 if toinenbactviral==1 & second_split==0
		tab multiinf_timedep
		*browse id ensibactviralpvm toinenbactviralpvm multiinf_timedep _t0 _t
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		
		
		*Dementia occurring from year 10 onwards
		local i 1
		foreach idg in ensibactviral toinenbactviral {
		display `i'
		replace _t0 = _t0 + 10 if multiinf_timedep==0 & `i'==1
		drop if _t0>=_t
		gen `idg'pvm10v = `idg'pvm + round(10*365.25) if `idg'pvm<entrypvm
		replace `idg'pvm10v = entrypvm if `idg'pvm10v<entrypvm
		replace _t0 = (`idg'pvm10v-syntpvm)/365.25 if `idg'pvm10v!=. & multiinf_timedep==`i'
		replace _t0 = _t0 + 10 if `idg'pvm10v==. & multiinf_timedep==`i'
		drop if _t0>=_t
		local ++i
		}
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		
		local irow 83
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Dementia occurring more than 10 years after the hospitalisation, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1  &_d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow 
		local ++irow 
		}
		
		*P for trend
		local ++prow
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)

log close



***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

capture log close
log using "Number of infections.log", replace

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of infections""'

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "Infection count.dta"
		tab anyinfcount
		replace anyinfcount=0 if anyanyinf==0
		recode anyinfcount 3=2 4=2
		drop if ses==.
		gen ensiinfcountpvm = ensianyinfpvm 
		
quietly: local irow = 23
quietly: local prow = `irow'
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = any`idg' if inlist(any`idg',1,2,3) & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab `idg'_timedep, matcell(x), if `idg'_split!=-1
		matrix list x
		tab `idg'_timedep ensidementia, matcell(y), if `idg'_split!=-1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)		
		
		local irow 24
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/3 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab `idg'_timedep, matcell(x), if `idg'_split!=-1
		matrix list x
		tab `idg'_timedep ensidementia, matcell(y), if `idg'_split!=-1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		local ++prow
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)		
		
		restore
		}	
log close		


*FIGURE sensitivity, CNS vs other infections
capture log close
log using "Sensitivity, CNS infection vs other infection.log", replace
local dg anycns noncns
local n_dg: word count `dg'
local dg_name `""Any CNS infection vs not" "Any extra-CNS infection vs not""'

		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		
quietly: local irow = 46
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
log close


*FIGURE sensitivity, infections predisposed to affect the CNS vs other infections
capture log close
log using "Sensitivity, predisposed towards CNS vs others.log", replace
local dg anycnspredisp noncnspredisp
local n_dg: word count `dg'
local dg_name `""Infection predisposed to affect the CNS vs no infection" "Infection not predisposed to affect the CNS vs no infection""'

		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.
		
quietly: local irow = 76
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
log close


*TYPE OF HERPESVIRUS

capture log close
log using "FIGURE herpestype.log", replace
local dg hsv12 nonhsv mildherpes
local n_dg: word count `dg'
local dg_name `""Herpes simplex virus 1 and 2 infection" "Other herpesvirusinfection" "Mild herpesvirusinfection""' 
quietly: local irow = 7

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "herpeskoodit.dta"
		drop _merge
		drop if ses==.
		tab anyhsv12
		gen anynonhsv = 0
		replace anynonhsv = 1 if anyherpes==1 & anyhsv12!=1
		
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		recode any`idg' .=0
		capture drop ensi`idg'pvm
		gen ensi`idg'pvm = ensiherpespvm if any`idg'==1
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE herpes sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE herpes sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		restore
		}

	log close

	

*MOST COMMON GRAM-POSITIVE INFECTIONS

capture log close
log using "FIGURE gramplus.log", replace
local dg streptoc cdiffic othergramplus
local n_dg: word count `dg'
local dg_name `""Streptococcal infection" "Enterocolitis due to Clostridium difficile" "Other Gram-positive infection""'
quietly: local irow = 7

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "grampluskoodit.dta"
		drop _merge
		drop if ses==.
		
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		recode any`idg' .=0
		gen ensi`idg'pvm = ensigrampluspvm if any`idg'==1
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE gramplus") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE gramplus") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		restore
		}

	log close

*MOST COMMON GRAM-NEGATIVE INFECTIONS

capture log close
log using "FIGURE gramminus.log", replace
local dg salmonella lyme othergramminus
local n_dg: word count `dg'
local dg_name `""Salmonella infection" "Lyme disease" "Other Gram-negative infection""'
quietly: local irow = 7

		quietly: use "Pooled infections_trimmattu", clear
		merge 1:1 id using "gramminuskoodit.dta"
		drop _merge
		drop if ses==.
		
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		recode any`idg' .=0
		replace anyothergramminus=1 if anyothergramminussepsis==1
		gen ensi`idg'pvm = ensigramminuspvm if any`idg'==1
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE gramminus") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE gramminus") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		restore
		}

	log close

	
	
*** FINE AND GRAY MODEL ***

capture log close
log using "FIGURE fine gray.log", replace
local dg anyinf bactinf viralinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Any bacterial infection" "Any viral infection""'
quietly: local irow = 7

		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.

		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab ensidementia
		*Coding of detailedfail: 1 early-onset dementia, 2 late-onset dementia, 3 death
		gen byte detailedfail = 1 if _d==1 & _t<65
		replace detailedfail = 2 if _d==1 & _t>=65
		replace detailedfail = 3 if failure==2
		tab ensidementia
		tab failure
		tab detailedfail
		
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		forvalues j=1/3 {
		preserve
		if `j'==1 | `j'== 2 {
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25)
		}
		else if `j'== 3 {
		stset exitpvm, id(id) failure(detailedfail==2) origin(time syntpvm) enter(time entrypvm) scale(365.25)
		}
		
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses detailedfail
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		if `j'==1 {
		stcox `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		}
		else if `j'== 2 {
		stcrreg `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort i.cohort, cl(cohort) compete(detailedfail==3)
		}
		else if `j'== 3 {
		stcrreg `idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort i.cohort, cl(cohort) compete(detailedfail==1 3)
		}
		
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE fine gray (appendix)") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			if `j'==1 {
			quietly: putexcel D`irow'=("`idg_name', Cox model")
			}
			else if `j'== 2 {
			quietly: putexcel D`irow'=("`idg_name', Fine & Gray model, competing outcome death + early-onset dementia")
			}
			else if `j'== 3 {
			quietly: putexcel D`irow'=("`idg_name', Fine & Gray model, competing outcome death")
			}
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE fine gray (appendix)") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		restore
		}
		local ++irow
		}
	log close	
	
	
	
*** TIME INTERACTION ***

capture log close
log using "Time interaction.log", replace

		quietly: use "Pooled infections_trimmattu", clear
		drop if ses==.

local dg anyinf bactinf viralinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Bacterial infection" "Viral infection""'
quietly: local irow = 7
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm ses
		egen entrypvm`idg' = rowmax(ensi`idg'pvm entrypvm)
		stset exitpvm, id(id) failure(ensidementia) origin(time entrypvm`idg') scale(365.25) 
		gen age`idg' = (entrypvm`idg'-syntpvm)/365.25
		gen age`idg'2 = age`idg'^2 
		
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.age`idg'#cohort c.age`idg'2#cohort, strata(cohort)
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.age`idg'#cohort c.age`idg'2#cohort, strata(cohort) tvc(any`idg') texp(ln(_t))
		restore
		}
log close		





****** REPLICATION ANALYSIS IN THE UK BIOBANK ******

capture log close
log using "Sipila, infections and dementia, UKB replication $S_DATE.log", replace

		
/*Generate a wide dataset (one row per participant) with the following variables:

variable name:			variable type				description

id						integer						participant id
dementia				binary (1/0)				Dementia during follow-up. Please use the established definition in the cohort.
dementia_AD_nonAD		binary (1/2)				Type of dementia: 1) Alzheimer's disease 2) other dementia	
dementiatype			categorical	(1,2,3,4,5,6)	Type of first dementia diagnosis: 1) Alzheimer's disease, 2) Frontotemporal dementia, 3) Parkinson's disease dementia, 4) Vascular dementia 5) Other specified dementia, 6) Unspecified dementia
exitdate				date						First of the following: date of first dementia diagnosis, date of death, end of follow-up (dementia diagnoses that occur before study entry should be included if such data are available)
birthday				date						Date of Birth
entrydate				date						entry
sex						binary (1/2)				Sex	
ses						categorical	(1,2,3)			Socioeconomic status or education (low, intermediate, high)
alcocl					categorical	(0,1,2,3)		non-drinkers, moderate drinkers, intermediate drinkers, and heavy drinkers
smoking					categorical	(1,2,3)			never smokers, ex-smokers, current smokers
bmi_who					categorical	(0,25,30)		normal (≤24.9 kg/m2), overweight (25.0-29.9 kg/m2), obese (≥30.0 kg/m2).
hypertensiondate		date						date of hypertension diagnosis
diabetesdate			date						date of diabetes diagnosis

*/

*Store the data to Stata format (e.g. SipilaELSA.dta)
*Please record all definitions of variables if they differ from those suggested above.



*Trim the data
use "SipilaELSA.dta", clear // Note that, despite the funny name, this file contains the UK Biobank data.
count
drop if exitdate<=entrydate
count
drop if ses==.
count

tab dementia
recode dementia .=0
tab dementia 

foreach idg in anyinf bactinf viralinf hiv {
replace ensi`idg'pvm=. if ensi`idg'pvm>=exitdate
replace any`idg' = 0 if ensi`idg'pvm==.
sum ensi`idg'pvm
tab any`idg'
}

*Create variables that have indentical names with those used in the primary analysis

gen diabetes_combpvm = diabetesdate
gen hypert_combpvm = hypertensiondate
gen tupakka = smoking
*alcocl does not change
*ses does not change
gen entrypvm = entrydate
gen exitpvm = exitdate
gen ensidementia = dementia
recode ensidementia 0=.
gen syntpvm = birthday
gen supu = sex
*dementiatype ei muutu

save "SipilaELSA.dta", replace


**** Analysis *****

use "SipilaELSA.dta", clear


local dg anyinf bactinf viralinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Bacterial infection" "Viral infection""'
quietly: local irow = 7

	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' , adjusted for sex and education and other variables"
		display ""
		display "UKB replication analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""

		/* Models
		MODEL 1: adjusted for age (as the time scale), sex and ses
		MODEL 2 = model 1 + those with HIV excluded
		MODEL 3 = model 2 + those with missing data on alcohol, smoking, and bmi excluded 
		MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes
		MODEL 5 = model 4 + those with data on apoe4 missing dropped
		MODEL 6 = model 5 + adjusted for apoe4 genotype
		*/
		
	forvalues j=1/6 {
	
		preserve
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25)
		tab any`idg'
		tab any`idg' _d
		sum diabetes_combpvm
		replace diabetes_combpvm = . if diabetes_combpvm>=exitpvm
		sum diabetes_combpvm
		gen dm = 1 if diabetes_combpvm < exitpvm
		tab dm
		sum hypert_combpvm
		replace hypert_combpvm = . if hypert_combpvm>=exitpvm
		sum hypert_combpvm
		gen htn = 1 if hypert_combpvm < exitpvm
		tab htn
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stsplit dm_split, at(0) after(time=diabetes_combpvm)
		gen dm_timedep =0
		replace dm_timedep = 1 if dm==1 & dm_split==0
		tab dm_timedep
		tab dm_timedep if dm_split!=-1 & `idg'_split!=-1
		
		stsplit htn_split, at(0) after(time=hypert_combpvm)
		gen htn_timedep =0
		replace htn_timedep = 1 if htn==1 & htn_split==0
		tab htn_timedep
		tab htn_timedep if htn_split!=-1 & `idg'_split!=-1 & dm_split!=-1
		
		if `j'==1 {
		di ""
		di "MODEL 1: adjusted for age (as the time scale), sex and ses"
		di ""
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==2 {
		di ""
		di "MODEL 2 = model 1 + those with HIV excluded"
		di ""
		drop if anyhiv==1
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==3 {
		di ""
		di "MODEL 3 = model 2 + those with missing data on alcohol, smoking, and bmi excluded"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==4 {
		di ""
		di "MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes (please note that this is model 3 in the published manuscript)"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep
		}
		else if `j'==5 {
		di ""
		di "MODEL 5 = model 4 + those with data on apoe4 missing dropped"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		drop if apoe==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep
		}
		else if `j'==6 {
		di ""
		di "MODEL 6 = model 5 + adjusted for apoe4 genotype (please note that this is model 4 in the published manuscript)"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		drop if apoe==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep i.apoe
		}
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE UKB") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' model `j'")
		tab any`idg', matcell(x), if `idg'_split!=-1 & htn_split!=-1 & dm_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1 & htn_split!=-1 & dm_split!=-1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE UKB") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		restore
		} 
		local ++irow 
		}



	
**** Analysis by dementia subtypes *****

use "SipilaELSA.dta", clear

local dg anyinf bactinf viralinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Bacterial infection" "Viral infection""'

foreach outcome in AD VD FTD PD {
	if "`outcome'" == "AD" {
	local type=1 
	}
	else if "`outcome'" == "VD" {
	local type=4
	}
	else if "`outcome'" == "FTD" {
	local type=2
	}
	else if "`outcome'" == "PD" {
	local type=3
	}

	quietly: local irow = 7
	
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' , adjusted for sex and education and other variables"
		display ""
		display "UKB replication analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""

		/* Models
		MODEL 1: adjusted for age (as the time scale), sex and ses
		MODEL 2 = model 1 + those with HIV excluded
		MODEL 3 = model 2 + those with missing data on alcohol, smoking, and bmi excluded 
		MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes (please note that this is model 3 in the published manuscript)
		MODEL 5 = model 4 + those with data on apoe4 missing dropped
		MODEL 6 = model 5 + adjusted for apoe4 genotype (please note that this is model 4 in the published manuscript)
		*/
		
	forvalues j=1/6 {
	
		preserve
		
		stset exitpvm, id(id) failure(dementiatype==`type') origin(time syntpvm) enter(time entrypvm) scale(365.25)
		tab any`idg'
		tab any`idg' _d
		sum diabetes_combpvm
		replace diabetes_combpvm = . if diabetes_combpvm>=exitpvm
		sum diabetes_combpvm
		gen dm = 1 if diabetes_combpvm < exitpvm
		tab dm
		sum hypert_combpvm
		replace hypert_combpvm = . if hypert_combpvm>=exitpvm
		sum hypert_combpvm
		gen htn = 1 if hypert_combpvm < exitpvm
		tab htn
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stsplit dm_split, at(0) after(time=diabetes_combpvm)
		gen dm_timedep =0
		replace dm_timedep = 1 if dm==1 & dm_split==0
		tab dm_timedep
		tab dm_timedep if dm_split!=-1 & `idg'_split!=-1
		
		stsplit htn_split, at(0) after(time=hypert_combpvm)
		gen htn_timedep =0
		replace htn_timedep = 1 if htn==1 & htn_split==0
		tab htn_timedep
		tab htn_timedep if htn_split!=-1 & `idg'_split!=-1 & dm_split!=-1
		
		if `j'==1 {
		di ""
		di "MODEL 1: adjusted for age (as the time scale), sex and ses"
		di ""
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==2 {
		di ""
		di "MODEL 2 = model 1 + those with HIV excluded"
		di ""
		drop if anyhiv==1
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==3 {
		di ""
		di "MODEL 3 = model 2 + those with missing data on alcohol, smoking, and bmi excluded"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		}
		else if `j'==4 {
		di ""
		di "MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes (please note that this is model 3 in the published manuscript)"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep
		}
		else if `j'==5 {
		di ""
		di "MODEL 5 = model 4 + those with data on apoe4 missing dropped"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		drop if apoe==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep
		}
		else if `j'==6 {
		di ""
		di "MODEL 6 = model 5 + adjusted for apoe4 genotype (please note that this is model 4 in the published manuscript)"
		di ""
		drop if anyhiv==1
		drop if tupakka==.
		drop if alcocl==.
		drop if bmi_who==.
		drop if apoe==.
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn_timedep i.dm_timedep i.apoe
		}
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE UKB, `outcome'") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' model `j'")
		tab any`idg', matcell(x), if `idg'_split!=-1 & htn_split!=-1 & dm_split!=-1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1 & htn_split!=-1 & dm_split!=-1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE UKB, `outcome'") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		restore
		} 
		local ++irow 
		}
}

log close


capture log close
log using "Sipila, UKB infection burden $S_DATE.log", replace

*Trim the data
use "SipilaELSA.dta", clear
count
drop if exitdate<=entrydate
count
drop if ses==.
count

tab dementia
recode dementia .=0
tab dementia 

foreach idg in anyinf bactinf viralinf hiv {
replace ensi`idg'pvm=. if ensi`idg'pvm>=exitdate
replace any`idg' = 0 if ensi`idg'pvm==.
sum ensi`idg'pvm
tab any`idg'
}

*Create variables that have identical names with those used in the primary analysis

capture gen diabetes_combpvm = diabetesdate
capture gen hypert_combpvm = hypertensiondate
capture gen tupakka = smoking
*alcocl does not change
*ses does not change
capture gen entrypvm = entrydate
capture gen exitpvm = exitdate
capture gen ensidementia = dementia
capture recode ensidementia 0=.
capture gen syntpvm = birthday
capture gen supu = sex
*dementiatype does not change

capture drop ensiinfcountpvm
gen ensiinfcountpvm = ensianyinfpvm 
		

save "SipilaELSA.dta", replace


**** Analysis *****


***** SINGLE VS MULTIPLE INFECTIONS *******		

count
tab toinenanyinf, mis
tab kolmasanyinf, mis

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm & ensitoinenanyinfpvm>=exitpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm  & ensikolmasanyinfpvm>=exitpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab ensidementia if ensianyinfpvm<exitpvm & ensitoinenanyinfpvm>=exitpvm
tab ensidementia if ensitoinenanyinfpvm<exitpvm & ensikolmasanyinfpvm>=exitpvm
tab ensidementia if ensikolmasanyinfpvm<exitpvm

sum ensianyinfpvm
replace ensianyinfpvm = . if ensianyinfpvm>=exitpvm
sum ensianyinfpvm		

sum ensitoinenanyinfpvm
replace ensitoinenanyinfpvm = . if ensitoinenanyinfpvm>=exitpvm
sum ensitoinenanyinfpvm		

sum ensikolmasanyinfpvm
replace ensikolmasanyinfpvm = . if ensikolmasanyinfpvm>=exitpvm
sum ensikolmasanyinfpvm		

replace toinenanyinf = . if ensitoinenanyinfpvm==.
replace kolmasanyinf = . if ensikolmasanyinfpvm==. 

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'
		preserve
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		stsplit first_split, at(0) after(time=ensianyinfpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anyanyinf==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=ensitoinenanyinfpvm)
		replace multiinf_timedep = 2 if toinenanyinf==1 & second_split==0
		tab multiinf_timedep
		stsplit third_split, at(0) after(time=ensikolmasanyinfpvm)
		replace multiinf_timedep = 3 if kolmasanyinf==1 & third_split==0
		tab multiinf_timedep
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/3 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 & third_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 & third_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)
		
		
		restore

		
***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local irow 37
local prow = `irow'

		preserve
		drop if ses==.
		egen ensibactviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
		gen toinenbactviralpvm = ensibactviralpvm if ensibactinfpvm==ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensibactinfpvm if ensibactinfpvm>ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensiviralinfpvm if ensiviralinfpvm>ensibactinfpvm & ensiviralinfpvm!=.
		gen anybactviral = 1 if ensibactviralpvm!=.
		gen toinenbactviral = 1 if toinenbactviralpvm!=.
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral
		tab anybactviral if toinenbactviral!=1
		tab toinenbactviral
		tab anybactviral ensidementia
		tab anybactviral ensidementia if toinenbactviral!=1
		tab toinenbactviral ensidementia
		
		stsplit first_split, at(0) after(time=ensibactviralpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anybactviral==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=toinenbactviralpvm)
		replace multiinf_timedep = 2 if toinenbactviral==1 & second_split==0
		tab multiinf_timedep
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		
		restore
		
***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

count
tab anyinfcount, mis
tab anyinfcount ensidementia, mis
replace anyinfcount = . if ensianyinfpvm>=exitpvm
tab anyinfcount, mis
tab anyinfcount ensidementia, mis

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of simultaneous infections""'
		
quietly: local irow = 23
local prow = `irow'
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		replace any`idg' = 2 if any`idg'>2 & any`idg'<.
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = any`idg' if inlist(any`idg',1,2) & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab `idg'_timedep, matcell(x), if `idg'_split!=-1
		matrix list x
		tab `idg'_timedep ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections
count
tab anyanycns
recode anyanycns .=0
tab anyanycns
tab anynoncns
recode anynoncns .=0
tab anynoncns
replace anyanycns=0 if ensianycnspvm>=exitpvm
replace anynoncns=0 if ensinoncnspvm>=exitpvm
replace ensianycnspvm=. if ensianycnspvm>=exitpvm
replace ensinoncnspvm=. if ensinoncnspvm>=exitpvm
tab anyanycns
tab anynoncns
tab anyanycns ensidementia
tab anynoncns ensidementia


local dg anycns noncns
local n_dg: word count `dg'
local dg_name `""Any CNS infection" "Any extra-CNS infection""'
quietly: local irow = 28
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "UKB replication analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		*Dementia occurring from year 5 onwards
		gen ensi`idg'pvm5v = ensi`idg'pvm + round(5*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm5v = entrypvm if ensi`idg'pvm5v<entrypvm
		replace _t0 = (ensi`idg'pvm5v-syntpvm)/365.25 if ensi`idg'pvm5v!=.
		replace _t0 = _t0 + 5 if ensi`idg'pvm5v==.
		drop if _t0>=_t
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_5y_exclusion")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("   Dementia occurring more than 5 years after the hospitalisation")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
		
*Difference between CNS infections and extra-CNS infections

	*CNS infection vs non cns infection 
	capture gen byte anycnsnotcns = 1 if ensianycnspvm==ensianyinfpvm & anyanycns==1 & ensianyinfpvm<exitpvm
	replace anycnsnotcns = 0 if ensinoncnspvm==ensianyinfpvm & anycnsnotcns!=1 & anynoncns==1 & ensianyinfpvm<exitpvm
	capture gen ensicnsnotcnspvm = ensianyinfpvm if ensianyinfpvm < exitpvm
	count
	sum ensicnsnotcnspvm
	tab anycnsnotcns
	tab anyanycns
	tab anynoncns
	tab anycnsnotcns ensidementia
	tab anyanycns ensidementia
	tab anynoncns ensidementia
	
	
local dg cnsnotcns
local n_dg: word count `dg'
local dg_name `""CNS vs extra-CNS infection""'
quietly: local irow = 34
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		gen byte anyeither = ensi`idg'pvm < exitpvm
		tab any`idg'
		tab anyeither
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if anyeither==1 & `idg'_split==0
		tab `idg'_timedep
		tab anyeither if `idg'_split!=-1

		stdescribe
		quietly: local failures = `r(N_fail)'

		tab anyeither if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`idg'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep supu ses
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		local ++irow		
		restore
		}
		





**** Analysis ALZHEIMER'S DISEASE *****



***** SINGLE VS MULTIPLE INFECTIONS *******		

count
tab toinenanyinf, mis
tab kolmasanyinf, mis

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm & ensitoinenanyinfpvm>=exitpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm  & ensikolmasanyinfpvm>=exitpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab dementiatype if ensianyinfpvm<exitpvm & ensitoinenanyinfpvm>=exitpvm
tab dementiatype if ensitoinenanyinfpvm<exitpvm & ensikolmasanyinfpvm>=exitpvm
tab dementiatype if ensikolmasanyinfpvm<exitpvm

sum ensianyinfpvm
replace ensianyinfpvm = . if ensianyinfpvm>=exitpvm
sum ensianyinfpvm		

sum ensitoinenanyinfpvm
replace ensitoinenanyinfpvm = . if ensitoinenanyinfpvm>=exitpvm
sum ensitoinenanyinfpvm		

sum ensikolmasanyinfpvm
replace ensikolmasanyinfpvm = . if ensikolmasanyinfpvm>=exitpvm
sum ensikolmasanyinfpvm		

replace toinenanyinf = . if ensitoinenanyinfpvm==.
replace kolmasanyinf = . if ensikolmasanyinfpvm==. 

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'
		preserve
		stset exitpvm, id(id) failure(dementiatype==1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		stsplit first_split, at(0) after(time=ensianyinfpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anyanyinf==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=ensitoinenanyinfpvm)
		replace multiinf_timedep = 2 if toinenanyinf==1 & second_split==0
		tab multiinf_timedep
		stsplit third_split, at(0) after(time=ensikolmasanyinfpvm)
		replace multiinf_timedep = 3 if kolmasanyinf==1 & third_split==0
		tab multiinf_timedep
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/3 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 & third_split!=-1
		matrix list x
		tab multiinf_timedep dementiatype, matcell(y), if first_split!=-1 & second_split!=-1 & third_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)
		
		restore

		
***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local irow 37
local prow = `irow'

		preserve
		drop if ses==.
		egen ensibactviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
		gen toinenbactviralpvm = ensibactviralpvm if ensibactinfpvm==ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensibactinfpvm if ensibactinfpvm>ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensiviralinfpvm if ensiviralinfpvm>ensibactinfpvm & ensiviralinfpvm!=.
		gen anybactviral = 1 if ensibactviralpvm!=.
		gen toinenbactviral = 1 if toinenbactviralpvm!=.
		tab anybactviral
		tab anybactviral if toinenbactviral!=1
		tab toinenbactviral
		tab anybactviral dementiatype
		tab anybactviral dementiatype if toinenbactviral!=1
		tab toinenbactviral dementiatype
		
		quietly: stset exitpvm, id(id) failure(dementiatype==1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral
		tab toinenbactviral
		
		stsplit first_split, at(0) after(time=ensibactviralpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anybactviral==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=toinenbactviralpvm)
		replace multiinf_timedep = 2 if toinenbactviral==1 & second_split==0
		tab multiinf_timedep
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		
		restore
		
		

***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

count
tab anyinfcount, mis
tab anyinfcount dementiatype, mis
replace anyinfcount = . if ensianyinfpvm>=exitpvm
tab anyinfcount, mis
tab anyinfcount dementiatype, mis

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of simultaneous infections""'
		
quietly: local irow = 23
local prow = `irow'
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(dementiatype==1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		replace any`idg' = 2 if any`idg'>2 & any`idg'<.
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = any`idg' if inlist(any`idg',1,2) & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab `idg'_timedep, matcell(x), if `idg'_split!=-1
		matrix list x
		tab `idg'_timedep dementiatype, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)		
		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections
count
tab anyanycns
recode anyanycns .=0
tab anyanycns
tab anynoncns
recode anynoncns .=0
tab anynoncns
replace anyanycns=0 if ensianycnspvm>=exitpvm
replace anynoncns=0 if ensinoncnspvm>=exitpvm
replace ensianycnspvm=. if ensianycnspvm>=exitpvm
replace ensinoncnspvm=. if ensinoncnspvm>=exitpvm
tab anyanycns
tab anynoncns
tab anyanycns dementiatype
tab anynoncns dementiatype


local dg anycns noncns
local n_dg: word count `dg'
local dg_name `""Any CNS infection" "Any extra-CNS infection""'
quietly: local irow = 28
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "UKB replication analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype==1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		tab any`idg' dementiatype, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		
		*Dementia occurring from year 5 onwards
		gen ensi`idg'pvm5v = ensi`idg'pvm + round(5*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm5v = entrypvm if ensi`idg'pvm5v<entrypvm
		replace _t0 = (ensi`idg'pvm5v-syntpvm)/365.25 if ensi`idg'pvm5v!=.
		replace _t0 = _t0 + 5 if ensi`idg'pvm5v==.
		drop if _t0>=_t
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_5y_exclusion")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("   Dementia occurring more than 5 years after the hospitalisation")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		tab any`idg' dementiatype, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
		
*Difference between CNS infections and extra-CNS infections

	*CNS infection vs non cns infection 
	capture gen byte anycnsnotcns = 1 if ensianycnspvm==ensianyinfpvm & anyanycns==1 & ensianyinfpvm<exitpvm
	replace anycnsnotcns = 0 if ensinoncnspvm==ensianyinfpvm & anycnsnotcns!=1 & anynoncns==1 & ensianyinfpvm<exitpvm
	capture gen ensicnsnotcnspvm = ensianyinfpvm if ensianyinfpvm < exitpvm
	count
	sum ensicnsnotcnspvm
	tab anycnsnotcns
	tab anyanycns
	tab anynoncns
	tab anycnsnotcns dementiatype
	tab anyanycns dementiatype
	tab anynoncns dementiatype
	
	
local dg cnsnotcns
local n_dg: word count `dg'
local dg_name `""CNS vs extra-CNS infection""'
quietly: local irow = 34
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype==1) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		gen byte anyeither = ensi`idg'pvm < exitpvm
		tab any`idg'
		tab anyeither
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if anyeither==1 & `idg'_split==0
		tab `idg'_timedep
		tab anyeither if `idg'_split!=-1

		stdescribe
		quietly: local failures = `r(N_fail)'

		tab anyeither if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`idg'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep supu ses
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		local ++irow		
		restore
		}
		
		

		

**** Analysis VASCULAR DEMENTIA *****



***** SINGLE VS MULTIPLE INFECTIONS *******		

count
tab toinenanyinf, mis
tab kolmasanyinf, mis

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab anyanyinf if ensianyinfpvm<exitpvm & ensianyinfpvm < ensitoinenanyinfpvm & ensitoinenanyinfpvm>=exitpvm 
tab toinenanyinf if ensitoinenanyinfpvm<exitpvm & ensitoinenanyinfpvm < ensikolmasanyinfpvm  & ensikolmasanyinfpvm>=exitpvm 
tab kolmasanyinf if ensikolmasanyinfpvm<exitpvm

tab dementiatype if ensianyinfpvm<exitpvm & ensitoinenanyinfpvm>=exitpvm
tab dementiatype if ensitoinenanyinfpvm<exitpvm & ensikolmasanyinfpvm>=exitpvm
tab dementiatype if ensikolmasanyinfpvm<exitpvm

sum ensianyinfpvm
replace ensianyinfpvm = . if ensianyinfpvm>=exitpvm
sum ensianyinfpvm		

sum ensitoinenanyinfpvm
replace ensitoinenanyinfpvm = . if ensitoinenanyinfpvm>=exitpvm
sum ensitoinenanyinfpvm		

sum ensikolmasanyinfpvm
replace ensikolmasanyinfpvm = . if ensikolmasanyinfpvm>=exitpvm
sum ensikolmasanyinfpvm		

replace toinenanyinf = . if ensitoinenanyinfpvm==.
replace kolmasanyinf = . if ensikolmasanyinfpvm==. 

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'
		preserve
		stset exitpvm, id(id) failure(dementiatype==4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		stsplit first_split, at(0) after(time=ensianyinfpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anyanyinf==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=ensitoinenanyinfpvm)
		replace multiinf_timedep = 2 if toinenanyinf==1 & second_split==0
		tab multiinf_timedep
		stsplit third_split, at(0) after(time=ensikolmasanyinfpvm)
		replace multiinf_timedep = 3 if kolmasanyinf==1 & third_split==0
		tab multiinf_timedep
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/3 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1 & third_split!=-1
		matrix list x
		tab multiinf_timedep dementiatype, matcell(y), if first_split!=-1 & second_split!=-1 & third_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)
		
		restore

***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local irow 37
local prow = `irow'

		preserve
		drop if ses==.
		egen ensibactviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
		gen toinenbactviralpvm = ensibactviralpvm if ensibactinfpvm==ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensibactinfpvm if ensibactinfpvm>ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensiviralinfpvm if ensiviralinfpvm>ensibactinfpvm & ensiviralinfpvm!=.
		gen anybactviral = 1 if ensibactviralpvm!=.
		gen toinenbactviral = 1 if toinenbactviralpvm!=.
		tab anybactviral
		tab anybactviral if toinenbactviral!=1
		tab toinenbactviral
		tab anybactviral dementiatype
		tab anybactviral dementiatype if toinenbactviral!=1
		tab toinenbactviral dementiatype

		quietly: stset exitpvm, id(id) failure(dementiatype==4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral
		tab toinenbactviral
		
		stsplit first_split, at(0) after(time=ensibactviralpvm)
		gen multiinf_timedep =0
		tab multiinf_timedep
		replace multiinf_timedep = 1 if anybactviral==1 & first_split==0
		tab multiinf_timedep
		stsplit second_split, at(0) after(time=toinenbactviralpvm)
		replace multiinf_timedep = 2 if toinenbactviral==1 & second_split==0
		tab multiinf_timedep
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf_timedep, matcell(x), if first_split!=-1 & second_split!=-1
		matrix list x
		tab multiinf_timedep ensidementia, matcell(y), if first_split!=-1 & second_split!=-1 &_d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		
		restore
		

***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

count
tab anyinfcount, mis
tab anyinfcount dementiatype, mis
replace anyinfcount = . if ensianyinfpvm>=exitpvm
tab anyinfcount, mis
tab anyinfcount dementiatype, mis

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of simultaneous infections""'
		
quietly: local irow = 23
local prow = `irow'
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(dementiatype==4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		replace any`idg' = 2 if any`idg'>2 & any`idg'<.
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = any`idg' if inlist(any`idg',1,2) & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `hr2' = m[1,3]
		scalar `ll2' = m[5,3]
		scalar `ul2' = m[6,3]
		scalar `hr3' = m[1,4]
		scalar `ll3' = m[5,4]
		scalar `ul3' = m[6,4]
		
		forvalues i =1/2 {
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab `idg'_timedep, matcell(x), if `idg'_split!=-1
		matrix list x
		tab `idg'_timedep dementiatype, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)		
		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections
count
tab anyanycns
recode anyanycns .=0
tab anyanycns
tab anynoncns
recode anynoncns .=0
tab anynoncns
replace anyanycns=0 if ensianycnspvm>=exitpvm
replace anynoncns=0 if ensinoncnspvm>=exitpvm
replace ensianycnspvm=. if ensianycnspvm>=exitpvm
replace ensinoncnspvm=. if ensinoncnspvm>=exitpvm
tab anyanycns
tab anynoncns
tab anyanycns dementiatype
tab anynoncns dementiatype


local dg anycns noncns
local n_dg: word count `dg'
local dg_name `""Any CNS infection" "Any extra-CNS infection""'
quietly: local irow = 28
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "UKB replication analysis, normal model and 5-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype==4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if any`idg'==1 & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox `idg'_timedep i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if `idg'_split!=-1
		matrix list x
		tab any`idg' dementiatype, matcell(y), if `idg'_split!=-1 & _d==1
		matrix list y
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		local ++irow 
		restore
		}
		
*Difference between CNS infections and extra-CNS infections

	*CNS infection vs non cns infection 
	capture gen byte anycnsnotcns = 1 if ensianycnspvm==ensianyinfpvm & anyanycns==1 & ensianyinfpvm<exitpvm
	replace anycnsnotcns = 0 if ensinoncnspvm==ensianyinfpvm & anycnsnotcns!=1 & anynoncns==1 & ensianyinfpvm<exitpvm
	capture gen ensicnsnotcnspvm = ensianyinfpvm if ensianyinfpvm < exitpvm
	count
	sum ensicnsnotcnspvm
	tab anycnsnotcns
	tab anyanycns
	tab anynoncns
	tab anycnsnotcns dementiatype
	tab anyanycns dementiatype
	tab anynoncns dementiatype
	
	
local dg cnsnotcns
local n_dg: word count `dg'
local dg_name `""CNS vs extra-CNS infection""'
quietly: local irow = 34
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype==4) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		
		gen byte anyeither = ensi`idg'pvm < exitpvm
		tab any`idg'
		tab anyeither
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if anyeither==1 & `idg'_split==0
		tab `idg'_timedep
		tab anyeither if `idg'_split!=-1

		stdescribe
		quietly: local failures = `r(N_fail)'

		tab anyeither if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`idg'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep supu ses
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		local ++irow		
		restore
		}
		
		
log close
