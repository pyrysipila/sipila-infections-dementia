version 16

**** PRIMARY ANALYSIS (FINNISH MULTICOHORT SAMPLE) ********



capture log close
log using "Jakauma-datan luominen v2.log", append

*Create the indicator variables for start of follow up for each infectious disease

use "Pooled infections_trimmattu", clear

*Random order
*seed from random.org from range 1 - 1 000 000 000
set rng mt64s
set rngstream 1
set seed 141476928 

gen double shuffle1 = runiform()
gen double shuffle2 = runiform()

sort shuffle1 shuffle2

*New id that is in the random order
gen shuffle_id = _n

*Difference between entry and time of infection
gen entry_to_inf = (ensianyinfpvm-entrypvm)/365.25
replace entry_to_inf = 0 if entry_to_inf < 0 
bysort cohort supu: sum entry_to_inf if entry_to_inf>0
*hist entry_to_inf
sum entry_to_inf
capture drop entry_to_inf_cut
egen entry_to_inf_cut = cut(entry_to_inf), at(0 0.0009765625 5(5)20 100)
tab entry_to_inf_cut
by cohort supu: tab entry_to_inf_cut
capture drop age_cut
egen age_cut = cut(age), at(18 30(10)60 100)
tab age_cut
bysort cohort supu age_cut: tab entry_to_inf_cut

levelsof cohort, local(cohorts)
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_cut, local(entry_lags)
local i 0
foreach cohort of local cohorts {
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0	
			foreach lag of local entry_lags {
				local ++i
			
				quietly: sum entry_to_inf if cohort==`cohort' & age_cut==`age' & supu==`sex' & entry_to_inf_cut==`lag'				
				
				if `i'==1 & `r(N)' == 0 {
				matrix mean_lag_`cohort'_`age'_`sex' = 0
				}
				else if `i'==1 {
				matrix mean_lag_`cohort'_`age'_`sex' = `r(mean)'
				}
				else if `r(N)' == 0 {
				matrix mean_lag_`cohort'_`age'_`sex' = mean_lag_`cohort'_`age'_`sex' \ 0
				}
				else {
				matrix mean_lag_`cohort'_`age'_`sex' = mean_lag_`cohort'_`age'_`sex' \ `r(mean)'
				}
				}
				matrix list mean_lag_`cohort'_`age'_`sex'				
				}
				}
				}

*Start of follow up when the infection occurred more than 10 years before dementia
gen entry_to_inf_10y = 10 + (ensianyinfpvm-entrypvm)/365.25
replace entry_to_inf_10y = 0 if entry_to_inf_10y < 0 
*hist entry_to_inf_10y
sum entry_to_inf_10y
capture drop entry_to_inf_10y_cut
egen entry_to_inf_10y_cut = cut(entry_to_inf_10y), at(0 0.0009765625 5(5)30 100)
tab entry_to_inf_10y_cut
sum entry_to_inf_10y if entry_to_inf_10y_cut==0

levelsof cohort, local(cohorts)
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_10y_cut, local(entry_lags)
local i 0
foreach cohort of local cohorts {
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0	
			foreach lag of local entry_lags {
				local ++i
			
				quietly: sum entry_to_inf_10y if cohort==`cohort' & age_cut==`age' & supu==`sex' & entry_to_inf_10y_cut==`lag'				
				
				if `i'==1 & `r(N)' == 0 {
				matrix mean_lag_10y_`cohort'_`age'_`sex' = 0
				}
				else if `i'==1 {
				matrix mean_lag_10y_`cohort'_`age'_`sex' = `r(mean)'
				}
				else if `r(N)' == 0 {
				matrix mean_lag_10y_`cohort'_`age'_`sex' = mean_lag_10y_`cohort'_`age'_`sex' \ 0
				}
				else {
				matrix mean_lag_10y_`cohort'_`age'_`sex' = mean_lag_10y_`cohort'_`age'_`sex' \ `r(mean)'
				}
				}
				matrix list mean_lag_10y_`cohort'_`age'_`sex'				
				}
				}
				}

*Return random order
capture drop sub_order
bysort cohort age_cut supu anyanyinf (shuffle1 shuffle2): gen sub_order = _n

*Similar distribution of delays in start of follow up for those not exposed to infections
capture drop entry_lag
gen entry_lag = .
levelsof cohort, local(cohorts)
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_cut, local(entry_lags)
foreach cohort of local cohorts {
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0
			local cum = 0
			sum id if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==0				
			local N_control = `r(N)'
			sum entry_to_inf_cut if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==1				
			local N_`cohort'_`age'_`sex' = `r(N)'
			foreach lag of local entry_lags {
				di "cohort `cohort', age_group `age', sex `sex'"
				local ++i
				di `i'
				sum entry_to_inf_cut if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==1 & entry_to_inf_cut==`lag'			
				local cum = `cum' + `r(N)'
				di `cum'
				di `N_`cohort'_`age'_`sex''
				di `N_control'
				di mean_lag_`cohort'_`age'_`sex'[`i',1]
				replace entry_lag = mean_lag_`cohort'_`age'_`sex'[`i',1] if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==0 & entry_lag==. & sub_order <= round((`cum' / `N_`cohort'_`age'_`sex'')*`N_control')
				}
				replace entry_lag = 0 if `N_`cohort'_`age'_`sex''==0 & entry_lag==. & anyanyinf==0
				}
				}
				}

*Delays for the 10y exclusion analysis
capture drop entry_lag_10y
gen entry_lag_10y = .
levelsof cohort, local(cohorts)
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_10y_cut, local(entry_lags)
foreach cohort of local cohorts {
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0
			local cum = 0
			sum id if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==0				
			local N_control = `r(N)'
			sum entry_to_inf_10y_cut if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==1				
			local N_`cohort'_`age'_`sex' = `r(N)'
			foreach lag of local entry_lags {
				di "cohort `cohort', age_group `age', sex `sex'"
				local ++i
				di `i'
				sum entry_to_inf_10y_cut if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==1 & entry_to_inf_10y_cut==`lag'			
				local cum = `cum' + `r(N)'
				di `cum'
				di `N_`cohort'_`age'_`sex''
				di `N_control'
				di mean_lag_10y_`cohort'_`age'_`sex'[`i',1]
				replace entry_lag_10y = mean_lag_10y_`cohort'_`age'_`sex'[`i',1] if cohort==`cohort' & age_cut==`age' & supu==`sex' & anyanyinf==0 & entry_lag_10y==. & sub_order <= round((`cum' / `N_`cohort'_`age'_`sex'')*`N_control')
				}
				replace entry_lag_10y = 0 if `N_`cohort'_`age'_`sex''==0 & entry_lag_10y==. & anyanyinf==0
				}
				}
				}
				
* Check  that the distribution are correct:				
sum entry_to_inf_10y
sum entry_lag_10y
levelsof cohort, local(cohorts)
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_10y_cut, local(entry_lags)
foreach cohort of local cohorts {
	foreach age of local age_groups {
		foreach sex of local sexes {
			sum entry_to_inf if cohort==`cohort' & age_cut==`age' & supu==`sex'
			sum entry_lag if cohort==`cohort' & age_cut==`age' & supu==`sex'
			sum entry_to_inf_10y if cohort==`cohort' & age_cut==`age' & supu==`sex'
			sum entry_lag_10y if cohort==`cohort' & age_cut==`age' & supu==`sex'
				}
				}
				}

				
save "jakauma v2.dta", replace

use "jakauma v2.dta", clear
*Update the infection variables so that they are coded missing if the participant does not have the index infection but has some other infection
local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral
foreach idg of local dg {
		tab any`idg'
		replace any`idg' = . if anyanyinf==1 & any`idg'==0 
		tab any`idg'
		
		egen lag_entry_`idg' = rowmax(entrypvm ensi`idg'pvm) if any`idg'==1
		replace lag_entry_`idg' = entrypvm + round(entry_lag*365.25) if any`idg'==0
		capture drop ensi`idg'pvm_10y
		gen ensi`idg'pvm_10y = ensi`idg'pvm + round(10*365.25)
		egen lag_entry_10y_`idg' = rowmax(entrypvm ensi`idg'pvm_10y) if any`idg'==1
		replace lag_entry_10y_`idg' = entrypvm + round(entry_lag_10y*365.25) if any`idg'==0
		}

save "Inf_dem_lag.dta", replace



****** Descriptive analysis and definition of variables *************
use "Inf_dem_lag.dta", replace

*Katsotaan eri infektioiden määrät
tab anyanyinf

tab anybactinf

tab anysysbact
tab anylocalbact

tab anybactsepsis
tab anybact_nonsepsi

tab anyextrac
tab anyintrac

tab anygramplus
tab anygramminus
tab anymycobact
tab anymycoplasma

tab anyviralinf
tab anyherpes
tab anypersviral
tab anyacuteviral

tab anyparasite
tab anymycose

*Piirretään histogrammi dementian iästä

stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 

gen dementia_age = int(_t)+0.5 if _d==1
histogram dementia_age, disc scheme(s1color) freq ///
	ysize(7.8) xsize(8.27) ///
	title(" ", size(medsmall) color(black)) ///
	graphregion(color(white)) ///
	xtitle("Age (years)", size(medsmall) color(black)) ///
	ytitle("Dementias, No.", size(medsmall) color(black))
graph export "eFigure4_new.emf", as(emf) replace
graph export "eFigure4_new.pdf", as(pdf) replace


		use "Inf_dem_lag.dta", replace
		merge 1:1 id using "Parkinson.dta"
		keep if _merge==3
		drop _merge


stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
	keep if _st==1
	count
	count if ensidementia==1
	
		label define supu 1"Men" 2"Women"
		label values supu supu
		label values diabetes yesno
		label values hypertensio yesno
		recode ensidementia .=0
		label values ensidementia yesno

		gen followup = _t-_t0
		
		drop age
		gen age = (lag_entry_anyinf - syntpvm)/365.25
		egen age_cat = cut(age), at(18,40,50,60,100)
		tab age_cat
				
		label define alcocl 0 "Non-drinker" 1 "Moderate" 2 "Intermediate" 3 "Heavy"
		label values alcocl alcocl
		
		foreach cm in hypert diabetes ihd cereb parkinson {
			gen `cm'_entry= `cm'_combpvm<=lag_entry_anyinf
			tab `cm'_entry
			}
		
		gen age_at_dementia = _t if _d==1
		
		keep if _st==1
		keep if anyanyinf!=.
		
		baselinetable ///
		age_cat(novarlabel afterhead("Age at entry, years")) ///
		supu(novarlabel afterhead("Sex")) ///
		ses(novarlabel afterhead("Education/Socioeconomic status")) ///	
		hypert_entry(novarlabel afterhead("Hypertension")) ///
		diabetes_entry(novarlabel afterhead("Diabetes mellitus")) ///
		ihd_entry(novarlabel afterhead("Ischaemic heart disease")) ///
		cereb_entry(novarlabel afterhead("Cerebrovascular disease")) ///
		parkinson_entry(novarlabel afterhead("Parkinson's disease")) ///
		tupakka(novarlabel afterhead("Smoking status")) ///
		alcocl(novarlabel afterhead("Alcohol drinking")) ///
		followup(cts novarlabel afterhead("Follow-up, median (IQR), years")) ///
		ensidementia (novarlabel afterhead("Dementia")) ///
		age_at_dementia(cts novarlabel afterhead("Age at dementia, median (IQR), years")) ///
		, by(cohort, total) ctsvartab(p50 (p25-p75)) missing ///
		exportexcel(eTable_3_25012021, cell(A6) replace)

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
tab bmi_who			//								body mass index (normal weight, overweight, obese)	
sum hypert_combpvm	//								date of hypertension diagnosis
sum diabetes_combpvm //								date of diabetes diagnosis
sum ihd_combpvm		//								date of ischaemic heart disease diagnosis
sum cereb_combpvm	//								date of cerebrovascular disease diagnosis
sum parkinson_combpvm	//								date of cerebrovascular disease diagnosis



*** PRIMARY ANALYSES USING THREE FINNISH COHORTS ***

*Incidence of infections

capture log close
log using "Incidence of infections and dementia.log", replace

use "Inf_dem_lag.dta", clear


capture log close
log using "FIGURE 2 dist.log", replace
stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
keep if _st==1
gen inf_at_entry = ensianyinfpvm<=entrypvm
tab inf_at_entry
gen byte inf_at_fu = ensianyinfpvm<exitpvm & inf_at_entry==0
tab inf_at_fu
egen infexit = rowmin(ensianyinfpvm exitpvm)
stset infexit, id(id) failure(inf_at_fu) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
stptime, per(100000)



*FIGURE 2.

use "Inf_dem_lag.dta", clear

local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral
local row 7 10 13 14 16 17 20 21 24 25 28 31 32 33
local n_dg: word count `dg'
local dg_name `""Any infectious disease vs no infection" "Any bacterial infection vs no infection" "Potentially invasive bacterial infection vs no infection"  "Localised bacterial infection vs no infection" "Bacterial infection with sepsis vs no infection" "Bacterial infection without sepsis vs no infection" "Extracellular bacterial infection vs no infection" "Intracellular bacterial infection vs no infection" "Gram-positive bacterial infection vs no infection" "Gram-negative bacterial infection vs no infection" "Any viral infection vs no infection" "Herpesvirus (persistent) infection vs no infection" "Other persistent viral infection vs no infection" "Acute viral infection vs no infection""'
forvalues i=1/`n_dg' {		
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local irow : word `i' of `row'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort, age as time-scale"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		preserve
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		stdescribe
		stdescribe if any`idg'==0
		stdescribe if any`idg'==1
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm lag_entry_`idg' lag_entry_10y_`idg' _st _d _origin _t _t0 ses
		keep if _st==1
		tab any`idg'
		tab any`idg' ensidementia
		
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local irow = `irow' + 29 
		
		*Dementia occurring from year 10 onwards
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_`idg') scale(365.25) 
		stdescribe
		stdescribe if any`idg'==0
		stdescribe if any`idg'==1

		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore	
		}
log close

***Age-standardised incidence of dementia

*()Information form FIGURE 2 that goes to appendix)
capture log close
log using "FIGURE 2 append.log", replace

use "Inf_dem_lag.dta", clear

local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral
local row 7 10 13 14 16 17 20 21 24 25 28 31 32 33
local n_dg: word count `dg'
local dg_name `""Any infectious disease vs no infection" "Any bacterial infection vs no infection" "Potentially invasive bacterial infection vs no infection"  "Localised bacterial infection vs no infection" "Bacterial infection with sepsis vs no infection" "Bacterial infection without sepsis vs no infection" "Extracellular bacterial infection vs no infection" "Intracellular bacterial infection vs no infection" "Gram-positive bacterial infection vs no infection" "Gram-negative bacterial infection vs no infection" "Any viral infection vs no infection" "Herpesvirus (persistent) infection vs no infection" "Other persistent viral infection vs no infection" "Acute viral infection vs no infection""'
forvalues i=1/`n_dg' {		
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local irow : word `i' of `row'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort, age as time-scale"
		display ""
		display "Pooled IPD-analysis, normal model and 10-year exclusion"
		display ""
		display "$S_TIME  $S_DATE"
		display ""

		if "`idg'" == "anyinf" {
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		keep if _st==1
		*define scalars
		tempname base_pt_18 base_pt_30 base_pt_35 base_pt_40 base_pt_45 base_pt_50 base_pt_55 base_pt_60 ///
		base_pt_65 base_pt_70 base_pt_75 base_pt_80 base_pt_85
		tempname base_rate_18_10y base_rate_30_10y base_rate_35_10y base_rate_40_10y base_rate_45_10y ///
		base_rate_50_10y base_rate_55_10y base_rate_60_10y base_rate_65_10y base_rate_70_10y ///
		base_rate_75_10y base_rate_80_10y base_rate_85_10y		
		tempname inf18 inf30 inf35 inf40 inf45 inf50 inf55 inf60 inf65 inf70 inf75 inf80 inf85
		tempname ptime_base incid_base stand_incid_base_10y stand_incid
		stptime, per(100000), if anyanyinf==0
		scalar `ptime_base' = r(ptime)
		scalar `incid_base' = r(rate)
		}
		
		preserve
		
		if "`idg'" != "anyinf" {
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		keep if _st==1
		}
		
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		sum _st if _st==1
		local N = `r(N)'
		sum _d if _d==1
		local failures = `r(N)'
		
		*Store person-times and rates to scalars
		stsplit age_split, at(18 30(5)85)
		tab age_split
		levelsof age_split, local(age_groups)
		di `age_groups'
		foreach age of local age_groups {
			capture noisily stptime, per(100000), if age_split==`age' & anyanyinf==0
			if _rc == 0 { // _rc is the return code of the command following capture. 0 means that there was no error (e.g., due to no observations)).
			scalar `base_pt_`age'' = r(ptime)
			}
			else {
			scalar `base_pt_`age'' = 0
			}
			capture noisily stptime, per(100000), if age_split==`age' & anyanyinf==1
			if _rc == 0  {
			scalar `inf`age'' = r(rate)
			}
			else {
			scalar `inf`age'' = 0
			}
		}
		*Compute age-standardised incidence rate (per 100 000 person years)
		scalar `stand_incid' = (`base_pt_18'*`inf18' + `base_pt_30'*`inf30' + `base_pt_35'*`inf35' + `base_pt_40'*`inf40' + `base_pt_45'*`inf45' ///
		+ `base_pt_50'*`inf50' + `base_pt_55'*`inf55' + `base_pt_60'*`inf60' + `base_pt_65'*`inf65' + `base_pt_70'*`inf70' + `base_pt_75'*`inf75' ///
		+ `base_pt_80'*`inf80' + `base_pt_85'*`inf85') / `ptime_base'
		di "Age-standardised incidence rate is " `stand_incid'
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 append") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
			quietly: putexcel K`irow'=(`N'), nformat(number)
			quietly: putexcel L`irow'=(`failures'), nformat(number)
			quietly: putexcel G`irow'=(`incid_base'), nformat(number_d2)
			quietly: putexcel H`irow'=(`stand_incid'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")

		restore
		
		*Dementia occurring from year 10 onwards
		local irow = `irow' + 29 
		
		preserve
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_`idg') scale(365.25) 
		keep if _st==1

		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		sum _st if _st==1
		local N = `r(N)'
		sum _d if _d==1
		local failures = `r(N)'
		
		*Store person-times and rates to scalars
		stsplit age_split, at(18 30(5)85)
		levelsof age_split, local(age_groups)
		foreach age of local age_groups {
			capture noisily stptime, per(100000), if age_split==`age' & anyanyinf==1
			if _rc == 0  {
			scalar `inf`age'' = r(rate)
			}
			else {
			scalar `inf`age'' = 0
			}
			}
			
		if "`idg'" == "anyinf" {
			levelsof age_split, local(age_groups)
			foreach age of local age_groups {
				capture noisily stptime, per(100000), if age_split==`age' & anyanyinf==0
				if _rc == 0  {
				scalar `base_rate_`age'_10y' = r(rate)
				}
				else {
				scalar `base_rate_`age'_10y' = 0
				}
			}
			}
		
		if "`idg'" == "anyinf" {
			scalar `stand_incid_base_10y' = (`base_pt_18'*`base_rate_18_10y' + `base_pt_30'*`base_rate_30_10y' + `base_pt_35'*`base_rate_35_10y' + ///
			`base_pt_40'*`base_rate_40_10y' + `base_pt_45'*`base_rate_45_10y' + `base_pt_50'*`base_rate_50_10y' ///
			+ `base_pt_55'*`base_rate_55_10y' + `base_pt_60'*`base_rate_60_10y' + `base_pt_65'*`base_rate_65_10y' + `base_pt_70'*`base_rate_70_10y' + ///
			`base_pt_75'*`base_rate_75_10y' + `base_pt_80'*`base_rate_80_10y' + `base_pt_85'*`base_rate_85_10y') / `ptime_base'
		
		di "Age-standardised incidence for comparison group in 10-year lag analysis is " `stand_incid_base_10y'
		
		}
		
		*Compute age-standardised incidence rate (per 100 000 person years)
		scalar `stand_incid' = (`base_pt_18'*`inf18' + `base_pt_30'*`inf30' + `base_pt_35'*`inf35' + `base_pt_40'*`inf40' + `base_pt_45'*`inf45' ///
		+ `base_pt_50'*`inf50' + `base_pt_55'*`inf55' + `base_pt_60'*`inf60' + `base_pt_65'*`inf65' + `base_pt_70'*`inf70' + `base_pt_75'*`inf75' ///
		+ `base_pt_80'*`inf80' + `base_pt_85'*`inf85') / `ptime_base'
		di "Age-standardised incidence rate is " `stand_incid' " in 10-year lag analysis"
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 append") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
			quietly: putexcel K`irow'=(`N'), nformat(number)
			quietly: putexcel L`irow'=(`failures'), nformat(number)
			quietly: putexcel G`irow'=(`stand_incid_base_10y'), nformat(number_d2)
			quietly: putexcel H`irow'=(`stand_incid'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")

		restore				
			}
log close

*** CUMULATIVE HAZARD *****

use "Inf_dem_lag.dta", clear
stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		sts graph, na by(anyanyinf) ///
		title("Any infectious disease" "vs no infection" " " ///
		, size(medsmall) color(black))  xti("Age (years)") yti("Nelson-Aalen cumulative hazard estimate", height(5)) ///
		legend(order(1 "No infection"  2 "Infection")) scheme(s1color) ///
		plot1opts(lcolor(blue)) plot2opts(lcolor(red)) xscale(range(18 100)) xlabel(20(20)100) ///
		risktable(20 40 60 80 100, size(small) failevents order(1 "No infection" 2 "Infection") ///
		title("N at risk (dementias)", justification(left) size(small)) rowtitle(, justification(left)))
		graph export "age_NA_anyinf.emf", as(emf) replace





***eFIGURE 8 ****
*(was previously figure 3)

capture log close
log using "FIGURE 3 dist.log", replace

use "Inf_dem_lag.dta", clear

*Use this first listed dementia diagnosis if there are several on the same hospitalisation (primary diagnosis instead of secondary and 1st secondary instead of 2nd secondary and so on)
drop dementiatyyppi
gen dementiatyyppi = 1 if dementiatyyppi_laaja=="AD"
replace dementiatyyppi = 2 if inlist(dementiatyyppi_laaja,"FTD","PD","VD","other","unspecified")

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
		preserve
		drop if ses==.
		quietly: stset exitpvm, id(id) failure(dementiatyyppi == 1) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses dementiatyyppi
		tab any`idg' if _st==1
		
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3 dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'_AD")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3 dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		restore
		
		*Non-Alzheimer's disease dementia
		preserve
		quietly: stset exitpvm, id(id) failure(dementiatyyppi == 2 3 4) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses dementiatyyppi
		tab any`idg' if _st==1
		
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3 dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'_nonAD")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 3 dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore	
		}
log close

	
*FIGURE sensitivity, acute and chronic infections
capture log close
log using "Sensitivity, acute, chronic, coinfectionsm, dist.log", replace
local dg nonpers anypers bactviralcoinf
local n_dg: word count `dg'
local dg_name `""Any acute infectious disease" "Any chronic infectious disease""'

		quietly: use "Inf_dem_lag.dta", clear
		merge 1:1 id using "Coinfections_trimmattu"
		drop _merge
		foreach idg of local dg {
			tab any`idg'
			replace any`idg' = . if any`idg'== 0 & anyanyinf == 1
			tab any`idg'
			egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
			
			gen ensi`idg'pvm_10y = ensi`idg'pvm + round(10*365.25)
			egen lag_entry_10y_`idg' = rowmax(lag_entry_10y_anyinf ensi`idg'pvm_10y)
		}
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
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		tab any`idg' ensidementia if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		restore
		
		*Dementia occurring from year 10 onwards
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		tab any`idg' ensidementia if _st==1

		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
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

		quietly: use "Inf_dem_lag.dta", clear
		merge 1:1 id using "Coinfections_trimmattu"
		drop _merge
		drop if ses==.
		egen lag_entry_multiinf = rowmax(lag_entry_anyinf ensitoinenanyinfpvm ensikolmasanyinfpvm)
		egen multiinfpvm = rowmax(ensianyinfpvm ensitoinenanyinfpvm ensikolmasanyinfpvm)
		gen multiinfpvm_10y = multiinfpvm +round(10*365.25)
		egen lag_entry_10y_multiinf = rowmax(lag_entry_10y_anyinf multiinfpvm_10y)
		
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_multiinf) scale(365.25) 
		keep if _st==1
		quietly: keep id anyanyinf supu cohort exitpvm ensidementia syntpvm entrypvm ensi*anyinfpvm *toinen* *kolmas*  _st _d _origin _t _t0 ses
		sum ensi*anyinfpvm
		
		gen byte multiinf = 0
		tab multiinf
		replace multiinf = 1 if anyanyinf==1
		tab multiinf
		replace multiinf = 2 if ensitoinenanyinfpvm!=.
		replace multiinf = 3 if ensikolmasanyinfpvm!=.
		tab multiinf
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3 N1 N2 N3 fail1 fail2 fail3
		forvalues i =1/3 {
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if inlist(multiinf,0,`i')
		quietly: matrix m = r(table)
		matrix list m
		scalar `N`i'' = e(N_sub)
		scalar `fail`i'' = e(N_fail)
		scalar `hr`i'' = m[1,2]
		scalar `ll`i'' = m[5,2]
		scalar `ul`i'' = m[6,2]
		}
		
		forvalues i =1/3 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf, matcell(x), if inlist(multiinf,0,`i') & _st==1
		matrix list x
		tab multiinf ensidementia, matcell(y), if inlist(multiinf,0,`i') & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		restore
		
		*Dementia occurring from year 10 onwards
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_multiinf) scale(365.25) 
		keep if _st==1
		quietly: keep id anyanyinf supu cohort exitpvm ensidementia syntpvm entrypvm ensi*anyinfpvm *toinen* *kolmas*  _st _d _origin _t _t0 ses
		sum ensi*anyinfpvm
		
		gen byte multiinf = 0
		tab multiinf
		replace multiinf = 1 if anyanyinf==1
		tab multiinf
		replace multiinf = 2 if ensitoinenanyinfpvm!=.
		replace multiinf = 3 if ensikolmasanyinfpvm!=.
		tab multiinf
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3 N1 N2 N3 fail1 fail2 fail3
		forvalues i =1/3 {
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if inlist(multiinf,0,`i')
		quietly: matrix m = r(table)
		matrix list m
		scalar `N`i'' = e(N_sub)
		scalar `fail`i'' = e(N_fail)
		scalar `hr`i'' = m[1,2]
		scalar `ll`i'' = m[5,2]
		scalar `ul`i'' = m[6,2]
		}
		
		local irow 17
		forvalues i =1/3 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("10-year lag, single vs multiple infections, infection No. `i'")
		tab multiinf, matcell(x), if inlist(multiinf,0,`i') & _st==1
		matrix list x
		tab multiinf ensidementia, matcell(y), if inlist(multiinf,0,`i') & _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		restore
				
log close

***** One or several different pathogens *******		

capture log close
log using "One vs at least two different pathogens.log", replace

local dg_name `""Multiple organisms""'
local irow 82
local prow = `irow'

		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.
		egen ensibactviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
		gen toinenbactviralpvm = ensibactviralpvm if ensibactinfpvm==ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensibactinfpvm if ensibactinfpvm>ensiviralinfpvm & ensibactinfpvm!=.
		replace toinenbactviralpvm = ensiviralinfpvm if ensiviralinfpvm>ensibactinfpvm & ensiviralinfpvm!=.
		gen anybactviral = 1 if ensibactviralpvm!=.
		gen toinenbactviral = 1 if toinenbactviralpvm!=.
		egen lag_entry_bactviral = rowmax(lag_entry_anyinf ensibactviralpvm toinenbactviralpvm)
		gen ensibactviralpvm_10y = ensibactviralpvm + round(10*365.25)
		gen toinenbactviralpvm_10y = toinenbactviralpvm + round(10*365.25)
		egen lag_entry_10y_bactviral = rowmax(lag_entry_10y_anyinf ensibactviralpvm_10y toinenbactviralpvm_10y)
		
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_bactviral) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses anyanyinf
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral if _st==1
		tab toinenbactviral if _st==1
		
		gen byte multiinf = 0 if anyanyinf==0
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 1 if anybactviral==1
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 2 if toinenbactviral==1
		tab multiinf
		tab multiinf if _st==1
		*browse id ensibactviralpvm toinenbactviralpvm multiinf _t0 _t
		
		quietly: stdescribe
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf, matcell(x), if _st==1 & inlist(multiinf,0,`i')
		matrix list x
		tab multiinf ensidementia, matcell(y), if _d==1 & inlist(multiinf,0,`i')
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		restore
		
		
		*Dementia occurring from year 10 onwards
		local i 1
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_bactviral) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses anyanyinf
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral if _st==1
		tab toinenbactviral if _st==1
		
		gen byte multiinf = 0 if anyanyinf==0
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 1 if anybactviral==1
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 2 if toinenbactviral==1
		tab multiinf
		tab multiinf if _st==1
		*browse id ensibactviralpvm toinenbactviralpvm multiinf _t0 _t
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Different micro-organisms after 10+ years, infection No. `i'")
		tab multiinf, matcell(x), if _st==1 & inlist(multiinf,0,`i')
		matrix list x
		tab multiinf ensidementia, matcell(y), if _d==1 & inlist(multiinf,0,`i')
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)
		restore
		

log close



***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

capture log close
log using "Number of infections.log", replace

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of infections""'

		quietly: use "Inf_dem_lag.dta", clear
		merge 1:1 id using "Infection count.dta"
		tab anyinfcount
		replace anyinfcount=0 if anyanyinf==0
		recode anyinfcount 3=2 4=2
		drop if ses==.
		gen ensiinfcountpvm = ensianyinfpvm 
		
quietly: local irow = 23
quietly: local prow = `irow'
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
		
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab any`idg', matcell(x), if _st==1 &  inlist(any`idg',0,`i')
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1 &  inlist(any`idg',0,`i')
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)		
		
		local irow 24
		restore
		
		*Dementia occurring from year 10 onwards
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_anyinf) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("10-year lag, No. of simultaneous infection diagnoses = `i'")
		tab any`idg', matcell(x), if _st==1 &  inlist(any`idg',0,`i')
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1 &  inlist(any`idg',0,`i')
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel Q`prow'=(`p'), nformat(number)		
		
		local irow 24
		restore		
		}	
log close		


*FIGURE sensitivity, CNS vs other infections
capture log close
log using "Sensitivity, CNS infection vs other infection.log", replace
local dg anycns noncns
local n_dg: word count `dg'
local dg_name `""Any CNS infection vs not" "Any extra-CNS infection vs not""'

		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.

foreach idg of local dg {
	tab any`idg'
	replace any`idg' = . if any`idg'!=1 & anyanyinf==1
	tab any`idg'
	egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
	gen ensi`idg'pvm_10y =  ensi`idg'pvm + round(10*365.25)
	egen lag_entry_10y_`idg' = rowmax(lag_entry_10y_anyinf ensi`idg'pvm_10y)
	}		
		
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
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		restore
		
		*Dementia occurring from year 10 onwards
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_`idg') scale(365.25) 
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow 
		restore
		}
log close


*FIGURE sensitivity, infections predisposed towards entering the CNS vs other infections
capture log close
log using "Sensitivity, predisposed towards CNS vs others.log", replace
local dg anycnspredisp noncnspredisp
local n_dg: word count `dg'
local dg_name `""Infection predisposed towards entering the CNS vs no infection" "Infection not predisposed towards entering the CNS vs no infection""'

		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.
	
		foreach idg of local dg {
			tab any`idg'
			replace any`idg' = . if any`idg'!=1 & anyanyinf==1
			tab any`idg'
			egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
			gen ensi`idg'pvm_10y =  ensi`idg'pvm + round(10*365.25)
			egen lag_entry_10y_`idg' = rowmax(lag_entry_10y_anyinf ensi`idg'pvm_10y)
			}		
		
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
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		restore
		
		*Dementia occurring from year 10 onwards
		preserve
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_`idg') scale(365.25) 
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("   `idg'_10y_exclusion")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("   Dementia occurring more than 10 years after the hospitalisation")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE sensitivity dist") modify
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

		quietly: use "Inf_dem_lag.dta", clear
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
		
		capture drop ensi`idg'pvm
		gen ensi`idg'pvm = ensiherpespvm if any`idg'==1
		
		tab any`idg'
		replace any`idg' = 0 if anyanyinf==0
		replace any`idg' = . if any`idg'!=1 & anyanyinf==1
		tab any`idg'
		egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		gen ensi`idg'pvm_10y =  ensi`idg'pvm + round(10*365.25)
		egen lag_entry_10y_`idg' = rowmax(lag_entry_10y_anyinf ensi`idg'pvm_10y)
				
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
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
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
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

		quietly: use "Inf_dem_lag.dta", clear
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
		
		recode any`idg' 0=.
		recode any`idg' .=0 if anyanyinf==0
		gen ensi`idg'pvm = ensigrampluspvm if any`idg'==1
		egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
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
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
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

		quietly: use "Inf_dem_lag.dta", clear
		merge 1:1 id using "gramminuskoodit.dta"
		tab anygramminus
		tab anysalmonella
		tab anylyme
		tab anyothergramminus
		tab anyothergramminussepsis
		replace anyothergramminus = 1 if anyothergramminussepsis==1
		tab anyothergramminus
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
		
		recode any`idg' 0=.
		recode any`idg' .=0 if anyanyinf==0
		gen ensi`idg'pvm = ensigramminuspvm if any`idg'==1
		egen lag_entry_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
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
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
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

		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.

		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
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
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25)
		}
		else if `j'== 3 {
		stset exitpvm, id(id) failure(detailedfail==2) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25)
		}
		
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses detailedfail lag_entry_`idg'
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		if `j'==1 {
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		}
		else if `j'== 2 {
		stcrreg any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort i.cohort, cl(cohort) compete(detailedfail==3)
		}
		else if `j'== 3 {
		stcrreg any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort i.cohort, cl(cohort) compete(detailedfail==1 3)
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
			quietly: putexcel D`irow'=("`idg_name', Fine & Gray model, competing outcome death")
			}
			else if `j'== 3 {
			quietly: putexcel D`irow'=("`idg_name', Fine & Gray model, competing outcome death + early-onset dementia")
			}
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if _d==1
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

		quietly: use "Inf_dem_lag.dta", clear
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
		
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm ses lag_entry_`idg'
		stset exitpvm, id(id) failure(ensidementia) origin(time lag_entry_`idg') scale(365.25) 
		gen age`idg' = (lag_entry_`idg'-syntpvm)/365.25
		gen age`idg'2 = age`idg'^2 
		
		tab any`idg' if _st==1
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.age`idg'#cohort c.age`idg'2#cohort, strata(cohort)
		matrix m = r(table)
		matrix list m
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.age`idg'#cohort c.age`idg'2#cohort, strata(cohort) tvc(any`idg') texp(ln(_t))
		matrix m2 = r(table)
		matrix list m2
		restore
		}
log close		

version 16

		use "Inf_dem_lag.dta", replace


*** eTABLE 4 ****
*(was previously TABLE 2) *
 
capture log close
log using "TABLE2.log", replace

use "Inf_dem_lag.dta", clear
quietly: drop if ses==.
quietly: gen hypertensio_combpvm = hypert_combpvm
merge 1:1 id using "Parkinson.dta"
keep if _merge==3
drop _merge

local dg anyinf
local eks `""no one" "hypertensio" "diabetes" "ihd" "cereb" "parkinson"'
local n_dg: word count `dg'
local n_eks: word count `eks' 
local dg_name `""Any infectious disease vs no infection""'
local eks_name `""no one" "hypertension" "diabetes" "ischaemic heart disease" "cerebrovascular disease" "Parkinson´s disease"' 
quietly: local irow = 8
	forvalues i=1/`n_dg' {
	forvalues j=1/`n_eks' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local ieks : word `j' of `eks'
		quietly: local ieks_name : word  `j' of `eks_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""	
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		if inlist("`ieks'","hypertensio", "diabetes", "ihd", "cereb", "parkinson") {
		drop if `ieks'_combpvm <= lag_entry_anyinf
		}
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		
		if "`ieks'"=="hypertensio" & "`idg'"=="anyinf" {
		local irow = 25
		}
		
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses `ieks'_comb
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Those with `ieks_name' excluded")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		restore
}
}

*Early- vs late-onset dementia (Gram-negative)
		use "Inf_dem_lag.dta", clear
		drop if ses==.
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		local idg anyinf
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		
		stsplit age65, at(64.999)
		recode age65 64.999=1
 
		preserve
		keep if age65==0
		local irow 12
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Early-onset dementia (onset before age 65)")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore
		
		preserve
		keep if age65==1
		local irow 13
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Late-onset dementia (onset at or after age 65)")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore
 
			
*Type of dementia
		use "Inf_dem_lag.dta", clear
tab dementiatyyppi_laaja anyanyinf if ensidementia==1, col chi
	local irow 15
	local j 1
	foreach type in AD FTD PD VD other unspecified  {
		preserve
		encode dementiatyyppi_laaja, gen(dementiatyyppi_laaja_nro)

		display ""
		display "`type'"
		display ""
		local idg gramminus
		drop if ses==.
		
		stset exitpvm, id(id) failure(dementiatyyppi_laaja_nro==`j') origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		
		local idg anyinf
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`type'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++j
		restore
		}

*Adjustments

local dg anyinf
local eks alcocl tupakka
local n_dg: word count `dg'
local n_eks: word count `eks' 
local dg_name `""Any infectious disease""'
local eks_name `""heavy drinking" "smoking""'
quietly: local irow = 31
	forvalues i=1/`n_dg' {
	forvalues j=1/`n_eks' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local ieks : word `j' of `eks'
		quietly: local ieks_name : word  `j' of `eks_name'
		display ""
		display `"`idg' followed by broad dementia, adjusted for sex and education and `eks_name' and stratified by cohort"'
		display ""	
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		use "Inf_dem_lag.dta", clear
		drop if ses==.
		
		if "`ieks'"=="alcocl" {
		replace lag_entry_anyinf = alkopvm if cohort==1 & alkopvm>=lag_entry_anyinf // 1 = FPS, in FPS, the questionnaire data was sometimes collected only after study entry. To avoid immortal time bias, the follow up is started only after the date of the questionnaire in the analyses using information from that questionnaire.
		tab `ieks'
		drop if `ieks' ==.
		}
		if "`ieks'"=="tupakka" {
		replace lag_entry_anyinf = tupakkapvm if cohort==1 & tupakkapvm>=lag_entry_anyinf  // 1 = FPS, in FPS, the questionnaire data was sometimes collected only after study entry. To avoid immortal time bias, the follow up is started only after the date of the questionnaire in the analyses using information from that questionnaire.
		tab `ieks'
		drop if `ieks' ==.
		}
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm lag_entry_anyinf ensi`idg'pvm  _st _d _origin _t _t0 ses `ieks'_comb
		tab any`idg'
		
		*Those with data on the covariate available but without adjustment
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Data available for `ieks_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		
		*Those with data on the covariate available and adjusted
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		if "`ieks'" == "alcocl" {
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort  1.`ieks'#cohort 2.`ieks'#cohort 3.`ieks'#cohort, strata(cohort)
		}
		else if "`ieks'" == "tupakka" {
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort  2.`ieks'#cohort 3.`ieks'#cohort, strata(cohort)
		}
		else if "`ieks'" == "bmi_who" {
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort  25.`ieks'#cohort 30.`ieks'#cohort, strata(cohort)
		}
		else if "`ieks'" == "bmi" {
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort  c.`ieks'#cohort c.`ieks'2#cohort, strata(cohort)
		}
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Additionally adjusted for `ieks_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		local ++irow 		
}
quietly: local irow = 59
}

*Potential period effect

local dg anyinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease""'
local irow 37
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, potential period effect"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		use "Inf_dem_lag.dta", clear
		drop if ses==.
		
		quietly: stset exitpvm, id(id) failure(failure ==1) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		keep if _st==1
		quietly: keep id any`idg' supu cohort exitpvm ensidementia failure syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses syntvuosi
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.syntvuosi#cohort, strata(cohort)
		
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("Adjusted for period effect")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "Pooled results.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local irow 65
		}

***** TEST FOR DIFFERENCE BETWEEN INFECTIONS *****

*Test for difference.

capture log close
log using "Test for difference v7.log", replace


	use "Inf_dem_lag.dta", replace
	drop if ses==.

	*Intra vs extracellular
	egen ensiintranotextracpvm = rowmin(ensiintracpvm ensiextracpvm)
	replace ensiintranotextracpvm = . if ensiintracpvm==ensiextracpvm
 	gen byte anyintranotextrac = 1 if ensiintracpvm==ensiintranotextracpvm & ensiintranotextracpvm!=. 
	replace anyintranotextrac = 0 if ensiextracpvm==ensiintranotextracpvm & ensiintranotextracpvm!=. & anyintranotextrac!=1
	sum ensiintranotextracpvm
	tab anyintranotextrac
	tab anyintrac
	tab anyextrac
	

	*Gram+ vs gram-
	egen ensigramplusnotminuspvm = rowmin(ensigrampluspvm ensigramminuspvm)
	replace ensigramplusnotminuspvm = . if ensigrampluspvm==ensigramminuspvm
 	gen byte anygramplusnotminus = 1 if ensigrampluspvm==ensigramplusnotminuspvm & ensigramplusnotminuspvm!=. 
	replace anygramplusnotminus = 0 if ensigramminuspvm==ensigramplusnotminuspvm & ensigramplusnotminuspvm!=. & anygramplusnotminus!=1
	sum ensigramplusnotminuspvm
	tab anygramplusnotminus
	tab anygramplus
	tab anygramminus
	

	*Bacterial infection with sepsis vs without sepsis 
	gen byte anybactsepsisnotsepsis = 1 if ensibactsepsispvm==ensibactinfpvm & anybactsepsis==1
	replace anybactsepsisnotsepsis = 0 if ensibact_nonsepsispvm==ensibactinfpvm & anybactsepsisnotsepsis!=1 & anybact_nonsepsis==1
	gen ensibactsepsisnotsepsispvm = ensibactinfpvm
	sum ensibactsepsisnotsepsispvm
	tab anybactsepsisnotsepsis
	tab anybact_nonsepsis
	tab anybactsepsis
	
	*CNS infection vs non cns infection 
	gen byte anycnsnotcns = 1 if ensianycnspvm==ensianyinfpvm & anyanycns==1
	replace anycnsnotcns = 0 if ensinoncnspvm==ensianyinfpvm & anycnsnotcns!=1 & anynoncns==1
	gen ensicnsnotcnspvm = ensianyinfpvm
	sum ensicnsnotcnspvm
	tab anycnsnotcns
	tab anyanycns
	tab anynoncns
	
	*CNS predisposed vs non cns predisposed 
	gen byte anycnspredispnotpredisp = 1 if ensianycnspredisppvm==ensianyinfpvm & anyanycnspredisp==1
	replace anycnspredispnotpredisp = 0 if ensinoncnspredisppvm==ensianyinfpvm & anycnspredispnotpredisp!=1 & anynoncnspredisp==1
	gen ensicnspredispnotpredisppvm = ensianyinfpvm
	sum ensicnspredispnotpredisppvm
	tab anycnspredispnotpredisp
	tab anyanycnspredisp
	tab anynoncnspredisp
	
	*Acute vs chronic
	gen byte anyacutenotchronic = 0 if ensianyperspvm==ensianyinfpvm & anyanypers==1
	replace anyacutenotchronic = 1 if ensinonperspvm==ensianyinfpvm & anyacutenotchronic!=1 & anynonpers==1
	gen ensiacutenotchronicpvm = ensianyinfpvm
	sum ensiacutenotchronicpvm
	tab anyacutenotchronic
	tab anyanypers
	tab anynonpers
	
	*Herpesvirus vs persistent viral infection
	egen ensiherpesnotpersviralpvm = rowmin(ensiherpespvm ensipersviralpvm)
	gen byte anyherpesnotpersviral = 1 if ensiherpespvm==ensiherpesnotpersviralpvm & anyherpes==1
	replace anyherpesnotpersviral = 0 if ensipersviralpvm==ensiherpesnotpersviralpvm & anyherpesnotpersviral!=1 & anypersviral==1
	sum ensiherpesnotpersviralpvm
	tab anyherpesnotpersviral
	tab anyherpes
	tab anypersviral
	
	*Herpesvirus vs acute viral infection
	egen ensiherpesnotacuteviralpvm = rowmin(ensiherpespvm ensiacuteviralpvm)
	gen byte anyherpesnotacuteviral = 1 if ensiherpespvm==ensiherpesnotacuteviralpvm & anyherpes==1
	replace anyherpesnotacuteviral = 0 if ensiacuteviralpvm==ensiherpesnotacuteviralpvm & anyherpesnotacuteviral!=1 & anyacuteviral==1
	sum ensiherpesnotacuteviralpvm
	tab anyherpesnotacuteviral
	tab anyherpes
	tab anyacuteviral
	
	*Non herpes persistent vs acute viral infection
	egen ensipersnotacuteviralpvm = rowmin(ensipersviralpvm ensiacuteviralpvm)
	gen byte anypersnotacuteviral = 1 if ensipersviralpvm==ensipersnotacuteviralpvm & anypersviral==1
	replace anypersnotacuteviral = 0 if ensiacuteviralpvm==ensipersnotacuteviralpvm & anypersnotacuteviral!=1 & anyacuteviral==1
	sum ensipersnotacuteviralpvm
	tab anypersnotacuteviral
	tab anypersviral
	tab anyacuteviral
	
	*Bacteria vs viruses
	egen ensibactnotviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
	gen byte anybactnotviral = 1 if ensibactinfpvm==ensibactnotviralpvm & anybactinf==1
	replace anybactnotviral = 0 if ensiviralinfpvm==ensibactnotviralpvm & anybactnotviral!=1 & anyviralinf==1
	sum ensibactnotviralpvm
	tab anybactnotviral
	tab anybactinf
	tab anyviralinf

	*Type of viral infection
	egen ensivirustypepvm = rowmin(ensiherpespvm ensipersviralpvm ensiacuteviralpvm)
	gen byte anyvirustype = 1 if ensiherpespvm==ensivirustypepvm & anyviralinf==1 & ensivirustypepvm!=.
	replace anyvirustype = 2 if ensipersviralpvm==ensivirustypepvm & anyviralinf==1 & anyvirustype==. & ensivirustypepvm!=.
	replace anyvirustype = 3 if ensiacuteviralpvm==ensivirustypepvm & anyviralinf==1 & anyvirustype==. & ensivirustypepvm!=.
	tab anyvirustype
	tab anyviralinf
	sum ensivirustypepvm ensiviralinfpvm
	*viisi tapausta, joilla on virusinfektio, mutta jonka tyyppi ei ole tiedossa.

	*Any persistent viral infection vs acute viral infection
	gen ensiacutenotanypersviralpvm = ensivirustypepvm
	recode anyvirustype 3=1 1=0 2=0, gen(anyacutenotanypersviral)
	
	
	*Invasive vs localised
	egen ensisysnotlocalpvm = rowmin(ensisysbactpvm ensilocalbactpvm)
	gen byte anysysnotlocal = 1 if ensisysbactpvm==ensisysnotlocalpvm & anysysbact==1
	replace anysysnotlocal = 0 if ensilocalbactpvm==ensisysnotlocalpvm & anysysnotlocal!=1 & anylocalbact==1
	sum ensisysnotlocalpvm
	tab anysysnotlocal
	tab anysysbact
	tab anylocalbact
	
	
	
	
	
local dg sysnotlocal intranotextrac gramplusnotminus bactnotviral acutenotanypersviral cnsnotcns acutenotchronic bactsepsisnotsepsis cnspredispnotpredisp
local n_dg: word count `dg'
local dg_name `""Invasive vs localised bacterial infection" "Intra vs extracellular bacterial infection" "Gram-positive vs gram-negative bacterial infection" "Bacterial vs viral infection" "Acute vs any persisten viral infection" "CNS vs extra-CNS infection" "Acute vs chronic infection" "Septic vs not septic bacterial infection" "Infection predisposed towards entering the CNS vs not predisposed""'
quietly: local irow = 7
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		egen le_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		gen ensi`idg'pvm_10y = ensi`idg'pvm + round(10*365.25)
		egen le_10y_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm_10y)
		
		gen byte anyeither = ensi`idg'pvm !=.
		replace anyeither = . if anyeither==0 & anyanyinf==1
		tab anyeither
				
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_`idg') scale(365.25) 
		quietly: keep id any`idg' ensi`idg'pvm supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm_10y anyeither ///
		le_`idg' le_10y_`idg' _st _d _origin _t _t0 ses 

		keep if _st==1
		gen `idg'_compare =0
		replace `idg'_compare = 1 if anyeither==1
		tab `idg'_compare
		recode `idg'_compare 1=2 if any`idg'==0
		tab `idg'_compare
		tab any`idg'
		tab anyeither
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 

		*Dementia occurring from year 10 onwards
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_10y_`idg') scale(365.25) 
		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
}

	log close

*Comparing viral infections
capture log close
log using "herpes vs other persistent viral v7.log", replace	
local dg herpesnotpersviral
local n_dg: word count `dg'
local dg_name `""Herpesvirus infection vs other potentially persistent viral infection""'
quietly: local irow = 40
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		egen le_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		gen ensi`idg'pvm_10y = ensi`idg'pvm + round(10*365.25)
		egen le_10y_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm_10y)
		
		gen byte anyeither = ensi`idg'pvm !=.
		replace anyeither = . if anyeither==0 & anyanyinf==1
		tab anyeither
				
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_`idg') scale(365.25) 
		quietly: keep id any`idg' ensi`idg'pvm supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm_10y anyeither ///
		le_`idg' le_10y_`idg' _st _d _origin _t _t0 ses 

		keep if _st==1
		gen `idg'_compare =0
		replace `idg'_compare = 1 if anyeither==1
		tab `idg'_compare
		recode `idg'_compare 1=2 if any`idg'==0
		tab `idg'_compare
		tab any`idg'
		tab anyeither

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 

		*Dementia occurring from year 10 onwards
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_10y_`idg') scale(365.25) 
		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
		}
	log close
	

capture log close
log using "Types of viral infection v7.log", replace	
	
local dg virustype 
local n_dg: word count `dg'
local dg_name `""Types of viral infection""'
quietly: local irow = 34
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		egen le_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm)
		gen ensi`idg'pvm_10y = ensi`idg'pvm + round(10*365.25)
		egen le_10y_`idg' = rowmax(lag_entry_anyinf ensi`idg'pvm_10y)
		
		gen byte anyeither = ensi`idg'pvm !=.
		replace anyeither = . if anyeither==0 & anyanyinf==1
		tab anyeither
				
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_`idg') scale(365.25) 
		quietly: keep id any`idg' ensi`idg'pvm supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm_10y anyeither ///
		le_`idg' le_10y_`idg' _st _d _origin _t _t0 ses anyviralinf anyanyinf

		keep if _st==1
		gen  `idg'_compare = any`idg' if anyeither==1
		replace `idg'_compare = 0 if anyanyinf==0
		tab `idg'_compare
		tab any`idg'
		tab anyeither
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		
		if "`idg'"=="virustype" {
		local comp "viralinf"
		replace anyviralinf = . if anyvirustype==.
		}
		tab any`comp'
		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare = 3.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`comp'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time le_10y_`idg') scale(365.25) 

		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_compare = 2.`idg'_compare = 3.`idg'_compare
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
		}

	log close	


	
***** NUMBER OF INFECTIONS *******		

capture log close
log using "P for number of simultaneous infections v7.log", replace

local dg infcount
local n_dg: word count `dg'
local dg_name `""No. of infections""'

		use "Inf_dem_lag.dta", replace
		merge 1:1 id using "Infection count.dta"
		tab anyinfcount
		recode anyinfcount 3=2 4=2
		drop if ses==.
		gen ensiinfcountpvm = ensianyinfpvm 
		
quietly: local irow = 37
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
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses ///
		lag_entry_anyinf lag_entry_10y_anyinf
		keep if _st==1
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.any`idg' = 2.any`idg'
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow

		*Dementia occurring from year 10 onwards
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_10y_anyinf) scale(365.25) 

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.any`idg' = 2.any`idg'
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow
		local ++irow
		restore
		}
log close		
	

capture log close	
log using "HSV vs other herpes v7.log", replace
		use "Inf_dem_lag.dta", replace
		merge 1:1 id using "herpeskoodit.dta"
		drop _merge
		drop if ses==.
		tab anyhsv12
		gen anynonhsv = 0 if anyanyinf==0
		replace anynonhsv = 1 if anyherpes==1 & anyhsv12!=1
		
		gen herpes_compare = 1 if anyhsv12==1
		replace herpes_compare = 2 if anynonhsv==1
		replace herpes_compare = 0 if anyanyinf==0
		tab herpes_compare
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_herpes) scale(365.25) 
		quietly: keep id anyherpes anyhsv12 anynonhsv supu cohort exitpvm ensidementia syntpvm entrypvm ensiherpespvm ///
		_st _d _origin _t _t0 ses herpes_compare
		tab anyherpes
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.herpes_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.herpes_compare=2.herpes_compare
log close		

capture log close	
log using "Type of Gram-positive and Gram-negative infections v7.log", replace

		quietly: use "Inf_dem_lag.dta", clear
		merge 1:1 id using "grampluskoodit.dta"
		drop _merge
		merge 1:1 id using "gramminuskoodit.dta"
		drop _merge
		drop if ses==.
		
		anystreptoc anycdiffic anyothergramplus anysalmonella anylyme anyothergramminussepsis
		
		gen gramplus_compare = 0 if anyanyinf==0
		replace gramplus_compare = 1 if anystreptoc==1
		replace gramplus_compare = 2 if anycdiffic==1
		replace gramplus_compare = 3 if anyothergramplus==1
		tab anygramplus
		tab gramplus_compare
		
		gen gramminus_compare = 0 if anyanyinf==0
		replace gramminus_compare = 1 if anysalmonella==1
		replace gramminus_compare = 2 if anylyme==1
		replace gramminus_compare = 3 if anyothergramminussepsis==1 | anyothergramminus==1
		tab anygramminus
		tab gramminus_compare
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_gramplus) scale(365.25) 
		quietly: keep id *gram* supu cohort exitpvm ensidementia syntpvm entrypvm  _st _d _origin _t _t0 ses
		tab anygramplus if _st==1
		tab gramplus_compare if _st==1
		
		stcox i.gramplus_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.gramplus_compare = 2.gramplus_compare
		test 1.gramplus_compare = 3.gramplus_compare
		test 2.gramplus_compare = 3.gramplus_compare
		test 1.gramplus_compare = 2.gramplus_compare = 3.gramplus_compare

		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_gramminus) scale(365.25) 
		quietly: keep id *gram* supu cohort exitpvm ensidementia syntpvm entrypvm  _st _d _origin _t _t0 ses
		tab anygramminus if _st==1
		tab gramminus_compare if _st==1
		
		stcox i.gramminus_compare 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.gramminus_compare = 2.gramminus_compare
		test 1.gramminus_compare = 3.gramminus_compare
		test 2.gramminus_compare = 3.gramminus_compare
		test 1.gramminus_compare = 2.gramminus_compare = 3.gramminus_compare

	log close


********* PROPORTIONAL HAZARDS ASSUMPTIONS **********

capture log close
log using "Proportional hazards assumptions v3.log", replace

local dg anyinf bactinf sysbact localbact bactsepsis bact_nonsepsis extrac intrac gramplus gramminus viralinf herpes persviral acuteviral anycns noncns
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Any bacterial infection" "Invasive bacterial infection" "Localised bacterial infection" "Bacterial infection with sepsis" "Bacterial infection without sepsis" "Extracellular bacterial infection" "Intracellular bacterial infection" "Gram-positive bacterial infection" "Gram-negative bacterial infection" "Any viral infection" "Herpesvirus infection" "Other potentially persistent viral infection" "Acute viral infection" "CNS infection" "Extra-CNS infection""'

		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.
		egen lag_entry_anycns = rowmax(lag_entry_anyinf ensianycnspvm)
		egen lag_entry_noncns = rowmax(lag_entry_anyinf ensinoncnspvm)
		replace anyanycns = . if anyanycns==0 & anyanyinf==1
		replace anynoncns = . if anynoncns==0 & anyanyinf==1


quietly: local irow = 7
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""	
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		preserve
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		keep if _st==1
		quietly: keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses lag_entry_`idg'
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		matrix r = r(table)
		estat phtest, detail
		matrix m = r(phtest)
		quietly: putexcel set "Pooled results.xlsm", sheet("IPD, PH-assumptions") modify
        quietly: putexcel A`irow'="`idg'", nformat(number)
        quietly: putexcel B`irow'="`idg_name'", nformat(number)
		quietly: putexcel C`irow'=`r(p)', nformat(number_d2)
		quietly: putexcel D`irow'=matrix(m[1,4]), nformat(number_d2)
		quietly: putexcel E`irow'=matrix(m[2,4]), nformat(number_d2)
		quietly: putexcel F`irow'=matrix(m[3,4]), nformat(number_d2)
		quietly: putexcel G`irow'=matrix(m[4,4]), nformat(number_d2)
		quietly: putexcel H`irow'=matrix(m[5,4]), nformat(number_d2)
		quietly: putexcel I`irow'=matrix(m[6,4]), nformat(number_d2)
		quietly: putexcel J`irow'=matrix(m[7,4]), nformat(number_d2)
		quietly: putexcel K`irow'=matrix(m[8,4]), nformat(number_d2)
		quietly: putexcel L`irow'=matrix(m[9,4]), nformat(number_d2)
		quietly: putexcel M`irow'=matrix(m[10,4]), nformat(number_d2)
		quietly: putexcel N`irow'=matrix(r[1,1]), nformat(number_d2)
	local ++irow
	restore
}
quietly: local irow = 6
		quietly: putexcel set "Pooled results.xlsm", sheet("IPD, PH-assumptions") modify
		quietly: putexcel E`irow'="Sex_FPS", nformat(number_d2)
		quietly: putexcel F`irow'="Sex_HeSSup", nformat(number_d2)
		quietly: putexcel G`irow'="Sex_STW", nformat(number_d2)
		quietly: putexcel E`irow'="2.edu_FPS", nformat(number_d2)
		quietly: putexcel F`irow'="2.edu_HeSSup", nformat(number_d2)
		quietly: putexcel G`irow'="2.edu_STW", nformat(number_d2)
		quietly: putexcel E`irow'="3.edu_FPS", nformat(number_d2)
		quietly: putexcel F`irow'="3.edu_HeSSup", nformat(number_d2)
		quietly: putexcel G`irow'="3.edu_STW", nformat(number_d2)
log close


*any infectious disease, invasive bacterial infections, extracellular bacterial infections, intracellular bacterial infections, Gram-negative bacterial infections, any viral infection, herpesvirus infections, other potentially persistent viral infections and non-cns infections had p<0.05 which means evidence for the violation of the proportional hazards assumption. Let's draw figures from them:


	**********************Scaled Schoenfeld residuals**********************************


	capture log close
	log using "Scaled Schoenfeld residuals.log", replace


	*Scaled schoenfeld residuals -figures:
*The code below that draws the figures has been modified from:
*Royston and Lambert: Flexible Parametric Survival Analysis Using Stata: Beyond the Cox Model. Chapter 7.3 What do we mean by a TD effect?


	local dg anyinf sysbact extrac intrac gramminus viralinf herpes persviral noncns
	local n_dg: word count `dg'
	local dg_name `""Any infectious disease (vs no infection)" "Invasive bacterial infection (vs no infection)" "Extracellular bacterial infection (vs no infection)" "Intracellular bacterial infection (vs no infection)" "Gram-negative bacterial infection (vs no infection)" "Any viral infection (vs no infection)" "Herpesvirus infection (vs no infection)" "Other potentially persistent viral infection (vs no infection)" "Extra-CNS infection (vs no infection)"'
		quietly: use "Inf_dem_lag.dta", clear
		drop if ses==.
		egen lag_entry_anycns = rowmax(lag_entry_anyinf ensianycnspvm)
		egen lag_entry_noncns = rowmax(lag_entry_anyinf ensinoncnspvm)
		replace anyanycns = . if anyanycns==0 & anyanyinf==1
		replace anynoncns = . if anynoncns==0 & anyanyinf==1
		drop if ses==.
	
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		
		preserve

		local j = 1
		foreach cohort in FPS HeSSup STW {
			gen `cohort'supu=0
			replace `cohort'supu=1 if supu==2 & cohort==`j'
			gen `cohort'ses2=0
			replace `cohort'ses2=1 if ses==2 & cohort==`j'
			gen `cohort'ses3=0
			replace `cohort'ses3=1 if ses==3 & cohort==`j'
			local ++j
			}
		tabulate cohort, gen(dummycohort)
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
			quietly: keep id any`idg' *supu* cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 *ses*
			keep if _st==1
			tab any`idg'
			
		stcox any`idg' FPSsupu HeSSupsupu STWsupu FPSses2 HeSSupses2 STWses2 FPSses3 HeSSupses3 STWses3, strata(cohort)
		capture drop *sca*
		predict sca*, scaledsch
		running sca1 _t if _d==1, gen(smooth_sca) gense(smooth_sca_se) nodraw
		gen smooth_esca = exp(smooth_sca)
		gen smooth_esca_lci = exp(smooth_sca - 1.96*smooth_sca_se)
		gen smooth_esca_uci = exp(smooth_sca + 1.96*smooth_sca_se)
		local betaround: di %7.2f exp(_b[any`idg'])
		tokenize "`betaround'"
		di "`1'"
		twoway (rarea smooth_esca_lci smooth_esca_uci _t, pstyle(ci) sort yaxis(1 2)) ///
			(line smooth_esca _t, sort clpattern(solid)) ///
			(function y = 1, lpattern(shortdash) range(_t)) ///
			(function y = `1', lpattern(longdash) range(_t)), ///
			legend(order(2 "Exponentiated scaled Schoenfeld residuals" 1 "95% confidence interval" 4 "Overall hazard ratio   ") holes(2) size(small)) ///
			title("`idg_name'" "and dementia", size(medsmall)) ///
			ytitle("Exponentiated scaled Schoenfeld residuals", size(small)) ///
			xtitle("Age (years)", size(small)) ///
			yscale(log range(0.1 100)) ///
			scheme(sj) ///
			ylabel(0.1 1 10 100, labsize(medsmall) angle(0)) ///
			ylabel("`1'", angle(0) axis(2))
		graph save sca_`idg'.gph, replace
		graph export sca_`idg'.pdf, as(pdf) replace
		restore
		}


	*Combine the figures:


	*An empty figure that is needed if the number of figures is not divisible by three:
	scatter cohort supu if _n==1, ysize(3.9) xsize(4.135) scheme(s1color) plotregion(style(none)) yscale(off noline) xscale(off noline) mcolor(%0) graphregion(color(white))
	graph save "sca_.gph", replace


	local dg anyinf sysbact extrac intrac gramminus viralinf herpes persviral noncns
	local n_dg: word count `dg'
	local i=1
	while `i' <= 9 {
		forvalues j=1/3 {
		local k = `i'+`j'-1
		quietly: local dg`j' : word `k' of `dg'
		}
	graph combine sca_`dg1'.gph sca_`dg2'.gph sca_`dg3'.gph, col(1) altshrink  scheme(s1color) graphregion(color(white))  plotregion(color(gs15)) ysize(8.27) xsize(5.83) // ysize ja xsize A5:n mukaan (tuumina)
	graph export "IPD_sca_`i'.emf", as(emf) replace
	local i = `i'+3
		}

	log close

		

***** TIME-DEPENDENT ANALYSES ********

****eFIGURE 11 *******

*(time-dependent version of the analyses in FIGURE 2)
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
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 revised") modify
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
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 revised") modify
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
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 revised") modify
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
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE 2 revised") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		}
log close


		
***** TIME-DEPENDENT ANALYSES: SINGLE VS MULTIPLE INFECTIONS *******		

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

***** TIME-DEPENDENT ANALYSES: One or several different pathogens *******		

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



***** TIME-DEPENDENT ANALYSES: NUMBER OF SIMULTANEOUS INFECTIONS *******		

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


***** TIME-DEPENDENT ANALYSES: FIGURE sensitivity, CNS vs other infections
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

version 16

***** TIME-DEPENDENT ANALYSES: Test for difference between infections

capture log close
log using "Test for difference.log", replace


	quietly: use "Pooled infections_trimmattu", clear
	drop if ses==.

	*Intra vs extracellular
	egen ensiintranotextracpvm = rowmin(ensiintracpvm ensiextracpvm)
	replace ensiintranotextracpvm = . if ensiintracpvm==ensiextracpvm
 	gen byte anyintranotextrac = 1 if ensiintracpvm==ensiintranotextracpvm & ensiintranotextracpvm!=. 
	replace anyintranotextrac = 0 if ensiextracpvm==ensiintranotextracpvm & ensiintranotextracpvm!=. & anyintranotextrac!=1
	sum ensiintranotextracpvm
	tab anyintranotextrac
	tab anyintrac
	tab anyextrac
	

	*Gram+ vs gram-
	egen ensigramplusnotminuspvm = rowmin(ensigrampluspvm ensigramminuspvm)
	replace ensigramplusnotminuspvm = . if ensigrampluspvm==ensigramminuspvm
 	gen byte anygramplusnotminus = 1 if ensigrampluspvm==ensigramplusnotminuspvm & ensigramplusnotminuspvm!=. 
	replace anygramplusnotminus = 0 if ensigramminuspvm==ensigramplusnotminuspvm & ensigramplusnotminuspvm!=. & anygramplusnotminus!=1
	sum ensigramplusnotminuspvm
	tab anygramplusnotminus
	tab anygramplus
	tab anygramminus
	

	*Bacterial infection with sepsis vs without sepsis 
	gen byte anybactsepsisnotsepsis = 1 if ensibactsepsispvm==ensibactinfpvm & anybactsepsis==1
	replace anybactsepsisnotsepsis = 0 if ensibact_nonsepsispvm==ensibactinfpvm & anybactsepsisnotsepsis!=1 & anybact_nonsepsis==1
	gen ensibactsepsisnotsepsispvm = ensibactinfpvm
	sum ensibactsepsisnotsepsispvm
	tab anybactsepsisnotsepsis
	tab anybact_nonsepsis
	tab anybactsepsis
	
	*CNS infection vs non cns infection 
	gen byte anycnsnotcns = 1 if ensianycnspvm==ensianyinfpvm & anyanycns==1
	replace anycnsnotcns = 0 if ensinoncnspvm==ensianyinfpvm & anycnsnotcns!=1 & anynoncns==1
	gen ensicnsnotcnspvm = ensianyinfpvm
	sum ensicnsnotcnspvm
	tab anycnsnotcns
	tab anyanycns
	tab anynoncns
	
	*Herpesvirus vs persistent viral infection
	egen ensiherpesnotpersviralpvm = rowmin(ensiherpespvm ensipersviralpvm)
	gen byte anyherpesnotpersviral = 1 if ensiherpespvm==ensiherpesnotpersviralpvm & anyherpes==1
	replace anyherpesnotpersviral = 0 if ensipersviralpvm==ensiherpesnotpersviralpvm & anyherpesnotpersviral!=1 & anypersviral==1
	sum ensiherpesnotpersviralpvm
	tab anyherpesnotpersviral
	tab anyherpes
	tab anypersviral
	
	*Herpesvirus vs acute viral infection
	egen ensiherpesnotacuteviralpvm = rowmin(ensiherpespvm ensiacuteviralpvm)
	gen byte anyherpesnotacuteviral = 1 if ensiherpespvm==ensiherpesnotacuteviralpvm & anyherpes==1
	replace anyherpesnotacuteviral = 0 if ensiacuteviralpvm==ensiherpesnotacuteviralpvm & anyherpesnotacuteviral!=1 & anyacuteviral==1
	sum ensiherpesnotacuteviralpvm
	tab anyherpesnotacuteviral
	tab anyherpes
	tab anyacuteviral
	
	*Non herpes persistent vs acute viral infection
	egen ensipersnotacuteviralpvm = rowmin(ensipersviralpvm ensiacuteviralpvm)
	gen byte anypersnotacuteviral = 1 if ensipersviralpvm==ensipersnotacuteviralpvm & anypersviral==1
	replace anypersnotacuteviral = 0 if ensiacuteviralpvm==ensipersnotacuteviralpvm & anypersnotacuteviral!=1 & anyacuteviral==1
	sum ensipersnotacuteviralpvm
	tab anypersnotacuteviral
	tab anypersviral
	tab anyacuteviral
	
	*Bacteria vs viruses
	egen ensibactnotviralpvm = rowmin(ensibactinfpvm ensiviralinfpvm)
	gen byte anybactnotviral = 1 if ensibactinfpvm==ensibactnotviralpvm & anybactinf==1
	replace anybactnotviral = 0 if ensiviralinfpvm==ensibactnotviralpvm & anybactnotviral!=1 & anyviralinf==1
	sum ensibactnotviralpvm
	tab anybactnotviral
	tab anybactinf
	tab anyviralinf

	*Type of viral infection
	egen ensivirustypepvm = rowmin(ensiherpespvm ensipersviralpvm ensiacuteviralpvm)
	gen byte anyvirustype = 1 if ensiherpespvm==ensivirustypepvm & anyviralinf==1 & ensivirustypepvm!=.
	replace anyvirustype = 2 if ensipersviralpvm==ensivirustypepvm & anyviralinf==1 & anyvirustype==. & ensivirustypepvm!=.
	replace anyvirustype = 3 if ensiacuteviralpvm==ensivirustypepvm & anyviralinf==1 & anyvirustype==. & ensivirustypepvm!=.
	tab anyvirustype
	tab anyviralinf
	sum ensivirustypepvm ensiviralinfpvm
	*viisi tapausta, joilla on virusinfektio, mutta jonka tyyppi ei ole tiedossa.

	*Any persistent viral infection vs acute viral infection
	gen ensiacutenotanypersviralpvm = ensivirustypepvm
	recode anyvirustype 3=1 1=0 2=0, gen(anyacutenotanypersviral)
	
	
	*Invasive vs localised
	egen ensisysnotlocalpvm = rowmin(ensisysbactpvm ensilocalbactpvm)
	gen byte anysysnotlocal = 1 if ensisysbactpvm==ensisysnotlocalpvm & anysysbact==1
	replace anysysnotlocal = 0 if ensilocalbactpvm==ensisysnotlocalpvm & anysysnotlocal!=1 & anylocalbact==1
	sum ensisysnotlocalpvm
	tab anysysnotlocal
	tab anysysbact
	tab anylocalbact
	
	
	
	
	
local dg sysnotlocal intranotextrac gramplusnotminus bactnotviral acutenotanypersviral cnsnotcns bactsepsisnotsepsis
local n_dg: word count `dg'
local dg_name `""Invasive vs localised bacterial infection" "Intra vs extracellular bacterial infection" "Gram-positive vs gram-negative bacterial infection" "Bacterial vs viral infection" "Acute vs any persisten viral infection" "CNS vs extra-CNS infection" "Septic vs not septic bacterial infection""'
quietly: local irow = 7
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' ensi`idg'pvm supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		
		gen byte anyeither = ensi`idg'pvm !=.
		tab any`idg'
		tab anyeither
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if anyeither==1 & `idg'_split==0
		tab `idg'_timedep
		tab anyeither if `idg'_split!=-1

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		tab anyeither if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`idg'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
		}

	log close

*Comparing viral infections
capture log close
log using "herpes vs other persistent viral.log", replace	
local dg herpesnotpersviral
local n_dg: word count `dg'
local dg_name `""Herpesvirus infection vs other potentially persistent viral infection""'
quietly: local irow = 40
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' ensi`idg'pvm supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		
		gen byte anyeither = ensi`idg'pvm !=.
		tab any`idg'
		tab anyeither
		
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = 1 if anyeither==1 & `idg'_split==0
		tab `idg'_timedep
		tab anyeither if `idg'_split!=-1

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		tab anyeither if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`idg'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
		}

	log close
	

capture log close
log using "Types of viral infection.log", replace	
	
local dg virustype 
local n_dg: word count `dg'
local dg_name `""Types of viral infection""'
quietly: local irow = 34
	forvalues i=1/`n_dg' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time entrypvm) scale(365.25) 
		quietly: keep id any`idg' *viralinf* supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		stsplit `idg'_split, at(0) after(time=ensi`idg'pvm)
		tab `idg'_split
		gen `idg'_timedep =0
		replace `idg'_timedep = any`idg' if inlist(any`idg',1,2,3) & `idg'_split==0
		tab `idg'_timedep
		tab any`idg' if `idg'_split!=-1

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'

		
		if "`idg'"=="virustype" {
		local comp "viralinf"
		replace anyviralinf = . if anyvirustype==.
		}
		tab any`comp' if `idg'_split!=-1
		tab `idg'_timedep
		recode `idg'_timedep 1=2 if any`comp'==0
		tab `idg'_timedep
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep = 3.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`comp'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow 
		
		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		stdescribe
		quietly: local failures = `r(N_fail)'
		
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep = 3.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow 
		local ++irow		
		restore
		}

	log close	


	
***** NUMBER OF INFECTIONS *******		

capture log close
log using "P for number of simultaneous infections.log", replace

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
		
quietly: local irow = 37
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
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name'")
		local ++irow

		*Dementia occurring from year 10 onwards
		gen ensi`idg'pvm10v = ensi`idg'pvm + round(10*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm10v = entrypvm if ensi`idg'pvm10v<entrypvm
		replace _t0 = (ensi`idg'pvm10v-syntpvm)/365.25 if ensi`idg'pvm10v!=.
		replace _t0 = _t0 + 10 if ensi`idg'pvm10v==.
		drop if _t0>=_t
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.`idg'_timedep 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		test 1.`idg'_timedep = 2.`idg'_timedep
		
		quietly: putexcel set "Pooled results.xlsm", sheet("FIGURE compare") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`r(p)'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("pooled")
			quietly: putexcel D`irow'=("`idg_name' from year 10 onwards")
		local ++irow
		local ++irow
		restore
		}
log close		












********** REPLICATION ANALYSIS (UK BIOBANK) ***************

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

foreach idg in anyinf bactinf viralinf hiv anycns anynoncns {
replace ensi`idg'pvm=. if ensi`idg'pvm>=exitdate
replace any`idg' = 0 if ensi`idg'pvm==.
sum ensi`idg'pvm
tab any`idg'
}

*Create variables that have indentical names with those used in the primary analysis

capture gen diabetes_combpvm = diabetesdate
capture gen hypert_combpvm = hypertensiondate
capture gen ihd_combpvm = ihddate
capture gen cereb_combpvm = cerebdate
capture gen parkinson_combpvm = parkinsondate
capture gen tupakka = smoking
*alcocl does not change
*ses does not change
capture gen entrypvm = entrydate
capture gen exitpvm = exitdate
capture gen ensidementia = dementia
recode ensidementia 0=.
capture gen syntpvm = birthday
capture gen supu = sex
capture gen dementiatyyppi_laaja = dementiatype
capture gen byte cohort==1

save "SipilaELSA.dta", replace

*** Lag periods ***

* Random lag periods for those without an infection drawn from the distribution of the lag periods of those with an infection.

use "SipilaELSA.dta", clear

version 16

*Random order
*seed from random.org from range 1 - 1 000 000 000
set rng mt64s
set rngstream 1
set seed 421100945

gen double shuffle1 = runiform()
gen double shuffle2 = runiform()

sort shuffle1 shuffle2

*New pseudo-id-variable that preserves the random order
gen shuffle_id = _n

*Delay between study entry and time of infection
gen entry_to_inf = (ensianyinfpvm-entrypvm)/365.25
replace entry_to_inf = 0 if entry_to_inf < 0 
bysort supu: sum entry_to_inf if entry_to_inf>0
sum entry_to_inf
capture drop entry_to_inf_cut
egen entry_to_inf_cut = cut(entry_to_inf), at(0 0.0009765625 2.5 5 7.5 100)
tab entry_to_inf_cut
by supu: tab entry_to_inf_cut
capture drop age
gen age = (entrypvm-syntpvm)/365.25
capture drop age_cut
egen age_cut = cut(age), at(18 50 60 100)
tab age_cut
bysort supu age_cut: tab entry_to_inf_cut

levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_cut, local(entry_lags)
local i 0
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0	
			foreach lag of local entry_lags {
				local ++i
			
				quietly: sum entry_to_inf if age_cut==`age' & supu==`sex' & entry_to_inf_cut==`lag'				
				
				if `i'==1 & `r(N)' == 0 {
				matrix mean_lag_`age'_`sex' = 0
				}
				else if `i'==1 {
				matrix mean_lag_`age'_`sex' = `r(mean)'
				}
				else if `r(N)' == 0 {
				matrix mean_lag_`age'_`sex' = mean_lag_`age'_`sex' \ 0
				}
				else {
				matrix mean_lag_`age'_`sex' = mean_lag_`age'_`sex' \ `r(mean)'
				}
				}
				matrix list mean_lag_`age'_`sex'				
				}
				}
				
*Return random order
capture drop sub_order
bysort age_cut supu anyanyinf (shuffle1 shuffle2): gen sub_order = _n

*Create a similar distribution of delays between study entry and start of follow up for those without infections as with those with infection
capture drop entry_lag
gen entry_lag = .
levelsof age_cut, local(age_groups)
levelsof supu, local(sexes)
levelsof entry_to_inf_cut, local(entry_lags)
	foreach age of local age_groups {
		foreach sex of local sexes {
			local i 0
			local cum = 0
			sum id if age_cut==`age' & supu==`sex' & anyanyinf==0				
			local N_control = `r(N)'
			sum entry_to_inf_cut if age_cut==`age' & supu==`sex' & anyanyinf==1				
			local N_`age'_`sex' = `r(N)'
			foreach lag of local entry_lags {
				di "age_group `age', sex `sex'"
				local ++i
				di `i'
				sum entry_to_inf_cut if age_cut==`age' & supu==`sex' & anyanyinf==1 & entry_to_inf_cut==`lag'			
				local cum = `cum' + `r(N)'
				di `cum'
				di `N_`age'_`sex''
				di `N_control'
				di mean_lag_`age'_`sex'[`i',1]
				replace entry_lag = mean_lag_`age'_`sex'[`i',1] if ///
				age_cut==`age' & supu==`sex' & anyanyinf==0 & entry_lag==. & sub_order <= round((`cum' / `N_`age'_`sex'')*`N_control')
				}
				replace entry_lag = 0 if `N_`age'_`sex''==0 & entry_lag==. & anyanyinf==0
				}
				}
		
*For each infection variable, code it to be missing if the person has some other infection
local dg anyinf
foreach idg of local dg {
		tab any`idg'
		egen lag_entry_`idg' = rowmax(entrypvm ensi`idg'pvm) if any`idg'==1
		replace lag_entry_`idg' = entrypvm + round(entry_lag*365.25) if any`idg'==0
		sum lag_entry_`idg'
		}
		
local dg bactinf viralinf anycns noncns
foreach idg of local dg {
		tab any`idg'
		replace any`idg' = . if any`idg'!=1 
		replace any`idg' = 0 if anyanyinf==0 & any`idg'!=1 
		tab any`idg'
		egen lag_entry_`idg' = rowmax(entrypvm ensi`idg'pvm) if any`idg'==1
		replace lag_entry_`idg' = entrypvm + round(entry_lag*365.25) if any`idg'==0
		sum lag_entry_`idg'
		}		
		

sum entry_lag
sum entry_to_inf
bysort supu age_cut: sum entry_lag entry_to_inf

save "UKB_lag.dta", replace




**** Analysis *****

use "UKB_lag.dta", clear


*Table 1


*Instal baselinetable package
capture net install st0524_1

preserve
stset exitdate, id(id) failure(dementia) origin(time birthday) enter(time lag_entry_anyinf) scale(365.25) 
	keep if _st==1 & ses!=.
	keep if anyanyinf!=.
	count
	count if ensidementia==1
	
		capture label define sex 1"Men" 2"Women"
		label values sex sex

		capture gen followup = _t-_t0
		
		capture gen age_entry = (lag_entry_anyinf - syntpvm)/365.25
		capture egen age_entry_cat = cut(age_entry), at(18,40,50,60,100)
		tab age_entry_cat
		sum age_entry, de
		sum age_entry, de, if supu==1
		sum age_entry, de, if supu==2
				
		label define alcocl 0 "Non-drinker" 1 "Moderate" 2 "Intermediate" 3 "Heavy"
		label values alcocl alcocl
		
		foreach cm in hypert diabetes ihd cereb parkinson {
			capture gen `cm'_entry= `cm'_combpvm<=lag_entry_anyinf
			capture tab `cm'_entry
			}
		
		capture gen age_at_dementia = _t if _d==1
		
		stdescribe
		sum _t, de
		sum _t if _d==1, de
		

		baselinetable ///
		age_entry_cat(novarlabel afterhead("Age at entry, years")) ///
		age_entry(cts novarlabel afterhead("Age at entry (IQR), years")) ///
		sex(novarlabel afterhead("Sex")) ///
		ses(novarlabel afterhead("Education/Socioeconomic status")) ///	
		hypert_entry(novarlabel afterhead("hypertensionn")) ///
		diabetes_entry(novarlabel afterhead("Diabetes mellitus")) ///
		ihd_entry(novarlabel afterhead("Ischaemic heart disease")) ///
		cereb_entry(novarlabel afterhead("Cerebrovascular disease")) ///
		parkinson_entry(novarlabel afterhead("Parkinson disease")) ///
		smoking(varlabel afterhead("Smoking status")) ///
		alcocl(novarlabel afterhead("Alcohol drinking")) ///
		bmi_who(novarlabel afterhead("Body mass index")) ///
		apoe(novarlabel afterhead("apolipoprotein E genotype")) ///
		followup(cts novarlabel afterhead("Follow-up, median (IQR), years")) ///
		dementia(novarlabel afterhead("Dementia")) ///
		age_at_dementia(cts novarlabel afterhead("Age at dementia, median (IQR), years")) ///
		, ctsvartab(p50 (p25-p75)) missing ///
		exportexcel(SipilaUKBtable1, cell(A6) replace)
restore

*FIGURE 3
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
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25)
		tab any`idg'
		tab any`idg' _d
		sum diabetes_combpvm
		replace diabetes_combpvm = . if diabetes_combpvm>=exitpvm
		sum diabetes_combpvm
		gen byte dm=0
		replace dm = 1 if diabetes_combpvm <= lag_entry_anyinf
		tab dm
		sum hypert_combpvm
		replace hypert_combpvm = . if hypert_combpvm>=exitpvm
		sum hypert_combpvm
		gen byte htn = 0
		replace htn = 1 if hypert_combpvm <= lag_entry_anyinf
		tab htn
		
		
		if `j'==1 {
		di ""
		di "MODEL 1: adjusted for age (as the time scale), sex and ses"
		di ""
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
		}
		else if `j'==2 {
		di ""
		di "MODEL 2 = model 1 + those with HIV excluded"
		di ""
		drop if anyhiv==1
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
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
		stcox any`idg' i.supu i.ses
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm i.apoe
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
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name' model `j'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
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
		
		stset exitpvm, id(id) failure(dementiatype==`type') origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25)
		tab any`idg'
		tab any`idg' _d
		sum diabetes_combpvm
		replace diabetes_combpvm = . if diabetes_combpvm>=exitpvm
		sum diabetes_combpvm
		gen byte dm = 0
		replace dm = 1 if diabetes_combpvm <= lag_entry_anyinf
		tab dm
		sum hypert_combpvm
		replace hypert_combpvm = . if hypert_combpvm>=exitpvm
		sum hypert_combpvm
		gen byte htn = 0
		replace htn = 1 if hypert_combpvm <= lag_entry_anyinf
		tab htn
			
		if `j'==1 {
		di ""
		di "MODEL 1: adjusted for age (as the time scale), sex and ses"
		di ""
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
		}
		else if `j'==2 {
		di ""
		di "MODEL 2 = model 1 + those with HIV excluded"
		di ""
		drop if anyhiv==1
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
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
		stcox any`idg' i.supu i.ses
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm
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
		stcox any`idg' i.supu i.ses i.tupakka i.alcocl i.bmi_who i.htn i.dm i.apoe
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
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("`idg_name' model `j'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
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
use "UKB_lag.dta", clear
count
drop if exitdate<=entrydate
count
drop if ses==.
count

tab dementia
recode dementia .=0
tab dementia 

foreach idg in anyinf bactinf viralinf hiv anycns noncns {
replace ensi`idg'pvm=. if ensi`idg'pvm>=exitdate
replace any`idg' = 0 if ensi`idg'pvm==.
replace any`idg' = . if ensi`idg'pvm==. & anyanyinf==1
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
		

save "UKB_lag.dta", replace


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
		
		egen lag_entry_multiinf = rowmax(lag_entry_anyinf ensitoinenanyinfpvm ensikolmasanyinfpvm)
		egen multiinfpvm = rowmax(ensianyinfpvm ensitoinenanyinfpvm ensikolmasanyinfpvm)

		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_multiinf) scale(365.25) 
		
		gen byte multiinf = 0
		tab multiinf
		replace multiinf = 1 if anyanyinf==1
		tab multiinf
		replace multiinf = 2 if ensitoinenanyinfpvm!=.
		replace multiinf = 3 if ensikolmasanyinfpvm!=.
		tab multiinf
		
		stdescribe
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3 N1 N2 N3 fail1 fail2 fail3
		forvalues i =1/3 {
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if inlist(multiinf,0,`i')
		quietly: matrix m = r(table)
		matrix list m
		scalar `N`i'' = e(N_sub)
		scalar `fail`i'' = e(N_fail)
		scalar `hr`i'' = m[1,2]
		scalar `ll`i'' = m[5,2]
		scalar `ul`i'' = m[6,2]
		}
		
		forvalues i =1/3 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf, matcell(x), if inlist(multiinf,0,`i') & _st==1
		matrix list x
		tab multiinf ensidementia, matcell(y), if inlist(multiinf,0,`i') & _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		
		restore

		
***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
		egen lag_entry_bactviral = rowmax(lag_entry_anyinf ensibactviralpvm toinenbactviralpvm)
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_bactviral) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses anyanyinf
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral if _st==1
		tab toinenbactviral if _st==1
		
		gen byte multiinf = 0 if anyanyinf==0
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 1 if anybactviral==1
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 2 if toinenbactviral==1
		tab multiinf
		tab multiinf if _st==1
		
		quietly: stdescribe
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf, matcell(x), if _st==1 & inlist(multiinf,0,`i')
		matrix list x
		tab multiinf ensidementia, matcell(y), if _d==1 & inlist(multiinf,0,`i')
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		restore
		
***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

count
tab anyinfcount, mis
tab anyinfcount ensidementia, mis
replace anyinfcount = . if anyanyinf!=0 & ensianyinfpvm>=exitpvm
replace anyinfcount = 0 if anyanyinf==0
replace anyinfcount = 2 if anyinfcount > 2 & anyinfcount < .
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
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		tab any`idg'
		tab any`idg' _d

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox any`idg' i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections
count
tab anyanycns
tab anynoncns
replace anyanycns=0 if ensianycnspvm>=exitpvm
replace anynoncns=0 if ensinoncnspvm>=exitpvm
replace ensianycnspvm=. if ensianycnspvm>=exitpvm
replace ensinoncnspvm=. if ensinoncnspvm>=exitpvm

replace anyanycns = . if anyanycns==0 & anyanyinf==1
tab anyanycns

replace anynoncns = . if anynoncns==0 & anyanyinf==1
tab anynoncns

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
		display "UKB replication analysis"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		tab any`idg'
		tab any`idg' _d
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
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
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		gen byte anyeither = 1 if ensi`idg'pvm < exitpvm
		replace anyeither = 0 if anyanyinf==0
		tab any`idg'
		tab any`idg' _d
		tab anyeither
		tab anyeither _d
	
		stdescribe
		
		gen `idg'_compare = 0 if anyanyinf==0
		replace `idg'_compare = 1 if any`idg'==1
		replace `idg'_compare = 2 if any`idg'==0
		tab `idg'_compare
		tab `idg'_compare _d
		stcox i.`idg'_compare supu ses
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sensitivity dist") modify
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

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'
		preserve
		
		egen lag_entry_multiinf = rowmax(lag_entry_anyinf ensitoinenanyinfpvm ensikolmasanyinfpvm)
		egen multiinfpvm = rowmax(ensianyinfpvm ensitoinenanyinfpvm ensikolmasanyinfpvm)

		stset exitpvm, id(id) failure(dementiatype = 1) origin(time syntpvm) enter(time lag_entry_multiinf) scale(365.25) 
		
		gen byte multiinf = 0
		tab multiinf
		replace multiinf = 1 if anyanyinf==1
		tab multiinf
		replace multiinf = 2 if ensitoinenanyinfpvm!=.
		replace multiinf = 3 if ensikolmasanyinfpvm!=.
		tab multiinf
		
		stdescribe
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3 N1 N2 N3 fail1 fail2 fail3
		forvalues i =1/3 {
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if inlist(multiinf,0,`i')
		quietly: matrix m = r(table)
		matrix list m
		scalar `N`i'' = e(N_sub)
		scalar `fail`i'' = e(N_fail)
		scalar `hr`i'' = m[1,2]
		scalar `ll`i'' = m[5,2]
		scalar `ul`i'' = m[6,2]
		}
		
		forvalues i =1/3 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf, matcell(x), if inlist(multiinf,0,`i') & _st==1
		matrix list x
		tab multiinf ensidementia, matcell(y), if inlist(multiinf,0,`i') & _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		
		restore

		
***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
		egen lag_entry_bactviral = rowmax(lag_entry_anyinf ensibactviralpvm toinenbactviralpvm)
		
		stset exitpvm, id(id) failure(dementiatype = 1) origin(time syntpvm) enter(time lag_entry_bactviral) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses anyanyinf
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral if _st==1
		tab toinenbactviral if _st==1
		
		gen byte multiinf = 0 if anyanyinf==0
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 1 if anybactviral==1
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 2 if toinenbactviral==1
		tab multiinf
		tab multiinf if _st==1
		
		quietly: stdescribe
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf, matcell(x), if _st==1 & inlist(multiinf,0,`i')
		matrix list x
		tab multiinf ensidementia, matcell(y), if _d==1 & inlist(multiinf,0,`i')
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		restore
		
***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

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
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype = 1) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		tab any`idg'
		tab any`idg' _d

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox any`idg' i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections

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
		
		stset exitpvm, id(id) failure(dementiatype = 1) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		tab any`idg'
		tab any`idg' _d
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		local ++irow 
		restore
		}
		
*Difference between CNS infections and extra-CNS infections

	*CNS infection vs non cns infection 
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
		
		stset exitpvm, id(id) failure(dementiatype = 1) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		gen byte anyeither = 1 if ensi`idg'pvm < exitpvm
		replace anyeither = 0 if anyanyinf==0
		tab any`idg'
		tab any`idg' _d
		tab anyeither
		tab anyeither _d
	
		stdescribe
		
		gen `idg'_compare = 0 if anyanyinf==0
		replace `idg'_compare = 1 if any`idg'==1
		replace `idg'_compare = 2 if any`idg'==0
		tab `idg'_compare
		tab `idg'_compare _d
		stcox i.`idg'_compare supu ses
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist AD") modify
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

local dg_name `""Multiple infections""'
local irow 16
local prow = `irow'
		preserve
		
		egen lag_entry_multiinf = rowmax(lag_entry_anyinf ensitoinenanyinfpvm ensikolmasanyinfpvm)
		egen multiinfpvm = rowmax(ensianyinfpvm ensitoinenanyinfpvm ensikolmasanyinfpvm)

		stset exitpvm, id(id) failure(dementiatype = 4) origin(time syntpvm) enter(time lag_entry_multiinf) scale(365.25) 
		
		gen byte multiinf = 0
		tab multiinf
		replace multiinf = 1 if anyanyinf==1
		tab multiinf
		replace multiinf = 2 if ensitoinenanyinfpvm!=.
		replace multiinf = 3 if ensikolmasanyinfpvm!=.
		tab multiinf
		
		stdescribe
		tempname hr1 ll1 ul1 hr2 ll2 ul2 hr3 ll3 ul3 N1 N2 N3 fail1 fail2 fail3
		forvalues i =1/3 {
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if inlist(multiinf,0,`i')
		quietly: matrix m = r(table)
		matrix list m
		scalar `N`i'' = e(N_sub)
		scalar `fail`i'' = e(N_fail)
		scalar `hr`i'' = m[1,2]
		scalar `ll`i'' = m[5,2]
		scalar `ul`i'' = m[6,2]
		}
		
		forvalues i =1/3 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Single vs multiple infections, infection No. `i'")
		tab multiinf, matcell(x), if inlist(multiinf,0,`i') & _st==1
		matrix list x
		tab multiinf ensidementia, matcell(y), if inlist(multiinf,0,`i') & _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		
		local ++irow
		local ++irow
		}
		
		*P for trend
		stcox multiinf i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		
		restore

		
***** One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
		egen lag_entry_bactviral = rowmax(lag_entry_anyinf ensibactviralpvm toinenbactviralpvm)
		
		stset exitpvm, id(id) failure(dementiatype = 4) origin(time syntpvm) enter(time lag_entry_bactviral) scale(365.25) 
		quietly: keep id anybactviral supu cohort exitpvm ensidementia syntpvm entrypvm ensi*pvm *toinen* _st _d _origin _t _t0 ses anyanyinf
		sum ensibactviralpvm toinenbactviralpvm
		tab anybactviral if _st==1
		tab toinenbactviral if _st==1
		
		gen byte multiinf = 0 if anyanyinf==0
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 1 if anybactviral==1
		tab multiinf
		tab multiinf if _st==1
		replace multiinf = 2 if toinenbactviral==1
		tab multiinf
		tab multiinf if _st==1
		
		quietly: stdescribe
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if multiinf!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg_name'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Different micro-organisms, infection No. `i'")
		tab multiinf, matcell(x), if _st==1 & inlist(multiinf,0,`i')
		matrix list x
		tab multiinf ensidementia, matcell(y), if _d==1 & inlist(multiinf,0,`i')
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++irow
		}
		*P for trend
		stcox multiinf 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)
		
		restore
		
***** NUMBER OF SIMULTANEOUS INFECTIONS *******		

count
tab anyinfcount, mis
tab anyinfcount dementiatype, mis
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
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		stset exitpvm, id(id) failure(dementiatype = 4) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		tab any`idg'
		tab any`idg' _d

		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=2
		quietly: matrix m = r(table)
		matrix list m
		tempname hr1 ll1 ul1 hr2 ll2 ul2 N1 fail1 N2 fail2
		scalar `hr1' = m[1,2]
		scalar `ll1' = m[5,2]
		scalar `ul1' = m[6,2]
		scalar `N1' = e(N_sub)
		scalar `fail1' = e(N_fail)
		stcox i.any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort), if any`idg'!=1
		quietly: matrix m = r(table)
		matrix list m		
		scalar `hr2' = m[1,2]
		scalar `ll2' = m[5,2]
		scalar `ul2' = m[6,2]
		scalar `N2' = e(N_sub)
		scalar `fail2' = e(N_fail)
		
		forvalues i =1/2 {
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel K`irow'=(`N`i''), nformat(number)
			quietly: putexcel L`irow'=(`fail`i''), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Simultaneous infection diagnoses")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("No. of simultaneous infection diagnoses = `i'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			local j = `i' + 1
			quietly: putexcel M`irow'=(x[`j',1]), nformat(number)
			quietly: putexcel N`irow'=(y[`j',1]), nformat(number)
		local ++irow
		local ++irow
		}

		*P for trend
		stcox any`idg' i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number_d4)		
		restore
		}	


*FIGURE sensitivity, CNS vs other infections
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
		
		stset exitpvm, id(id) failure(dementiatype = 4) origin(time syntpvm) enter(time lag_entry_`idg') scale(365.25) 
		tab any`idg'
		tab any`idg' _d
		
		stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' i.supu i.ses
		quietly: matrix m = r(table)
		matrix list m
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("`idg_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		local ++irow 
		restore
		}
		
*Difference between CNS infections and extra-CNS infections

	*CNS infection vs non cns infection 
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
		
		stset exitpvm, id(id) failure(dementiatype = 4) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		
		gen byte anyeither = 1 if ensi`idg'pvm < exitpvm
		replace anyeither = 0 if anyanyinf==0
		tab any`idg'
		tab any`idg' _d
		tab anyeither
		tab anyeither _d
	
		stdescribe
		
		gen `idg'_compare = 0 if anyanyinf==0
		replace `idg'_compare = 1 if any`idg'==1
		replace `idg'_compare = 2 if any`idg'==0
		tab `idg'_compare
		tab `idg'_compare _d
		stcox i.`idg'_compare supu ses
		test 1.`idg'_compare = 2.`idg'_compare
		
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("FIGURE sens dist VaD") modify
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

*** eTABLE 4 (former TABLE 2 *** 
capture log close
log using "TABLE2.log", replace

capture gen hypertensio_combpvm = hypert_combpvm

local dg anyinf
local eks `""no one" "hypertensio" "diabetes" "ihd" "cereb" "parkinson"'
local n_dg: word count `dg'
local n_eks: word count `eks' 
local dg_name `""Gram negative bacteria" "Herpesviruses""'
local eks_name `""no one" "hypertension" "diabetes" "ischaemic heart disease" "cerebrovascular disease" "Parkinson´s disease"' 
quietly: local irow = 8
	forvalues i=1/`n_dg' {
	forvalues j=1/`n_eks' {
		preserve
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local ieks : word `j' of `eks'
		quietly: local ieks_name : word  `j' of `eks_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""	
		display "Pooled IPD-analysis, normal model"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		if inlist("`ieks'","hypertensio", "diabetes", "ihd", "cereb", "parkinson") {
		drop if `ieks'_combpvm <= lag_entry_anyinf
		}
		
		stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
			
		if "`ieks'"=="hypertensio" & "`idg'"=="anyinf" {
		local irow = 25
		}
		
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses `ieks'_comb
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Those with `ieks_name' excluded")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		restore
}
}

*Early- vs late-onset dementia
		local idg anyinf
		quietly: stset exitpvm, id(id) failure(ensidementia) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		local idg anyinf
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		
		stsplit age65, at(64.999)
		recode age65 64.999=1
 
		preserve
		keep if age65==0
		local irow 12
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Early-onset dementia (onset before age 65)")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore
		
		preserve
		keep if age65==1
		local irow 13
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Late-onset dementia (onset at or after age 65)")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		restore
 
			
*Type of dementia
	use "UKB_lag.dta", clear
	tab dementiatyyppi_laaja anyanyinf if ensidementia==1, col chi
	local idg anyinf
	local irow 15
	local j 1
	foreach type in AD FTD PD VD other unspecified  {
		preserve
		
		display ""
		display "`type'"
		display ""
		local idg gramminus
		drop if ses==.
		
		stset exitpvm, id(id) failure(dementiatype==`j') origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25)
		keep if _st==1
		
		local idg anyinf
		quietly: capture keep id any`idg' supu cohort exitpvm ensidementia syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort, strata(cohort)
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("`type'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' _d, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		local ++j
		restore
		}

*Adjustments

*Potential confounders	
local dg anyinf
local eks alcocl smoking bmi_who all
local n_dg: word count `dg'
local n_eks: word count `eks' 
local dg_name `""Any infectious disease""'
local eks_name `""heavy drinking" "smoking" "body mass index" "all three""'
quietly: local irow = 31
	forvalues i=1/`n_dg' {
	forvalues j=1/`n_eks' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		quietly: local ieks : word `j' of `eks'
		quietly: local ieks_name : word  `j' of `eks_name'
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		preserve
		drop if ses==.

		if "`ieks'" == "all" {
		drop if alcocl==. | smoking==. | bmi_who==.
		count
		}
		else {
		tab `ieks'
		drop if `ieks' ==.
		}
		stset exitdate, id(id) failure(dementia) origin(time birthday) enter(time lag_entry_`idg') scale(365.25)
		keep if _st==1 & ses!=.
		tab any`idg'
		tab any`idg' _d
		
		*Those with data on the covariate available but without adjustment
		stdescribe
		stcox any`idg' i.sex i.ses 
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel D`irow'=("Data available for `ieks_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow
		
		*Those with data on the covariate available and adjusted
		stdescribe
		quietly: local failures = `r(N_fail)'
		if "`ieks'" == "all" {
		stcox any`idg' sex i.ses i.alcocl i.smoking i.bmi_who
		}
		else  {
		stcox any`idg' sex i.ses i.`ieks'
		}
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel D`irow'=("Additionally adjusted for `ieks_name'")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local ++irow 
		local ++irow
		restore
}
}

*Potential period effect

local dg anyinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease""'
local irow 37
	forvalues i=1/`n_dg' {
		quietly: local idg : word `i' of `dg'
		quietly: local idg_name : word  `i' of `dg_name'
		display ""
		display "`idg' followed by broad dementia, adjusted for sex and education and stratified by cohort"
		display ""
		display "Pooled IPD-analysis, potential period effect"
		display ""
		display "$S_TIME  $S_DATE"
		display ""
		
		drop if ses==.
		
		quietly: stset exitpvm, id(id) failure(failure ==1) origin(time syntpvm) enter(time lag_entry_anyinf) scale(365.25) 
		keep if _st==1
		quietly: keep id any`idg' supu cohort exitpvm ensidementia failure syntpvm entrypvm ensi`idg'pvm  _st _d _origin _t _t0 ses syntvuosi
		tab any`idg'
		
		quietly: stdescribe
		quietly: local failures = `r(N_fail)'
		stcox any`idg' 2.supu#cohort 2.ses#cohort 3.ses#cohort c.syntvuosi#cohort, strata(cohort)
		
		quietly: matrix m = r(table)
		tempname hr ll ul
		scalar `hr' = m[1,1]
		scalar `ll' = m[5,1]
		scalar `ul' = m[6,1]
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr'), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll'), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul'), nformat(number_d2)
			quietly: putexcel B`irow'=("`idg'")
			quietly: putexcel C`irow'="UKB"
			quietly: putexcel D`irow'=("Adjusted for period effect")
		tab any`idg', matcell(x), if _st==1
		matrix list x
		quietly: tab any`idg' ensidementia, matcell(y), if _d==1
		matrix list y
		quietly: putexcel set "SipilaUKB_$S_DATE.xlsm", sheet("TABLE2") modify
			quietly: putexcel M`irow'=(x[2,1]), nformat(number)
			quietly: putexcel N`irow'=(y[2,1]), nformat(number)
		local irow 65
		}
log close




******** TIME-DEPENDENT ANALYSIS **************

capture log close
log using "Sipila, infections and dementia, UKB replication $S_DATE, time-dependent.log", replace

***** eFIGURE 12 (TIME-DEPENDENT VERSION OF FIGURE 3) ********

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
		di "MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes"
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
		di "MODEL 6 = model 5 + adjusted for apoe4 genotype"
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
			quietly: putexcel C`irow'=("UKB")
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

	
**** By type of dementia *****

use "SipilaELSA.dta", clear

local dg anyinf bactinf viralinf
local n_dg: word count `dg'
local dg_name `""Any infectious disease" "Bacterial infection" "Viral infection""'

*Dementiatyypit 1 = AD, 2 = V
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
		MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes
		MODEL 5 = model 4 + those with data on apoe4 missing dropped
		MODEL 6 = model 5 + adjusted for apoe4 genotype
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
		di "MODEL 4 = model 3 + adjusted for alcohol, smoking, bmim hypertension, diabetes"
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
		di "MODEL 6 = model 5 + adjusted for apoe4 genotype"
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
			quietly: putexcel C`irow'=("UKB")
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
log using "Sipila, UKB infection burden $S_DATE, time-dependent.log", replace


***** TIME-DEPENDENT: SINGLE VS MULTIPLE INFECTIONS *******		

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

		
*****TIME-DEPENDENT:  One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
			quietly: putexcel C`irow'=("UKB")
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
		
*****TIME-DEPENDENT:  NUMBER OF SIMULTANEOUS INFECTIONS *******		

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
		display "Pooled IPD-analysis, normal model"
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
		
		local irow 24
		
		restore
		}	


*TIME-DEPENDENT: FIGURE sensitivity, CNS vs other infections
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
		local ++irow 
		restore
		}
		
*TIME-DEPENDENT: Difference between CNS infections and extra-CNS infections

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
		






















****TIME-DEPENDENT:  Analysis ALZHEIMER'S DISEASE *****



*****TIME-DEPENDENT:  SINGLE VS MULTIPLE INFECTIONS *******		

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

		
*****TIME-DEPENDENT:  One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
			quietly: putexcel C`irow'=("UKB")
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
		
		

*****TIME-DEPENDENT:  NUMBER OF SIMULTANEOUS INFECTIONS *******		

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
		display "Pooled IPD-analysis, normal model"
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
		
		local irow 24
		restore
		}	


*TIME-DEPENDENT: FIGURE sensitivity, CNS vs other infections
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
		local ++irow 
		restore
		}
		
*TIME-DEPENDENT: Difference between CNS infections and extra-CNS infections

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
		
		













**** TIME-DEPENDENT: Analysis VASCULAR DEMENTIA *****



*****TIME-DEPENDENT:  SINGLE VS MULTIPLE INFECTIONS *******		

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
		
		
		*Dementia occurring from year 5 onwards
		di "Dementia occurring from year 5 onwards"
		local i 1
		foreach idg in anyinf toinenanyinf kolmasanyinf {
		display `i'
		replace _t0 = _t0 + 5 if multiinf_timedep==0 & `i'==1
		drop if _t0>=_t
		gen ensi`idg'pvm5v = ensi`idg'pvm + round(5*365.25) if ensi`idg'pvm<entrypvm
		replace ensi`idg'pvm5v = entrypvm if ensi`idg'pvm5v<entrypvm
		replace _t0 = (ensi`idg'pvm5v-syntpvm)/365.25 if ensi`idg'pvm5v!=. & multiinf_timedep==`i'
		replace _t0 = _t0 + 5 if ensi`idg'pvm5v==. & multiinf_timedep==`i'
		drop if _t0>=_t
		local ++i
		}
		
		stdescribe
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
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel K`irow'=(`e(N_sub)'), nformat(number)
			quietly: putexcel L`irow'=(`e(N_fail)'), nformat(number)
			quietly: putexcel E`irow'=(`hr`i''), nformat(number_d2)
			quietly: putexcel F`irow'=(`ll`i''), nformat(number_d2)
			quietly: putexcel G`irow'=(`ul`i''), nformat(number_d2)
			quietly: putexcel B`irow'=("Multiinf")
			quietly: putexcel C`irow'=("UKB")
			quietly: putexcel D`irow'=("Dementia occurring more than 5 years after the hospitalisation, infection No. `i'")
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
		local ++prow
		stcox multiinf_timedep i.supu i.ses
		quietly: matrix n = r(table)
		matrix list n
		tempname p
		scalar `p' = n[4,1]
		quietly: putexcel set "UKB infection burden $S_DATE.xlsm", sheet("FIGURE sensitivity VaD") modify
			quietly: putexcel O`prow'=(`p'), nformat(number)
		restore

*****TIME-DEPENDENT:  One or several different pathogens *******		

local dg_name `""Multiple organisms""'
local idg_name `""Multiple_organisms""'
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
			quietly: putexcel C`irow'=("UKB")
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
		

*****TIME-DEPENDENT:  NUMBER OF SIMULTANEOUS INFECTIONS *******		

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
		display "Pooled IPD-analysis, normal model"
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


*TIME-DEPENDENT: FIGURE sensitivity, CNS vs other infections
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
		
*TIME-DEPENDENT: Difference between CNS infections and extra-CNS infections

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
		
		restore
		}
		
		
log close
