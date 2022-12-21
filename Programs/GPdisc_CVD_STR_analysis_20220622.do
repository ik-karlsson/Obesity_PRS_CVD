////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: BMI_dem_analysis_GOSH_survival
Purpose: Survival analysis of BMI on dementia, with BMI as time-dependent exposure
Study: FORTE postdoc aim 1: Study how genetic factors influence the association between 
longitudinal changes in BMI (cognitive decline across domains) and dementia.

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20190325 by Ida Karlsson (IK)
	Updated: 	20190411 by IK - Updated data, optimized code for time-varying effects
				20190604 by IK - Updated data, adjust for smoking. PRSs already adjusted for PCs.
				20190912 by IK - Modify to look at BMI in 10-year age bands, to correspond to HRS
								 Better format of output.
				20190918 by IK - Go back to 15 year age spans. 
								 Use mean cetered BMI instead of standardized, presenting by 5 unit differennce.
				20191004 by IK - Add p-value for the interaction; drop those with entry age = exit age early on
								 Correct how ID is created
				20191025 by IK - Updated data after fixing bug in study sample

Input: P:\Dementia_IK\FORTE_postdoc\Aim_1\Data\Derived_data\BMI_dem_adata_GOSH_20190604.dta

STATA v15.1		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										// DATA MANAGEMENT //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log using session, replace
set more off

use "C:\Users\idamat\Desktop\Temp\GPdisc_CVD_STR_adata_20220622_l.dta", clear

// Drop those missing for dementia and covariates..
drop if reg_CVDany==.
//n=0
drop if smoke==.
//=28
drop if educ==.
//n=6

///// Generating age categories
sum BMI_age

///// Generating age categories
gen agecat=1 if BMI_age >=20 & BMI_age <35
replace agecat=2 if BMI_age >=35 & BMI_age <50
replace agecat=3 if BMI_age >=50 & BMI_age <65
replace agecat=4 if BMI_age >=65 & BMI_age <80
replace agecat=5 if BMI_age >=80

label define agecatl 1 "20-34" 2 "35-49" 3 "50-64" 4 "65-79" 5 ">80" 
label values agecat agecatl  

tab agecat

/*   agecat |      Freq.     Percent        Cum.
------------+-----------------------------------
      20-34 |      7,311       17.05       17.05
      35-49 |      8,099       18.89       35.93
      50-64 |     16,381       38.20       74.13
      65-79 |      9,230       21.52       95.65
        >80 |      1,864        4.35      100.00
------------+-----------------------------------
      Total |     42,885      100.00
*/

// Look at mean age and SD
sort agecat

by agecat: egen age_mean= mean(BMI_age)
by agecat: egen age_sd  = sd(BMI_age)

tab age_mean agecat
tab age_sd agecat

/// Mean center BMI and divide by 5, within each age category
sort agecat

by agecat: egen BMI_mean= mean(BMI)
by agecat: gen  z_bmi = (BMI-BMI_mean)/5

tab BMI_mean agecat

// Look at SD
by agecat: egen BMI_sd  = sd(BMI)
tab BMI_sd agecat

// ID within each age category
gen double ID=twinnr*10+agecat

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								///// SETTING UP THE SURVIVAL DATA /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate age at end of follow-up
replace lastage =reg_CVDany_age if reg_CVDany==1
drop if BMI_age==lastage
//466 deleted
drop if BMI_age > lastage
// 3 deleted

// Assign age_exit as next BMI assessment, except last assessment within each age category
sort ID BMI_age
by ID: gen age_exit=BMI_age[_n+1]

list twinnr reg_CVDany agecat BMI_age reg_CVDany_age lastage age_exit  in 1/20, noobs

// Assign age_exit as last_age for the last BMI measure (or only) within the same age category
replace age_exit=lastage if age_exit==.

// Recoding so that only CVD within each follow-up period is coded as 1
gen cvd_fup=0 if reg_CVDany==0
replace cvd_fup=0 if (reg_CVDany==1 & reg_CVDany_age > age_exit) 
replace cvd_fup=1 if (reg_CVDany==1 & reg_CVDany_age > BMI_age & reg_CVDany_age <= age_exit) 

tab reg_CVDany cvd_fup
// n=X events

list twinnr ID reg_CVDany agecat BMI_age reg_CVDany_age lastage age_exit cvd_fup if reg_CVDany==1 in 1/100, noobs

// Setting up the survival data with dementia as outcome
stset age_exit, failure(cvd_fup==1) id(ID) enter(BMI_age)

list twinnr agecat BMI_age reg_CVDany_age age_exit _t0 _t _d reg_CVDany cvd_fup _st in 1/20, noobs

stdes


//list twinnr agecat BMI_age reg_CVDany_age age_exit _t0 _t _d reg_CVDany cvd_fup _st if _st!=1, noobs
drop if _st!=1
//n=1138; all under age 20 at BMI measurement

// Number of individuals
egen uniqind = nvals(twinnr)
tab uniqind
// NOTE! Not working - complains about nvals?

// Number of individuals per age category
egen uniqind_age = nvals(twinnr), by(agecat) 
tab uniqind_age agecat
/*
uniqind_ag |                         agecat
         e |     20-35      35-50      50-65      65-80        >80 |     Total
-----------+-------------------------------------------------------+----------
      2569 |         0          0          0          0      4,396 |     4,396 
      8101 |     8,101          0          0          0          0 |     8,101 
     10425 |         0     15,911          0          0          0 |    15,911 
     11941 |         0          0     18,753          0          0 |    18,753 
     13227 |         0          0          0     17,658          0 |    17,658 
-----------+-------------------------------------------------------+----------
     Total |     8,101     15,911     18,753     17,658      4,396 |    64,819 
*/

// Number of events and follow-up time by age category 
stdes if agecat==1
stdes if agecat==2
stdes if agecat==3
stdes if agecat==4
stdes if agecat==5

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// SURVIVAL ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////// Any CVD

cap postclose coxreg_dem
	postfile coxreg_dem str20  (agecat base PRSbmi PRSbmi_c1 PRSbmi_c2 PRSbmi_c3 PRSbmi_pval) ///
								using C:\Users\idamat\Desktop\Temp\GPdisc_CVD_STR_analysis_20220622.dta, replace
foreach i of num 1/5 {
	set more off
	// Base model
	stcox z_bmi i.sex educ smoke if(agecat==`i'), vce(cluster pairid) strata(study)
	lincom z_bmi 
		local m_base = r(estimate)
		local se_base = r(se)
		local Lci_base = exp(`m_base'-1.96*`se_base')
		local Uci_base = exp(`m_base'+1.96*`se_base')
		local base= string(exp(`m_base'),"%9.2f")+" ("+string(`Lci_base',"%9.2f")+"-"+string(`Uci_base',"%9.2f")+")"  

	// PRSbmi model
	stcox z_bmi i.sex educ smoke PRS_BMI_PCadj if(agecat==`i'), vce(cluster pairid) strata(study)
	lincom z_bmi 
		local m_PRSbmi = r(estimate)
		local se_PRSbmi = r(se)
		local Lci_PRSbmi = exp(`m_PRSbmi'-1.96*`se_PRSbmi')
		local Uci_PRSbmi = exp(`m_PRSbmi'+1.96*`se_PRSbmi')
		local PRSbmi= string(exp(`m_PRSbmi'),"%9.2f")+" ("+string(`Lci_PRSbmi',"%9.2f")+"-"+string(`Uci_PRSbmi',"%9.2f")+")"  

	/// By tertiles of PRSbmi
	stcox tert_PRS_BMI i.tert_PRS_BMI#c.z_bmi i.sex educ smoke if(agecat==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#c.z_bmi
		local m_PRSbmi_c1 = r(estimate)
		local se_PRSbmi_c1 = r(se)
		local Lci_PRSbmi_c1 = exp(`m_PRSbmi_c1'-1.96*`se_PRSbmi_c1')
		local Uci_PRSbmi_c1 = exp(`m_PRSbmi_c1'+1.96*`se_PRSbmi_c1')
		local PRSbmi_c1= string(exp(`m_PRSbmi_c1'),"%9.2f")+" ("+string(`Lci_PRSbmi_c1',"%9.2f")+"-"+string(`Uci_PRSbmi_c1',"%9.2f")+")"  

	/// 2nd quartile
	lincom 2.tert_PRS_BMI#c.z_bmi
		local m_PRSbmi_c2 = r(estimate)
		local se_PRSbmi_c2 = r(se)
		local Lci_PRSbmi_c2 = exp(`m_PRSbmi_c2'-1.96*`se_PRSbmi_c2')
		local Uci_PRSbmi_c2 = exp(`m_PRSbmi_c2'+1.96*`se_PRSbmi_c2')
		local PRSbmi_c2 = string(exp(`m_PRSbmi_c2'),"%9.2f")+" ("+string(`Lci_PRSbmi_c2',"%9.2f")+"-"+string(`Uci_PRSbmi_c2',"%9.2f")+")"  

	/// 3rd quartile
	lincom 3.tert_PRS_BMI#c.z_bmi
		local m_PRSbmi_c3 = r(estimate)
		local se_PRSbmi_c3 = r(se)
		local Lci_PRSbmi_c3 = exp(`m_PRSbmi_c3'-1.96*`se_PRSbmi_c3')
		local Uci_PRSbmi_c3 = exp(`m_PRSbmi_c3'+1.96*`se_PRSbmi_c3')
		local PRSbmi_c3= string(exp(`m_PRSbmi_c3'),"%9.2f")+" ("+string(`Lci_PRSbmi_c3',"%9.2f")+"-"+string(`Uci_PRSbmi_c3',"%9.2f")+")"  
	
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##c.z_bmi i.sex educ smoke if(agecat==`i'), vce(cluster pairid) strata(study)
	lincom c.tert_PRS_BMI#c.z_bmi
	local PRSbmi_pval = string(r(p),"%9.2f")
	
	post coxreg_dem ("`i'") ("`base'") ("`PRSbmi'") ("`PRSbmi_c1'") ("`PRSbmi_c2'") ("`PRSbmi_c3'") ("`PRSbmi_pval'") 
}
postclose coxreg_dem


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_1\Logs\Analyses\GOSH\BMI_dem_analysis_GOSH_survival_20191025.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

