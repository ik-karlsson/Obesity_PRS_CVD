////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_analyses_20220214_late
Purpose: Survival analysis of BMI on CVD, and the effect of PRS for BMI
Study: SfoEpi grant 

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20220214 by Ida Karlsson (IK)
 	Updated: 	20221021 by IK and Elsa Ojalehto - use updated data after removing those <40 (only affecting midlife analyses)
								 and refining the studypop 	
				20221212 by IK - use updated data including those <60 at end of follow-up (based on rev comments)
								 Affects only midlife analyses, but also changing the descriptives to only
								 include those in ST
	
STATA v17	
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log using session, replace
set more off

use "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Data\GOSH\GPdisc_CVD_STR_adata_20221212.dta", clear

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										// DATA MANAGEMENT //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Dropping those missing latelife BMI
drop if BMIlate==.
tab sex, mi
// remaining n=6592

///// Look at distributions
summarize BMIlate
graph box BMIlate

/// CAtegorize
gen BMIlate_cat=1 if BMIlate>=18.5 & BMIlate <25
replace BMIlate_cat=2 if BMIlate>=25 & BMIlate <30
replace BMIlate_cat=3 if BMIlate!=. & BMIlate >=30

tab BMIlate_cat, mi

/// Categorize the PRSs
graph box PRS_BMI_PCadj, by(study)
// SALT-Y looks odd.. 
tab study
// Only 16 in SALT-Y - will add them to TwinGene for standardization and study adjustment as they are both from SALT, and the PRSs are close to identical in the full samples
replace study="TwinGene" if study=="SALT-Y"

egen PRS_BMI_grpm= mean(PRS_BMI_PCadj), by(study)
egen PRS_BMI_grpsd= sd(PRS_BMI_PCadj), by(study)
gen z_PRS_BMI_PCadj= (PRS_BMI_PCadj-PRS_BMI_grpm)/PRS_BMI_grpsd

// Tertiles
xtile tert_PRS_BMI = z_PRS_BMI_PCadj, nquantiles(3)
graph box z_PRS_BMI_PCadj, by(tert_PRS_BMI)

// Regress the PRS on BMI
regress BMIlate z_PRS_BMI_PCadj BMIlate_age sex,vce(cluster pairid)
/*----------------------------------------------------------------------------
             |               Robust
     BMIlate | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
z_PRS_BMI_~j |   1.165453   .0505943    23.04   0.000     1.066262    1.264644
------------------------------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							///// SETTING UP THE SURVIVAL DATA, ANY CVD /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate age at end of follow-up
gen age_exit=lastage if reg_CVDany==0
replace age_exit=reg_CVDany_age if reg_CVDany==1

// Setting up the survival data with cvd as outcome
stset age_exit, failure(reg_CVDany==1) id(twinnr) enter(BMIlate_age)
/*------------------------------------------------------------------------
      6,592  total observations
      1,104  observations end on or before enter()
--------------------------------------------------------------------------
      5,488  observations remaining, representing
      5,488  subjects
      1,681  failures in single-failure-per-subject data
     50,614  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        65
                                          Last observed exit t =       104 */
stdes
/*
                                   |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects          5488   
Number of records           5488           1           1          1          1

Entry time (first)                  71.75711          65         71         94
Exit time (final)                   80.97977          66         80        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk               50614    9.222668           1         10         36

Failures                    1681    .3063047           0          0          1
------------------------------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									/////  DESCRIPTIVE STATISTICS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Convert education
gen educ_r=1 if educ==0
replace educ_r=0 if educ==1

cap postclose descriptive
	postfile descriptive str20 (exp_var all all_ps ctrl ctrl_ps cvd_case cvd_case_ps cvd_p_value) ///
	using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_late_descriptive.dta, replace

	foreach exp_var in sex educ_r smoke {
			set more off
			// All
			estpost tab `exp_var' if(_st==1)
				mat b = e(b)
				mat pct = e(pct)
				local all = string(b[1,2],"%9.0f")  
				local all_ps = "("+string(pct[1,2],"%9.2f")+")" 
			// CVD
			estpost tab `exp_var' reg_CVDany if(_st==1), chi2 nototal
				mat b = e(b)
				mat colpct = e(colpct)
				local ctrl = string(b[1,2],"%9.0f")  
				local ctrl_ps = "("+string(colpct[1,2],"%9.2f")+")"    
				local cvd_case = string(b[1,4],"%9.0f")  
				local cvd_case_ps = "("+string(colpct[1,4],"%9.2f")+")"    
				local cvd_p_value = string(e(p))  
				
         post descriptive 	("`exp_var'") ("`all'") ("`all_ps'") ///
							("`ctrl'") ("`ctrl_ps'") ("`cvd_case'") ("`cvd_case_ps'") ("`cvd_p_value'") 
        }
	
		foreach exp_var in BMIlate_age lastage age_death BMIlate PRS_BMI_PCadj {
			set more off
			// All
			mean `exp_var' if(_st==1)
				mat b = e(b)
				mat sd  = e(sd)
				local all = string(b[1,1],"%9.2f")  
				local all_ps = "("+string(sd[1,1],"%9.2f")+")"    
			// Dementia
			ttest `exp_var' if(_st==1), by(reg_CVDany)
				local ctrl = string(r(mu_1),"%9.2f")  
				local ctrl_ps = "("+string(r(sd_1),"%9.2f")+")"    
				local cvd_case = string(r(mu_2),"%9.2f")  
				local cvd_case_ps = "("+string(r(sd_2),"%9.2f")+")"    
				local cvd_p_value = string(r(p))  

		post descriptive ("`exp_var'") ("`all'") ("`all_ps'") ("`ctrl'") ("`ctrl_ps'") ///
						 ("`cvd_case'") ("`cvd_case_ps'") ("`cvd_p_value'") 
		} 
postclose descriptive

// Separate look at age of onset (only for those in _st)
summarize reg_CVDany_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_CVDany~e |      1,681    80.24807    6.739788         66        102 */

tab BMIlate_cat reg_CVDany if(_st==1), column mi
/*

BMIlate_cat|      reg_CVDany
           |         0          1 |     Total
-----------+----------------------+----------
         1 |     1,737        685 |     2,422 
           |     45.63      40.75 |     44.13 
-----------+----------------------+----------
         2 |     1,601        756 |     2,357 
           |     42.05      44.97 |     42.95 
-----------+----------------------+----------
         3 |       469        240 |       709 
           |     12.32      14.28 |     12.92 
-----------+----------------------+----------
     Total |     3,807      1,681 |     5,488 
           |    100.00     100.00 |    100.00     */
		   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// SURVIVAL ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// test ph assumption
// K-M curves and log-log plots
// By BMI
sts graph, by(BMIlate_cat) 
stphplot, by(BMIlate_cat) 
// Questionable..

// Check with ph-test
stcox i.BMIlate_cat i.sex educ smoke, vce(cluster pairid)
estat phtest, detail
// Clear issue with sex here too, but smoking looks fine! Will stratify on sex anyway, and will include both in "strata"

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c1_ob PRS_c2_ow PRS_c2_ob PRS_c3_ow PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_late_survCVDany.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
	}		
postclose coxreg_cvd 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								///// SURVIVAL ANALYSES, non-stroke CVD /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

stset, clear
drop age_exit

// Generate age at end of follow-up
gen age_exit=lastage if reg_CVDns==0
replace age_exit=reg_CVDns_age if reg_CVDns==1
// NOTE! Censoring individuals at stroke
replace age_exit=reg_stroke_age if reg_stroke==1 & reg_stroke_age < age_exit

// Setting up the survival data with cvd as outcome
stset age_exit, failure(reg_CVDns==1) id(twinnr) enter(BMIlate_age)
/*--------------------------------------------------------------------------
       6,592  total observations
      1,104  observations end on or before enter()
--------------------------------------------------------------------------
      5,488  observations remaining, representing
      5,488  subjects
      1,246  failures in single-failure-per-subject data
     50,614  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        65
                                          Last observed exit t =       104 */
stdes
/*                                   |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects          5488   
Number of records           5488           1           1          1          1

Entry time (first)                  71.75711          65         71         94
Exit time (final)                   80.97977          66         80        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk               50614    9.222668           1         10         36

Failures                    1246    .2270408           0          0          1
------------------------------------------------------------------------------
*/

// Separate look at age of onset (only for those in _st)
summarize reg_CVDns_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_CVDns_~e |      1,246    80.64526    6.902372         66        102 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c1_ob PRS_c2_ow PRS_c2_ob PRS_c3_ow PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_late_survCVDns.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
	}		
postclose coxreg_cvd 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								///// SURVIVAL ANALYSES, STROKE /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

stset, clear
drop age_exit

// Generate age at end of follow-up
gen age_exit=lastage if reg_stroke==0
replace age_exit=reg_stroke_age if reg_stroke==1

// Setting up the survival data with cvd as outcome
stset age_exit, failure(reg_stroke==1) id(twinnr) enter(BMIlate_age)
/*------------------------------------------------------------------------
      6,592  total observations
        353  observations end on or before enter()
--------------------------------------------------------------------------
      6,239  observations remaining, representing
      6,239  subjects
        831  failures in single-failure-per-subject data
     61,236  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        65
                                          Last observed exit t =       104 */
										  
stdes
/*                                 |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects          6239   
Number of records           6239           1           1          1          1

Entry time (first)                  71.98349          65         72         94
Exit time (final)                   81.79853          66         81        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk               61236    9.815034           1         11         36

Failures                     831    .1331944           0          0          1
------------------------------------------------------------------------------ */

// Separate look at age of onset (only for those in _st)
summarize reg_stroke_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_stroke~e |        831    80.79783    6.432976         66         98 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c1_ob PRS_c2_ow PRS_c2_ob PRS_c3_ow PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_late_survSTROKE.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"  
		
	// PRSbmi model
	stcox z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"  

	// BMI + PRS model
	stcox i.BMIlate_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"  
		
	// Interaction model
	stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ow'") ("`PRS_c2_ob'")  ("`PRS_c3_ow'") ("`PRS_c3_ob'") 
	}		
postclose coxreg_cvd 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20221212_late.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

