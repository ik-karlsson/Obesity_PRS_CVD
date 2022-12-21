////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_analyses_20220214_mid
Purpose: Survival analysis of BMI on CVD, and the effect of PRS for BMI
Study: SfoEpi grant 

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20220214 by Ida Karlsson (IK)
	Updated: 	20221021 by IK - use updated data after removing those <40 and refining the studypop
				20221212 by IK - use updated data including those <60 at end of follow-up (based on rev comments)
								 Refine descriptives to only include those in ST
	
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

/// Dropping those missing midlife BMI
drop if BMImid==.
tab sex, mi
// remaining n=16123

///// Look at distributions
summarize BMImid
graph box BMImid

/// CAtegorize
gen BMImid_cat=1 if BMImid>=18.5 & BMImid <25
replace BMImid_cat=2 if BMImid>=25 & BMImid <30
replace BMImid_cat=3 if BMImid!=. & BMImid >=30

tab BMImid_cat, mi

/// Categorize the PRSs
graph box PRS_BMI_PCadj, by(study)
// Looks good, but will standardize within study before categorizing

egen PRS_BMI_grpm= mean(PRS_BMI_PCadj), by(study)
egen PRS_BMI_grpsd= sd(PRS_BMI_PCadj), by(study)
gen z_PRS_BMI_PCadj= (PRS_BMI_PCadj-PRS_BMI_grpm)/PRS_BMI_grpsd

// Tertiles
xtile tert_PRS_BMI = z_PRS_BMI_PCadj, nquantiles(3)
graph box z_PRS_BMI_PCadj, by(tert_PRS_BMI)

// Regress the PRS on BMI
regress BMImid z_PRS_BMI_PCadj BMImid_age sex,vce(cluster pairid)
/*----------------------------------------------------------------------------
             |               Robust
      BMImid | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
z_PRS_BMI_~j |   1.123931   .0281897    39.87   0.000     1.068674    1.179188
 ------------------------------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							///// SETTING UP THE SURVIVAL DATA, BMI MID /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate age at end of follow-up
gen age_exit=lastage if reg_CVDany==0
replace age_exit=reg_CVDany_age if reg_CVDany==1

// Setting up the survival data with cvd as outcome
stset age_exit, failure(reg_CVDany==1) id(twinnr) enter(BMImid_age)
/*------------------------------------------------------------------------
     16,123  total observations
        337  observations end on or before enter()
--------------------------------------------------------------------------
     15,786  observations remaining, representing
     15,786  subjects
      3,123  failures in single-failure-per-subject data
    278,705  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        40
                                          Last observed exit t =       104	*/
										  
stdes
/*
                                   |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects         15786   
Number of records          15786           1           1          1          1

Entry time (first)                  52.08419          40         52         64
Exit time (final)                   69.73939          45         69        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk              278705     17.6552           1         16         62

Failures                    3123    .1978335           0          0          1
------------------------------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									/////  DESCRIPTIVE STATISTICS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Convert education
gen educ_r=1 if educ==0
replace educ_r=0 if educ==1

cap postclose descriptive
	postfile descriptive str20 (exp_var all all_ps ctrl ctrl_ps cvd_case cvd_case_ps cvd_p_value) ///
	using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_mid_descriptive.dta, replace

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
	
		foreach exp_var in BMImid_age lastage age_death BMImid PRS_BMI_PCadj {
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

// Separate look at CVD onset
summarize reg_CVDany_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_CVDany~e |      3,123    70.03138    10.29194         45        102 */

tab BMImid_cat reg_CVDany if(_st==1), column mi
/*
           |      reg_CVDany
BMImid_cat |         0          1 |     Total
-----------+----------------------+----------
         1 |     7,361      1,570 |     8,931 
           |     58.13      50.27 |     56.58 
-----------+----------------------+----------
         2 |     4,433      1,274 |     5,707 
           |     35.01      40.79 |     36.15 
-----------+----------------------+----------
         3 |       869        279 |     1,148 
           |      6.86       8.93 |      7.27 
-----------+----------------------+----------
     Total |    12,663      3,123 |    15,786 
           |    100.00     100.00 |    100.00    */
		   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// SURVIVAL ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// test ph assumption
// K-M curves and log-log plots
// By BMI
sts graph, by(BMImid_cat) 
stphplot, by(BMImid_cat) 
// Questionable..

// Check with ph-test
stcox i.BMImid_cat i.sex educ smoke, vce(cluster pairid)
estat phtest, detail
// Clear issue with sex and smoking! Will stratify on sex anyway, and will include both in "strata"

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_pval_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob PRS_pval_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_mid_survCVDany.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
		
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'") 
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
stset age_exit, failure(reg_CVDns==1) id(twinnr) enter(BMImid_age)
/*--------------------------------------------------------------------------
     16,123  total observations
        337  observations end on or before enter()
--------------------------------------------------------------------------
     15,786  observations remaining, representing
     15,786  subjects
      2,414  failures in single-failure-per-subject data
    278,705  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        40
                                          Last observed exit t =       104
*/

stdes
/*                                   |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects         15786   
Number of records          15786           1           1          1          1

Entry time (first)                  52.08419          40         52         64
Exit time (final)                   69.73939          45         69        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk              278705     17.6552           1         16         62

Failures                    2414    .1529203           0          0          1
------------------------------------------------------------------------------
*/

// age at onset
summarize reg_CVDns_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_CVDns_~e |      2,414    70.13836    10.55371         45        102 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_pval_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob PRS_pval_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_mid_survCVDns.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
			
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
			
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")
	
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'")
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
stset age_exit, failure(reg_stroke==1) id(twinnr) enter(BMImid_age)
/*------------------------------------------------------------------------
     16,123  total observations
         89  observations end on or before enter()
--------------------------------------------------------------------------
     16,034  observations remaining, representing
     16,034  subjects
      1,156  failures in single-failure-per-subject data
    298,485  total analysis time at risk and under observation
                                                At risk from t =         0
                                     Earliest observed entry t =        40
                                          Last observed exit t =       104
 */
stdes
/*                                 |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects         16034   
Number of records          16034           1           1          1          1

Entry time (first)                  52.11669          40         52         64
Exit time (final)                   70.73244          47         70        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk              298485    18.61575           1         16         62

Failures                    1156    .0720968           0          0          1
------------------------------------------------------------------------------ */

// age at onset
summarize reg_stroke_age if(_st==1)
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
reg_stroke~e |      1,156    72.60986    10.16236         47         98 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cap postclose coxreg_cvd
	postfile coxreg_cvd str20  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_pval_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob PRS_pval_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20221212_mid_survSTROKE.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
			
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ, vce(cluster pairid) strata(study sex smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"   
	
	/// By tertiles of PRSbmi
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"  
	/// 2nd quartile
lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"  
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"  
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"  
	lincom 3.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"  
			
	/// P-value for the interaction
	stcox c.tert_PRS_BMI##i.BMImid_cat educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom c.tert_PRS_BMI#2.BMImid_cat
		local PRS_pval_ow = string(r(p),"%9.2f")
	lincom c.tert_PRS_BMI#3.BMImid_cat
		local PRS_pval_ob = string(r(p),"%9.2f")
	
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_pval_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") ("`PRS_pval_ob'")
	}		 

postclose coxreg_cvd

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20221212_mid.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////




