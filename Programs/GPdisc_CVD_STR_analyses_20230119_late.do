////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_analyses_20220214_late
Purpose: Survival analysis of BMI on CVD, and the effect of PRS for BMI
Study: SfoEpi grant 

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20220214 by Ida Karlsson (IK)
 	Updated: 	20221021 by IK - use updated data after removing those <40 (only affecting midlife analyses)
								 and refining the studypop 	
				20221212 by IK - use updated data including those <60 at end of follow-up (based on rev comments)
								 Affects only midlife analyses, but also changing the descriptives to only
								 include those in ST
				20230119 by IK - Add p-values to tables and adjust PRS -> BMI for smoking and education, based 
								 on reviewer comments. Refine PH test, and use bootstrap to adjust SEs.
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
/*-------------------------------------------------------------------------------
                |               Robust
        BMIlate | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
----------------+----------------------------------------------------------------
z_PRS_BMI_PCadj |   1.165453   .0505943    23.04   0.000     1.066262    1.264644
--------------------------------------------------------------------------------- */

// Adjusted for smoking and education
regress BMIlate z_PRS_BMI_PCadj BMIlate_age sex smoke educ,vce(cluster pairid)
/*-------------------------------------------------------------------------------
                |               Robust
        BMIlate | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
----------------+----------------------------------------------------------------
z_PRS_BMI_PCadj |   1.160436   .0508113    22.84   0.000      1.06082    1.260053
--------------------------------------------------------------------------------- */

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
	using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_late_descriptive.dta, replace

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
sts graph, by(BMIlate_cat) name(km_bmi, replace)
stphplot, by(BMIlate_cat) name(surv_bmi, replace)
// By PRS 
sts graph, by(tert_PRS_BMI) name(km_prs, replace)
stphplot, by(tert_PRS_BMI) name(surv_prs, replace)
// By sex
sts graph, by(sex) name(km_sex, replace)
stphplot, by(sex) name(surv_sex, replace)
// By smoking
sts graph, by(smoke) name(km_smoke, replace)
stphplot, by(smoke) name(surv_smoke, replace)
// By educ
sts graph, by(educ) name(km_educ, replace)
stphplot, by(educ) name(surv_educ, replace)
// Save the plots
graph combine km_bmi surv_bmi km_prs surv_prs km_sex surv_sex km_smoke surv_smoke km_educ surv_educ, iscale(.25)
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_late_PHplots.pdf", replace
// Check with ph-test
stcox i.BMIlate_cat z_PRS_BMI_PCadj i.sex educ smoke, vce(cluster pairid)
estat phtest, detail
/*
--------------------------------------------------------
             |        rho     chi2       df    Prob>chi2
-------------+------------------------------------------
1b.BMIlate~t |          .        .        1           .
2.BMIlate_~t |    0.00015     0.00        1       0.9950
3.BMIlate_~t |   -0.04178     3.11        1       0.0776
z_PRS_BMI_~j |    0.03082     1.77        1       0.1828
      1b.sex |          .        .        1           .
       2.sex |    0.09931    16.26        1       0.0001
        educ |    0.02199     0.83        1       0.3628
       smoke |   -0.02204     0.82        1       0.3639
-------------+------------------------------------------
 Global test |               24.98        6       0.0003
-------------------------------------------------------- */
// Clear issue with sex here too, but others look fine. Will stratify on sex anyway, and use bootstrapping here too to obtain less biased SEs

// Test for colinearity in joint model of BMI and the PRS, using strata for sex and smoke
// With continous BMI
stcox BMIlate z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
estat vce, corr
// Correlation -0.36

// Look at the variance inflation factor from logistic regression (does not translate to Cox)
// Countinous BMI
regress reg_CVDany BMIlate z_PRS_BMI_PCadj educ sex smoke, vce(cluster pairid)
estat vif
/*  Variable |       VIF       1/VIF  
-------------+----------------------

z_PRS_BMI_~j |      1.11    0.897724
     BMIlate |      1.11    0.901679
       smoke |      1.06    0.944433
         sex |      1.06    0.947778
        educ |      1.00    0.995167
-------------+----------------------
    Mean VIF |      1.07			*/
	
// Without the PRS
regress reg_CVDany BMIlate  educ sex smoke, vce(cluster pairid)
estat vif
/*
    Variable |       VIF       1/VIF  
-------------+----------------------
       smoke |      1.05    0.948537
         sex |      1.05    0.948619
        educ |      1.00    0.997113
     BMIlate |      1.00    0.997324
-------------+----------------------
    Mean VIF |      1.03			*/

// Only linear BMI and the PRS
regress reg_CVDany BMIlate z_PRS_BMI_PCadj, vce(cluster pairid)
estat vif
// VIF 1.11

cap postclose coxreg_cvd
	postfile coxreg_cvd str26  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_202230119_late_survCVDany.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'")
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
	postfile coxreg_cvd str26  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_late_survCVDns.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") 
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
	postfile coxreg_cvd str26  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_late_survSTROKE.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj educ, vce(cluster pairid) strata(study sex smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   

	post coxreg_cvd ("All") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'")
					
///// Sex-stratified
	foreach i of num 1/2 {
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMIlate_cat
		local m_BMIind_ob = r(estimate)
		local se_BMIind_ob = r(se)
		local Lci_BMIind_ob = exp(`m_BMIind_ob'-1.96*`se_BMIind_ob')
		local Uci_BMIind_ob = exp(`m_BMIind_ob'+1.96*`se_BMIind_ob')
		local p_BMIind_ob = r(p)
		local BMIind_ob= string(exp(`m_BMIind_ob'),"%9.2f")+" ("+string(`Lci_BMIind_ob',"%9.2f")+"-"+string(`Uci_BMIind_ob',"%9.2f")+")"+", p="+string(`p_BMIind_ob',"%9.3f")
		
	// PRSbmi model
	bootstrap, reps(100) seed(7982): stcox z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat z_PRS_BMI_PCadj educ if(sex==`i'), vce(cluster pairid) strata(study smoke)
	lincom 2.BMIlate_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMIlate_cat
		local m_BMIjnt_ob = r(estimate)
		local se_BMIjnt_ob = r(se)
		local Lci_BMIjnt_ob = exp(`m_BMIjnt_ob'-1.96*`se_BMIjnt_ob')
		local Uci_BMIjnt_ob = exp(`m_BMIjnt_ob'+1.96*`se_BMIjnt_ob')
		local p_BMIjnt_ob = r(p)
		local BMIjnt_ob= string(exp(`m_BMIjnt_ob'),"%9.2f")+" ("+string(`Lci_BMIjnt_ob',"%9.2f")+"-"+string(`Uci_BMIjnt_ob',"%9.2f")+")"+", p="+string(`p_BMIjnt_ob',"%9.3f") 
	lincom z_PRS_BMI_PCadj 
		local m_PRSjnt = r(estimate)
		local se_PRSjnt = r(se)
		local Lci_PRSjnt = exp(`m_PRSjnt'-1.96*`se_PRSjnt')
		local Uci_PRSjnt = exp(`m_PRSjnt'+1.96*`se_PRSjnt')
		local p_PRSjnt = r(p)
		local PRSjnt= string(exp(`m_PRSjnt'),"%9.2f")+" ("+string(`Lci_PRSjnt',"%9.2f")+"-"+string(`Uci_PRSjnt',"%9.2f")+")"+", p="+string(`p_PRSjnt',"%9.3f") 
		
	// Interaction model
	bootstrap, reps(100) seed(7982): stcox i.BMIlate_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMIlate_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMIlate_cat
		local m_BMIint_ob = r(estimate)
		local se_BMIint_ob = r(se)
		local Lci_BMIint_ob = exp(`m_BMIint_ob'-1.96*`se_BMIint_ob')
		local Uci_BMIint_ob = exp(`m_BMIint_ob'+1.96*`se_BMIint_ob')
		local p_BMIint_ob = r(p)
		local BMIint_ob= string(exp(`m_BMIint_ob'),"%9.2f")+" ("+string(`Lci_BMIint_ob',"%9.2f")+"-"+string(`Uci_BMIint_ob',"%9.2f")+")"+", p="+string(`p_BMIint_ob',"%9.3f")  
	lincom z_PRS_BMI_PCadj 
		local m_PRSint = r(estimate)
		local se_PRSint = r(se)
		local Lci_PRSint = exp(`m_PRSint'-1.96*`se_PRSint')
		local Uci_PRSint = exp(`m_PRSint'+1.96*`se_PRSint')
		local p_PRSint = r(p)
		local PRSint= string(exp(`m_PRSint'),"%9.2f")+" ("+string(`Lci_PRSint',"%9.2f")+"-"+string(`Uci_PRSint',"%9.2f")+")"+", p="+string(`p_PRSint',"%9.3f")   
	lincom 2.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMIlate_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMIlate_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMIlate_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMIlate_cat
		local m_PRS_c3_ob = r(estimate)
		local se_PRS_c3_ob = r(se)
		local Lci_PRS_c3_ob = exp(`m_PRS_c3_ob'-1.96*`se_PRS_c3_ob')
		local Uci_PRS_c3_ob = exp(`m_PRS_c3_ob'+1.96*`se_PRS_c3_ob')
		local p_PRS_c3_ob = r(p)
		local PRS_c3_ob= string(exp(`m_PRS_c3_ob'),"%9.2f")+" ("+string(`Lci_PRS_c3_ob',"%9.2f")+"-"+string(`Uci_PRS_c3_ob',"%9.2f")+")"+", p="+string(`p_PRS_c3_ob',"%9.3f")   
		
	post coxreg_cvd ("`i'") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'")
	}		 

postclose coxreg_cvd

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20230119_late.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

