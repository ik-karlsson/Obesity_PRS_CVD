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

// Adjusted for smoking and education
regress BMImid z_PRS_BMI_PCadj BMImid_age sex smoke educ,vce(cluster pairid)
/*----------------------------------------------------------------------------
             |               Robust
      BMImid | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
z_PRS_BMI_P~j |   1.123334   .0282912    39.71   0.000     1.067878     1.17879
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
	using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_mid_descriptive.dta, replace

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
sts graph, by(BMImid_cat) name(km_bmi, replace)
stphplot, by(BMImid_cat) name(surv_bmi, replace)
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
graph combine km_bmi surv_bmi km_prs surv_prs km_sex surv_sex km_smoke surv_smoke km_educ surv_educ, iscale(.3)
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_mid_PHplots.pdf", replace
// Check with ph-test
stcox i.BMImid_cat z_PRS_BMI_PCadj i.sex educ smoke, vce(cluster pairid)
estat phtest, detail
/*------------------------------------------------------
             |        rho     chi2       df    Prob>chi2
-------------+------------------------------------------
1b.BMImid_~t |          .        .        1           .
2.BMImid_cat |   -0.04920     7.72        1       0.0055
3.BMImid_cat |   -0.05050     9.20        1       0.0024
z_PRS_BMI_~j |   -0.00251     0.02        1       0.8836
      1b.sex |          .        .        1           .
       2.sex |    0.06688    14.64        1       0.0001
        educ |   -0.01865     1.13        1       0.2883
       smoke |   -0.07109    16.30        1       0.0001
-------------+------------------------------------------
 Global test |               59.86        6       0.0000
-------------------------------------------------------- */
// Clear issue with BMI, sex, and smoking. Will stratify on sex anyway, but use bootstrapping to obtain less biased SEs
// Questionable, but will keep BMI category as is (ref Hernan)

// Test for colinearity in joint model of BMI and the PRS, using strata for sex and smoke
// With continous BMI
stcox BMImid z_PRS_BMI_PCadj sex smoke educ, vce(cluster pairid) strata(study)
estat vce, corr
// Correlation -0.35

// Look at the variance inflation factor from logistic regression (does not translate to Cox)
// Countinous BMI
regress reg_CVDany BMImid z_PRS_BMI_PCadj educ sex smoke, vce(cluster pairid)
estat vif
/*  Variable |       VIF       1/VIF  
-------------+----------------------
      BMImid |      1.16    0.861729
z_PRS_BMI_~j |      1.14    0.877795
         sex |      1.04    0.962067
       smoke |      1.01    0.986039
        educ |      1.01    0.993486
-------------+----------------------
    Mean VIF |      1.07			*/

// Without the PRS
regress reg_CVDany BMImid educ sex smoke, vce(cluster pairid)
estat vif
/*  Variable |       VIF       1/VIF  
-------------+----------------------
         sex |      1.03    0.968900
      BMImid |      1.03    0.973129
       smoke |      1.01    0.991999
        educ |      1.01    0.994990
-------------+----------------------
    Mean VIF |      1.02			*/
	
// Only linear BMI and the PRS
regress reg_CVDany BMImid z_PRS_BMI_PCadj, vce(cluster pairid)
estat vif
// VIF 1.12

cap postclose coxreg_cvd
	postfile coxreg_cvd str26 (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_mid_survCVDany.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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
	postfile coxreg_cvd str26  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_mid_survCVDns.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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
	postfile coxreg_cvd str26  (sex BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_mid_survSTROKE.dta, replace
	set more off
	// BMI, categorical
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p) 
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj i.sex smoke educ, vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat i.sex smoke educ, vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_c3_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIind_ow = r(estimate)
		local se_BMIind_ow = r(se)
		local Lci_BMIind_ow = exp(`m_BMIind_ow'-1.96*`se_BMIind_ow')
		local Uci_BMIind_ow = exp(`m_BMIind_ow'+1.96*`se_BMIind_ow')
		local p_BMIind_ow = r(p)
		local BMIind_ow= string(exp(`m_BMIind_ow'),"%9.2f")+" ("+string(`Lci_BMIind_ow',"%9.2f")+"-"+string(`Uci_BMIind_ow',"%9.2f")+")"+", p="+string(`p_BMIind_ow',"%9.3f")
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIjnt_ow = r(estimate)
		local se_BMIjnt_ow = r(se)
		local Lci_BMIjnt_ow = exp(`m_BMIjnt_ow'-1.96*`se_BMIjnt_ow')
		local Uci_BMIjnt_ow = exp(`m_BMIjnt_ow'+1.96*`se_BMIjnt_ow')
		local p_BMIjnt_ow = r(p)
		local BMIjnt_ow= string(exp(`m_BMIjnt_ow'),"%9.2f")+" ("+string(`Lci_BMIjnt_ow',"%9.2f")+"-"+string(`Uci_BMIjnt_ow',"%9.2f")+")"+", p="+string(`p_BMIjnt_ow',"%9.3f") 
	lincom 3.BMImid_cat
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
	bootstrap, reps(100) seed(7982): stcox i.BMImid_cat##c.z_PRS_BMI_PCadj smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	lincom 2.BMImid_cat
		local m_BMIint_ow = r(estimate)
		local se_BMIint_ow = r(se)
		local Lci_BMIint_ow = exp(`m_BMIint_ow'-1.96*`se_BMIint_ow')
		local Uci_BMIint_ow = exp(`m_BMIint_ow'+1.96*`se_BMIint_ow')
		local p_BMIint_ow = r(p)
		local BMIint_ow= string(exp(`m_BMIint_ow'),"%9.2f")+" ("+string(`Lci_BMIint_ow',"%9.2f")+"-"+string(`Uci_BMIint_ow',"%9.2f")+")"+", p="+string(`p_BMIint_ow',"%9.3f")  
	lincom 3.BMImid_cat
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
	lincom 2.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ow = r(estimate)
		local se_BMI_PRSint_ow = r(se)
		local Lci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'-1.96*`se_BMI_PRSint_ow')
		local Uci_BMI_PRSint_ow = exp(`m_BMI_PRSint_ow'+1.96*`se_BMI_PRSint_ow')
		local p_BMI_PRSint_ow = r(p)
		local BMI_PRSint_ow= string(exp(`m_BMI_PRSint_ow'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ow',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ow',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ow',"%9.3f")   
	lincom 3.BMImid_cat#z_PRS_BMI_PCadj
		local m_BMI_PRSint_ob = r(estimate)
		local se_BMI_PRSint_ob = r(se)
		local Lci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'-1.96*`se_BMI_PRSint_ob')
		local Uci_BMI_PRSint_ob = exp(`m_BMI_PRSint_ob'+1.96*`se_BMI_PRSint_ob')
		local p_BMI_PRSint_ob = r(p)
		local BMI_PRSint_ob= string(exp(`m_BMI_PRSint_ob'),"%9.2f")+" ("+string(`Lci_BMI_PRSint_ob',"%9.2f")+"-"+string(`Uci_BMI_PRSint_ob',"%9.2f")+")"+", p="+string(`p_BMI_PRSint_ob',"%9.3f")   
	
	/// By tertiles of PRSbmi
	bootstrap, reps(100) seed(7982): stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat smoke educ if(sex==`i'), vce(cluster pairid) strata(study)
	/// 1st quartile
	lincom 1.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c1_ow = r(estimate)
		local se_PRS_c1_ow = r(se)
		local Lci_PRS_c1_ow = exp(`m_PRS_c1_ow'-1.96*`se_PRS_c1_ow')
		local Uci_PRS_c1_ow = exp(`m_PRS_c1_ow'+1.96*`se_PRS_c1_ow')
		local p_PRS_c1_ow = r(p)
		local PRS_c1_ow= string(exp(`m_PRS_c1_ow'),"%9.2f")+" ("+string(`Lci_PRS_c1_ow',"%9.2f")+"-"+string(`Uci_PRS_c1_ow',"%9.2f")+")"+", p="+string(`p_PRS_c1_ow',"%9.3f")  
	lincom 1.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c1_ob = r(estimate)
		local se_PRS_c1_ob = r(se)
		local Lci_PRS_c1_ob = exp(`m_PRS_c1_ob'-1.96*`se_PRS_c1_ob')
		local Uci_PRS_c1_ob = exp(`m_PRS_c1_ob'+1.96*`se_PRS_c1_ob')
		local p_PRS_c1_ob = r(p)
		local PRS_c1_ob= string(exp(`m_PRS_c1_ob'),"%9.2f")+" ("+string(`Lci_PRS_c1_ob',"%9.2f")+"-"+string(`Uci_PRS_c1_ob',"%9.2f")+")"+", p="+string(`p_PRS_c1_ob',"%9.3f")    
	/// 2nd quartile
	lincom 2.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c2_ow = r(estimate)
		local se_PRS_c2_ow = r(se)
		local Lci_PRS_c2_ow = exp(`m_PRS_c2_ow'-1.96*`se_PRS_c2_ow')
		local Uci_PRS_c2_ow = exp(`m_PRS_c2_ow'+1.96*`se_PRS_c2_ow')
		local p_PRS_c2_ow = r(p)
		local PRS_c2_ow= string(exp(`m_PRS_c2_ow'),"%9.2f")+" ("+string(`Lci_PRS_c2_ow',"%9.2f")+"-"+string(`Uci_PRS_c2_ow',"%9.2f")+")"+", p="+string(`p_PRS_c2_ow',"%9.3f")    
	lincom 2.tert_PRS_BMI#3.BMImid_cat
		local m_PRS_c2_ob = r(estimate)
		local se_PRS_c2_ob = r(se)
		local Lci_PRS_c2_ob = exp(`m_PRS_c2_ob'-1.96*`se_PRS_c2_ob')
		local Uci_PRS_c2_ob = exp(`m_PRS_c2_ob'+1.96*`se_PRS_c2_ob')
		local p_PRS_c2_ob = r(p)
		local PRS_c2_ob= string(exp(`m_PRS_c2_ob'),"%9.2f")+" ("+string(`Lci_PRS_c2_ob',"%9.2f")+"-"+string(`Uci_PRS_c2_ob',"%9.2f")+")"+", p="+string(`p_PRS_c2_ob',"%9.3f")    
	/// 3rd quartile
	lincom 3.tert_PRS_BMI#2.BMImid_cat
		local m_PRS_c3_ow = r(estimate)
		local se_PRS_c3_ow = r(se)
		local Lci_PRS_c3_ow = exp(`m_PRS_c3_ow'-1.96*`se_PRS_c3_ow')
		local Uci_PRS_c3_ow = exp(`m_PRS_c3_ow'+1.96*`se_PRS_c3_ow')
		local p_PRS_cc_ow = r(p)
		local PRS_c3_ow= string(exp(`m_PRS_c3_ow'),"%9.2f")+" ("+string(`Lci_PRS_c3_ow',"%9.2f")+"-"+string(`Uci_PRS_c3_ow',"%9.2f")+")"+", p="+string(`p_PRS_c3_ow',"%9.3f")    
	lincom 3.tert_PRS_BMI#3.BMImid_cat
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

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20230119_mid.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////




