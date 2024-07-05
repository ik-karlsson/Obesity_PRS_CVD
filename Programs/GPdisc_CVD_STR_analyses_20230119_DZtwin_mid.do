////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_analyses_20220629_twin_mid
Purpose: Competing risk regression as sensitivity analysis of BMI on CVD, and the effect of PRS for BMI
Study: SfoEpi grant 

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20220629 by Ida Karlsson (IK) (based on main analysis code)
	Updated: 	20221021 by IK - use updated data after removing those <40 and refining the studypop
								 Add DZ only analyses; drop sex-stratified analyses (low power)
				20221102 by IK - run the main model, including the interaction, in DZ only
				20221212 by IK - use updated data including those <60 at end of follow-up (based on rev comments)
				20230119 by IK - Add p-values to tables based on reviewer comments.
	
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

/// CAtegorize BMI
gen BMImid_cat=1 if BMImid>=18.5 & BMImid <25
replace BMImid_cat=2 if BMImid>=25 & BMImid <30
replace BMImid_cat=3 if BMImid!=. & BMImid >=30

tab BMImid_cat, mi

/// Categorize the PRSs
egen PRS_BMI_grpm= mean(PRS_BMI_PCadj), by(study)
egen PRS_BMI_grpsd= sd(PRS_BMI_PCadj), by(study)
gen z_PRS_BMI_PCadj= (PRS_BMI_PCadj-PRS_BMI_grpm)/PRS_BMI_grpsd

// Tertiles
xtile tert_PRS_BMI = z_PRS_BMI_PCadj, nquantiles(3)
	 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							///// SETTING UP THE SURVIVAL DATA, BMI MID /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate age at end of follow-up
gen age_exit=lastage if reg_CVDany==0
replace age_exit=reg_CVDany_age if reg_CVDany==1

// Setting up the survival data with cvd as outcome
stset age_exit, failure(reg_CVDany==1) id(twinnr) enter(BMImid_age)

///// Pair descriptives done in the code for MZ twins
drop if _st!=1
// n=337 deleted

// Keep only complete pairs!
sort pairid 
by pairid: gen pair = _N
drop if pair!=2

stdes

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// SURVIVAL ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Numeric studyn
encode study, gen(studyn)
tab study studyn, mi

cap postclose coxreg_cvd
	postfile coxreg_cvd str26  (Model BMIind_ow BMIind_ob PRSind BMIjnt_ow BMIjnt_ob PRSjnt ///
								BMIint_ow BMIint_ob PRSint BMI_PRSint_ow BMI_PRSint_ob ///
								PRS_c1_ow PRS_c2_ow PRS_c3_ow PRS_c1_ob PRS_c2_ob PRS_c3_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_DZtwin_mid_CVDany.dta, replace
	set more off
	// BMI, categorical
	stcox i.BMImid_cat educ sex smoke if(zygosity!=1), strata(pairid)
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
	stcox z_PRS_BMI_PCadj educ sex smoke if(zygosity!=1), strata(pairid)
	lincom z_PRS_BMI_PCadj 
		local m_PRSind = r(estimate)
		local se_PRSind = r(se)
		local Lci_PRSind = exp(`m_PRSind'-1.96*`se_PRSind')
		local Uci_PRSind = exp(`m_PRSind'+1.96*`se_PRSind')
		local p_PRSind = r(p)
		local PRSind= string(exp(`m_PRSind'),"%9.2f")+" ("+string(`Lci_PRSind',"%9.2f")+"-"+string(`Uci_PRSind',"%9.2f")+")"+", p="+string(`p_PRSind',"%9.3f")

	// BMI + PRS model
	stcox i.BMImid_cat z_PRS_BMI_PCadj educ sex smoke if(zygosity!=1), strata(pairid)
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
	stcox i.BMImid_cat##c.z_PRS_BMI_PCadj educ sex smoke if(zygosity!=1), strata(pairid)
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
	stcox i.tert_PRS_BMI i.tert_PRS_BMI#i.BMImid_cat educ sex smoke if(zygosity!=1), strata(pairid)
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
		
	post coxreg_cvd ("Within DZ") ("`BMIind_ow'") ("`BMIind_ob'") ("`PRSind'") ("`BMIjnt_ow'") ("`BMIjnt_ob'") ("`PRSjnt'") ///
					("`BMIint_ow'") ("`BMIint_ob'") ("`PRSint'") ("`BMI_PRSint_ow'") ("`BMI_PRSint_ob'") ///
					("`PRS_c1_ow'") ("`PRS_c2_ow'")  ("`PRS_c3_ow'") ("`PRS_c1_ob'") ("`PRS_c2_ob'") ("`PRS_c3_ob'") 

postclose coxreg_cvd

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20230119_DZtwin_mid.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

