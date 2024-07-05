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

///// Pair descriptives
drop if _st!=1
// n=337 deleted

// Keep only complete pairs!
sort pairid 
by pairid: gen pair = _N
// CVD by pair
bysort pairid: egen reg_CVDany_pair = total(reg_CVDany)

drop if pair!=2
tab zygosity reg_CVDany_pair
/*         |         reg_CVDany_pair
  zygosity |         0          1          2 |     Total
-----------+---------------------------------+----------
         1 |     2,780        886        374 |     4,040 
         2 |     2,316        952        380 |     3,648 
         4 |     1,904        586        110 |     2,600 
-----------+---------------------------------+----------
     Total |     7,000      2,424        864 |    10,288 

--> 5144 complete pairs; 1212 are discordant for CVD
--> 2020 complete MZ pairs; 443 are discordant for CVD
--> 3124 complete DZ pairs; 769 are discordant for CVD */
	 
stdes
/*                                 |-------------- Per subject --------------|
Category                   Total        Mean         Min     Median        Max
------------------------------------------------------------------------------
Number of subjects         10288   
Number of records          10288           1           1          1          1

Entry time (first)                  52.25855          40         53         64
Exit time (final)                   70.26516          45         70        104

Subjects with gap              0   
Time on gap                    0           .           .          .          .
Time at risk              185252    18.00661           1         16         62

Failures                    2076    .2017885           0          0          1
------------------------------------------------------------------------------
 */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// SURVIVAL ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Numeric studyn
encode study, gen(studyn)
tab study studyn, mi

cap postclose coxreg_cvd
	postfile coxreg_cvd str26  (PRScat BMIall_ow BMItwin_ow BMIdz_ow BMImz_ow BMIall_ob BMItwin_ob BMIdz_ob BMImz_ob) ///
								using P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\GPdisc_CVD_STR_analyses_20230119_twin_mid_CVDany.dta, replace
	set more off
			
	// BMI, categorical, co-twin control
	stcox i.BMImid_cat educ sex smoke, strata(pairid)
		lincom 2.BMImid_cat
		local m_BMItwin_ow = r(estimate)
		local se_BMItwin_ow = r(se)
		local Lci_BMItwin_ow = exp(`m_BMItwin_ow'-1.96*`se_BMItwin_ow')
		local Uci_BMItwin_ow = exp(`m_BMItwin_ow'+1.96*`se_BMItwin_ow')
		local p_BMItwin_ow = r(p)
		local BMItwin_ow= string(exp(`m_BMItwin_ow'),"%9.2f")+" ("+string(`Lci_BMItwin_ow',"%9.2f")+"-"+string(`Uci_BMItwin_ow',"%9.2f")+")"+", p="+string(`p_BMItwin_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMItwin_ob = r(estimate)
		local se_BMItwin_ob = r(se)
		local Lci_BMItwin_ob = exp(`m_BMItwin_ob'-1.96*`se_BMItwin_ob')
		local Uci_BMItwin_ob = exp(`m_BMItwin_ob'+1.96*`se_BMItwin_ob')
		local p_BMItwin_ob = r(p)
		local BMItwin_ob= string(exp(`m_BMItwin_ob'),"%9.2f")+" ("+string(`Lci_BMItwin_ob',"%9.2f")+"-"+string(`Uci_BMItwin_ob',"%9.2f")+")"+", p="+string(`p_BMItwin_ob',"%9.3f")

	// BMI, categorical, co-twin control, DZ only
	stcox i.BMImid_cat educ sex smoke if(zygosity!=1), strata(pairid)
	lincom 2.BMImid_cat
		local m_BMIdz_ow = r(estimate)
		local se_BMIdz_ow = r(se)
		local Lci_BMIdz_ow = exp(`m_BMIdz_ow'-1.96*`se_BMIdz_ow')
		local Uci_BMIdz_ow = exp(`m_BMIdz_ow'+1.96*`se_BMIdz_ow')
		local p_BMIdz_ow = r(p)
		local BMIdz_ow= string(exp(`m_BMIdz_ow'),"%9.2f")+" ("+string(`Lci_BMIdz_ow',"%9.2f")+"-"+string(`Uci_BMIdz_ow',"%9.2f")+")"+", p="+string(`p_BMIdz_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMIdz_ob = r(estimate)
		local se_BMIdz_ob = r(se)
		local Lci_BMIdz_ob = exp(`m_BMIdz_ob'-1.96*`se_BMIdz_ob')
		local Uci_BMIdz_ob = exp(`m_BMIdz_ob'+1.96*`se_BMIdz_ob')
		local p_BMIdz_ob = r(p)
		local BMIdz_ob= string(exp(`m_BMIdz_ob'),"%9.2f")+" ("+string(`Lci_BMIdz_ob',"%9.2f")+"-"+string(`Uci_BMIdz_ob',"%9.2f")+")"+", p="+string(`p_BMIdz_ob',"%9.3f")

	// BMI, categorical, co-twin control, MZ only
	stcox i.BMImid_cat educ sex smoke if(zygosity==1), strata(pairid)
	lincom 2.BMImid_cat
		local m_BMImz_ow = r(estimate)
		local se_BMImz_ow = r(se)
		local Lci_BMImz_ow = exp(`m_BMImz_ow'-1.96*`se_BMImz_ow')
		local Uci_BMImz_ow = exp(`m_BMImz_ow'+1.96*`se_BMImz_ow')
		local p_BMImz_ow = r(p)
		local BMImz_ow= string(exp(`m_BMImz_ow'),"%9.2f")+" ("+string(`Lci_BMImz_ow',"%9.2f")+"-"+string(`Uci_BMImz_ow',"%9.2f")+")"+", p="+string(`p_BMImz_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMImz_ob = r(estimate)
		local se_BMImz_ob = r(se)
		local Lci_BMImz_ob = exp(`m_BMImz_ob'-1.96*`se_BMImz_ob')
		local Uci_BMImz_ob = exp(`m_BMImz_ob'+1.96*`se_BMImz_ob')
		local p_BMImz_ob = r(p)
		local BMImz_ob= string(exp(`m_BMImz_ob'),"%9.2f")+" ("+string(`Lci_BMImz_ob',"%9.2f")+"-"+string(`Uci_BMImz_ob',"%9.2f")+")"+", p="+string(`p_BMImz_ob',"%9.3f")

		post coxreg_cvd ("All") ("`BMIall_ow'") ("`BMItwin_ow'") ("`BMIdz_ow'") ("`BMImz_ow'") ("`BMIall_ob'") ("`BMItwin_ob'") ("`BMIdz_ob'") ("`BMImz_ob'") 
	
	/// By tertiles of PRSbmi
	foreach j of num 1/3 {
	set more off
	// BMI, categorical, co-twin control
	stcox i.BMImid_cat educ sex smoke if(tert_PRS_BMI==`j'), strata(pairid)
		lincom 2.BMImid_cat
		local m_BMItwin_ow = r(estimate)
		local se_BMItwin_ow = r(se)
		local Lci_BMItwin_ow = exp(`m_BMItwin_ow'-1.96*`se_BMItwin_ow')
		local Uci_BMItwin_ow = exp(`m_BMItwin_ow'+1.96*`se_BMItwin_ow')
		local p_BMItwin_ow = r(p)
		local BMItwin_ow= string(exp(`m_BMItwin_ow'),"%9.2f")+" ("+string(`Lci_BMItwin_ow',"%9.2f")+"-"+string(`Uci_BMItwin_ow',"%9.2f")+")"+", p="+string(`p_BMItwin_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMItwin_ob = r(estimate)
		local se_BMItwin_ob = r(se)
		local Lci_BMItwin_ob = exp(`m_BMItwin_ob'-1.96*`se_BMItwin_ob')
		local Uci_BMItwin_ob = exp(`m_BMItwin_ob'+1.96*`se_BMItwin_ob')
		local p_BMItwin_ob = r(p)
		local BMItwin_ob= string(exp(`m_BMItwin_ob'),"%9.2f")+" ("+string(`Lci_BMItwin_ob',"%9.2f")+"-"+string(`Uci_BMItwin_ob',"%9.2f")+")"+", p="+string(`p_BMItwin_ob',"%9.3f")

	// BMI, categorical, co-twin control, DZ only
	stcox i.BMImid_cat educ sex smoke if(tert_PRS_BMI==`j' & zygosity!=1), strata(pairid)
	lincom 2.BMImid_cat
		local m_BMIdz_ow = r(estimate)
		local se_BMIdz_ow = r(se)
		local Lci_BMIdz_ow = exp(`m_BMIdz_ow'-1.96*`se_BMIdz_ow')
		local Uci_BMIdz_ow = exp(`m_BMIdz_ow'+1.96*`se_BMIdz_ow')
		local p_BMIdz_ow = r(p)
		local BMIdz_ow= string(exp(`m_BMIdz_ow'),"%9.2f")+" ("+string(`Lci_BMIdz_ow',"%9.2f")+"-"+string(`Uci_BMIdz_ow',"%9.2f")+")"+", p="+string(`p_BMIdz_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMIdz_ob = r(estimate)
		local se_BMIdz_ob = r(se)
		local Lci_BMIdz_ob = exp(`m_BMIdz_ob'-1.96*`se_BMIdz_ob')
		local Uci_BMIdz_ob = exp(`m_BMIdz_ob'+1.96*`se_BMIdz_ob')
		local p_BMIdz_ob = r(p)
		local BMIdz_ob= string(exp(`m_BMIdz_ob'),"%9.2f")+" ("+string(`Lci_BMIdz_ob',"%9.2f")+"-"+string(`Uci_BMIdz_ob',"%9.2f")+")"+", p="+string(`p_BMIdz_ob',"%9.3f")

	// BMI, categorical, co-twin control, MZ only
	stcox i.BMImid_cat educ sex smoke if(tert_PRS_BMI==`j' & zygosity==1), strata(pairid)
	lincom 2.BMImid_cat
		local m_BMImz_ow = r(estimate)
		local se_BMImz_ow = r(se)
		local Lci_BMImz_ow = exp(`m_BMImz_ow'-1.96*`se_BMImz_ow')
		local Uci_BMImz_ow = exp(`m_BMImz_ow'+1.96*`se_BMImz_ow')
		local p_BMImz_ow = r(p)
		local BMImz_ow= string(exp(`m_BMImz_ow'),"%9.2f")+" ("+string(`Lci_BMImz_ow',"%9.2f")+"-"+string(`Uci_BMImz_ow',"%9.2f")+")"+", p="+string(`p_BMImz_ow',"%9.3f")
	lincom 3.BMImid_cat
		local m_BMImz_ob = r(estimate)
		local se_BMImz_ob = r(se)
		local Lci_BMImz_ob = exp(`m_BMImz_ob'-1.96*`se_BMImz_ob')
		local Uci_BMImz_ob = exp(`m_BMImz_ob'+1.96*`se_BMImz_ob')
		local p_BMImz_ob = r(p)
		local BMImz_ob= string(exp(`m_BMImz_ob'),"%9.2f")+" ("+string(`Lci_BMImz_ob',"%9.2f")+"-"+string(`Uci_BMImz_ob',"%9.2f")+")"+", p="+string(`p_BMImz_ob',"%9.3f")

		post coxreg_cvd ("`j'") ("`BMIall_ow'") ("`BMItwin_ow'") ("`BMIdz_ow'") ("`BMImz_ow'") ("`BMIall_ob'") ("`BMItwin_ob'") ("`BMIdz_ob'") ("`BMImz_ob'") 
	}
postclose coxreg_cvd

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_analyses_20230119_twin_mid.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

