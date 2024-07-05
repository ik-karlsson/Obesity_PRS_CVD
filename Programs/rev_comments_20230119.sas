/* Look at Ns etc, reviewer comments */
/* Ida Karlsson, 20230119 */

/* 1. How many are in SATSA or GENDER, but not SALT? */

data gs;
	merge	salt.v_salt_admin (in=a) 
			gender.INCIDENT_COMPLETE_080812 (in=b);
	by twinnr;
	if b and not a;
run;
*n=63 in GENDER but not SALT;

data satsa;
	set satsa.RESP0517;
	if ipt1=1 or ipt2=1 or ipt3=1 or ipt5=1 or ipt6=1;
	keep twinnr;
proc sort;
	by twinnr;
run;

data ss;
	merge	salt.v_salt_admin (in=a) 
			satsa (in=b);
	by twinnr;
	if b and not a;
run;
*n=351 in SATSA but not SALT;

/* Base sample, SALT, SATSA, or GENDER */
data base;
	merge	salt.v_salt_admin (in=a) 
			gender.INCIDENT_COMPLETE_080812 (in=b)
			satsa (in=c);
	by twinnr;
	if a then salt=1;
	if b then gender=1;
	if c then satsa=1;
	keep twinnr salt gender satsa;
run;
*n=45337;

/* PRS in respective sample */
PROC IMPORT OUT= WORK.PRS
            DATAFILE= "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Data\GOSH\PRS_BMI_WHR_adjPC_20220213.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;
*n=19361;

proc sort data=prs;
	by twinnr;
run;

data base_prs;
	merge	base (in=a) 
			prs (in=b keep=twinnr study);
	by twinnr;
	if study ne ' ' then prs=1;
		else prs=0;
run;

proc freq data=base_prs;
	tables study study*(salt gender satsa)/missprint norow nocol nopercent;
run;
/* 	SATSA: 622
	GENDER: 377
	HARMONY: 2052-622-377 = 1017
	TwinGene: 10906
	SALTY: 6403 */

/***********************************************************************************************/
/* Sample characteristics */
/***********************************************************************************************/
/* NOTE! Using same codes as in main analysisdata code */

libname GPdisc 'P:\Dementia_IK\SfoEpi\GPdisc_CVD\Data\GOSH';
libname aim1 'P:\Dementia_IK\FORTE_postdoc\Aim_1\Data\Derived_data\Aim_1a\GHOST';

/* Birthyear, age at death, sex, from STR */
data cov1;
	merge 	stradmin.person_info (in=a keep=twinnr pairid tvab sex birthyr birthdate vitalstatus deathdate deathyr deathmon deathday)
			base_prs (in=c);
	by twinnr;
	if a and c;

	/***** Calculate age at death *****/
	/* NOTE! Recode all who died >2016 as alive */
	if deathyr gt 2016 then do;
		deathdate=' ';
		deathyr=' ';
		deathmon=' ';
		deathday=' ';
		vitalstatus=0;
	end;
	/* Convert birthdate to date format */
	/* NOTE! Gives error message for three twins with birthdate=1000000 */
	/* They are kept as we have GOSH info (this is dealt with in the dementia file), but birthdate left as is */
	bdate=input(put(birthdate,8.),yymmdd8.);
	/* 1 twin is dead, but without a date. Assigning date of death from CDR. */
	/* NOTE! Twinnr and birthdate removed prior to posting to GitHub/IK */
    if twinnr=XXXXXXX then do;
        deathdate='XXXXXXXX';
        deathyr='XXXX';
        deathmon='XX';
        deathday='XX';
	end;
	/* Deal with faulty deathdates, imputing to last day of month/year if missing */
	if deathdate ne ' ' then do;

		if deathmon in(' ', '00') then deathmon='12';
		if deathday in (' ', '00') then do;
			if deathmon in ('01', '03', '05', '07', '08', '10', '12') then deathday='31';
			else if deathmon in ('04', '06', '09', '11') then deathday='30';
			else if deathmon='02' then deathday='28';
		end;
		/* Create deathdate in date format */
		ddate=MDY(deathmon, deathday, deathyr);
		/* Calculate age at death */
		age_death=int(YRDIF(bdate, ddate, 'actual'));
	end;
	/* Create last age as age at death or end of 2016 (register linkage) */
	if age_death ne . then lastage=age_death;
	else lastage=int(YRDIF(bdate, MDY(12,31,2016), 'actual'));

	birthyear=1*birthyr;
	drop bdate ddate deathyr deathmon deathday;
	
run;

/* Education and smoking */
data cov3;
	merge 	cov1 (in=a)
			GPdisc.educ_all_20220214 (in=b keep=twinnr educ)
			aim1.smoking_GOSH_20190603 (in=c keep=twinnr STR_smoke);
	by twinnr;
	if a;
	rename STR_smoke=smoke;
run;

/* Age at baseline from SALT, SATSA, or GENDER */
data s_age;
	set satsa.RESP0517;
	if ipt1=1 or ipt2=1 or ipt3=1 or ipt5=1 or ipt6=1;
	s_age=min(iage1, iage2, iage3, iage5, iage6);
	keep twinnr s_age;
proc sort;
	by twinnr;
run;
data check;
	set s_age;
	if s_age=.;
run;

data cov4;
	merge	cov3 (in=a) 
			salt.v_salt_admin (in=b keep=twinnr interview_age_value)
			gender.INCIDENT_COMPLETE_080812 (in=c keep=twinnr g_age1)
			s_age (in=d);
	by twinnr;
	if a;
	age_b=min(interview_age_value, g_age1, s_age);
run;

/*Female sex, N (%)
Low education, N (%)
Smokers, N (%)
Age at baseline, M (SD)
Age at last follow up, M (SD)
Age at death, M (SD)
BMI at baseline, M (SD)
*/

/* Frequencies for all, and those with/without genotype data */
proc freq data=cov4;
	tables smoke educ/missprint norow;
run;

proc freq data=cov4;
	tables smoke educ prs*(smoke educ)/missprint norow nopercent;
run;

proc sort data=cov4;
	by prs;
run;

proc means data=cov4;
	var birthyear age_b age_death;
run;

proc means data=cov4;
	var birthyear age_b age_death;
	class prs;
run;
