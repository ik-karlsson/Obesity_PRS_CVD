/***********************************************************************************************
Name: GPdisc_CVD_STR_adata_20220214
Purpose: Combine data for analysis of the association between BMI and CVD, by PRS, in the STR

Study: SfoEpi project: Study how genetic predisposition to higher BMI/WHR influence the association
between BMI/WHR and CVD 

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20220214 by Ida Karlsson (IK)
	Updated: 	20221021 by IK - Add first available measure (regardless of mid/late-life);
								 Include only those 40 or older at BMI measure; drop WHR
				20221212 by IK - Include also those <60 at end of follow-up (based on rev. comments)

	Steps:		1. Libnames 
				2. PRS data - defines the study population
				3. CVD
				4. Covariates and additional variables
				5. BMI data (long format)
				6. Save the data


  	SAS v9.4

************************************************************************************************/
/**************************************** 1. Libnames ******************************************/
/***********************************************************************************************/

libname GPdisc 'P:\Dementia_IK\SfoEpi\GPdisc_CVD\Data\GOSH';
libname bmi 'P:\Dementia_IK\BMI\BMI_data\Data\Derived_data';
libname s2 'P:\Dementia_IK\Archived\PhD project\Study II - CHD_genes_dem\Data\Derived_data';
libname aim1 'P:\Dementia_IK\FORTE_postdoc\Aim_1\Data\Derived_data\Aim_1a\GHOST';

/***********************************************************************************************/
/********************** 2. Read in PRS data to define the study population *********************/
/***********************************************************************************************/

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

proc freq data=prs;
	table study;
run;
/*                                                    Cumulative    Cumulative
                 study       Frequency     Percent     Frequency      Percent
                 ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒE
                 GOSH            2052       10.60          2052        10.60
                 SALT-Y          6403       33.07          8455        43.67
                 TwinGene       10906       56.33         19361       100.00 */

/***********************************************************************************************/
/**************************************** 3. CVD data ******************************************/
/***********************************************************************************************/
PROC IMPORT OUT= WORK.cvd
            DATAFILE= "P:\Dementia_IK\Various data across and outside studies\Register_diagnoses\Data\CVD_REGs_STRlink_20200422.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;

data prs_cvd;
	merge	cvd (in=a keep= twinnr reg_cvd reg_cvd_age reg_stroke reg_stroke_age) 
			prs (in=b)
			stradmin.V_REQUESTED_PERSONS_TWINNR_2018 (in=c keep=twinnr) ;
	by twinnr;
	if b and c;
	if reg_cvd=. then reg_cvd=0;
	if reg_stroke=. then reg_stroke=0;

	reg_cvd_age=int(reg_cvd_Age);
	reg_stroke_age=int(reg_stroke_age);

	reg_CVDany=max(reg_CVD, reg_stroke);
	if reg_cvd=1 or reg_stroke=1 then reg_CVDany_age=min(reg_cvd_age, reg_stroke_age);
	rename 	reg_CVD=reg_CVDns
			reg_CVD_age=reg_CVDns_age;
run;
*n=19311 (n=50 missing register data);
 
proc means data=prs_Cvd;
	var reg_CVDany_age reg_CVDns_age reg_stroke_age;
run;

proc freq data=prs_cvd;
	tables reg_CVDany*(reg_CVDns reg_stroke)/missprint norow nocol nopercent;
run;

/***********************************************************************************************/
/**************************** 4. Covariates and additional variables ***************************/
/***********************************************************************************************/

/* Birthyear, age at death, sex, and zygosity from STR */
data cov1;
	merge 	stradmin.person_info (in=a keep=twinnr pairid tvab sex birthdate vitalstatus deathdate deathyr deathmon deathday)
			stradmin.zygosity (in=b keep=twinnr zygosity)
			prs_cvd (in=c);
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
	if twinnr=[edited out for GitHub/IK] then do;
		deathdate='20060825';
		deathyr='2006';
		deathmon='08';
		deathday='25';
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

	drop bdate ddate deathyr deathmon deathday;

run;
*n=19292 (n=19 not in person info. NOTE! Included as register info in flow-chart, same interpretation);

/* DROPPED 20221212, after reviewer comments on 1st submission (Circulation)
Drop those <60 at lastage or <40 at CVD 
data cov2;	
	set cov1;
	if lastage ge 60;
run;
*n=18084 (n=1208 <60);

/* Education and smoking */
data cov3;
	merge 	cov1 (in=a)
			GPdisc.educ_all_20220214 (in=b keep=twinnr educ)
			aim1.smoking_GOSH_20190603 (in=c keep=twinnr STR_smoke);
	by twinnr;
	if a;
	rename STR_smoke=smoke;
run;
proc freq data=cov3;
	table educ smoke /missprint norow nocol nopercent;
run;
*n=198 misssing education, n=249 missing smoking;

/***********************************************************************************************/
/**************************************** 6. BMI data ******************************************/
/***********************************************************************************************/
/* Select first measure, and the measure taken closest to age 50 and 75 */
/* NOTE! Adding age of CVD onset */
/* BMI */
data bmi_40 bmi_mid bmi_late;
	merge	 bmi.Bmi_data_long_20210510 (in=a) 
			 cov3 (in=b keep=twinnr reg_CVDany_age reg_CVDany_age);
	by twinnr;
	if a and b and bmi ne .;
	age=int(age);
	if reg_CVDany_age ne . and (reg_CVDany_age le age) then prev=1;
	ag50=abs((age+0.1)-50); /* Adding 0.1 to select the youngest, if two of equal length before/after */
	ag75=abs((age-0.1)-75); /* Subtracting 0.1 to select the oldest, if two of equal length before/after */
	if age ge 40 then output bmi_40;
	if age ge 40 and age lt 65 then output bmi_mid;
	else if age ge 65 then output bmi_late;
run;

proc sort data=bmi_40;
	by twinnr age;
proc sort data=bmi_40 nodupkey;
	by twinnr;
run;
*n=18309;

proc sort data=bmi_mid;
	by twinnr prev ag50;
proc sort data=bmi_mid nodupkey;
	by twinnr;
run;
*n=16447;

proc sort data=bmi_late;
	by twinnr ag75;
proc sort data=bmi_late nodupkey;
	by twinnr;
run;
*n=6891;

/* Clean up */
data bmi_402;
	set bmi_40;
	BMI_age=int(age);
	if bmi lt 18.5 then delete;
	if bmi gt 55 then delete;
	keep twinnr study bmi BMI_age;
	rename 	bmi=BMI
			study=BMI_study;
run;
*n=18155;

data bmi_mid2;
	set bmi_mid;
	BMImid_age=int(age);
	if bmi lt 18.5 then delete;
	if bmi gt 55 then delete;
	keep twinnr study bmi BMImid_age;
	rename 	bmi=BMImid
			study=BMImid_study;
run;
*n=16301;

data bmi_late2;
	set bmi_late;
	BMIlate_age=int(age);
	if bmi lt 18.5 then delete;
	if bmi gt 55 then delete;
	keep twinnr study bmi BMIlate_age;
	rename 	bmi=BMIlate
			study=BMIlate_study;
run;
*n=6788;

/* Add to the data */
data bmi_merge;
	merge 	cov3 (in=a)
			BMI_402 (in=b)
			BMI_mid2 (in=c)
			BMI_late2 (in=d);
	by twinnr;
	if a and (b or c or d);
run;
*n=18191 (n=1101 with no measures taken, no measures at age >40, or BMI <18.5 or >55);

proc means data=bmi_merge;
	var BMI_age BMImid_age BMIlate_age;
run;

/*     Variable           N            Mean         Std Dev         Minimum         Maximum
     ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒ
     BMI_age        18155      53.2916001       8.5575972      40.0000000      93.0000000
     BMImid_age     16301      52.2126250       5.5101303      40.0000000      64.0000000
     BMIlate_age     6788      72.0583382       5.1893695      65.0000000      94.0000000
     ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒ
*/

/* Drop those prevalent at first measrue */
data bmi_inc;
	set bmi_merge;
	if reg_CVDany_age= . or (reg_CVDany_age gt BMI_age) or (reg_CVDany_age gt BMImid_age) or (reg_CVDany_age gt BMIlate_age);
run;
*n=18022 (n=169 prevalent); 

/* Drop those missing covariates */
data bmi_nomiss;
	set bmi_inc;
	if educ ne . and smoke ne .;
run;
*n=17988 (n=34 missing covariates);

/***********************************************************************************************/
/************************************* 7. Save the data ****************************************/
/***********************************************************************************************/

data GPdisc.GPdisc_CVD_STR_adata_20221212;
	retain	twinnr
			pairid
			tvab
			sex 
			zygosity
			birthdate
			deathdate
			age_death
			vitalstatus
			lastage
			educ
			smoke

			reg_CVDany
			reg_CVDns
			reg_stroke
			reg_CVDany_age
			reg_CVDns_age
			reg_stroke_age
			
			BMI
			BMI_age
			BMI_study
			BMImid
			BMIlate
			BMImid_age
			BMIlate_age
			BMImid_study
			BMIlate_study
			
			PRS_BMI_PCadj
			PRS_WHR_PCadj
			PRS_WHRadjBMI_PCadj
			;
	set bmi_nomiss;
run;

/* STATA format */
PROC EXPORT DATA= GPdisc.GPdisc_CVD_STR_adata_20221212
	OUTFILE= "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Data\GOSH\GPdisc_CVD_STR_adata_20221212.dta" 
    DBMS=STATA REPLACE;
RUN;

/***********************************************************************************************/
/**************************************** END OF FILE ******************************************/
/***********************************************************************************************/
