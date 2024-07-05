////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_forestplots_mz
Purpose: Create forest plots for the publication on BMI on CVD, and the effect of PRS for BMI; MZ twins only
Study: SfoEpi grant 

Created by Ida Karlsson (IK) and Elsa Ojalehto (EO)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm
	
NOTE! Finishing touches for revison of publication; plots named accordingly
STATA v17	
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log using session, replace
set more off

////////////////////////////////////////////////////////////////////////////////
////////// Midlife
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 grp
10 . . . .
9 . . .	0
8 2.04 1.43 2.91 0
7 1.21 0.95 2.90 1
6 . . . 0
5 1.81 1.41 2.32 0
4 1.02 0.95 2.58 1
3 . . . 0
2 1.54 1.27 1.87 0
1 1.29 0.95 3.22 1
end

label define jj      10 "                         " ///
					 9 "                         " ///
					 8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. & grp==0, yaxis(1) msize(vlarge) msymbol(d) mcolor(purple))  ///
           (rcap min95 max95 order if grp==0, yaxis(1) horizontal lc(purple) title(""))  ///
		   (scatter order estimate if order !=. & grp==1, yaxis(1) msize(vlarge) msymbol(d) mcolor(purple*0.6))  ///
           (rcap min95 max95 order if grp==1, yaxis(1) horizontal lc(purple*0.6) title("")),  /// 
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.95 3.0))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(9 0.8 "{bf:Low PGS{subscript:BMI}}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(8 0.8 "Total sample", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(7 0.8 "MZ twin pairs", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(7 0.95 "[0.3]", yaxis(1) place(w) j(right) size(small))  ///
				text(7 3.3 "[4.8]", yaxis(1) place(w) j(right) size(small))  ///
				text(6 0.8 "{bf:Medium PGS{subscript:BMI}}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.8 "Total sample", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.8 "MZ twin pairs", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.95 "[0.4]", yaxis(1) place(w) j(right) size(small))  ///
				text(3 0.8 "{bf:High PGS{subscript:BMI}}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.8 "Total sample", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.8 "MZ twin pairs", yaxis(1) place(w) j(right) size(medsmall)) ///
				text(1 0.95 "[0.5]", yaxis(1) place(w) j(right) size(small)) 
				
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_3_r2.tif", replace 
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_3_r2.eps", fontface("Times New Roman") replace 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_forestplots_mz_20230214.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////