////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: GPdisc_CVD_STR_forestplots
Purpose: Create forest plots for the publication on BMI on CVD, and the effect of PRS for BMI
Study: SfoEpi grant 

Created by Ida Karlsson (IK) and Elsa Ojalehto (EO)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm
	
NOTE! Finishing touches for revison of publication; plots named accordingly
STATA v17	
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
log using session, replace
set more off

graph set window fontface "Times New Roman"

////////////////////////////////////////////////////////////////////////////////
////////// Midlife
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 
8 . . .	
7 1.24 1.09 1.41
6 1.24 1.08 1.42
5 1.33 1.17 1.52
4 . . .
3 2.08 1.42 3.05
2 1.82 1.42 2.34
1 1.55 1.25 1.92



end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(purple))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(purple)  ///
		   title("{bf:Total sample}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(med))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(med))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(total, replace)

////////////////////////////////////////////////////////////////////////////////
// Men
clear
input order estimate min95 max95 
8 . . .	
7 1.28 1.06 1.53
6 1.18 0.99 1.40
5 1.37 1.10 1.70
4 . . .	
3 1.96 1.15 3.27
2 1.72 1.07 2.40
1 1.44 0.92 1.94
end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(ebblue))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(ebblue)  ///
		   title("{bf:       Men}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(med))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(med))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(men, replace)

////////////////////////////////////////////////////////////////////////////////
// Women
clear
input order estimate min95 max95
8 . . .	
7 1.18 0.94 1.48
6 1.33 1.11 1.60
5 1.26 1.02 1.56
4 . . .	
3 2.18 1.27 3.74
2 1.93 1.37 2.71
1 1.66 1.23 2.22
end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(maroon))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(maroon)  ///
		   title("{bf:     Women}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(med))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(med))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(women, replace)
				
////////////////////////////////////////////////////////////////////////////////
// Combine and export
graph combine total men women, cols(3) title("{bf:a) Midlife}", ///
position(11) justification(left)) ///
scheme(s1mono) plotregion(style(none)) ysize(5cm) xsize(10.7cm)

graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_2a_r2.eps", fontface("Times New Roman") replace
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_2a_r2.tif", replace

////////////////////////////////////////////////////////////////////////////////
////////// Late-life
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 
8 . . .	
7 1.24 1.02 1.51
6 1.13 0.93 1.36
5 1.22 1.01 1.46
4 . . .
3 1.80 1.29 2.51
2 1.15 0.87 1.52
1 1.33 1.05 1.68



end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(purple))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(purple)  ///
		   title("{bf:Total sample}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(med))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(med))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(total, replace)
		
////////////////////////////////////////////////////////////////////////////////
// Men
clear
input order estimate min95 max95 
8 . . .	
7 1.33 1.03 1.71
6 1.18 0.93 1.50
5 1.30 0.97 1.74
4 . . .	
3 1.65 0.79 3.42
2 1.26 0.86 1.84
1 1.31 0.92 1.87
end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(ebblue))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(ebblue)  ///
		   title("{bf:       Men}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(men, replace)

////////////////////////////////////////////////////////////////////////////////
// Women
clear
input order estimate min95 max95 
8 . . .	
7 1.14 0.89 1.46
6 1.03 0.81 1.29
5 1.11 0.85 1.47
4 . . .	
3 1.84 1.22 2.79
2 1.03 0.67 1.56
1 1.34 1.05 1.71
end

label define jj      8 "                         " ///
					 7 "                         " ///
					 6 "                         " ///
                     5 "                         " ///
					 4 "                         " ///
					 3 "                         " ///
                     2 "                         " ///
					 1 "                         " ///
                     0 "                         ", replace
label value order jj 

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(maroon))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(maroon)  ///
		   title("{bf:     Women}", position(11) justification(center))),  /// /* No title */
		   xsize(12) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(1 2 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.7 3.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(off)    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(8 0.6 "{bf:Overweight}", yaxis(1) place(w) j(right) size(med))  ///
				text(7 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.6 "{bf:Obesity}", yaxis(1) place(w) j(right) size(med))  ///
				text(3 0.6 "Low PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.6 "Medium PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall)) ///
		name(women, replace)	
	
////////////////////////////////////////////////////////////////////////////////
// Combine and export
graph combine total men women, cols(3) title("{bf:b) Late-life}", ///
position(11) justification(left)) ///
scheme(s1mono) plotregion(style(none)) ysize(5cm) xsize(10.7cm)

graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_2b_r2.eps", fontface("Times New Roman") replace
graph export "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Output\GOSH\Forestplots_revision\Figure_2b_r2.tif", replace
				
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\SfoEpi\GPdisc_CVD\Log\GOSH\GPdisc_CVD_STR_forestplots_20230214.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
