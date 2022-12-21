
////////////////////////////////////////////////////////////////////////////////
////////// Midlife
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 
9 . . .	
8 2.04 1.43 2.91
7 1.21 0.95 2.90
6 . . . 
5 1.81 1.41 2.32
4 1.02 0.95 2.58
3 . . .
2 1.54 1.27 1.87
1 1.29 0.95 3.22
end

label define jj      9 "                         " ///
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

twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(gs5))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(gs5)  ///
		   title("")),  /// /* No title */
		   /*title("Total sample", box bexpand)),  /// If we want to, we could use a title for each, but may be prettier to add later.. */ ///
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
				
graph export "C:\Users\idamat\OneDrive - Karolinska Institutet\Elsa\GPdisc_CVD\Documents\Results\Fortestplot_mz_up.tif", replace 
