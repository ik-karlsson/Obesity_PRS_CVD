
////////////////////////////////////////////////////////////////////////////////
////////// Midlife
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 
8 . . .	
7 1.18 1.03 1.36
6 1.16 1.03 1.31
5 1.19 1.05 1.35
4 . . .
3 2.08 1.46 2.96
2 1.61 1.26 2.06
1 1.39 1.15 1.68



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
		   title("")),  /// /* No title */
		   /*title("Total sample", box bexpand)),  /// If we want to, we could use a title for each, but may be prettier to add later.. */ ///
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))

////////////////////////////////////////////////////////////////////////////////
// Men
clear
input order estimate min95 max95 
8 . . .	
7 1.20 1.01 1.43
6 1.11 0.96 1.29
5 1.20 1.02 1.41
4 . . .	
3 1.89 1.15 3.11
2 1.47 1.07 2.04
1 1.20 0.92 1.57
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
		   title("")),  /// /* Title */
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))

////////////////////////////////////////////////////////////////////////////////
// Women
clear
input order estimate min95 max95 
8 . . .	
7 1.13 0.89 1.44
6 1.25 1.02 1.53
5 1.14 0.94 1.38
4 . . .	
3 2.26 1.35 3.78
2 1.79 1.25 2.56
1 1.59 1.23 2.06
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
		   title("")),  /// /* Title */
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))
				

////////////////////////////////////////////////////////////////////////////////
////////// Late-life
////////////////////////////////////////////////////////////////////////////////
// Total sample
clear
input order estimate min95 max95 
8 . . .	
7 1.23 1.03 1.47
6 1.12 0.94 1.33
5 1.24 1.02 1.51
4 . . .
3 1.74 1.22 2.49
2 1.17 0.88 1.55
1 1.35 1.07 1.69



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
		   title("")),  /// /* No title */
		   /*title("Total sample", box bexpand)),  /// If we want to, we could use a title for each, but may be prettier to add later.. */ ///
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))

////////////////////////////////////////////////////////////////////////////////
// Men
clear
input order estimate min95 max95 
8 . . .	
7 1.30 1.02 1.66
6 1.19 0.94 1.50
5 1.33 1.02 1.73
4 . . .	
3 1.54 0.83 2.84
2 1.34 0.90 1.99
1 1.30 0.92 1.82
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
		   title("")),  /// /* Title */
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))


////////////////////////////////////////////////////////////////////////////////
// Women
clear
input order estimate min95 max95 
8 . . .	
7 1.16 0.88 1.51
6 1.00 0.77 1.29
5 1.14 0.86 1.51
4 . . .	
3 1.82 1.16 2.87
2 1.02 0.68 1.51
1 1.38 1.02 1.88
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
		   title("")),  /// /* Title */
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
				text(1 0.6 "High PGS{subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))
								
				
				
				




////////////////////////////////////////////////////////////////////////////////
// Co-twin control, midlife
// NOTE! Skippar nog den här plotten och har tabell - blir inte snyggt när det 
// är så breda CIs..
clear
input order estimate min95 max95 
12 . . .
11 2.09 1.35 3.25
10 1.55 1.15 2.08
9 1.45 1.16 1.80
8
7 1.21 0.43 3.39
6 1.40 0.66 2.95
5 1.23 0.73 2.10
4
3 1.20 0.30 4.82
2 1.19 0.46 3.03
1 1.24 0.53 2.89



end

// Here I just define an empty label, don't really know why, haha.. 
label define jj      12 "                         " ///
					 11 "                         " ///
					 10 "                         " ///
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

// Next comest the actual plot
// You can modify the size, shape, and color of the symols
twoway     (scatter order estimate if order !=. , yaxis(1) msize(vlarge) msymbol(d) mcolor(purple))  ///
           (rcap min95 max95 order, yaxis(1) horizontal lc(purple)  ///
		   title("")),  /// /* Title */
		   xsize(18) ysize(12) /// /* The size of the plot - often needs modifications to get a good ratio */
           xline(1, lc(black) lw(thin) lpattern(line))   /// /* This just gives you a vertical line at HR=1 */
		   xlabel(0.5 1 1.5 3, angle(horiz) format(%2.1fc) labsize(medsmall))    /// /* This specifies where numbers are given on the x-axis */
           xscale(log r(0.3 4.8))           /// /* This specifies the range of the x-axis */
           ylabel(,valuelabel angle(horiz) noticks nogrid)      yscale(noline r(0.5 3.2)) /// /* And the scale of the y-axis */
        legend(order(1 "Hazard ratios (HR)" 2 "95% Confidence interval (CI)") row(2) pos(7) size(small) symxsize(medium))    /// /* Legend, text size etc */
           scheme(s1mono) plotregion(style(none))    ///
        xtitle(" ") ytitle(" ") ///
				text(12 0.35 "{bf:Total sample}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(11 0.3 "Medium PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(10 0.3 "Low PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(9 0.3 "High PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(8 0.35 "{bf:Co-twin control, all twins}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(7 0.3 "Low PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(6 0.3 "Medium PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(5 0.3 "High PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(4 0.35 "{bf:Co-twin control, MZ twins}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(3 0.3 "Low PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(2 0.3 "Medium PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))  ///
				text(1 0.3 "High PGS {subscript:BMI}", yaxis(1) place(w) j(right) size(medsmall))
				

