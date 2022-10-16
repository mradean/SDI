
***************************************************************
* This do-file replicates the examples from 'sdii_manual.pdf' *
* ('sdii.ado' must be previously installed)					  *	
***************************************************************

* Example 1: Compute SDIs to compare independent distributions
sdii, sample1(10 2) sample2(5 1)

* Example 2: Compute SDIs to indicate substantive significance
sdii, sample1(10 2) sample2(5 1) mvalue(1) difference

* Example 3: Compute SDIs to compare correlated distributions
sdii, sample1(10 2) sample2(5 1) m(1) corr(.5) diff

* Example 4: Compute SDIs to compare paired samples
sdii, sample1(40 10 2) sample2(40 5 4) level(90)

* Examples 5: Compute SDIs to compare unpaired samples
sdii, sample1(60 10 2) sample2(40 5 4) unpaired

* Examples 6: Compute SDIs to compare unpaired samples with unequal population variances
sdii, sample1(60 10 2) sample2(40 5 4) unpaired unequal welch

* Example 7: Compute SDIs to compare paired variable samples
use https://www.stata-press.com/data/r17/fuel, clear
sdii mpg1 mpg2

* Example 8: Compute percentile SDIs
sdii mpg1 mpg2, percentile median

* Example 9: Compute SDIs to compare unpaired variable samples
sdii mpg1 mpg2, unpaired

* Example 10: Compute SDIs to compare group samples within the same variable
stack mpg1 mpg2, into(mpg) clear
sdii mpg, by(_stack) unpaired

*Using sdi interactively
* Examples 11: Compute SDIs to compare estimated statistics
webuse nhanes2f, clear
logit diabetes age i.race, nolog
nlcom _b[2.race]-_b[3.race]
	
matrix eb = e(b)
matrix eV = e(V)
mat list eb
mat list eV
tempname m1 s1 m2 s2 cor12
scalar `m1' = el(eb,1,3)
scalar `s1' = sqrt(el(eV,3,3))
scalar `m2' = el(eb,1,4)
scalar `s2' = sqrt(el(eV,4,4))
scalar `cor12' = el(eV,4,3)/(sqrt(el(eV,3,3))*sqrt(el(eV,4,4)))
sdii, sample1(`=`m1'' `=`s1'') sample2(`=`m2'' `=`s2'') corr(`=`cor12'') diff

// Figures
matrix rsdi = r(sdii)
gen		pr = rsdi[1,1] in 1
gen 	ll = rsdi[1,7] in 1
gen 	ul = rsdi[1,8] in 1
gen 	lc = rsdi[1,1]- rsdi[1,2]*invnormal(1-.5*(1-.95)) in 1
gen 	uc = rsdi[1,1]+ rsdi[1,2]*invnormal(1-.5*(1-.95)) in 1
replace pr = rsdi[2,1] in 2
replace ll = rsdi[2,7] in 2
replace ul = rsdi[2,8] in 2
replace lc = rsdi[2,1]- rsdi[2,2]*invnormal(1-.5*(1-.95)) in 2
replace uc = rsdi[2,1]+ rsdi[2,2]*invnormal(1-.5*(1-.95)) in 2
replace pr = rsdi[3,1] in 3
replace ll = rsdi[3,7] in 3
replace ul = rsdi[3,8] in 3

gen xval = _n in 1/3
gen yline = 0 in 1/2
gen ylval = .5 in 1
replace ylval = 2.5 in 2

#delimit ;
// The difference with the standard CI;
graph tw scatter pr xval in 3, msymbol(S) msize(medlarge) mcolor(black)
    || rcap ll ul xval in 3, lcolor(black) lwidth(medthick) lpattern(solid)
	|| line yline ylval, lcolor(black) lwidth(medthin) lpattern(solid)
	|| ,
		xlabel(3 `" "Blacks {&minus}" "Other Minorities" "', labsize(3.4))
		ylabel(-.4(.4)1.2, axis(1) labsize(3.4) nogrid angle(0))
		legend(off)
		ytitle("Difference in Coefficients", axis(1) size(3.6))
		xsca(titlegap(3) ra(2.5 5.5))
		ysca(titlegap(2.5))
		yline(0, lcolor(black) lwidth(medthin) lpattern(solid))
		scheme(s1mono) graphregion(fcolor(white)) aspectratio(1)
		note("Note: 95% CIs.", size(3.5));

// The estimates with the standard CI;
graph tw scatter pr xval in 1/2, msymbol(S) msize(medlarge) mcolor(black)
    || rcap lc uc xval in 1/2, lcolor(black) lwidth(medthick) lpattern(solid)
	|| line yline ylval, lcolor(black) lwidth(medthin) lpattern(solid)
	|| ,
		xlabel(1 "Blacks" 2 "Other Minorities", labsize(3.4))
		ylabel(-.6(.4)1, axis(1) labsize(3.4) nogrid angle(0))
		legend(off)
		//xtitle(War Performance, size(3.6))
		ytitle("Estimated Coefficients", axis(1) size(3.6))
		xsca(titlegap(3))
		ysca(titlegap(2.5))
		yline(0, lcolor(black) lwidth(medthin) lpattern(solid))
		scheme(s1mono) graphregion(fcolor(white)) aspectratio(1)
		note("Note: 95% CIs.", size(3.5));
		
// The estimates with SDIs;
graph tw scatter pr xval in 1/2, msymbol(S) msize(medlarge) mcolor(black)
    || rcap ll ul xval in 1, lcolor(black) lwidth(medthick) lpattern(solid)
    || rcap ll ul xval in 2, lcolor(black) lwidth(medthick) lpattern(dash)
	|| line yline ylval, lcolor(black) lwidth(medthin) lpattern(solid)
	|| ,
		xlabel(1 "Blacks" 2 "Other Minorities", labsize(3.4))
		ylabel(-.6(.4)1, axis(1) labsize(3.4) nogrid angle(0))
		legend(off)
		ytitle("Estimated Coefficients", axis(1) size(3.6))
		xsca(titlegap(3))
		ysca(titlegap(2.5))
		yline(0, lcolor(black) lwidth(medthin) lpattern(solid))
		scheme(s1mono) graphregion(fcolor(white)) aspectratio(1)
		note("Note: `=string(rsdi[1,9],"%10.1f")'% SDI.", size(3.5));
#delimit cr

* Examples 12: Compute SDIs to compare 'margins' results
quietly sum age, detail
margins, at(age=(`r(p25)' `r(p75)')) post
nlcom _b[1._at]-_b[2._at]
		
quietly logit diabetes age i.race
quietly sum age, detail
quietly margins, at(age=(`r(p25)' `r(p75)'))
tempname m1 s1 m2 s2 cor12
scalar `m1' = el(r(b),1,1)
scalar `s1' = sqrt(el(r(V),1,1))
scalar `m2' = el(r(b),1,2)
scalar `s2' = sqrt(el(r(V),2,2))
scalar `cor12' = el(r(V),2,1)/(sqrt(el(r(V),1,1))*sqrt(el(r(V),2,2)))
sdii, sample1(`=`m1'' `=`s1'') sample2(`=`m2'' `=`s2'') corr(`=`cor12'') diff

exit

