
clear
version 14

*****************************
* create the simulated data *
*****************************

set seed 32307
set obs 1000
set more off

local b1 = -.45

gen x = exp(runiform())
gen e  = rlogistic()+2.5
gen y = cond(invlogit(`b1'*x + e) > .5, 1, 0)
keep y x

*****************
* run the model *
*****************
logit y x

exit


*********
* Fig 1 *
*********

// min vs max
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
margins, at(x=(`r1' `r2')) post noatlegend
nlcom _b[2._at]-_b[1._at]
qui logit y x
sdi, xofi(x=(`r1' `r2')) range diff

// 4 out of 51 estimates
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u4 = (`r2'-`r1')/3 // unit increase
margins, at(x=(`r1'(`u4')`r2')) noatlegend
sdi, xofi(x=(`r1'(`u4')`r2')) range

// Full set of resuts
* Fig 1a*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
margins, at(x=(`r1'(`u')`r2')) noatlegend
mat res = r(table)
mat resq = res[1,1...]'
mat resl = res[5,1...]'
mat resu = res[6,1...]'
svmat resq , name(q)
svmat resl , name(l1_ci)
svmat resu , name(h1_ci)
rename (l1_ci1 h1_ci1) (l1_ci h1_ci)

gen xaxis = _n in 1/51
sum xaxis, meanonly
#del ;
tw line q1 xaxis, lwidth(medium)
|| line l1_ci xaxis, lpattern(dash) lwidth(medium)
|| line h1_ci xaxis, lpattern(dash) lwidth(medium)
|| ,
	yline(`=l1_ci[1]', lpattern(shortdash) lwidth(thin) lcolor(gs10))
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.01) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(.7(.05).9, labsize(3.4) nogrid angle(0))
	ytitle(Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5))
	legend(off)
	note("Note: 95% CIs.", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
#del cr
exit

* Fig 1b*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
sdi, xofi(x=(`r1'(`u')`r2')) range nolegend nonotes

mat res = r(sdi)
mat resq = res[1...,1]
mat resl = res[1..2,6]
mat resu = res[1..2,7]
local sdi = res[1,8]
svmat resq , name(q)
svmat resl , name(lp_sdi)
svmat resu , name(hp_sdi)
rename (lp_sdi1 hp_sdi1) (lp_sdi hp_sdi)

gen xaxis = 1 in 1
replace xaxis = 51 in 2
replace xaxis = _n-1 in 3/51
sort xaxis
sum xaxis, meanonly
#del ;
tw line q1 xaxis, lwidth(medium) lcolor(black)
|| scatter q1 xaxis if lp_sdi!=., msymbol(S i i) msize(small) mcolor(black) clwidth(medium) clcolor(black)
|| rcap lp_sdi hp_sdi xaxis, lcolor(black) lwidth(medthick) clwidth(thin) clcolor(black)
|| ,
	yline(`=lp_sdi[1]', lpattern(shortdash) lwidth(thin) lcolor(gs10))
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.01) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(.7(.05).9, labsize(3.4) nogrid angle(0))
	ytitle(Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5))
	legend(off)
	note("Note: `=string(round(`sdi',.1))'% SDIs.", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
#del cr
exit


*****************
* Fig 2a and 2b *
*****************

// 4 out of 51 estimates
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u4 = (`r2'-`r1')/3 // unit increase
margins, at(x=(`r1'(`u4')`r2')) at(x=(`=`r1'+1'(`u4')`=`r2'+1')) noatlegend
sdi, xofi(x=(`r1'(`u4')`r2')) nunit(1)

// Full set of resuts
* Fig 2a*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
margins, at(x=(`r1'(`u')`r2')) at(x=(`=`r1'+1'(`u')`=`r2'+1')) noatlegend
mat res = r(table)
mat resq1 = res[1,1..51]'
mat resl1 = res[5,1..51]'
mat resu1 = res[6,1..51]'
mat resq2 = res[1,52..102]'
mat resl2 = res[5,52..102]'
mat resu2 = res[6,52..102]'
svmat resq1 , name(q1)
svmat resl1 , name(l1_ci)
svmat resu1 , name(h1_ci)
svmat resq2 , name(q2)
svmat resl2 , name(l2_ci)
svmat resu2 , name(h2_ci)
rename (q11 l1_ci1 h1_ci1 q21 l2_ci1 h2_ci1) (q1 l1_ci h1_ci q2 l2_ci h2_ci)

gen xaxis = _n in 1/51
sum xaxis, meanonly
#del ;
tw line q1 xaxis, lpattern(solid) lwidth(medium) lcolor(red)
|| line q2 xaxis, lpattern(solid) lwidth(medium) lcolor(blue)
|| line l1_ci xaxis, lpattern(dash) lwidth(medium) lcolor(red)
|| line h1_ci xaxis, lpattern(dash) lwidth(medium) lcolor(red)
|| line l2_ci xaxis, lpattern(shortdash) lwidth(medium) lcolor(blue)
|| line h2_ci xaxis, lpattern(shortdash) lwidth(medium) lcolor(blue)
|| ,
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.01) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(.5(.1).9, labsize(3.4) nogrid angle(0))
	ytitle(Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5))
	legend(off)
	note("Note: 95% CIs.", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
#del cr
exit

* Fig 2b*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
sdi, xofi(x=(`r1'(`u')`r2')) nunit(1) many nolegend nonotes

mata st_matrix("res",st_matrix("r(sdi)")[range(1,101,2),]\st_matrix("r(sdi)")[range(2,102,2),])
mat resq1 = res[1..51,1]
mat resl1 = res[1..51,6]
mat resu1 = res[1..51,7]
mat resq2 = res[52..102,1]
mat resl2 = res[52..102,6]
mat resu2 = res[52..102,7]
mat ress = res[1...,8]
svmat resq1 , name(q1)
svmat resl1 , name(l1_sdi)
svmat resu1 , name(h1_sdi)
svmat resq2 , name(q2)
svmat resl2 , name(l2_sdi)
svmat resu2 , name(h2_sdi)
svmat ress , name(rsdi)
rename (q11 l1_sdi1 h1_sdi1 q21 l2_sdi1 h2_sdi1) (q1 l1_sdi h1_sdi q2 l2_sdi h2_sdi)

sum rsdi1, meanonly
local mi = r(min)
local ma = r(max)
gen xaxis = _n in 1/51
sum xaxis, meanonly
#del ;
tw line q1 xaxis, lpattern(solid) lwidth(medium) lcolor(red)
|| line q2 xaxis, lpattern(solid) lwidth(medium) lcolor(blue)
|| line l1_sdi xaxis, lpattern(dash) lwidth(medium) lcolor(red)
|| line h1_sdi xaxis, lpattern(dash) lwidth(medium) lcolor(red)
|| line l2_sdi xaxis, lpattern(shortdash) lwidth(medium) lcolor(blue)
|| line h2_sdi xaxis, lpattern(shortdash) lwidth(medium) lcolor(blue)
|| ,
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.01) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(.6(.07).88, labsize(3.4) nogrid angle(0))
	ytitle(Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5))
	legend(off)
	note("Note: SDI Range [`=string(round(`mi',.1))'%, `=string(round(`ma',.1))'%].", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
#del cr
exit


*****************
* Fig 2c and 2d *
*****************

// 4 out of 51 estimates
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u4 = (`r2'-`r1')/3 // unit increase
margins, at(x=(`r1'(`u4')`r2')) at(x=(`=`r1'+1'(`u4')`=`r2'+1')) post noatlegend
forval i = 1/4 {
	nlcom _b[`=4+`i''._at]-_b[`i'._at]
}
qui logit y x
sdi, xofi(x=(`r1'(`u4')`r2')) nunit(1) firstdiff range

// Full set of resuts
* Fig 2c*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
margins, at(x=(`r1'(`u')`r2')) at(x=(`=`r1'+1'(`u')`=`r2'+1')) post noatlegend
mat resq = J(51,1,.) 
mat resv = J(51,1,.) 
forval i = 1/51 {
	nlcom _b[`=51+`i''._at]-_b[`i'._at]
	mat resq[`i',1] = r(b)
	mat resv[`i',1] = r(V)
}
mata st_matrix("resv",sqrt(st_matrix("resv"))*invnormal(.95+(1-.95)/2))
qui logit y x
svmat resq , name(qd)
svmat double resv , name(pm) // plus-minus
rename (qd1 pm1) (qd pm)
gen ld_ci = qd-pm
gen hd_ci = qd+pm

gen xaxis = _n in 1/51
sum xaxis, meanonly
#del ;
tw line qd xaxis, lpattern(solid) lwidth(medium) lcolor(black)
|| line ld_ci xaxis, lpattern(dash) lwidth(medium) lcolor(black)
|| line hd_ci xaxis, lpattern(dash) lwidth(medium) lcolor(black)
|| ,
	yline(0, lpattern(shortdash) lwidth(medthin) lcolor(black))
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.01) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(-.18(.05).02 0, labsize(3.4) nogrid angle(0))
	ytitle(Change in Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5))
	legend(off)
	note("Note: 95% CIs.", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
	//graph export diff_ci.pdf, replace;
#del cr
exit

* Fig 2d*
keep y x
sum x, meanonly
local r1 = r(min)
local r2 = r(max)
local u = (`r2'-`r1')/50 // unit increase
sdi, xofi(x=(`r1'(`u')`r2')) nunit(1) firstdiff range nolegend nonote

mat res = r(sdi)
mat resq = res[1...,1]
mat resl = res[1..2,6]
mat resu = res[1..2,7]
local sdi = res[1,8]
svmat resq , name(qd)
svmat resl , name(ld_sdi)
svmat resu , name(hd_sdi)
rename (qd1 ld_sdi1 hd_sdi1) (qd ld_sdi hd_sdi)

gen xaxis = 1 in 1
replace xaxis = 51 in 2
replace xaxis = _n-1 in 3/51
sort xaxis

qui margins, at(x=(`r1'(`u')`r2')) at(x=(`=`r1'+1'(`u')`=`r2'+1')) post noatlegend
mat resv = J(51,1,.) 
forval i = 1/51 {
	nlcom _b[`=51+`i''._at]-_b[`i'._at]
	mat resv[`i',1] = r(V)
}
mata st_matrix("resv",sqrt(st_matrix("resv"))*invnormal(.95+(1-.95)/2))
qui logit y x
svmat double resv , name(pm) // plus-minus
rename (pm1) (pm)
gen hd_ci = qd+pm
gen sigd = (qd<0 & hd_ci<0) // significance of diff

sum xaxis, meanonly
#del ;
tw line qd xaxis if sigd==1, lwidth(medium) lcolor(black)
|| line qd xaxis if sigd==0, lpattern(dash) lwidth(medium) lcolor(black)
|| scatter qd xaxis if (ld_sdi!=. & sigd==1), msymbol(S i i) msize(small) mcolor(black) clwidth(medium) clcolor(black)
|| scatter qd xaxis if (ld_sdi!=. & sigd==0), msymbol(Sh i i) msize(small) mcolor(black) clpattern(dash) clwidth(medium) clcolor(black)
|| rcap ld_sdi hd_sdi xaxis if sigd==1, lcolor(black) lwidth(medthick) clwidth(thin) clcolor(black)
|| rcap ld_sdi hd_sdi xaxis if sigd==0, lpattern(dash) lcolor(black) lwidth(medthick) clpattern(dash) clwidth(thin) clcolor(black)
|| ,
	yline(0, lpattern(shortdash) lwidth(medthin) lcolor(black))
	xtitle({it:x}, size(3.6))
	xlabel(`=round(`r(min)',.01)' (`=round((round(`r(max)',.1) - round(`r(min)',.01))/5,.01)') `=round(`r(max)',.01)', labsize(3.4))
	ylabel(-.09(.02)-.03, labsize(3.4) nogrid angle(0))
	ytitle(Change in Pr({it:y}), size(3.6))
	xsca(titlegap(3))
	ysca(titlegap(2.5) ra(-.1 -.03))
	legend(off)
	note("Note: `=string(round(`sdi',.1))'% SDIs.", size(3.5))
	scheme(s2mono) graphregion(fcolor(white)) bgcolor(white);
#del cr
exit
