{smcl}
{* *! version 1 7Dec2022}{...}
{viewerjumpto "Syntax" "sdii##syntax"}{...}
{viewerjumpto "Description" "sdii##description"}{...}
{viewerjumpto "Options" "sdii##options"}{...}
{viewerjumpto "Examples" "sdii##examples"}{...}
{viewerjumpto "The User's Manual" "sdii##manual"}{...}
{viewerjumpto "Stored results" "sdii##results"}{...}
{title:Title}

{phang}
{cmd:sdii} {hline 2} compute significance of differences intervals (SDIs)


{marker syntax}{...}
{title:Syntax}

{p 4 17 2}
{it:Compare two distributions or paired samples with set characteristics}

{p 8 17 2}
{cmdab:sdii}
{cmd: ,}
sample1([{it:#obs1}] {it:#mean1 #sd1}) sample2([{it:#obs2}] {it:#mean2 #sd2}) [{it:options options1}]

{p 4 17 2}
{it:Compare two unpaired samples with set characteristics}

{p 8 17 2}
{cmdab:sdii}
{cmd: ,}
sample1({it:#obs1 #mean1 #sd1}) sample2({it:#obs2 #mean2 #sd2}) {ul:unp}aired [{it:options options2}]

{p 4 17 2}
{it:Compare two paired variable samples}

{p 8 17 2}
{cmdab:sdii}
{it:varname1 varname2}
{ifin}
{weight}
[{cmd:,}
{it:options options1}]

{p 6 17 2}
or

{p 8 17 2}
{cmdab:sdii}
{it:varname1}
{ifin}
{weight}
{cmd: ,}
varname2({it:exp}) [{it:options options1}]

{p 4 17 2}
{it:Compare two unpaired variable samples}

{p 8 17 2}
{cmdab:sdii}
{it:varname1 varname2}
{ifin}
{weight}
{cmd: ,}
{ul:unp}aired [{it:options options2}]

{p 6 17 2}
or

{p 8 17 2}
{cmdab:sdii}
{it:varname1}
{ifin}
{weight}
{cmd: ,}
varname2({it:exp}) {ul:unp}aired [{it:options options2}]

{p 6 17 2}
or

{p 8 17 2}
{cmdab:sdii}
{it:varname1}
{ifin}
{weight}
{cmd: ,}
by({it:groupvar}) [{it:options options2}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt diff:erence}} report the difference in estimates{p_end}
{synopt:{opt l:evel(#)}} set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt m:value(#)}} set meaningful value {it:m}; default is 0{p_end}
{synopt:{opt nol:egend}} suppress output legend{p_end}
{synopt:{cmd:power}[(#)]} power; default is {cmd:power(0.8)}{p_end}
{synopt:{opt prec:ision(#)}} set level of precision for SDIs; default is 1{p_end}
{synopt:{opt rev:erse}} reverse order for the difference in estimates computation{p_end}
{synopt:{cmd:sample1(}[{it:#obs1}] {it:#mean1 #sd1}{cmd:)}} input the number of observations, mean, and standard deviation of the first sample{p_end}
{synopt:{cmd:s1(}[{it:#obs1}] {it:#mean1 #sd1}{cmd:)}} shorthand for {cmd:sample1(}[{it:#obs1}] {it:#mean1 #sd1}{cmd:)}{p_end}
{synopt:{cmd:sample2(}[{it:#obs2}] {it:#mean2 #sd2}{cmd:)}} input the number of observations, mean, and standard deviation of the second sample{p_end}
{synopt:{cmd:s2(}[{it:#obs2}] {it:#mean2 #sd2}{cmd:)}} shorthand for {cmd:sample2(}[{it:#obs2}] {it:#mean2 #sd2}{cmd:)}{p_end}
{synopt:{opt unp:aired}} treat data as unpaired{p_end}
{synopt:{cmd:variable2(}{it:exp}{cmd:)}} generates the second comparison sample as defined by the expression {it:exp}{p_end}
{synopt:{cmd:var2(}{it:exp}{cmd:)}} shorthand for {cmd:variable2(}{it:exp}{cmd:)}{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 30 tabbed}{...}
{synopthdr:options1 (with paired data)}
{synoptline}
{synopt:{opt corr:elation}(#)} set correlation level; default is 0{p_end}
{synopt:{opt dist:ribution}(#)} treat data as sampling distributions{p_end}
{synopt:{cmd:percentile}[({it:{ul:alt}def)}]} compute percentile SDIs instead of the default, standard error-based SDIs; with suboption {it:altdef} an alternative formula is used to calculate percentiles{p_end}
{synopt:{bf:pctile}[({it:{ul:alt}def)}]} shorthand for {cmd:percentile}[({it:{ul:alt}def)}]{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 30 tabbed}{...}
{synopthdr:options2 (with unpaired data)}
{synoptline}
{synopt:{cmd:by(}{it:groupvar}{cmd:)}} variable defining the groups{p_end}
{synopt:{opt une:qual}} unpaired data have unequal variances{p_end}
{synopt:{opt w:elch}} use Welch's approximation{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
Note: Syntax elements within square brackets [] are optional. Underlining indicates minimal abbreviation.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:sdii} computes SDIs to indicate whether two points estimates (typically the mean of distributions or of samples from the study population) are statistically or substantively distinct		


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt diff:erence} reports the point estimate for the difference in estimates (typically the difference in means) with its standard CI. The default is the 95% CI or as set by {cmd:set level}.

{phang}
{opt l:evel(#)} specifies the confidence level, as a percentage, for confidence intervals. The default is {cmd:level(95)} or as set by {cmd:set level}.

{phang}
{opt m:value(#)} computes SDIs that indicate significance of differences other than zero. By typing {cmd:mvalue(1.3)}  we test whether the difference between the compared estimates is greater than 1.3 rather than 0. # can be either positive or negative; by default # = 0.

{phang}
{opt nolegend} specifies that the legend detailing the output items and the explanatory notes be suppressed.

{phang}
{cmd:power}[(#)] computes the minimum detectable value of the effect size δ (the difference in estimates), given power, significance level α, sample characteristics (size and standard deviation), and, for paired data only, the correlation level ρ. Power is technically defined as (1-β), where β is the probability of type II error. # sets the power within the (0, 1) range; by default # = 0.8. The significance level is calculated from the specified level for confidence intervals, α = (1-{cmd:level}/100).

{phang}
{opt prec:ision(#)} specifies the decimal precision of reported SDIs, ranging from no decimals to up to six decimals (# ∈ {0,1,2,3,4,5,6}). The difference between typing {cmd:precision(0)} or {cmd:precision(3)} is that an 85-ish% SDI will be rounded to either 85% or 85.xxx%, respectively. The default is to report one digit precision SDIs, i.e., # = 1.

{phang}
{opt reverse} reverses the order of the two samples when calculating the difference in estimates. By default {cmd:sdii} observes the order in which the samples are listed in {it:varlist} or {cmd:sample*()}. Specifically, {it:var2} or {cmd:sample2()} is subtracted from {it:var1} or {cmd:sample1()}, respectively. With {cmd:by()}, the group corresponding to the largest value in the variable in {cmd:by()} is subtracted from the group with the smallest value in {cmd:by()}. {cmd:reverse} reverses this behavior and the order in which variables appear in the table. {cmd:reverse} also reverses the sign of the meaningful value {it:m} (i.e., {it:m} = {cmd:mvalue(-#)}).

{phang}
{cmd:sample1(}[{it:#obs1}] {it:#mean1 #sd1}{cmd:)} specifies the number of observations, mean, and standard deviation, respectively, of the first sample. When {cmd:sample1()} receives only two arguments, these numbers are assumed to represent the mean and standard deviation of a normal distribution rather than of a sample (e.g., {cmd:sample1(}{it:#mean1 #sd1}{cmd:)}).

{phang}
{cmd:sample2(}[{it:#obs2}] {it:#mean2 #sd2}{cmd:)} same as above for the second sample. Specifying samples with different number of observations (i.e., {it:#obs1} != {it:#obs2}), implies {cmd:unpaired}. 

{phang}
{opt unpaired} specifies that the data be treated as unpaired.

{phang}
{cmd:variable2(}{it:exp}{cmd:)} generates the second sample as defined by the expression in the suboption {it:exp}. {it:var2} can be a transformation of the compared sample {it:var1} (e.g., {cmd:variable2(var1^2)}), or of another exiting variable (e.g., {cmd:variable2(var3+1)}).

{dlgtab:With paired data}

{phang}
{opt corr:elation}(#) indicates the correlation level, ρ, for paired data. It cannot be combined with {it:varlist} as in this case ρ is the observed correlation level between the two variables. By default ρ = 0, and its range is [-1,1].

{phang}
{opt dist:ribution}(#) specifies that the compared variables are sampling distributions rather than sample groups. Sampling distributions contain a number of realizations of the quantity of interest (e.g., sample mean), and are typically computed via simulations using either canned Stata commands (e.g., {cmd:bootstrap}) or user-written programs (e.g., Clarify). The relevant difference between sampling distributions and samples from the population, is that the standard error for the latter is a function of the sample size, whereas for the former is simply the standard deviation. {cmd:distribution} does not affect percentile-based SDIs, since they are not a function of the standard error.

{phang}
{opt median} reports the median as the summary statistic value, whereas mean is the default. {cmd:median} can be used only in combination with the {cmd:percentile} option.

{phang}
{cmd:percentile}[({it:{ul:alt}def)}] specifies that percentile SDIs be calculated instead of standard error-based SDIs, which is the default. The default method for calculating percentiles is to invert the empirical distribution function by using averages, (x_i + x_{i+1})/2, where the function is flat. When the suboption {it:altdef} is specified, an alternative formula that uses an interpolation method is employed. Weights cannot be used when {it:altdef} is specified. 

{dlgtab:With unpaired data}

{phang}
{cmd:by(}{it:groupvar}{cmd:)} specifies the {it:groupvar} that defines the two group samples to be compared. Specifying {cmd:by()} implies {cmd:unpaired}.

{phang}
{opt unequal} specifies that the unpaired data not to be assumed to have equal variances.

{phang}
{opt welch} specifies that the approximate degrees of freedom for the test be obtained from Welch's (1947) formula rather than from the Satterthwaite's (1946) approximation formula, which is the default when {cmd:unequal} is specified. Specifying {cmd:welch} implies {cmd:unequal}  and {cmd:unpaired}.


{marker examples}{...}
{title:Examples}

{pstd}
These examples are intended for quick reference. For a more detailed overview of {cmd:sdii} and examples with discussion, see {browse "sdii_manual.sdii":The {cmd:sdii} User's Manual}.  

{pstd}
Example 1: Compute SDIs to compare independent distributions

{phang2}{cmd: . sdii, sample1(10 2) sample2(5 1)}{p_end}

{pstd}
Example 2: Compute SDIs to indicate substantive significance

{phang2}{cmd: . sdii, sample1(10 2) sample2(5 1) mvalue(1) difference}{p_end}

{pstd}
Example 3: Compute SDIs to compare correlated distributions

{phang2}{cmd: . sample1(10 2) sample2(5 1) m(1) corr(.5) diff}{p_end}

{pstd}
Example 4: Compute SDIs to compare paired samples

{phang2}{cmd: . sdii, sample1(40 10 2) sample2(40 5 4) level(90)}{p_end}

{pstd}
Example 5: Compute SDIs to compare unpaired samples

{phang2}{cmd: . sdii, sample1(60 10 2) sample2(40 5 4) unpaired}{p_end}

{pstd}
Example 6: Compute SDIs to compare unpaired samples with unequal population variances

{phang2}{cmd: . sdii, sample1(60 10 2) sample2(40 5 4) unpaired unequal welch}{p_end}

{pstd}
Example 7: Compute SDIs to compare paired variable samples

{phang2}{cmd: . use https://www.stata-press.com/data/r17/fuel, clear}{p_end}
{phang2}{cmd: . sdii mpg1 mpg2}{p_end}

{pstd}
Example 8: Computing percentile SDIs

{phang2}{cmd: . sdii mpg1 mpg2, percentile median}{p_end}

{pstd}
Example 9: Compute SDIs to compare unpaired variable samples

{phang2}{cmd: . sdii mpg1 mpg2, unpaired}{p_end}

{pstd}
Example 10: Compute SDIs to compare group samples within the same variable

{phang2}{cmd: . stack mpg1 mpg2, into(mpg) clear}{p_end}
{phang2}{cmd: . sdii mpg, by(_stack) unpaired}{p_end}

{pstd}
Example 11: Compute SDIs to compare estimated statistics

{phang2}{cmd: . webuse nhanes2f, clear}{p_end}
{phang2}{cmd: . logit diabetes age i.race, nolog}{p_end}
{phang2}{cmd: . matrix eb = e(b)}{p_end}
{phang2}{cmd: . matrix eV = e(V)}{p_end}
{phang2}{cmd: . tempname m1 s1 m2 s2 cor12}{p_end}
{phang2}{cmd: . scalar `m1' = el(eb,1,3)}{p_end}
{phang2}{cmd: . scalar `s1' = sqrt(el(eV,3,3))}{p_end}
{phang2}{cmd: . scalar `m2' = el(eb,1,4)}{p_end}
{phang2}{cmd: . scalar `s2' = sqrt(el(eV,4,4))}{p_end}
{phang2}{cmd: . scalar `cor12' = el(eV,4,3)/(sqrt(el(eV,3,3))*sqrt(el(eV,4,4)))}{p_end}
{phang2}{cmd: . sdii, sample1(`=`m1'' `=`s1'') sample2(`=`m2'' `=`s2'') corr(`=`cor12'') diff}{p_end}

{pstd}
Example 12: Compute SDIs to compare {cmd:margins} results

{phang2}{cmd: . quietly sum age, detail}{p_end}
{phang2}{cmd: . margins, at(age=(`r(p25)' `r(p75)'))}{p_end}
{phang2}{cmd: . tempname m1 s1 m2 s2 cor12}{p_end}
{phang2}{cmd: . scalar `m1' = el(r(b),1,1)}{p_end}
{phang2}{cmd: . scalar `s1' = sqrt(el(r(V),1,1))}{p_end}
{phang2}{cmd: . scalar `m2' = el(r(b),1,2)}{p_end}
{phang2}{cmd: . scalar `s2' = sqrt(el(r(V),2,2))}{p_end}
{phang2}{cmd: . scalar `cor12' = el(r(V),2,1)/(sqrt(el(r(V),1,1))*sqrt(el(r(V),2,2)))}{p_end}
{phang2}{cmd: . sdii, sample1(`=`m1'' `=`s1'') sample2(`=`m2'' `=`s2'') corr(`=`cor12'') diff}{p_end}


{marker manual}{...}
{title:The User's Manual}

{pstd}
{browse "sdii_manual.pdf":The {cmd:sdii} User's Manual}{p_end}

	   
{marker results}{...}
{title:Stored results}

{pstd}
{cmd:sdii} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(level)}}confidence level of confidence intervals{p_end}
{synopt:{cmd:r(mvalue)}}meaningful value {it:m}{p_end}
{synopt:{cmd:r(sd_1)}}standard deviation for distribution or sample 1{p_end}
{synopt:{cmd:r(sd_2)}}standard deviation for distribution or sample 2{p_end}
{synopt:{cmd:r(sd_d)}}with {cmd:difference}; (combined) standard deviation for the difference in estimates{p_end}

{synopt:{it:For samples}}{p_end}
{synopt:{cmd:r(alpha)}}with {cmd:power}; significance level{p_end}
{synopt:{cmd:r(delta)}}with {cmd:power}; minimum effect size{p_end}
{synopt:{cmd:r(N_1)}}sample size for sample 1{p_end}
{synopt:{cmd:r(N_2)}}sample size for sample 2{p_end}
{synopt:{cmd:r(N_d)}}with {cmd:difference}; (combined) sample size for the difference in estimates{p_end}
{synopt:{cmd:r(power)}}with {cmd:power}; power{p_end}

{synopt:{it:For paired data}}{p_end}
{synopt:{cmd:r(rho)}}correlation level{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:sdii}{p_end}
{synopt:{cmd:r(cmdline)}}command as typed{p_end}
{synopt:{cmd:r(type)}}type of the uncertainty interval{p_end}

{p2col 5 20 24 2:Matrices}{p_end}
{synopt:{cmd:r(sdii)}}matrix containing the compared statistics with their standard errors, test statistics, p-values, degrees of freedom, critical values, lower and upper confidence limits, and confidence level of the uncertainty{p_end}
{p2colreset}{...}
