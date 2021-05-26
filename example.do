********************************************************************************
* Example with simulated data
********************************************************************************

use ado/data/df_hom.dta, clear

do ado/did2s.ado

********************************************************************************
* Static
********************************************************************************

gen u = runiform()
gen pw = 1/u

did2s dep_var, first_stage(i.state i.year) treat_formula(i.treat) treat_var(treat) vce(cluster state)

* Example esttab
esttab, nobaselevels se


********************************************************************************
* Event Study
********************************************************************************

* factors can't be negative
gen rel_year_shift = rel_year + 20

did2s dep_var, first_stage(i.state i.year) treat_formula(i.rel_year_shift) treat_var(treat)


********************************************************************************
* Castle Doctrine
********************************************************************************

use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear


* Covariates
global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44
global spending l_exp_subsidy l_exp_pubwelfare
global xvar l_police unemployrt poverty l_income l_prisoner l_lagprisoner $demo $spending

* No Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) treat_formula(i.post) treat_var(post) vce(cluster sid)

* Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $xvar) treat_formula(i.post) treat_var(post) vce(cluster sid)

