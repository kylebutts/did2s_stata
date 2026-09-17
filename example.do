********************************************************************************
* Example with simulated data
********************************************************************************

use https://github.com/kylebutts/did2s_stata/raw/main/data/df_hom.dta, clear

* net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/") replace
do ado/did2s.ado

********************************************************************************
* Static
********************************************************************************

did2s dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(state)

* Example esttab
esttab, nobaselevels se


********************************************************************************
* Event Study
********************************************************************************

* factors can't be negative
gen rel_year_shift = rel_year + 20
replace rel_year_shift = 100 if rel_year_shift == .

did2s dep_var, absorb(unit year) second_stage(b100.rel_year_shift) treatment(treat) cluster(state)


********************************************************************************
* Castle Doctrine
********************************************************************************

use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear


* Covariates
global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44
global spending l_exp_subsidy l_exp_pubwelfare
global xvar l_police unemployrt poverty l_income l_prisoner l_lagprisoner $demo $spending

* No Covariates
did2s l_homicide [aweight=popwt], absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)

* Covariates
did2s l_homicide [aweight=popwt], first_stage($xvar) absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)
