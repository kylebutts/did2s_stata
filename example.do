********************************************************************************
* Example with simulated data
********************************************************************************

use https://github.com/kylebutts/did2s_stata/raw/main/data/df_hom.dta, clear

* net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")
do ado/did2s.ado

********************************************************************************
* Static
********************************************************************************

** Manual 2SDiD (with incorrect standard errors)
* Step 1: Manually
reg dep_var i.unit i.year if treat == 0, nocons

* Step 2: Regress transformed outcome onto treatment status for all units
predict adj, residuals
reg adj i.treat, vce(cluster state) nocons

** 2SDiD with correct se's
did2s dep_var, first_stage(i.unit i.year) treat_formula(i.treat) treat_var(treat) cluster(state)

* Example esttab
esttab, nobaselevels se


********************************************************************************
* Event Study
********************************************************************************

* factors can't be negative
gen rel_year_shift = rel_year + 20
replace rel_year_shift = 100 if rel_year_shift == .

did2s dep_var, first_stage(i.unit i.year) treat_formula(b100.rel_year_shift) treat_var(treat) cluster(state)


********************************************************************************
* Castle Doctrine
********************************************************************************

use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear


* Covariates
global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44
global spending l_exp_subsidy l_exp_pubwelfare
global xvar l_police unemployrt poverty l_income l_prisoner l_lagprisoner $demo $spending

* No Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) treat_formula(i.post) treat_var(post) cluster(sid)

* Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $xvar) treat_formula(i.post) treat_var(post) cluster(sid)

