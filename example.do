********************************************************************************
* Example with simulated data
********************************************************************************

use data/df_hom.dta, clear

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

