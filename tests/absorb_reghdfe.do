*! Regression tests for the reghdfe-based did2s estimator

clear all
set more off
set seed 20260917

capture which reghdfe
if (_rc != 0) {
    display as error "tests/absorb_reghdfe.do requires reghdfe"
    exit 199
}

run "ado/did2s.ado"
use "data/df_het.dta", clear

gen double x = rnormal()
gen double wt = 1 + mod(unit, 7)
gen int rel = rel_year + 20
replace rel = 100 if missing(rel)

* Static reference values documented by the package example.
quietly did2s dep_var, absorb(state year) second_stage(i.treat) ///
    treatment(treat) cluster(state)
matrix b = e(b)
matrix V = e(V)
local k = colsof(b)
assert abs(b[1, `k'] - 2.380156) < 1e-5
assert abs(sqrt(V[`k', `k']) - .0614383) < 1e-5
assert e(N) == _N

* Event study.
quietly did2s dep_var, absorb(state year) second_stage(ib100.rel) ///
    treatment(treat) cluster(state)
assert e(N) == _N
matrix b = e(b)
assert colsof(b) > 1

* Weights and a continuous first-stage covariate.
quietly did2s dep_var [aw=wt], first_stage(x) absorb(unit year) ///
    second_stage(i.treat) treatment(treat) cluster(state)
assert e(N) == _N
matrix V = e(V)
assert mreldif(V, V') < 1e-12

* if restrictions flow through both stages and e(sample).
quietly count if mod(unit, 3)
local restricted_N = r(N)
quietly did2s dep_var if mod(unit, 3), absorb(unit year) ///
    second_stage(i.treat) treatment(treat) cluster(state)
assert e(N) == `restricted_N'
gen byte restricted_sample = e(sample)
assert restricted_sample == 0 if !mod(unit, 3)
drop restricted_sample

* Unsupported fixed-effect levels are dropped, exactly as if excluded by hand.
preserve
replace treat = 1 if unit <= 10

quietly did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(unit)
matrix b_auto = e(b)
matrix V_auto = e(V)
scalar N_auto = e(N)
gen byte absorb_sample = e(sample)

quietly did2s dep_var if unit > 10, absorb(unit year) ///
    second_stage(i.treat) treatment(treat) cluster(unit)
matrix b_manual = e(b)
matrix V_manual = e(V)

assert mreldif(b_auto, b_manual) < 1e-8
assert mreldif(V_auto, V_manual) < 1e-6
assert absorb_sample == 0 if unit <= 10
assert N_auto == _N - 310
restore

* absorb() is required.
capture noisily did2s dep_var, second_stage(i.treat) ///
    treatment(treat) cluster(state)
assert _rc == 198

preserve
replace treat = 2 in 1
capture noisily did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(state)
assert _rc == 459
restore

capture noisily did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(state) tolerance(0)
assert _rc == 198

capture noisily did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(state) maxiter(0)
assert _rc == 198

display as result "All did2s regression tests passed."
