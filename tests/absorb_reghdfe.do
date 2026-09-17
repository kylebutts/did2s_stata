*! Regression tests for the reghdfe-based absorb() path

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

capture program drop compare_paths
program define compare_paths
    args label

    quietly did2s $classic
    matrix b_classic = e(b)
    matrix V_classic = e(V)
    scalar N_classic = e(N)

    quietly did2s $absorbed
    matrix b_absorbed = e(b)
    matrix V_absorbed = e(V)
    scalar N_absorbed = e(N)

    mata: st_numscalar("db", max(abs(st_matrix("b_classic") :- st_matrix("b_absorbed"))))
    mata: st_numscalar("dV", max(abs(st_matrix("V_classic") :- st_matrix("V_absorbed"))))

    assert scalar(db) < 1e-8
    assert scalar(dV) < 1e-6
    assert N_classic == N_absorbed
    display as result "PASS: `label'"
end

global classic  "dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global absorbed "dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(state)"
compare_paths "static"

global classic  "dep_var, first_stage(i.unit i.year) second_stage(ib100.rel) treatment(treat) cluster(state)"
global absorbed "dep_var, absorb(unit year) second_stage(ib100.rel) treatment(treat) cluster(state)"
compare_paths "event study"

global classic  "dep_var [aw=wt], first_stage(i.unit i.year x) second_stage(i.treat) treatment(treat) cluster(state)"
global absorbed "dep_var [aw=wt], first_stage(x) absorb(unit year) second_stage(i.treat) treatment(treat) cluster(state)"
compare_paths "weights and covariate"

global classic  "dep_var if mod(unit, 3), first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global absorbed "dep_var if mod(unit, 3), absorb(unit year) second_stage(i.treat) treatment(treat) cluster(state)"
compare_paths "if restriction"

preserve
replace treat = 1 if unit <= 10

quietly did2s dep_var if unit > 10, first_stage(i.unit i.year) ///
    second_stage(i.treat) treatment(treat) cluster(unit)
matrix b_screened = e(b)
matrix V_screened = e(V)

quietly did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(unit)
matrix b_absorbed = e(b)
matrix V_absorbed = e(V)
gen byte absorb_sample = e(sample)

mata: st_numscalar("db", max(abs(st_matrix("b_screened") :- st_matrix("b_absorbed"))))
mata: st_numscalar("dV", max(abs(st_matrix("V_screened") :- st_matrix("V_absorbed"))))
assert scalar(db) < 1e-8
assert scalar(dV) < 1e-6
assert absorb_sample == 0 if unit <= 10
assert e(N) == _N - 310
restore

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

display as result "All absorb() regression tests passed."
