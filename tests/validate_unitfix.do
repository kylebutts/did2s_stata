*------------------------------------------------------------------------------
* VALIDATION for the unit()+weights fix (did2s_unitfix.ado)
*
* Two claims, both must hold:
*   A. WITH weights, unit() now reproduces the equivalent full-dummy
*      specification (the FWL ground truth) to machine precision, at every
*      weight dispersion. Before the fix the error grew with dispersion.
*   B. WITHOUT weights (and with constant within-unit weights), behavior is
*      UNCHANGED relative to the current package version.
*
* SETUP: did2s_unitfix.ado (program name did2s) and did2s_orig.ado in the
* working directory.
*------------------------------------------------------------------------------

clear all
set more off
run "did2s_orig.ado"
run "did2s_unitfix.ado"

* ---- data (same DGP as the bug report) ----------------------------------
set seed 4242
set obs 60
gen int state = _n
gen int g = cond(state <= 30, 8, 9999)
expand 25
bysort state: gen int year = _n
gen byte treat = year >= g
gen double state_fe = rnormal()
bysort state: replace state_fe = state_fe[1]
gen double year_fe = year / 10
gen double y = state_fe + year_fe + 2*treat + rnormal()

* ---- A: weighted unit() vs full dummies, sweeping dispersion ------------
display as txt _n "{hline 72}"
display as txt "A. fixed unit() vs full dummies under weights (rel. error must be ~1e-15)"
display as txt "{hline 72}"
display as txt "  spread     dummies (truth)      unit(), fixed    rel. error"
display as txt "{hline 72}"

foreach spread in 1 2 5 20 100 {
    capture drop wt
    gen double wt = 1 + (`spread' - 1) * ((state - 1) / 59)

    qui did2s y [aw = wt], first_stage(i.state i.year) ///
        second_stage(i.treat) treatment(treat) cluster(state)
    scalar b_truth = _b[1.treat]
    scalar se_truth = _se[1.treat]

    qui did2s y [aw = wt], first_stage(i.year) unit(state) ///
        second_stage(i.treat) treatment(treat) cluster(state)
    scalar b_fix = _b[1.treat]
    scalar se_fix = _se[1.treat]

    scalar relb  = abs(b_fix - b_truth) / abs(b_truth)
    scalar relse = abs(se_fix - se_truth) / abs(se_truth)

    display as txt "  " %6.0f `spread' "   " as res %14.8f b_truth ///
        "  " %14.8f b_fix as txt "   b: " as res %9.1e relb ///
        as txt "  se: " as res %9.1e relse
}

* ---- B: unweighted, fixed version vs current package --------------------
display as txt _n "{hline 72}"
display as txt "B. no weights: fixed version must match the current package exactly"
display as txt "{hline 72}"

qui did2s_orig y, first_stage(i.year) unit(state) ///
    second_stage(i.treat) treatment(treat) cluster(state)
scalar b_o  = _b[1.treat]
scalar se_o = _se[1.treat]
qui did2s y, first_stage(i.year) unit(state) ///
    second_stage(i.treat) treatment(treat) cluster(state)
scalar b_f  = _b[1.treat]
scalar se_f = _se[1.treat]
display as txt "  rel diff b  = " as res %10.2e abs(b_f - b_o)/abs(b_o)
display as txt "  rel diff se = " as res %10.2e abs(se_f - se_o)/abs(se_o)

* constant within-unit weights: also unchanged (weighted mean == mean)
capture drop wt
gen double wt = 1 + mod(state, 7)
qui did2s_orig y [aw = wt], first_stage(i.year) unit(state) ///
    second_stage(i.treat) treatment(treat) cluster(state)
scalar b_o  = _b[1.treat]
scalar se_o = _se[1.treat]
qui did2s y [aw = wt], first_stage(i.year) unit(state) ///
    second_stage(i.treat) treatment(treat) cluster(state)
scalar b_f  = _b[1.treat]
scalar se_f = _se[1.treat]
display as txt "  constant-within-unit weights:"
display as txt "  rel diff b  = " as res %10.2e abs(b_f - b_o)/abs(b_o)
display as txt "  rel diff se = " as res %10.2e abs(se_f - se_o)/abs(se_o)

display as txt _n "{hline 72}"
display as txt "PASS = all of A at ~1e-15 or smaller, all of B at ~1e-15 or smaller."
display as txt "{hline 72}"
