*! test_did2s.do -- rev 3 (2026-10-21): adds e(converged), e(depvar), c.x##i.g routing, and the two
*!   non-convergence cases the help file describes (E13-E14)
*------------------------------------------------------------------------------
* did2s v0.7: TEST SUITE
*
* v0.7 solves every first stage by preconditioned conjugate gradients. Plain
* factor main effects in first_stage() (i.x, ib#.x, ibn.x, ...) and unit()
* are absorbed rather than entered as dummies, so there is one code path.
*
* What is checked:
*   A-C  v0.7 reproduces v0.6 (tests/reference/did2s_ref.ado, the frozen
*        upstream main) wherever every fixed-effect level has untreated
*        observations: same sample, b to 1e-8, V to 1e-6 relative.
*   D    Routing: which terms are absorbed, and that every way of writing the
*        same model gives the same answer.
*   E    Fixed-effect levels with no untreated observations: v0.7 drops those
*        rows and must then equal v0.6 run on the remaining rows.
*   F    Errors.
*   G    Solver tolerance: tightening past the default must move nothing.
*
* Thresholds. Point estimates agree to ~1e-12. V cannot agree as tightly:
* v0.6 takes a generalized inverse of a SINGULAR dummy cross-product with
* invsym(), whose rank determination is tolerance-based, and CG's minimum-norm
* inverse is a different one. The observed floor is 4e-10 to 4e-8 and does not
* move with the solver tolerance (tolerance_sweep.do), hence 1e-6 for V.
*
* SETUP: run from the repository root or from tests/. Needs ado/did2s.ado,
* tests/reference/did2s_ref.ado and data/df_het.dta. About 5 minutes, most
* of it v0.6 on the 1,000-unit specifications.
*
* Every line must read PASS. The last line counts failures.
*------------------------------------------------------------------------------

clear all
set more off

capture confirm file "ado/did2s.ado"
if (_rc == 0) global ROOT "."
else global ROOT ".."
run "$ROOT/ado/did2s.ado"
run "$ROOT/tests/reference/did2s_ref.ado"

global NFAIL 0
global NTEST 0

display as txt _n "{hline 92}"
display as txt "test_did2s.do  [rev 3]"
display as txt "{hline 92}"

*-- helpers -------------------------------------------------------------------
capture mata mata drop _t_rep()
mata:
/* Compares b_a/V_a with b_b/V_b. Writes local "pass" back to Stata. */
void _t_rep(string scalar lbl, real scalar tolb, real scalar tolv, real scalar sameN)
{
    real matrix ba, bb, Va, Vb
    real scalar rb, rv, sb, sv, ok

    ba = st_matrix("b_a"); bb = st_matrix("b_b")
    Va = st_matrix("V_a"); Vb = st_matrix("V_b")
    if (cols(ba) != cols(bb)) {
        printf("{err}%-52s  different number of coefficients   **FAIL**\n", lbl)
        st_local("pass", "0")
        return
    }
    sb = max(abs(ba)); if (sb <= 0) sb = 1
    sv = max(abs(Va)); if (sv <= 0) sv = 1
    rb = max(abs(ba - bb)) / sb
    rv = max(abs(Va - Vb)) / sv
    ok = (rb <= tolb) & (rv <= tolv) & sameN
    printf("%-52s  b %8.1e  V %8.1e  %s  %s\n", lbl, rb, rv,
        (sameN ? "     " : "N DIFF"), (ok ? "PASS" : "{err}**FAIL**{txt}"))
    st_local("pass", strofreal(ok))
}
end

capture program drop _t_tally
program define _t_tally
    args pass
    global NTEST = $NTEST + 1
    if ("`pass'" != "1") global NFAIL = $NFAIL + 1
end

* Runs $A and $B (full commands), compares. Optional thresholds.
capture program drop cmp
program define cmp
    args lbl tolb tolv
    if ("`tolb'" == "") local tolb 1e-8
    if ("`tolv'" == "") local tolv 1e-6

    capture qui $A
    if (_rc != 0) {
        display as error %-52s "`lbl'" "  first command errored, rc = " _rc "   **FAIL**"
        _t_tally 0
        exit
    }
    matrix b_a = e(b)
    matrix V_a = e(V)
    local N_a = e(N)

    capture qui $B
    if (_rc != 0) {
        display as error %-52s "`lbl'" "  second command errored, rc = " _rc "   **FAIL**"
        _t_tally 0
        exit
    }
    matrix b_b = e(b)
    matrix V_b = e(V)
    local N_b = e(N)

    local pass 0
    mata: _t_rep("`lbl'", `tolb', `tolv', `N_a' == `N_b')
    _t_tally `pass'
end

* Checks e(absorb) from the last estimation against an expected list.
capture program drop chkabs
program define chkabs
    args lbl expected
    local got "`e(absorb)'"
    local g : list sort got
    local x : list sort expected
    local ok : list g == x
    if (`ok') display as txt %-52s "`lbl'" "  e(absorb) = `got'" _col(98) "PASS"
    else display as error %-52s "`lbl'" "  e(absorb) = `got', expected `expected'   **FAIL**"
    _t_tally `ok'
end

* Checks a scalar condition.
capture program drop chk
program define chk
    args lbl cond detail
    local ok = (`cond')
    if (`ok') display as txt %-52s "`lbl'" "  `detail'" _col(98) "PASS"
    else display as error %-52s "`lbl'" "  `detail'   **FAIL**"
    _t_tally `ok'
end

* Checks that a command fails with the expected return code.
capture program drop chkrc
program define chkrc
    args lbl rc
    capture qui $A
    local got = _rc
    local ok = (`got' == `rc')
    if (`ok') display as txt %-52s "`lbl'" "  rc = `got'" _col(98) "PASS"
    else display as error %-52s "`lbl'" "  rc = `got', expected `rc'   **FAIL**"
    _t_tally `ok'
end

*==============================================================================
* DATA
*==============================================================================
use "$ROOT/data/df_het.dta", clear

set seed 90210
gen double wt      = 1 + runiform()
gen double wt2     = wt*2
gen int    wt_int  = 1 + int(runiform()*4)
gen byte   keeprow = runiform() > 0.10
gen long   region  = ceil(state/4)
gen double xcov    = rnormal()
gen byte   treat2  = treat
gen int    cohort  = mod(unit, 7) + 1
gen rel_year_shift = rel_year + 20
qui replace rel_year_shift = 100 if missing(rel_year_shift)
gen str8   unit_s  = "u" + string(unit, "%04.0f")

global OPT "treatment(treat) cluster(state)"
global ES  "second_stage(ib100.rel_year_shift)"

*==============================================================================
display as txt _n "A. v0.7 vs v0.6, same command"
display as txt "{hline 92}"
*==============================================================================

global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
cmp "A01 static, state+year"

global A "did2s     dep_var, first_stage(i.state i.year) $ES $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year) $ES $OPT"
cmp "A02 event study, state+year"

global A "did2s     dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s_ref dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "A03 static, unit+year (1,000 units)"

global A "did2s     dep_var, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
global B "did2s_ref dep_var, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
cmp "A04 event study, unit+year"

global A "did2s     dep_var, first_stage(i.year) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.year) second_stage(i.treat) $OPT"
cmp "A05 one dimension"

global A "did2s     dep_var, first_stage(i.state i.year i.cohort) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year i.cohort) second_stage(i.treat) $OPT"
cmp "A06 three dimensions"

global A "did2s     dep_var, first_stage(i.state i.year xcov) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year xcov) second_stage(i.treat) $OPT"
cmp "A07 FEs plus a covariate"

global A "did2s     dep_var, first_stage(xcov) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(xcov) second_stage(i.treat) $OPT"
cmp "A08 covariate only (constant absorbed)"

global A "did2s     dep_var, first_stage(i.state i.year i.state#c.year) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year i.state#c.year) second_stage(i.treat) $OPT"
cmp "A09 state-specific trends"

global A "did2s     dep_var, first_stage(i.state i.year c.xcov##i.region) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year c.xcov##i.region) second_stage(i.treat) $OPT"
cmp "A10 c.x##i.g (nested in state)"

global A "did2s     dep_var, first_stage(i.state i.year i(1/3).cohort) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year i(1/3).cohort) second_stage(i.treat) $OPT"
cmp "A11 level-restricted factor i(1/3).x"

global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat#c.xcov) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat#c.xcov) $OPT"
cmp "A12 interaction in stage 2"

global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat2) $OPT"
global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat2) $OPT"
cmp "A13 collinear stage-2 term (omitted)"

global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "A14 cluster finer than FE"

global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(region)"
global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(region)"
cmp "A15 cluster coarser than FE"

global A "did2s     dep_var if keeprow, first_stage(i.state i.year) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var if keeprow, first_stage(i.state i.year) second_stage(i.treat) $OPT"
cmp "A16 if restriction"

preserve
    qui replace dep_var = . if mod(_n, 37) == 0
    global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
    global B "did2s_ref dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
    cmp "A17 missing outcome values"
restore

preserve
    set seed 4711
    qui drop if runiform() < 0.15
    global A "did2s     dep_var, first_stage(i.state i.year) $ES $OPT"
    global B "did2s_ref dep_var, first_stage(i.state i.year) $ES $OPT"
    cmp "A18 unbalanced panel"
restore

*==============================================================================
display as txt _n "B. v0.7 vs v0.6, weights"
display as txt "{hline 92}"
*==============================================================================

foreach wt in "aw = wt" "fw = wt_int" "pw = wt" "iw = wt" {
    global A "did2s     dep_var [`wt'], first_stage(i.state i.year) second_stage(i.treat) $OPT"
    global B "did2s_ref dep_var [`wt'], first_stage(i.state i.year) second_stage(i.treat) $OPT"
    cmp "B [`wt']"
}

global A "did2s     dep_var [aw = wt], first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
global B "did2s_ref dep_var [aw = wt], first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
cmp "B aweight, unit+year event study"

* v0.6 fails on a weight EXPRESSION with first_stage() (it passes the raw
* expression to st_view()). v0.7 materializes it; check against the
* precomputed variable.
global A "did2s     dep_var [aw = wt*2], first_stage(i.state i.year) second_stage(i.treat) $OPT"
global B "did2s_ref dep_var [aw = wt2],  first_stage(i.state i.year) second_stage(i.treat) $OPT"
cmp "B weight expression (v0.6 needs a variable)"

*==============================================================================
display as txt _n "C. unit()"
display as txt "{hline 92}"
*==============================================================================

global A "did2s     dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s_ref dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "C01 unit() vs v0.6 unit()"

global A "did2s     dep_var [aw = wt], first_stage(i.year) unit(unit) $ES treatment(treat) cluster(unit)"
global B "did2s_ref dep_var [aw = wt], first_stage(i.year) unit(unit) $ES treatment(treat) cluster(unit)"
cmp "C02 weighted unit() vs v0.6 unit()"

global A "did2s     dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s     dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "C03 unit(x) same as absorb(x)" 1e-11 1e-9
qui did2s dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)
chkabs "C04 unit() reported in e(absorb)" "year unit"

global A "did2s     dep_var, first_stage(i.year) unit(unit_s) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s     dep_var, first_stage(i.year) unit(unit)   second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "C05 string unit()" 1e-11 1e-9

*==============================================================================
display as txt _n "D. routing"
display as txt "{hline 92}"
*==============================================================================

qui did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT
chkabs "D01 i.state i.year" "state year"

global A "did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
global B "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT"
cmp "D02 first_stage(i.x i.y) identical to absorb(x y)" 1e-14 1e-14

qui did2s dep_var, first_stage(ib(first).state ib1995.year ibn.cohort) second_stage(i.treat) $OPT
chkabs "D03 ib(first). ib#. ibn. are routed" "state year cohort"
global A "did2s dep_var, first_stage(ib(first).state ib1995.year ibn.cohort) second_stage(i.treat) $OPT"
global B "did2s dep_var, absorb(state year cohort) second_stage(i.treat) $OPT"
cmp "D04   ... same results as absorb()" 1e-11 1e-9

* fvunab expands i.(a b); if a Stata version did not, the term would stay a
* regressor and D06 would still pass (correct, just slower).
qui did2s dep_var, first_stage(i.(state year)) second_stage(i.treat) $OPT
display as txt %-52s "D05 i.(state year)" "  e(absorb) = `e(absorb)'" _col(98) "(info)"
global A "did2s dep_var, first_stage(i.(state year)) second_stage(i.treat) $OPT"
global B "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT"
cmp "D06   ... same results as absorb()" 1e-11 1e-9

qui did2s dep_var, first_stage(i.state i.year i.state#c.year) second_stage(i.treat) $OPT
chkabs "D07 interaction i.x#c.y stays a regressor" "state year"

qui did2s dep_var, first_stage(i.state i.year i(1/3).cohort) second_stage(i.treat) $OPT
chkabs "D08 i(1/3).x stays a regressor" "state year"

qui did2s dep_var, first_stage(xcov) second_stage(i.treat) $OPT
chkabs "D09 no FE: nothing reported" ""

global A "did2s dep_var, first_stage(i.state i.year xcov) second_stage(i.treat) $OPT"
global B "did2s dep_var, first_stage(xcov) absorb(state year) second_stage(i.treat) $OPT"
cmp "D10 mixed: covariate + routed FEs = absorb()" 1e-14 1e-14

global A "did2s dep_var, first_stage(i.state i.year) absorb(state) second_stage(i.treat) $OPT"
global B "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT"
cmp "D11 duplicate absorb(x) + i.x" 1e-14 1e-14
qui $A
chkabs "D12   ... reported once" "state year"

global A "did2s dep_var, first_stage(i.unit i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
global B "did2s dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmp "D13 duplicate unit(x) + i.x" 1e-14 1e-14

global A "did2s dep_var, first_stage(i.year i.state) second_stage(i.treat) $OPT"
global B "did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
cmp "D14 term order does not matter" 1e-11 1e-9

qui did2s dep_var, first_stage(i.state i.year c.xcov##i.region) second_stage(i.treat) $OPT
display as txt %-52s "D15 c.x##i.g" "  e(absorb) = `e(absorb)'" _col(98) "(info)"

qui did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT
chk "D16 e(converged) on a regular run" "e(converged) == 1" "e(converged) = `e(converged)'"
local dv "`e(depvar)'"
local okdv = ("`dv'" == "dep_var")
chk "D17 e(depvar)" "`okdv'" "e(depvar) = `dv'"

*==============================================================================
display as txt _n "E. fixed-effect levels with no untreated observations"
display as txt "{hline 92}"
*==============================================================================

preserve
    * ten always-treated units
    qui replace treat = 1 if inrange(unit, 1, 10)
    qui replace rel_year_shift = 20 if inrange(unit, 1, 10)
    gen byte noid = inrange(unit, 1, 10)
    qui count if noid
    local nnoid = r(N)

    qui did2s dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)
    chk "E01 always-treated units dropped" "e(N_noid) == `nnoid'" "e(N_noid) = `e(N_noid)', expected `nnoid'"
    chk "E02   ... e(N) reduced" "e(N) == _N - `nnoid'" "e(N) = `e(N)'"

    global A "did2s     dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
    global B "did2s_ref dep_var if !noid, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
    cmp "E03   ... equals v0.6 on the remaining rows"

    global A "did2s     dep_var, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
    global B "did2s_ref dep_var if !noid, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
    cmp "E04   ... event study"

    * v0.6's unit() already dropped these rows (the demeaned variables are
    * missing there), so v0.7 must match it without the if.
    global A "did2s     dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
    global B "did2s_ref dep_var, first_stage(i.year) unit(unit) second_stage(i.treat) treatment(treat) cluster(unit)"
    cmp "E05   ... unit() equals v0.6 unit() unchanged"

    global A "did2s     dep_var, first_stage(i.unit i.year xcov) second_stage(i.treat) treatment(treat) cluster(unit)"
    global B "did2s_ref dep_var if !noid, first_stage(i.unit i.year xcov) second_stage(i.treat) treatment(treat) cluster(unit)"
    cmp "E06   ... with a covariate"
restore

preserve
    * a final year in which every unit is treated
    qui replace treat = 1 if year == 2020
    qui replace rel_year_shift = 20 if year == 2020 & rel_year_shift == 100
    gen byte noid = year == 2020
    qui count if noid
    local nnoid = r(N)

    qui did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT
    chk "E07 fully treated year dropped" "e(N_noid) == `nnoid'" "e(N_noid) = `e(N_noid)', expected `nnoid'"

    global A "did2s     dep_var, first_stage(i.state i.year) second_stage(i.treat) $OPT"
    global B "did2s_ref dep_var if !noid, first_stage(i.state i.year) second_stage(i.treat) $OPT"
    cmp "E08   ... equals v0.6 on the remaining rows"
restore

preserve
    * 20 units untreated only in the first year: thin but identified, no drop
    qui replace treat = 1 if unit <= 20 & year > 1990
    global A "did2s     dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(state)"
    global B "did2s_ref dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(state)"
    cmp "E09 units untreated in one year only"
    qui $A
    chk "E10   ... nothing dropped" "e(N_noid) == 0" "e(N_noid) = `e(N_noid)'"
restore

preserve
    * No never-treated units: every unit is treated by the last periods, so
    * those periods have no untreated support. (Case 11 of the 2026-10-19
    * comparison, where the published reghdfe branch failed with r(503).)
    clear
    set seed 101
    qui set obs 500
    gen long u = _n
    gen int  g = floor(runiform() * 14) + 4
    gen double a = rnormal()
    gen int  cl = mod(u, 200) + 1
    qui expand 20
    bysort u: gen int t = _n
    gen byte treat = t >= g
    gen double y = a + 0.1 * t + rnormal() + cond(treat, 1 + 0.05 * (t - g), 0)
    gen int rel = t - g + 20
    bysort t: egen byte anyu = max(treat == 0)
    qui count if !anyu
    local nnoid = r(N)

    global A "did2s     y, first_stage(i.u i.t) second_stage(ib19.rel) treatment(treat) cluster(cl)"
    global B "did2s_ref y if anyu, first_stage(i.u i.t) second_stage(ib19.rel) treatment(treat) cluster(cl)"
    cmp "E11 no never-treated units, event study"
    qui $A
    chk "E12   ... fully treated periods dropped" "e(N_noid) == `nnoid'" "e(N_noid) = `e(N_noid)', expected `nnoid'"
restore

preserve
    * Disconnected untreated support. Units 1-50 are untreated in years 1-10,
    * units 51-100 in years 11-20. Every unit and every year has untreated
    * observations, so nothing is dropped, but units 1-50 in years 11-20 have
    * no identified counterfactual: their unit and year effects are estimated
    * in separate components with separate normalizations. The variance solve
    * cannot converge, and did2s must say so.
    clear
    set seed 202
    qui set obs 100
    gen long u = _n
    gen double a = rnormal()
    qui expand 20
    bysort u: gen int t = _n
    gen byte treat = cond(u <= 50, t > 10, t <= 10)
    gen double y = a + 0.1 * t + treat + rnormal()
    qui did2s y, first_stage(i.u i.t) second_stage(i.treat) treatment(treat) cluster(u)
    chk "E13 disconnected untreated support: nothing dropped" "e(N_noid) == 0" "e(N_noid) = `e(N_noid)'"
    chk "E14   ... and flagged as not converged" "e(converged) == 0" "e(converged) = `e(converged)'"
restore

preserve
    * A covariate that is zero on every untreated row but not on treated
    * rows: its coefficient is not identified from untreated data.
    gen double xz = xcov * treat
    qui did2s dep_var, first_stage(i.state i.year xz) second_stage(i.treat) $OPT
    chk "E15 covariate with no untreated variation: not converged" "e(converged) == 0" "e(converged) = `e(converged)'"
restore

*==============================================================================
display as txt _n "F. errors"
display as txt "{hline 92}"
*==============================================================================

global A "did2s dep_var, second_stage(i.treat) $OPT"
chkrc "F01 neither first_stage() nor absorb()" 198
global A "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT tolerance(0)"
chkrc "F02 tolerance(0)" 198
global A "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT maxiter(0)"
chkrc "F03 maxiter(0)" 198
preserve
    qui replace treat = 2 in 1
    global A "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT"
    chkrc "F04 non-0/1 treatment" 450
restore
preserve
    qui replace treat = 1
    global A "did2s dep_var, absorb(state year) second_stage(i.treat) $OPT"
    chkrc "F05 everything treated" 2000
restore

*==============================================================================
display as txt _n "G. solver tolerance: default 1e-12 vs 1e-14"
display as txt "{hline 92}"
*==============================================================================

global A "did2s dep_var, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit)"
global B "did2s dep_var, first_stage(i.unit i.year) $ES treatment(treat) cluster(unit) tolerance(1e-14) maxiter(50000)"
cmp "G01 event study, unit+year"

global A "did2s dep_var [aw = wt], first_stage(i.state i.year i.cohort) second_stage(i.treat) $OPT"
global B "did2s dep_var [aw = wt], first_stage(i.state i.year i.cohort) second_stage(i.treat) $OPT tolerance(1e-14) maxiter(50000)"
cmp "G02 weighted, three dimensions"

*==============================================================================
display as txt _n "{hline 92}"
if ($NFAIL == 0) display as txt "ALL $NTEST CHECKS PASSED"
else display as error "$NFAIL OF $NTEST CHECKS FAILED"
display as txt "{hline 92}"
