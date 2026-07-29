*! equivalence_test_absorb.do -- rev 5
*! section C targets the shipping default (1e-12 vs 1e-14); section D documents keepnoid as a known difference, not an equivalence claim
*! Expects did2s.ado v0.7 with default tolerance(1e-12) maxiter(2000).

*------------------------------------------------------------------------------
* did2s: EQUIVALENCE TEST -- absorb() vs the classic dummy first stage
*
* absorb() is a different ALGORITHM for the same estimator, not a different
* estimator. Wherever every fixed-effect level has untreated support, the two
* paths use the same estimation sample and must agree to solver tolerance.
*
* This complements equivalence_test.do, which checks the classic path against
* the pre-PR package. Here both sides are the current did2s.ado.
*
* SETUP: did2s.ado and df_het.dta in the working directory.
*
* PASS = relative diff < 1e-8 for b, < 1e-6 for V. The two differ; see the
* note in the reporting helper below.
* Section D is a deliberate non-equivalence and is labelled as such.
*------------------------------------------------------------------------------

display as txt _n "{hline 84}"
display as txt "equivalence_test_absorb.do  [rev 5]"
display as txt "{hline 84}"

clear all
set more off
run "did2s.ado"

*-- reporting helper ----------------------------------------------------------
capture mata mata drop _d2s_rep()
mata:
void _d2s_rep(string scalar lbl)
{
    real scalar db, dv, sb, sv, rb, rv
    string scalar vb, vv, verdict

    db = max(abs(st_matrix("D_b")))
    dv = max(abs(st_matrix("D_V")))
    sb = max(abs(st_matrix("b_a")))
    sv = max(abs(st_matrix("V_a")))
    if (sb <= 0) sb = 1
    if (sv <= 0) sv = 1
    rb = db / sb
    rv = dv / sv

    /* Separate thresholds, because the two columns have different floors.
       Point estimates agree with the dummy path to ~1e-13 once the solver
       tolerance is tight. The variance matrix does not: the dummy path takes
       a generalized inverse of a SINGULAR k1 x k1 matrix with invsym(), whose
       rank determination is tolerance-based, and CG's minimum-norm inverse is
       a different one. The two agree in exact arithmetic; in floating point
       they differ by the conditioning, worst under dispersed weights. The
       observed floor across real datasets is 4e-10 to 4e-08 and does NOT move
       when the solver tolerance is tightened -- see tolerance_sweep.do. */
    vb = (rb < 1e-8 ? "ok" : "DIFF")
    vv = (rv < 1e-6 ? "ok" : "DIFF")
    verdict = ((rb < 1e-8) & (rv < 1e-6) ? "PASS" : "**FAIL**")

    printf("%-40s  b rel=%8.1e %-4s   V rel=%8.1e %-4s   %s\n", lbl, rb, vb, rv, vv, verdict)
}
end

*-- run classic then absorb and compare ---------------------------------------
capture program drop cmppair
program define cmppair
    args lbl

    capture qui did2s $SPEC_CLASSIC
    if (_rc != 0) {
        display as error "  `lbl'   CLASSIC errored, rc = " _rc
        exit
    }
    matrix b_a = e(b)
    matrix V_a = e(V)
    scalar N_a = e(N)

    capture qui did2s $SPEC_ABSORB
    if (_rc != 0) {
        display as error "  `lbl'   ABSORB errored, rc = " _rc
        exit
    }
    matrix b_b = e(b)
    matrix V_b = e(V)
    scalar N_b = e(N)

    if (colsof(b_a) != colsof(b_b)) {
        display as error "  `lbl'   DIFFERENT # OF COEFFICIENTS"
        exit
    }
    if (N_a != N_b) {
        display as error "  `lbl'   e(N) DIFFERS: " N_a " vs " N_b ///
            "  (samples not comparable -- check the identification screen)"
    }

    matrix D_b = b_a - b_b
    matrix D_V = V_a - V_b
    mata: _d2s_rep(st_local("lbl"))
end

*==============================================================================
* DATA
*==============================================================================
use "df_het.dta", clear

set seed 90210
gen double wt      = 1 + runiform()
gen int    wt_int  = 1 + int(runiform()*4)
gen byte   keeprow = runiform() > 0.10
gen long   region  = ceil(state/4)
gen double xcov    = rnormal()
gen byte   treat2  = treat
gen rel_year_shift = rel_year + 20
qui replace rel_year_shift = 100 if missing(rel_year_shift)

* a third categorical dimension, crossed with the others
gen int cohort = mod(unit, 7) + 1

display as txt _n "{hline 84}"
display as txt "A. CORE EQUIVALENCE -- absorb() must reproduce the dummy path"
display as txt "{hline 84}"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "A01 static, state+year FE"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
cmppair "A02 event study, state+year FE"

global SPEC_CLASSIC "dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
global SPEC_ABSORB  "dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmppair "A03 static, UNIT+year FE"

global SPEC_CLASSIC "dep_var, first_stage(i.unit i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)"
global SPEC_ABSORB  "dep_var, absorb(unit year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)"
cmppair "A04 event study, UNIT+year FE"

global SPEC_CLASSIC "dep_var, first_stage(i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "A05 one dimension only"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year i.cohort) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(state year cohort) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "A06 three dimensions"

display as txt _n "{hline 84}"
display as txt "B. WEIGHTS, RESTRICTIONS, COVARIATES"
display as txt "{hline 84}"

global SPEC_CLASSIC "dep_var [aw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var [aw = wt], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B01 aweight"

global SPEC_CLASSIC "dep_var [fw = wt_int], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var [fw = wt_int], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B02 fweight"

global SPEC_CLASSIC "dep_var [pw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var [pw = wt], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B03 pweight"

global SPEC_CLASSIC "dep_var [iw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var [iw = wt], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B04 iweight"

* Weight EXPRESSIONS. absorb() materializes the expression into a variable
* because Mata's st_data() needs a name. first_stage() still hands the raw
* expression to st_view() and errors -- pre-existing, untouched by this PR,
* and worth reporting separately. Checked here against the equivalent
* pre-computed weight rather than against the classic path.
capture drop wt2
gen double wt2 = wt*2
global SPEC_CLASSIC "dep_var [aw = wt2], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var [aw = wt*2], absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B05 weight expression (absorb only)"

global SPEC_CLASSIC "dep_var if keeprow, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var if keeprow, absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B06 if restriction"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year xcov) second_stage(i.treat) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, first_stage(xcov) absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
cmppair "B07 continuous covariate in stage 1"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat#c.xcov) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat i.treat#c.xcov) treatment(treat) cluster(state)"
cmppair "B08 interaction in stage 2"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat2) treatment(treat) cluster(state)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat i.treat2) treatment(treat) cluster(state)"
cmppair "B09 collinear stage-2 term (omitted)"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(unit)"
cmppair "B10 cluster finer than FE"

global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(region)"
global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(region)"
cmppair "B11 cluster coarser than FE"

preserve
    replace dep_var = . if mod(_n, 37) == 0
    global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
    global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
    cmppair "B12 missing outcome values"
restore

preserve
    drop if runiform() < 0.15
    global SPEC_CLASSIC "dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
    global SPEC_ABSORB  "dep_var, absorb(state year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
    cmppair "B13 unbalanced panel"
restore

*==============================================================================
* C. TOLERANCE SENSITIVITY
*
* The absorb path is iterative. Tightening past the default should move nothing
* at reported precision; if it does, the default is too loose. Compared against
* 1e-14 rather than a looser setting, so this tests the shipping default rather
* than documenting an old one.
*==============================================================================
display as txt _n "{hline 84}"
display as txt "C. TOLERANCE SENSITIVITY -- 1e-12 (default) vs 1e-14"
display as txt "{hline 84}"

capture program drop tolpair
program define tolpair
    args lbl
    qui did2s $SPEC_ABSORB tolerance(1e-12)
    matrix b_a = e(b)
    matrix V_a = e(V)
    qui did2s $SPEC_ABSORB tolerance(1e-14) maxiter(50000)
    matrix D_b = b_a - e(b)
    matrix D_V = V_a - e(V)
    mata: _d2s_rep(st_local("lbl"))
end

global SPEC_ABSORB "dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)"
tolpair "C01 static, state+year"

global SPEC_ABSORB "dep_var, absorb(unit year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)"
tolpair "C02 event study, unit+year"

global SPEC_ABSORB "dep_var [aw = wt], absorb(state year cohort) second_stage(i.treat) treatment(treat) cluster(state)"
tolpair "C03 weighted, three dimensions"

*==============================================================================
* D. NON-IDENTIFIED FIXED-EFFECT LEVELS
*
* Add units that are treated in EVERY period. Their unit effect has no
* untreated observation to identify it.
*
* Expected:
*   - default absorb() drops those rows and reports the count
*   - keepnoid retains them; the first-stage system is then inconsistent, so
*     the solver should report non-convergence
*   - keepnoid does NOT reproduce the classic dummy path. This was tested and
*     is false: point estimates differ by ~1e-3 relative. The two paths zero
*     the unidentified level inside different parameterizations -- the classic
*     path runs reg with an explicit constant column and reg's own omission
*     rules, absorb maps unsupported levels to zero with no separate constant.
*     They agree wherever every level has support and have no reason to agree
*     where one does not.
*
* D01 is therefore EXPECTED TO FAIL in both columns. It is kept as a regression
* test on the size of the disagreement, not as an equivalence claim.
*==============================================================================
display as txt _n "{hline 84}"
display as txt "D. NON-IDENTIFIED LEVELS -- expected differences, read the notes"
display as txt "{hline 84}"

preserve
    * ten always-treated units
    qui replace treat = 1 if inrange(unit, 1, 10)
    qui replace rel_year_shift = 20 if inrange(unit, 1, 10)

    display as txt _n "  default absorb() -- expect a drop message:"
    qui did2s dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(unit)
    display as txt "    e(N) with screen = " e(N)

    display as txt _n "  keepnoid vs classic dummies (expected to differ):"
    global SPEC_CLASSIC "dep_var, first_stage(i.unit i.year) second_stage(i.treat) treatment(treat) cluster(unit)"
    global SPEC_ABSORB  "dep_var, absorb(unit year) second_stage(i.treat) treatment(treat) cluster(unit) keepnoid"
    cmppair "D01 keepnoid vs dummies [EXPECTED FAIL]"
restore

display as txt _n "{hline 84}"
display as txt "Sections A, B and C must all read PASS."
display as txt "Section D: D01 is expected to FAIL in both columns -- see the notes."
display as txt "Run it without capture/qui once to confirm the non-convergence warning fires."
display as txt "{hline 84}"
