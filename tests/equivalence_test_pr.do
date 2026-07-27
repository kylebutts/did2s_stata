*------------------------------------------------------------------------------
* did2s: EQUIVALENCE TEST -- did2s_pr.ado vs the ORIGINAL package (version 0.5)
*
* Runs the original and the proposed (performance-only) code side by side on
* identical data across a battery of specifications, comparing e(b), e(V)
* and e(N).
*
* The PR claims "identical results, faster". Accordingly, EVERY test here
* must PASS, including test 19: the PR deliberately changes no behavior, so
* it reproduces the original's unit()+weights output exactly -- even though
* that output is wrong (see the note at the bottom and the separate bug
* report). Test 20 errors identically in both versions.
*
* SETUP: put did2s_pr.ado, did2s_orig.ado and df_het.dta in the working
* directory (or edit the paths below).
*
* PASS = relative diff < 1e-8 for both b and V, for all of tests 01-19.
*------------------------------------------------------------------------------

clear all
set more off

run "did2s_orig.ado"
run "did2s_pr.ado"

*-- reporting helper -----------------------------------------------------------
capture mata mata drop _d2s_report()
mata:
void _d2s_report(string scalar lbl)
{
    real scalar db, dv, sb, sv, rb, rv
    string scalar verdict

    db = max(abs(st_matrix("D_b")))
    dv = max(abs(st_matrix("D_V")))
    sb = max(abs(st_matrix("b_o")))
    sv = max(abs(st_matrix("V_o")))
    if (sb <= 0) {
        sb = 1
    }
    if (sv <= 0) {
        sv = 1
    }
    rb = db / sb
    rv = dv / sv
    verdict = "**FAIL**"
    if (rb < 1e-8) {
        if (rv < 1e-8) {
            verdict = "PASS"
        }
    }
    printf("%-38s  db=%9.2e  rel=%8.1e   dV=%9.2e  rel=%8.1e   %s\n", lbl, db, rb, dv, rv, verdict)
}
end

capture program drop runboth
program define runboth
    args lbl

    capture qui did2s_orig $SPEC
    local rc_o = _rc
    if (`rc_o' != 0) {
        display as error "  " abbrev("`lbl'", 36) "   ORIGINAL errored, rc = `rc_o'"
        exit
    }
    matrix b_o = e(b)
    matrix V_o = e(V)
    scalar N_o = e(N)

    capture qui did2s $SPEC
    local rc_n = _rc
    if (`rc_n' != 0) {
        display as error "  " abbrev("`lbl'", 36) "   OPTIMIZED errored, rc = `rc_n'"
        exit
    }
    matrix b_n = e(b)
    matrix V_n = e(V)
    scalar N_n = e(N)

    if (colsof(b_o) != colsof(b_n)) {
        display as error "  `lbl'   DIFFERENT # OF COEFFICIENTS"
        exit
    }
    if (N_o != N_n) {
        display as error "  `lbl'   e(N) DIFFERS: " N_o " vs " N_n
    }

    matrix D_b = b_o - b_n
    matrix D_V = V_o - V_n
    mata: _d2s_report(st_local("lbl"))
end

*==============================================================================
* DATA
*==============================================================================
use "df_het.dta", clear

* sanity check that the expected variables are present
foreach v in dep_var state year unit treat rel_year {
    capture confirm variable `v'
    if _rc {
        display as error "df_het.dta has no variable `v' -- edit the data prep block"
        exit 111
    }
}

* a weight variable, an unbalanced indicator, and a coarser cluster
set seed 90210
gen double wt      = 1 + runiform()
gen int    wt_int  = 1 + int(runiform()*4)
gen byte   keeprow = runiform() > 0.10
gen long   region  = ceil(state/4)          // coarser than state, nests it
gen double xcov    = rnormal()
gen rel_year_shift = rel_year + 20
qui replace rel_year_shift = 100 if missing(rel_year_shift)

* a copy of treat, to force a collinear (omitted) second-stage term
gen byte treat2 = treat

display as txt _n "{hline 78}"
display as txt "CLASSIC PATH -- must match the original exactly"
display as txt "{hline 78}"

*-- basic ---------------------------------------------------------------------
global SPEC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "01 static, cluster(state)"

global SPEC "dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
runboth "02 event study, cluster(state)"

global SPEC "dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)"
runboth "03 event study, cluster(unit)"

global SPEC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(region)"
runboth "04 static, coarse cluster(region)"

*-- weights -------------------------------------------------------------------
global SPEC "dep_var [aw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "05 aweight"

global SPEC "dep_var [fw = wt_int], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "06 fweight"

global SPEC "dep_var [pw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "07 pweight"

global SPEC "dep_var [iw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "08 iweight"

*-- sample restrictions -------------------------------------------------------
global SPEC "dep_var if keeprow, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "09 if restriction"

global SPEC "dep_var if year >= 3, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "10 if on time window"

*-- specification variants ----------------------------------------------------
global SPEC "dep_var, first_stage(i.state i.year xcov) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "11 continuous covariate in stage 1"

global SPEC "dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat#c.xcov) treatment(treat) cluster(state)"
runboth "12 interaction in stage 2"

global SPEC "dep_var, first_stage(i.state i.year) second_stage(i.treat i.treat2) treatment(treat) cluster(state)"
runboth "13 collinear stage-2 term (omitted)"

global SPEC "dep_var, first_stage(i.state i.year i.state#c.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "14 state-specific trends (i.#c.)"

*-- unit() option -------------------------------------------------------------
global SPEC "dep_var, first_stage(i.year) unit(state) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "15 unit(), no weights"

global SPEC "dep_var, first_stage(i.year) unit(state) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
runboth "16 unit(), event study"

*==============================================================================
* MISSING DATA
*==============================================================================
preserve
    replace dep_var = . if mod(_n, 37) == 0
    global SPEC "dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
    runboth "17 missing outcome values"
restore

preserve
    drop if runiform() < 0.15
    global SPEC "dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)"
    runboth "18 unbalanced panel"
restore

*==============================================================================
* EDGE CASES: known quirks of the original, faithfully reproduced
*==============================================================================
display as txt _n "{hline 78}"
display as txt "EDGE CASES -- the PR must reproduce these exactly, quirks included"
display as txt "{hline 78}"

* The original's unit() demeaning uses UNWEIGHTED untreated-unit means even
* under weights. That is arguably a bug (see the separate issue), but this
* PR changes no behavior, so the two versions must still agree exactly.
global SPEC "dep_var [aw = wt], first_stage(i.year) unit(state) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "19 unit() + weights [MUST MATCH]"

* Both versions pass the raw weight EXPRESSION to Mata's st_view(), which
* only accepts a bare variable name, so a non-trivial expression errors
* inside Mata in BOTH. runboth reports the original's error and stops;
* confirm the PR version errors the same way:
global SPEC "dep_var [aw = wt*2], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
runboth "20 weight expression [BOTH ERROR]"
capture qui did2s $SPEC
display as txt "  20b PR version return code (expect 3598): " as res _rc

display as txt _n "{hline 78}"
display as txt "CONTEXT FOR TEST 19 (documented separately as a bug, NOT fixed here)"
display as txt "{hline 78}"
display as txt "By Frisch-Waugh, first_stage(i.year) unit(state) should reproduce"
display as txt "first_stage(i.state i.year) exactly. Without weights it does (test 15)."
display as txt "With weights, BOTH versions demean by unweighted means, so both deviate"
display as txt "from the full-dummy ground truth by the same amount -- the check below"
display as txt "is EXPECTED to fail, identically for either version:"
display as txt ""

* unit() vs full dummies under weights: deviates in both versions alike
global SPEC "dep_var [aw = wt], first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)"
qui did2s $SPEC
matrix b_o = e(b)
matrix V_o = e(V)
qui did2s dep_var [aw = wt], first_stage(i.year) unit(state) second_stage(i.treat) treatment(treat) cluster(state)
matrix b_n = e(b)
matrix V_n = e(V)
matrix D_b = b_o - b_n
matrix D_V = V_o - V_n
mata: _d2s_report("19b unit() vs dummies [KNOWN BUG]")

display as txt _n "{hline 78}"
display as txt "Tests 01-19 must all read PASS before opening the pull request."
display as txt "(19b documents a pre-existing bug, filed separately; its FAIL is expected.)"
display as txt "{hline 78}"
