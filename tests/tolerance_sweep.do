*! tolerance_sweep.do -- rev 4 (classic = frozen v0.6 reference, tests/reference/did2s_ref.ado)
*! four-column output (b/V vs classic AND vs absorb 1e-14); section 4 measures the cost of the tighter default
*! Expects did2s.ado v0.7 with default tolerance(1e-12) maxiter(2000).

*------------------------------------------------------------------------------
* did2s absorb(): WHERE IS THE ACCURACY FLOOR, AND WHOSE IS IT?
*
* The real-data equivalence checks fail at 1e-7 to 1e-8 on event-study
* specifications while passing on statics. Both paths carry error:
*
*   absorb()      stops conjugate gradients at a finite relative residual, so
*                 solution error is bounded by roughly tol * kappa, not tol
*   first_stage() inverts a SINGULAR k1 x k1 matrix with invsym(), which has
*                 its own floating-point error
*
* Comparing the two against each other cannot tell them apart. This does, by
* also comparing absorb() against ITSELF at a very tight tolerance:
*
*   "b/V vs a(1e-14)"  -> absorb's OWN convergence error at that tolerance.
*                         Must fall to zero as the tolerance tightens.
*   "b/V vs class"     -> the total gap against the dummy path.
*
* The decisive comparison is the V pair. If "V vs a(1e-14)" goes to zero while
* "V vs class" plateaus, then absorb has fully converged to its own answer and
* the remaining gap belongs to the dummy path -- invsym() taking a generalized
* inverse of a singular matrix by tolerance-based rank detection, which is a
* different generalized inverse than CG's minimum-norm one. If instead BOTH V
* columns plateau at the same value, the difference is ours and needs
* explaining before the PR goes out.
*
* Read the plateau as the achievable agreement between the two algorithms. It
* is the number the PR should quote, and it is what the equivalence threshold
* should be set against.
*
* SETUP: did2s.ado in the working directory; internet for the datasets.
*------------------------------------------------------------------------------

display as txt _n "{hline 84}"
display as txt "tolerance_sweep.do  [rev 3]"
display as txt "{hline 84}"

clear all
set more off
capture confirm file "ado/did2s.ado"
if (_rc == 0) global ROOT "."
else global ROOT ".."
run "$ROOT/ado/did2s.ado"
run "$ROOT/tests/reference/did2s_ref.ado"

capture mata mata drop _d2s_sweep()
mata:
void _d2s_sweep(string scalar tol)
{
    real scalar sb, sv
    sb = max(abs(st_matrix("b_cl"))); if (sb <= 0) sb = 1
    sv = max(abs(st_matrix("V_cl"))); if (sv <= 0) sv = 1
    printf("    %-8s   %10.2e   %10.2e   %10.2e   %10.2e\n", tol,
           max(abs(st_matrix("D1"))) / sb,
           max(abs(st_matrix("D2"))) / sv,
           max(abs(st_matrix("D3"))) / sb,
           max(abs(st_matrix("D4"))) / sv)
}
end

capture program drop tolsweep
program define tolsweep
    args lbl

    * reference: absorb at the tightest tolerance we will use
    qui did2s $SPEC_ABSORB tolerance(1e-14) maxiter(50000)
    matrix b_ref = e(b)
    matrix V_ref = e(V)

    * the classic path (v0.6 explicit dummies), which has no tolerance to set
    capture qui did2s_ref $SPEC_CLASSIC
    if (_rc != 0) {
        display as error "    classic path errored, rc = " _rc
        exit
    }
    matrix b_cl = e(b)
    matrix V_cl = e(V)

    display as txt _n "  `lbl'"
    display as txt "    tol        b vs class   V vs class   b vs a(1e-14)  V vs a(1e-14)"
    foreach t in 1e-6 1e-8 1e-10 1e-12 1e-14 {
        qui did2s $SPEC_ABSORB tolerance(`t') maxiter(50000)
        matrix D1 = b_cl - e(b)
        matrix D2 = V_cl - e(V)
        matrix D3 = b_ref - e(b)
        matrix D4 = V_ref - e(V)
        mata: _d2s_sweep("`t'")
    }
end


*==============================================================================
* 1. CASTLE -- the two specifications that behaved differently
*==============================================================================
display as txt _n "{hline 78}"
display as txt "CASTLE DOCTRINE"
display as txt "{hline 78}"

capture use "castle.dta", clear
if _rc {
    capture use "https://github.com/scunning1975/mixtape/raw/master/castle.dta", clear
}
if _rc {
    display as error "castle.dta unavailable -- skipping"
}
else {
    gen int rel = time_til + 20
    qui replace rel = 100 if missing(time_til)

    * How thin are the event-time cells? This is the conditioning story.
    display as txt _n "  observations per event-time cell (rel = time_til + 20):"
    tab rel

    global SPEC_CLASSIC "l_homicide, first_stage(i.sid i.year) second_stage(i.post) treatment(post) cluster(sid)"
    global SPEC_ABSORB  "l_homicide, absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)"
    tolsweep "1a static, unweighted  [passed at 1e-8]"

    global SPEC_CLASSIC "l_homicide [aw = popwt], first_stage(i.sid i.year) second_stage(ib100.rel) treatment(post) cluster(sid)"
    global SPEC_ABSORB  "l_homicide [aw = popwt], absorb(sid year) second_stage(ib100.rel) treatment(post) cluster(sid)"
    tolsweep "1c event study, weighted  [failed at 1e-8]"
}


*==============================================================================
* 2. MEDICAID EXPANSION
*==============================================================================
display as txt _n "{hline 78}"
display as txt "ACA MEDICAID EXPANSION"
display as txt "{hline 78}"

capture use "https://raw.githubusercontent.com/Mixtape-Sessions/Advanced-DID/main/Exercises/Data/ehec_data.dta", clear
if _rc {
    display as error "ehec_data.dta unavailable -- skipping"
}
else {
    gen byte treated = cond(missing(yexp2), 0, year >= yexp2)
    gen int  rel     = cond(missing(yexp2), 100, year - yexp2 + 20)

    display as txt _n "  observations per event-time cell:"
    tab rel

    global SPEC_CLASSIC "dins, first_stage(i.stfips i.year) second_stage(ib100.rel) treatment(treated) cluster(stfips)"
    global SPEC_ABSORB  "dins, absorb(stfips year) second_stage(ib100.rel) treatment(treated) cluster(stfips)"
    tolsweep "2b event study  [failed at 1e-8]"
}


*==============================================================================
* 3. NLSWORK -- the covariate case, which was the worst of the five
*==============================================================================
display as txt _n "{hline 78}"
display as txt "NLS YOUNG WOMEN"
display as txt "{hline 78}"

capture webuse nlswork, clear
if _rc {
    display as error "nlswork unavailable -- skipping"
}
else {
    gen byte u1 = (union == 1) if !missing(union)
    bysort idcode (year): egen int onset = min(cond(u1 == 1, year, .))
    gen byte treat = cond(missing(onset), 0, year >= onset)
    gen int  rel   = cond(missing(onset), 100, year - onset + 20)
    bysort idcode: egen byte has0 = max(treat == 0)

    keep if has0 == 1
    set seed 20260729
    bysort idcode: gen byte first = (_n == 1)
    gen double r = runiform() if first
    bysort idcode (first): replace r = r[_N]
    keep if r <= 0.09

    global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year) second_stage(ib100.rel) treatment(treat) cluster(idcode)"
    global SPEC_ABSORB  "ln_wage, absorb(idcode year) second_stage(ib100.rel) treatment(treat) cluster(idcode)"
    tolsweep "3b event study  [failed at 1e-8]"

    global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year hours) second_stage(i.treat) treatment(treat) cluster(idcode)"
    global SPEC_ABSORB  "ln_wage, first_stage(hours) absorb(idcode year) second_stage(i.treat) treatment(treat) cluster(idcode)"
    tolsweep "3c continuous covariate  [worst: 1.6e-7]"
}


*==============================================================================
* 4. WHAT DOES THE TIGHTER TOLERANCE COST?
*
* On a balanced, unweighted, two-dimension panel the tolerance is not binding:
* CG's residual drops past 1e-12 in the same iteration it drops past 1e-8, so
* every tolerance gives the same iterate and the same runtime. That is why
* benchmark_scaling.do is insensitive to it.
*
* The cost only appears where the operator is poorly conditioned. This is that
* case: dispersed weights, three absorbed dimensions, an event-study second
* stage (one solve per column of X2).
*==============================================================================
display as txt _n "{hline 78}"
display as txt "TOLERANCE COST -- ill-conditioned design, 4,000 units"
display as txt "{hline 78}"

clear
set seed 5150
set obs 4000
gen long unit = _n
gen int g = ceil(runiform()*8) + 4
qui replace g = 9999 if runiform() < 0.25
gen double u_fe = rnormal()
gen int cohort = mod(unit, 7) + 1
gen double wt = exp(3*runiform())          // weights spanning roughly 20x
expand 20
bysort unit: gen int year = _n
gen byte treat = year >= g
gen int rel_year_shift = cond(g == 9999, 100, year - g + 20)
gen double dep_var = u_fe + 0.02*year + cond(treat, 2, 0) + rnormal()

display as txt "    tol         absorb (s)   vs 1e-8"
local base = .
foreach t in 1e-8 1e-10 1e-12 1e-14 {
    timer clear 5
    timer on 5
    qui did2s dep_var [aw = wt], absorb(unit year cohort) ///
        second_stage(ib100.rel_year_shift) treatment(treat) ///
        cluster(unit) tolerance(`t') maxiter(50000)
    timer off 5
    qui timer list 5
    local tt = r(t5)
    if ("`t'" == "1e-8") local base = `tt'
    local ratio = `tt' / `base'
    display as txt "    " %-10s "`t'" as res %10.2f `tt' as txt "s    " ///
        as res %6.2f `ratio' as txt "x"
}

display as txt _n "If 1e-12 costs more than about 2x here, consider 1e-10 as the"
display as txt "default: absorb is far enough ahead that either is acceptable, but"
display as txt "the PR should quote whichever it ships."

display as txt _n "{hline 78}"
display as txt "HOW TO READ THIS"
display as txt "{hline 78}"
display as txt "Column 3 (b vs abs(1e-14)) must fall monotonically. If it does not,"
display as txt "the solver is not converging and maxiter or the preconditioner is at"
display as txt "fault, not the tolerance."
display as txt ""
display as txt "Column 1 (b vs class) should fall and then flatten. The flat value is"
display as txt "the agreement ceiling between the two algorithms, set by the classic"
display as txt "path's own invsym() error on a singular matrix. Quote it in the PR,"
display as txt "and set the default tolerance one or two orders tighter than the point"
display as txt "where the plateau is reached."
display as txt "{hline 78}"
