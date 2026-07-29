*! validate_real_data.do -- rev 3
*! split PASS thresholds (b < 1e-8, V < 1e-6); castle URL from did2s help examples
*! Expects did2s.ado v0.7 with default tolerance(1e-12) maxiter(2000).

*------------------------------------------------------------------------------
* did2s absorb(): VALIDATION ON REAL DATA
*
* Simulated panels are balanced, complete, integer-keyed and evenly weighted.
* Real panels are none of those things, and that is what breaks code. These
* three exercise the paths simulated data never reaches: string identifiers,
* gappy time variables, skewed population weights, missing outcomes, units
* observed only once, and fixed-effect levels with no untreated support.
*
* Note: the point is NUMERICAL agreement between the two algorithms, not the
* credibility of any research design. Dataset 3 in particular uses a
* constructed treatment; it is a stress test, not an estimate of anything.
*
* SETUP: did2s.ado in the working directory; internet access for datasets 2
* and 3. Dataset 1 ships with the R package (data-raw/castle.dta) and is also
* on the Mixtape Sessions repo.
*------------------------------------------------------------------------------

display as txt _n "{hline 84}"
display as txt "validate_real_data.do  [rev 3]"
display as txt "{hline 84}"

clear all
set more off
run "did2s.ado"

capture mata mata drop _d2s_rrep()
mata:
void _d2s_rrep(string scalar lbl)
{
    real scalar db, dv, sb, sv, rb, rv
    db = max(abs(st_matrix("D_b")))
    dv = max(abs(st_matrix("D_V")))
    sb = max(abs(st_matrix("b_a"))); if (sb <= 0) sb = 1
    sv = max(abs(st_matrix("V_a"))); if (sv <= 0) sv = 1
    rb = db / sb
    rv = dv / sv
    /* b threshold 1e-8, V threshold 1e-6: V has a floor set by the dummy
       path's invsym() on a singular matrix, not by the solver. */
    printf("%-42s  b rel=%8.1e   V rel=%8.1e   %s\n", lbl, rb, rv,
           ((rb < 1e-8) & (rv < 1e-6) ? "PASS" : "**CHECK**"))
}
end

capture program drop cmpreal
program define cmpreal
    args lbl
    timer clear 1
    timer on 1
    capture qui did2s $SPEC_CLASSIC
    timer off 1
    if (_rc != 0) {
        display as error "  `lbl'  CLASSIC errored, rc = " _rc
        exit
    }
    qui timer list 1
    local t_c = r(t1)
    matrix b_a = e(b)
    matrix V_a = e(V)
    local N_a = e(N)

    timer clear 2
    timer on 2
    capture qui did2s $SPEC_ABSORB
    timer off 2
    if (_rc != 0) {
        display as error "  `lbl'  ABSORB errored, rc = " _rc
        exit
    }
    qui timer list 2
    local t_b = r(t2)
    local N_b = e(N)

    if (`N_a' != `N_b') {
        display as error "  `lbl'  e(N) DIFFERS: `N_a' vs `N_b'"
    }
    matrix D_b = b_a - e(b)
    matrix D_V = V_a - e(V)
    mata: _d2s_rrep(st_local("lbl"))
    display as txt "      classic " as res %8.2f `t_c' as txt "s   absorb " ///
        as res %8.2f `t_b' as txt "s   N = `N_a'"
end


*==============================================================================
* 1. CASTLE DOCTRINE  --  Cheng and Hoekstra (2013)
*
* State-year panel, 2000-2010, 50 states, 550 observations. Twenty states pass
* "stand your ground" laws between 2005 and 2009; the rest never do.
*
* Why this one: it is the dataset the R package documents (?did2s, data(castle))
* and it carries POPULATION WEIGHTS (popwt) that span two orders of magnitude
* across states. Weight dispersion of that size is exactly what surfaced the
* unit() demeaning bug, so it is the natural weighted test case. It is also
* small enough to cross-check against the R implementation line by line.
*
* Watch out: `state' is a STRING here. Use `sid' for absorb() and cluster().
*==============================================================================
display as txt _n "{hline 84}"
display as txt "1. CASTLE DOCTRINE (Cheng & Hoekstra 2013) -- 550 obs, 50 states"
display as txt "{hline 84}"

capture use "castle.dta", clear
if _rc {
    * the URL used in did2s's own help file examples
    capture use "https://github.com/scunning1975/mixtape/raw/master/castle.dta", clear
}
if _rc {
    display as error "castle.dta not found -- copy it from the R package (data-raw/castle.dta)"
}
else {
    * event-study relative time, never-treated pooled into a reference bin
    capture confirm variable time_til
    if !_rc {
        gen int rel = time_til + 20
        qui replace rel = 100 if missing(time_til)
    }

    global SPEC_CLASSIC "l_homicide, first_stage(i.sid i.year) second_stage(i.post) treatment(post) cluster(sid)"
    global SPEC_ABSORB  "l_homicide, absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)"
    cmpreal "1a static, unweighted"

    global SPEC_CLASSIC "l_homicide [aw = popwt], first_stage(i.sid i.year) second_stage(i.post) treatment(post) cluster(sid)"
    global SPEC_ABSORB  "l_homicide [aw = popwt], absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)"
    cmpreal "1b static, POPULATION WEIGHTED"

    capture confirm variable rel
    if !_rc {
        global SPEC_CLASSIC "l_homicide [aw = popwt], first_stage(i.sid i.year) second_stage(ib100.rel) treatment(post) cluster(sid)"
        global SPEC_ABSORB  "l_homicide [aw = popwt], absorb(sid year) second_stage(ib100.rel) treatment(post) cluster(sid)"
        cmpreal "1c event study, weighted"
    }
}


*==============================================================================
* 2. ACA MEDICAID EXPANSION  --  ehec_data.dta
*
* State-level panel of insurance coverage among low-income childless adults.
* `dins' is the coverage share, `yexp2' the year the state expanded Medicaid
* (missing for never-expanders), `stfips' the state identifier.
*
* Why this one: it is the worked example in the R package README, so the Stata
* and R answers are directly comparable, and it has a genuine never-treated
* group with staggered timing among the rest -- the canonical case the
* estimator is designed for.
*==============================================================================
display as txt _n "{hline 84}"
display as txt "2. ACA MEDICAID EXPANSION (ehec_data)"
display as txt "{hline 84}"

capture use "https://raw.githubusercontent.com/Mixtape-Sessions/Advanced-DID/main/Exercises/Data/ehec_data.dta", clear
if _rc {
    display as error "could not download ehec_data.dta -- skipping section 2"
}
else {
    gen byte treated = cond(missing(yexp2), 0, year >= yexp2)
    gen int  rel     = cond(missing(yexp2), 100, year - yexp2 + 20)

    qui levelsof stfips, local(nst)
    display as txt "  states: " wordcount("`nst'")
    qui count if missing(yexp2)
    display as txt "  never-expander observations: " r(N)

    global SPEC_CLASSIC "dins, first_stage(i.stfips i.year) second_stage(i.treated) treatment(treated) cluster(stfips)"
    global SPEC_ABSORB  "dins, absorb(stfips year) second_stage(i.treated) treatment(treated) cluster(stfips)"
    cmpreal "2a static"

    global SPEC_CLASSIC "dins, first_stage(i.stfips i.year) second_stage(ib100.rel) treatment(treated) cluster(stfips)"
    global SPEC_ABSORB  "dins, absorb(stfips year) second_stage(ib100.rel) treatment(treated) cluster(stfips)"
    cmpreal "2b event study"
}


*==============================================================================
* 3. NLS YOUNG WOMEN  --  webuse nlswork
*
* 28,534 person-year observations on 4,711 women, 1968-1988. Ships with Stata,
* so anyone reviewing the PR can run this without downloading anything.
*
* Why this one: it is the only dataset here with a LARGE fixed-effect
* dimension, which is the whole point of absorb(). It is also badly behaved in
* useful ways -- the year variable has gaps, the panel is very unbalanced, many
* women appear only once or twice, and several variables have missing values.
*
* The treatment is CONSTRUCTED (first year observed in a union, absorbing) and
* has no causal interpretation. It is here to produce realistic staggered
* timing on a large real panel, nothing more.
*
* Note the screen will fire: women whose first observed year already has
* union == 1 have no untreated observation, so their person effect is not
* identified. That is a real occurrence, not a defect, and is why the
* equivalence check below restricts to women with untreated support.
*==============================================================================
display as txt _n "{hline 84}"
display as txt "3. NLS YOUNG WOMEN (webuse nlswork) -- 4,711 person effects"
display as txt "{hline 84}"

capture webuse nlswork, clear
if _rc {
    display as error "could not load nlswork -- skipping section 3"
}
else {
    * absorbing union onset
    gen byte u1 = (union == 1) if !missing(union)
    bysort idcode (year): egen int onset = min(cond(u1 == 1, year, .))
    gen byte treat = cond(missing(onset), 0, year >= onset)
    gen int  rel   = cond(missing(onset), 100, year - onset + 20)

    * untreated support at the person level
    bysort idcode: egen byte has0 = max(treat == 0)

    qui count if has0 == 0
    display as txt "  observations on women with NO untreated year: " r(N)
    display as txt "  (the identification screen should drop exactly these)"

    display as txt _n "  -- full sample, absorb() only, default screen --"
    timer clear 3
    timer on 3
    qui did2s ln_wage, absorb(idcode year) second_stage(i.treat) ///
        treatment(treat) cluster(idcode)
    timer off 3
    qui timer list 3
    display as txt "     absorb: " as res %8.2f r(t3) as txt "s   e(N) = " e(N)

    *-- equivalence on a subsample small enough for the dummy path -------------
    * The classic path forms and inverts a k1 x k1 matrix with k1 ~ 4,700 on
    * the full sample. Set `dofull' to 1 only if you have time and memory.
    local dofull 0

    preserve
        keep if has0 == 1
        * random ~400 women, so both paths run in seconds
        set seed 20260729
        bysort idcode: gen byte first = (_n == 1)
        gen double r = runiform() if first
        bysort idcode (first): replace r = r[_N]
        keep if r <= 0.09

        qui levelsof idcode, local(ids)
        display as txt _n "  -- equivalence subsample: " wordcount("`ids'") " women --"

        global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year) second_stage(i.treat) treatment(treat) cluster(idcode)"
        global SPEC_ABSORB  "ln_wage, absorb(idcode year) second_stage(i.treat) treatment(treat) cluster(idcode)"
        cmpreal "3a static, person+year FE"

        global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year) second_stage(ib100.rel) treatment(treat) cluster(idcode)"
        global SPEC_ABSORB  "ln_wage, absorb(idcode year) second_stage(ib100.rel) treatment(treat) cluster(idcode)"
        cmpreal "3b event study"

        global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year hours) second_stage(i.treat) treatment(treat) cluster(idcode)"
        global SPEC_ABSORB  "ln_wage, first_stage(hours) absorb(idcode year) second_stage(i.treat) treatment(treat) cluster(idcode)"
        cmpreal "3c covariate with missing values"
    restore

    if (`dofull' == 1) {
        preserve
            keep if has0 == 1
            display as txt _n "  -- FULL sample equivalence (slow: k1 ~ 4,700) --"
            global SPEC_CLASSIC "ln_wage, first_stage(i.idcode i.year) second_stage(i.treat) treatment(treat) cluster(idcode)"
            global SPEC_ABSORB  "ln_wage, absorb(idcode year) second_stage(i.treat) treatment(treat) cluster(idcode)"
            cmpreal "3d full sample, static"
        restore
    }
}

display as txt _n "{hline 84}"
display as txt "All rows should read PASS. Timing on section 3 is the headline"
display as txt "number for the PR: it is the only one with a large FE dimension."
display as txt "{hline 84}"
