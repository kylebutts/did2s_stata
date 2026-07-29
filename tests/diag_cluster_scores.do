*! diag_cluster_scores.do -- rev 1
*! Why does the identification screen cut the clustered SE by ~20% while moving
*! the point estimate only ~1.9%?
*! Expects did2s.ado v0.7 and df_het.dta in the working directory.
*------------------------------------------------------------------------------
* PART 1 is the decisive test: run first_stage() on the SCREENED sample. If it
* reproduces absorb()'s screened result, the whole difference is attributable to
* the 310 dropped observations and the screen is provably equivalent to "drop
* these rows, then proceed exactly as before."
*
* PART 2 explains the magnitude: it rebuilds construct_V()'s score matrix and
* compares the per-cluster contributions to the meat for the ten flagged units
* against the other 990.
*
* Hypothesis under test: in the explicit-dummy path the flagged units get no
* unit dummy, because fvrevar runs on the untreated subsample. Their residual
* therefore carries the entire unit effect rather than just noise, so ten of a
* thousand clusters contribute scores an order of magnitude too large, and the
* meat matrix -- a sum of squared cluster scores -- is inflated accordingly.
*------------------------------------------------------------------------------

clear all
set more off
run "did2s.ado"

display as txt _n "{hline 78}"
display as txt "diag_cluster_scores.do  [rev 1]"
display as txt "{hline 78}"

use "df_het.dta", clear
qui replace treat = 1 if inrange(unit, 1, 10)
gen byte flagged = inrange(unit, 1, 10)

qui count if flagged
local nflag = r(N)
display as txt "  flagged observations: `nflag'   (10 always-treated units)"

*==============================================================================
* PART 1 -- is the screen equivalent to dropping those rows by hand?
*==============================================================================
display as txt _n "{hline 78}"
display as txt "PART 1: absorb() screened  vs  first_stage() on the screened sample"
display as txt "{hline 78}"

qui did2s dep_var, first_stage(i.unit i.year) second_stage(i.treat) ///
    treatment(treat) cluster(unit)
scalar b_full  = _b[1.treat]
scalar se_full = _se[1.treat]
scalar n_full  = e(N)

qui did2s dep_var if !flagged, first_stage(i.unit i.year) second_stage(i.treat) ///
    treatment(treat) cluster(unit)
scalar b_drop  = _b[1.treat]
scalar se_drop = _se[1.treat]
scalar n_drop  = e(N)

qui did2s dep_var, absorb(unit year) second_stage(i.treat) ///
    treatment(treat) cluster(unit)
scalar b_abs  = _b[1.treat]
scalar se_abs = _se[1.treat]
scalar n_abs  = e(N)

display as txt _n "                                            b            SE         N"
display as txt "  first_stage(), full sample     " as res %14.8f b_full %14.8f se_full %10.0f n_full
display as txt "  first_stage(), flagged dropped " as res %14.8f b_drop %14.8f se_drop %10.0f n_drop
display as txt "  absorb(), screened (default)   " as res %14.8f b_abs  %14.8f se_abs  %10.0f n_abs

display as txt _n "  absorb() vs first_stage()-on-screened-sample:"
display as txt "    rel diff in b  = " as res %10.2e abs(b_abs - b_drop)/abs(b_drop)
display as txt "    rel diff in SE = " as res %10.2e abs(se_abs - se_drop)/abs(se_drop)
display as txt ""
display as txt "  If those are at float noise, the screen changes nothing except which"
display as txt "  observations are used, and the entire SE movement is caused by the 310"
display as txt "  dropped rows rather than by the absorption algorithm."

*==============================================================================
* PART 2 -- where in the meat matrix does the difference live?
*
* Rebuilds the classic path by hand so the score matrix from construct_V() can
* be inspected per cluster.
*==============================================================================
display as txt _n "{hline 78}"
display as txt "PART 2: per-cluster score contributions under first_stage()"
display as txt "{hline 78}"

capture mata mata drop _diag_scores()
mata:
void _diag_scores(string scalar treatname, string scalar clname,
                  string scalar u1name,    string scalar u2name,
                  string scalar tousename, string scalar x1names,
                  string scalar x2names,   string scalar flagname)
{
    real colvector treat, clv, u1, u2, w, w0, p, flag, fs, A2, isflag, hi, lo
    real matrix X1, X2, V, S, A, info
    real scalar G, g, mh, ml

    treat = st_data(., treatname, tousename)
    clv   = st_data(., clname,    tousename)
    u1    = st_data(., u1name,    tousename)
    u2    = st_data(., u2name,    tousename)
    flag  = st_data(., flagname,  tousename)
    X1    = st_data(., tokens(x1names), tousename)
    X2    = st_data(., tokens(x2names), tousename)

    w  = J(rows(X1), 1, 1)
    w0 = w :* (1 :- treat)

    /* exactly construct_V()'s algebra */
    V = cross(X2, w, X1) * invsym(cross(X1, w0, X1))
    S = (X2 :* (w :* u2)) - (X1 :* (w :* u1)) * V'

    p    = order(clv, 1)
    S    = S[p, .]
    fs   = flag[p]
    info = panelsetup(clv[p], 1)
    G    = rows(info)

    A      = J(G, cols(S), 0)
    isflag = J(G, 1, 0)
    for (g = 1; g <= G; g++) {
        A[g, .]   = colsum(panelsubmatrix(S, g, info))
        isflag[g] = max(panelsubmatrix(fs, g, info))
    }

    /* contribution of cluster g to the meat, summed over coefficients */
    A2 = rowsum(A :^ 2)
    hi = select(A2, isflag :== 1)
    lo = select(A2, isflag :== 0)
    mh = mean(hi)
    ml = mean(lo)

    printf("\n  clusters: %g   flagged: %g   unflagged: %g\n", G, rows(hi), rows(lo))
    printf("  mean squared cluster score, flagged   = %12.4e\n", mh)
    printf("  mean squared cluster score, unflagged = %12.4e\n", ml)
    printf("  ratio                                 = %12.1f\n", mh / ml)
    printf("  flagged share of total meat           = %11.2f%%\n",
           100 * sum(hi) / sum(A2))
    printf("\n  %g of %g clusters (%4.1f%% of them) contributing %4.1f%% of the meat\n",
           rows(hi), G, 100 * rows(hi) / G, 100 * sum(hi) / sum(A2))
}
end

*-- replicate the classic path so the intermediates are available -------------
tempvar touse ones adj first_u second_u
qui gen byte `touse' = 1
qui gen `ones' = 1

fvrevar i.unit i.year if (`touse' == 1) & (treat == 0)
local full_first_stage "`r(varlist)' `ones'"

qui reg dep_var `full_first_stage' if (`touse' == 1) & (treat == 0), nocons
qui predict double `adj' if `touse', residual

tempname b_first omit
matrix `b_first' = e(b)
_ms_omit_info `b_first'
matrix `omit' = r(omit)
local vars_first : colnames e(b)
local i = 1
local newlist
foreach var in `vars_first' {
    if `omit'[1, `i'] == 0 {
        if ("`var'" != "_cons") local newlist "`newlist' `var'"
    }
    local ++i
}
local vars_first "`newlist'"

qui gen double `first_u' = `adj' * (1 - treat) if `touse'

fvrevar i.treat if `touse'
local full_second_stage `r(varlist)'
qui reg `adj' `full_second_stage' if `touse', nocons

tempname b_second omit2
matrix `b_second' = e(b)
_ms_omit_info `b_second'
matrix `omit2' = r(omit)
local vars_second : colnames e(b)
local i = 1
local newlist
foreach var in `vars_second' {
    if `omit2'[1, `i'] == 0 {
        if ("`var'" != "_cons") local newlist "`newlist' `var'"
    }
    local ++i
}
local vars_second "`newlist'"

qui predict double `second_u' if `touse', residual

*-- residual magnitudes, flagged vs not --------------------------------------
display as txt _n "  second-stage residual, by group:"
tabstat `second_u', by(flagged) stats(n mean sd min max) format(%12.4f)

display as txt _n "  first-stage residual carried into the score (first_u):"
display as txt "  (zero by construction on treated rows, so identically zero for"
display as txt "   the flagged units -- their entire score is the X2 u2 term)"
tabstat `first_u', by(flagged) stats(n mean sd) format(%12.4f)

*-- per-cluster meat contributions -------------------------------------------
mata: _diag_scores("treat", "unit", "`first_u'", "`second_u'", "`touse'", ///
                   "`vars_first'", "`vars_second'", "flagged")

display as txt _n "{hline 78}"
display as txt "HOW TO READ THIS"
display as txt "{hline 78}"
display as txt "PART 1 settles whether the screen is equivalent to dropping the rows"
display as txt "by hand. If it is, the SE change belongs to the sample, not the method."
display as txt ""
display as txt "PART 2 gives the magnitude. If 1% of clusters carry a double-digit"
display as txt "share of the meat, the inflated-cluster explanation is confirmed and"
display as txt "can be stated in the PR rather than proposed."
display as txt "{hline 78}"
