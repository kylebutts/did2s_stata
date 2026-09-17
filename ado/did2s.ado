*! version 0.7  absorb(): first stage with absorbed fixed effects

* Main function
capture program drop did2s
program define did2s, eclass sortpreserve
  *-> Setup

    version 13
    syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw /], second_stage(varlist fv) TREATment(varname) cluster(varname) [ first_stage(varlist fv) absorb(varlist numeric) unit(varname) TOLerance(real 1e-10) MAXITer(integer 2000) ]

    if ("`absorb'" == "" & "`first_stage'" == "") {
      display as error "option {bf:first_stage()} required (or use {bf:absorb()})"
      exit(198)
    }
    if ("`absorb'" != "" & "`unit'" != "") {
      display as error "options {bf:absorb()} and {bf:unit()} may not be combined; put the unit fixed effect in {bf:absorb()}"
      exit(198)
    }
    if (`tolerance' <= 0) {
      display as error "option {bf:tolerance()} must be positive"
      exit(198)
    }
    if (`maxiter' <= 0) {
      display as error "option {bf:maxiter()} must be a positive integer"
      exit(198)
    }

    * to use
    tempvar touse
    mark `touse' `if' `in'

    * confirm cluster is a numeric variable
    capture confirm numeric variable `cluster'

    if _rc != 0 { 
      display as error "cluster variable {bf:`clustervar'} is not a numeric variable."
      exit(198)
    }

    if("`weight'" == "") {
      local weightexp = ""
    } 
    else {
      local weightexp "`weight'=`exp'"
    }

  *============================================================================
  * PATH A: absorb().  Path B below is unchanged.
  *============================================================================
  if ("`absorb'" != "") {

    capture which reghdfe
    if (_rc != 0) {
      display as error "option {bf:absorb()} requires {bf:reghdfe}"
      display as error "install it with {bf:ssc install reghdfe}"
      exit(199)
    }

    * Mata reads the weight with st_data(), which requires a variable name,
    * so an expression such as [aw = pop*2] must be materialized first.
    if ("`weight'" != "") {
      tempvar wvar
      qui gen double `wvar' = `exp'
      local hdfeweightexp "`weight'=`wvar'"
    }
    else {
      local wvar = ""
      local hdfeweightexp = ""
    }

    markout `touse' `varlist' `treatment' `cluster' `absorb'
    if ("`wvar'" != "") markout `touse' `wvar'
    if ("`first_stage'" != "") {
      fvrevar `first_stage', list
      markout `touse' `r(varlist)'
    }
    fvrevar `second_stage', list
    markout `touse' `r(varlist)'

    quietly count if `touse' & !inlist(`treatment', 0, 1)
    if (r(N) > 0) {
      display as error "option {bf:treatment()} must identify treatment with only 0 and 1"
      exit(459)
    }

    local Qvars ""
    if ("`first_stage'" != "") {
      fvrevar `first_stage' if `touse'
      local Qvars "`r(varlist)'"
    }

    fvrevar `second_stage' if `touse'
    local full_second_stage `r(varlist)'

    *-> Stage 1. Save each absorbed effect separately so it can be propagated
    *   to treated observations sharing a level observed in the untreated
    *   estimation sample.
    local absorb_spec ""
    local fehatvars ""
    foreach fe of varlist `absorb' {
      tempvar fehat
      local absorb_spec "`absorb_spec' `fehat'=`fe'"
      local fehatvars "`fehatvars' `fehat'"
    }

    quietly reghdfe `varlist' `first_stage' [`hdfeweightexp'] ///
      if `touse' & `treatment' == 0, absorb(`absorb_spec') ///
      noconstant keepsingletons verbose(-1)

    tempvar xb y0hat adj
    if ("`first_stage'" != "") {
      quietly predict double `xb' if `touse', xb
    }
    else {
      quietly gen double `xb' = 0 if `touse'
    }
    quietly gen double `y0hat' = `xb' if `touse'

    local i = 1
    foreach fe of varlist `absorb' {
      local fehat : word `i' of `fehatvars'
      gsort `fe' -`fehat'
      quietly by `fe': replace `fehat' = `fehat'[1] if missing(`fehat') & `touse'
      quietly replace `y0hat' = `y0hat' + `fehat' if `touse'
      local ++i
    }

    quietly gen double `adj' = `varlist' - `y0hat' if `touse'
    quietly count if `touse' & missing(`adj')
    local ndrop = r(N)
    if (`ndrop' > 0) {
      quietly replace `touse' = 0 if missing(`adj')
      display as txt "note: `ndrop' observation(s) dropped because an absorbed fixed-effect"
      display as txt "      level was not observed in the untreated estimation sample."
    }
    quietly count if `touse'
    if (r(N) == 0) {
      display as error "no observations have identified first-stage predictions"
      exit(2000)
    }

    tempvar first_u
    quietly gen double `first_u' = `adj' * (1 - `treatment') if `touse'

    capture mata: mata drop __D2S
    mata: __D2S = _d2s_setup("`Qvars'", "`absorb'", "`treatment'", "`touse'", "`wvar'", "`first_u'")

    *-> Second Stage
    quietly reg `adj' `full_second_stage' [`weightexp'] if `touse', nocons

      tempname b_second omit
      matrix `b_second' = e(b)
      _ms_omit_info `b_second'
      matrix `omit' = r(omit)
      local vars_second: colnames e(b)

      local i = 1
      local newlist
      foreach var in `vars_second' {
        if `omit'[1,`i'] == 0 {
          if ("`var'" != "_cons") local newlist "`newlist' `var'"
        }
        local ++i
      }
      local vars_second "`newlist'"

    tempvar second_u
    quietly predict double `second_u' if `touse', residual

    *-> SE correction
    mata: V = _d2s_hdfe_V(__D2S, "`vars_second'", "`second_u'", "`cluster'", "`touse'", `tolerance', `maxiter')

    if (scalar(__d2s_rmax) > `tolerance') {
      display as txt "warning: variance-correction solver did not converge (max relative residual "
      display as txt "         " %8.1e scalar(__d2s_rmax) ", tolerance `tolerance')."
      display as txt "         DO NOT report these standard errors. Try {bf:maxiter()}, a coarser"
      display as txt "         fixed-effect specification, or a cluster bootstrap."
    }
    capture scalar drop __d2s_rmax

    *-> Export
    quietly reg `adj' `second_stage' [`weightexp'] if `touse', nocons robust depname(`varlist')
    tempname b
    matrix `b' = e(b)
    local V_names: rownames e(V)
    local N = e(N)
    tempvar finalsample
    quietly gen byte `finalsample' = e(sample)

    tempname V_final
    mata: st_matrix(st_local("V_final"), construct_V_final(V))

    matrix rownames `V_final' = `V_names'
    matrix colnames `V_final' = `V_names'

    ereturn clear
    ereturn post `b' `V_final', esample(`finalsample')
    ereturn local cmdline `"`0'"'
    ereturn local cmd "did2s"
    ereturn local clustvar "`cluster'"
    ereturn local absorb "`absorb'"
    ereturn scalar N = `N'

    ereturn display

    capture mata: mata drop __D2S
    exit
  }


  *============================================================================
  * PATH B: classic (explicit dummies, optional unit() demeaning) -- unchanged
  *============================================================================

    if("`unit'" != "") {
      qui preserve
    }

  *-> First Stage 

    * Using `nocons` and adding a constant manually (for i.state i.year)
    tempvar ones
    gen `ones' = 1

    fvrevar `first_stage' if (`touse' == 1) & (`treatment' == 0)
    local full_first_stage "`r(varlist)' `ones'"
    
    * Manually demean all the first_stage vars by id (when treat == 0).
    * The first stage is estimated by WLS when weights are supplied, so the
    * demeaning must use the WEIGHTED mean over untreated observations within
    * unit; otherwise Frisch-Waugh fails and unit() stops reproducing the
    * equivalent full-dummy specification. Reduces to the simple mean when
    * weights are absent or constant within unit. (Issue #19)
    if ("`unit'" != "") {
      tempvar dmwt mean wsum xsum
      if ("`weight'" == "") {
        quietly gen double `dmwt' = 1
      }
      else {
        quietly gen double `dmwt' = `exp'
      }
      foreach var of varlist `full_first_stage' `varlist' {
        cap drop `wsum'
        cap drop `xsum'
        cap drop `mean'
        * weighted mean, ignoring rows where `var' is missing (as egen's
        * mean() does)
        quietly egen double `wsum' = total(cond(`touse' & `treatment' == 0 ///
            & !missing(`var'), `dmwt', .)), by(`unit')
        quietly egen double `xsum' = total(cond(`touse' & `treatment' == 0, ///
            `dmwt' * `var', .)), by(`unit')
        quietly gen double `mean' = `xsum' / `wsum'
        quietly replace `var' = `var' - `mean'
      }
    }

    * First stage regression (with clustering and weights)
    quietly reg `varlist' `full_first_stage' [`weightexp'] if (`touse' == 1) & (`treatment' == 0), nocons

    
    * Residualize outcome variable
    tempvar adj
    quietly predict double `adj' if `touse', residual

    **-> Get names of non-omitted variables
      * https://www.stata.com/support/faqs/programming/factor-variable-support/
      
      tempname b_first omit
      matrix `b_first' = e(b)
      _ms_omit_info `b_first'
      matrix `omit' = r(omit)
      local vars_first: colnames e(b)
      
      local i = 1
      local newlist
      foreach var in `vars_first' {
        if `omit'[1,`i'] == 0 {
          if ("`var'" != "_cons") {
            local newlist "`newlist' `var'"
          }
        }
        local ++i
      }

      local vars_first "`newlist'"
      * disp "`vars_first'"

    **-> Create first_u, with 0s in row where D_it = 1
    tempvar first_u
    quietly gen `first_u' = `adj' * (1 - `treatment') if `touse'

  *-> Second Stage

    fvrevar `second_stage' if `touse'
    local full_second_stage `r(varlist)'

    * Second stage regression
    quietly reg `adj' `full_second_stage' [`weightexp'] if `touse', nocons

    **-> Get names of non-omitted variables
      * https://www.stata.com/support/faqs/programming/factor-variable-support/
      
      tempname b_second omit 
      matrix `b_second' = e(b)
      _ms_omit_info `b_second'
      matrix `omit' = r(omit)
      local vars_second: colnames e(b)
      
      local i = 1
      local newlist
      foreach var in `vars_second' {
        if `omit'[1,`i'] == 0 {
          if ("`var'" != "_cons") {
            local newlist "`newlist' `var'"
          }
        }
        local ++i
      }
      local vars_second "`newlist'"
      * disp "`vars_second'"

      * get number of 2nd stage variables 
      local n_non_omit_second: word count `vars_second'

    **-> Create second_u
    tempvar second_u
    quietly predict double `second_u' if `touse', residual 
      

  *-> Standard Error Adjustment

    * Keep only esample for second_stage
    tempvar touse_x_esample
    qui gen `touse_x_esample' = (`touse' == 1) & (e(sample) == 1)
    mata: V = construct_V("`treatment'", "`cluster'", "`first_u'", "`second_u'", "`touse_x_esample'", "`vars_first'", "`vars_second'", "`exp'", `n_non_omit_second')

  *-> Export

    * Second stage regression (with pretty display)
    quietly reg `adj' `second_stage' [`weightexp'] if `touse', nocons robust depname(`varlist')
    tempname b
    matrix `b' = e(b)
    local V_names: rownames e(V)
    local N = e(N)

    * Fill in V for omitted variables
    tempname V_final
    mata: st_matrix(st_local("V_final"), construct_V_final(V)) 

    matrix rownames `V_final' = `V_names'
    matrix colnames `V_final' = `V_names'

    ereturn clear
    ereturn post `b' `V_final', esample(`touse')
    ereturn local cmdline `"`0'"'
    ereturn local cmd "did2s"
    ereturn local clustvar "`cluster'"
    ereturn scalar N = `N'

    ereturn display

    if("`unit'" != "") {
      qui restore
    }
end

* Mata Functions
version 13
capture mata mata drop construct_V()
capture mata mata drop construct_V_final()
capture mata mata drop _d2s_prob()
capture mata mata drop _d2s_grpsum()
capture mata mata drop _d2s_Xv()
capture mata mata drop _d2s_XtWv()
capture mata mata drop _d2s_M0v()
capture mata mata drop _d2s_pcg()
capture mata mata drop _d2s_setup()
capture mata mata drop _d2s_hdfe_V()
mata: 

  struct _d2s_prob {
    real scalar     N, q, D, k1, kfe, rmax
    real colvector  w, w0, treat, first_u, dprec, ends
    real matrix     Q, LEV, P2
    real rowvector  Ld, has0, eoff, elen, feoff
  }

  /* Group sums of the columns of x, given a sort order p and the group end
     positions (in sorted order). */
  real matrix _d2s_grpsum(real matrix x, real colvector p, real colvector ends)
  {
    real matrix out
    real colvector cs
    real scalar G, j

    G = rows(ends)
    out = J(G, cols(x), 0)
    for (j = 1; j <= cols(x); j++) {
      cs = quadrunningsum(x[p, j])
      if (G == 1) {
        out[., j] = cs[ends]
      }
      else {
        out[., j] = cs[ends] - (0 \ cs[ends[|1 \ G - 1|]])
      }
    }
    return(out)
  }

  /* X1 * v */
  real matrix _d2s_Xv(struct _d2s_prob scalar s, real matrix v)
  {
    real matrix xv, vd
    real scalar d, off, m

    m = cols(v)
    if (s.q > 0) {
      xv = s.Q * v[|1, 1 \ s.q, m|]
    }
    else {
      xv = J(s.N, m, 0)
    }
    off = s.q
    for (d = 1; d <= s.D; d++) {
      vd = v[|off + 1, 1 \ off + s.Ld[d], m|] \ J(1, m, 0)
      xv = xv + vd[s.LEV[., d], .]
      off = off + s.Ld[d]
    }
    return(xv)
  }

  /* X1' * t */
  real matrix _d2s_XtWv(struct _d2s_prob scalar s, real matrix t)
  {
    real matrix out, g
    real colvector ends
    real scalar d

    if (s.q > 0) {
      out = cross(s.Q, t)
    }
    else {
      out = J(0, cols(t), 0)
    }
    for (d = 1; d <= s.D; d++) {
      ends = s.ends[|s.eoff[d] + 1 \ s.eoff[d] + s.elen[d]|]
      g = _d2s_grpsum(t, s.P2[., d], ends)
      if (s.has0[d] == 1) {
        g = g[|2, 1 \ rows(g), cols(g)|]
      }
      out = out \ g
    }
    return(out)
  }

  /* (X10' W X10) v */
  real matrix _d2s_M0v(struct _d2s_prob scalar s, real matrix v)
  {
    return(_d2s_XtWv(s, s.w0 :* _d2s_Xv(s, v)))
  }

  /* Jacobi-preconditioned conjugate gradients, one solve per column of C.
     Any two solutions differ only in directions annihilated by X10, so they
     produce the same fitted-value correction. A divergence guard stops
     runaway iterations on inconsistent systems. */
  real matrix _d2s_pcg(struct _d2s_prob scalar s, real matrix C, real scalar tol, real scalar maxiter)
  {
    real matrix B
    real colvector c, b, r, z, pp, Ap
    real scalar j, it, normc, rz, rznew, alpha, beta, relres, pAp

    B = J(s.k1, cols(C), 0)
    for (j = 1; j <= cols(C); j++) {
      c = C[., j]
      normc = sqrt(cross(c, c))
      if (normc <= 0) {
        continue
      }
      b = J(s.k1, 1, 0)
      r = c
      z = r :/ s.dprec
      pp = z
      rz = cross(r, z)
      relres = 1
      for (it = 1; it <= maxiter; it++) {
        Ap = _d2s_M0v(s, pp)
        pAp = cross(pp, Ap)
        if (pAp <= 0) {
          break
        }
        alpha = rz / pAp
        b = b + alpha * pp
        r = r - alpha * Ap
        relres = sqrt(cross(r, r)) / normc
        if (relres < tol) {
          break
        }
        if (relres > 1e6) {
          break
        }
        z = r :/ s.dprec
        rznew = cross(r, z)
        beta = rznew / rz
        pp = z + beta * pp
        rz = rznew
      }
      if (relres > s.rmax) {
        s.rmax = relres
      }
      B[., j] = b
    }
    return(B)
  }

  struct _d2s_prob scalar _d2s_setup(string scalar qnames, string scalar fenames, string scalar treatname, string scalar tousename, string scalar wname, string scalar firstuname)
  {
    struct _d2s_prob scalar s
    string rowvector fes
    real colvector f, p, s0, sw0g, supp, gid, lev, ends, gg, dq
    real matrix inf
    real scalar d, g, Lall, L

    s.treat = st_data(., treatname, tousename)
    s.first_u = st_data(., firstuname, tousename)
    s.N = rows(s.treat)
    if (wname != "") {
      s.w = st_data(., wname, tousename)
    }
    else {
      s.w = J(s.N, 1, 1)
    }
    s.w0 = s.w :* (1 :- s.treat)

    if (qnames != "") {
      s.Q = st_data(., tokens(qnames), tousename)
      s.q = cols(s.Q)
    }
    else {
      s.Q = J(s.N, 0, 0)
      s.q = 0
    }

    /* s.w0 is the actual first-stage weight vector. The support indicator is
       used only to enumerate levels represented in the untreated sample. */
    s0 = (s.treat :== 0) :* (s.w :> 0)

    fes = tokens(fenames)
    s.D = cols(fes)
    s.LEV  = J(s.N, s.D, .)
    s.P2   = J(s.N, s.D, .)
    s.Ld   = J(1, s.D, 0)
    s.has0 = J(1, s.D, 0)
    s.eoff = J(1, s.D, 0)
    s.elen = J(1, s.D, 0)
    s.feoff = J(1, s.D, 0)
    s.ends = J(0, 1, 0)
    s.rmax = 0

    for (d = 1; d <= s.D; d++) {
      f = st_data(., fes[d], tousename)
      p = order(f, 1)
      inf = panelsetup(f[p], 1)
      Lall = rows(inf)
      sw0g = _d2s_grpsum(s0, p, inf[., 2])
      supp = sw0g :> 0
      L = sum(supp)
      if (L == 0) {
        _error("did2s: a fixed effect has no levels with untreated observations")
      }
      gid = J(Lall, 1, 0)
      gid[selectindex(supp)] = (1::L)
      lev = J(s.N, 1, 0)
      for (g = 1; g <= Lall; g++) {
        lev[p[|inf[g, 1] \ inf[g, 2]|]] = J(inf[g, 2] - inf[g, 1] + 1, 1, gid[g])
      }
      s.Ld[d] = L
      s.LEV[., d] = lev + (lev :== 0) :* (L + 1)
      p = order(lev, 1)
      s.P2[., d] = p
      inf = panelsetup(lev[p], 1)
      s.has0[d] = (lev[p[1]] == 0)
      s.eoff[d] = rows(s.ends)
      s.elen[d] = rows(inf)
      s.ends = s.ends \ inf[., 2]
    }

    s.feoff[1] = 0
    for (d = 2; d <= s.D; d++) {
      s.feoff[d] = s.feoff[d - 1] + s.Ld[d - 1]
    }
    s.kfe = s.feoff[s.D] + s.Ld[s.D]
    s.k1  = s.q + s.kfe

    /* Jacobi preconditioner: diag(X10' W X10) */
    s.dprec = J(0, 1, 0)
    if (s.q > 0) {
      dq = colsum(s.w0 :* (s.Q :^ 2))'
      s.dprec = dq
    }
    for (d = 1; d <= s.D; d++) {
      ends = s.ends[|s.eoff[d] + 1 \ s.eoff[d] + s.elen[d]|]
      gg = _d2s_grpsum(s.w0, s.P2[., d], ends)
      if (s.has0[d] == 1) {
        gg = gg[|2, 1 \ rows(gg), 1|]
      }
      s.dprec = s.dprec \ gg
    }
    s.dprec = s.dprec + (s.dprec :<= 0)

    return(s)
  }

  /* Corrected VCV, following the R implementation with the first-stage
     inverse applied iteratively:

       B = (X10' W X10)^- X1' W X2
       Z = X1 B
       IF_i = [X2_i w_i u2_i - Z_i w_i u1_i] (X2' W X2)^-1
       V = sum_g (sum_{i in g} IF_i)' (sum_{i in g} IF_i)                 */
  real matrix _d2s_hdfe_V(struct _d2s_prob scalar s, string scalar x2names, string scalar u2name, string scalar clname, string scalar tousename, real scalar tol, real scalar maxiter)
  {
    real matrix X2, first_stage_rhs, first_stage_coef, first_stage_fitted
    real matrix influence, cluster_influence, info, x2tx2_inv, cov
    real colvector u2, clv, p
    real scalar k2, G, g

    X2  = st_data(., tokens(x2names), tousename)
    u2  = st_data(., u2name, tousename)
    clv = st_data(., clname, tousename)
    k2  = cols(X2)

    /* Map estimation noise from the untreated first stage into each
       second-stage coefficient without constructing X10' W X10. */
    first_stage_rhs = _d2s_XtWv(s, s.w :* X2)
    first_stage_coef = _d2s_pcg(s, first_stage_rhs, tol, maxiter)
    first_stage_fitted = _d2s_Xv(s, first_stage_coef)

    /* This follows the R implementation directly. Each row is the
       observation-level influence function:

         (X2' W X2)^-1 [X2_i w_i u2_i - Z_i w_i u1_i]

       written transposed so observations are rows. */
    x2tx2_inv = invsym(cross(X2, s.w, X2))
    influence = ((X2 :* (s.w :* u2)) -
      (first_stage_fitted :* (s.w :* s.first_u))) * x2tx2_inv

    /* Sum influence functions within cluster, then take their cross-product. */
    p         = order(clv, 1)
    influence = influence[p, .]
    info = panelsetup(clv[p], 1)
    G    = rows(info)
    cluster_influence = J(G, k2, 0)
    for (g = 1; g <= G; g++) {
      cluster_influence[g, .] = colsum(panelsubmatrix(influence, g, info))
    }
    cov = cross(cluster_influence, cluster_influence)

    st_numscalar("__d2s_rmax", s.rmax)
    return(cov)
  }

  /* ---- classic path, unchanged ---- */

  matrix construct_V(string scalar treatment_str, string scalar cluster_str, string scalar first_u_str, string scalar second_u_str, string scalar touse_str, string scalar vars_first_str, string scalar vars_second_str, string scalar weights_str, real scalar k2) {

    real colvector treat, cluster_var, first_u, second_u, weights, w0, p
    real matrix X1, X2, first_stage_map, first_stage_score
    real matrix second_stage_score, influence, cluster_influence, info
    real matrix x2tx2_inv, cov
    real scalar G, g

    st_view(treat = ., ., treatment_str, touse_str)
    st_view(cluster_var = ., ., cluster_str, touse_str)
    st_view(first_u = ., ., first_u_str, touse_str)
    st_view(second_u = ., ., second_u_str, touse_str)

    st_view(X1 = ., ., vars_first_str, touse_str)
    st_view(X2 = ., ., vars_second_str, touse_str)

    if(weights_str != "") {
      st_view(weights = ., ., weights_str, touse_str)
    }
    else {
      weights = J(rows(X1), 1, 1)
    }

    /* (X10' W0 X10) is cross(X1, .) with the weight zeroed on treated rows,
       so no row-selected copy of X1 is needed. */
    w0 = weights :* (1 :- treat)

    /* Map estimation noise from the untreated first stage into each
       second-stage coefficient. */
    first_stage_map = cross(X2, weights, X1) *
      invsym(cross(X1, w0, X1))

    /* Follow the R implementation directly. The two score terms are N x k2;
       multiplying their difference by the second-stage bread gives the
       observation-level influence function with observations in rows. */
    second_stage_score = X2 :* (weights :* second_u)
    first_stage_score = (X1 :* (weights :* first_u)) * first_stage_map'
    x2tx2_inv = invsym(cross(X2, weights, X2))
    influence = (second_stage_score - first_stage_score) * x2tx2_inv

    /* The cluster-robust covariance is the cross-product of influence
       functions after summing them within cluster. */
    p         = order(cluster_var, 1)
    influence = influence[p, .]
    info = panelsetup(cluster_var[p], 1)
    G    = rows(info)

    cluster_influence = J(G, k2, 0)
    for(g = 1; g <= G; g++) {
      cluster_influence[g, .] = colsum(panelsubmatrix(influence, g, info))
    }

    cov = cross(cluster_influence, cluster_influence)

    return(cov)
  }

  matrix construct_V_final(numeric matrix V_adj){
    matrix V_final, omit, idx
    scalar i, j
 
    omit = st_matrix(st_local("omit"))
    V_final = J(length(omit), length(omit), 0)

    /* index of non-omitted variables */
    idx = select(1..length(omit), omit :== 0)

    for(i = 1; i <= length(idx); i++) {
      for(j = 1; j <= i; j++) {
        V_final[idx[i], idx[j]] = V_adj[i,j]
        V_final[idx[j], idx[i]] = V_adj[i,j]
      }
    }
    
    return(V_final)
  }
end 
