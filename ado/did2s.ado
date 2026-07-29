*! version 0.7  absorb(): first stage with absorbed fixed effects

* Main function
capture program drop did2s
program define did2s, eclass
  *-> Setup

    version 13
    syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw /], second_stage(varlist fv) TREATment(varname) cluster(varname) [ first_stage(varlist fv) absorb(varlist numeric) unit(varname) TOLerance(real 1e-12) MAXITer(integer 2000) KEEPNOid ]

    if ("`absorb'" == "" & "`first_stage'" == "") {
      display as error "option {bf:first_stage()} required (or use {bf:absorb()})"
      exit(198)
    }
    if ("`absorb'" != "" & "`unit'" != "") {
      display as error "options {bf:absorb()} and {bf:unit()} may not be combined; put the unit fixed effect in {bf:absorb()}"
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

    * Mata reads the weight with st_data(), which requires a variable name,
    * so an expression such as [aw = pop*2] must be materialized first.
    if ("`weight'" != "") {
      tempvar wvar
      qui gen double `wvar' = `exp'
    }
    else {
      local wvar = ""
    }

    markout `touse' `varlist' `treatment' `cluster' `absorb'
    if ("`wvar'" != "") markout `touse' `wvar'
    if ("`first_stage'" != "") {
      fvrevar `first_stage', list
      markout `touse' `r(varlist)'
    }
    fvrevar `second_stage', list
    markout `touse' `r(varlist)'

    *-> Identification screen: FE levels with no untreated support.
    * One pass suffices. Every row flagged here has zero untreated weight (it
    * is treated, or untreated with non-positive weight), and a level is kept
    * iff it contains at least one untreated row with positive weight, so
    * removing flagged rows leaves every other level's untreated support --- in
    * every dimension --- exactly as it was. There is nothing for a second pass
    * to find.
    tempvar noid
    qui gen byte `noid' = 0
    mata: st_numscalar("__d2s_nid", _d2s_flag_nosupp("`absorb'", "`treatment'", "`touse'", "`wvar'", "`noid'"))

    local ndrop = 0
    if (scalar(__d2s_nid) > 0 & "`keepnoid'" == "") {
      local ndrop = scalar(__d2s_nid)
      qui replace `touse' = 0 if `noid'
    }

    if (`ndrop' > 0) {
      display as txt "note: `ndrop' observation(s) dropped: their fixed-effect level had no untreated"
      display as txt "      observations, so the first-stage counterfactual is not identified there."
      display as txt "      (Option {bf:keepnoid} keeps them; the estimate is then normalization-dependent.)"
    }
    if (scalar(__d2s_nid) > 0 & "`keepnoid'" != "") {
      display as txt "warning: " scalar(__d2s_nid) " observation(s) have a fixed-effect level with no untreated"
      display as txt "         observations. {bf:keepnoid} keeps them, which sets those levels' effects to"
      display as txt "         zero. Point estimates are normalization-dependent and do NOT match the"
      display as txt "         equivalent {bf:first_stage()} specification. The variance calculation is"
      display as txt "         inconsistent under this option: DO NOT REPORT THESE STANDARD ERRORS."
    }
    capture scalar drop __d2s_nid

    qui count if `touse'
    if (r(N) == 0) {
      display as error "no observations left after the identification screen"
      exit(2000)
    }

    local Qvars ""
    if ("`first_stage'" != "") {
      fvrevar `first_stage' if `touse'
      local Qvars "`r(varlist)'"
    }

    fvrevar `second_stage' if `touse'
    local full_second_stage `r(varlist)'

    *-> Stage 1
    tempvar adj
    qui gen double `adj' = .
    capture mata: mata drop __D2S
    mata: __D2S = _d2s_setup("`varlist'", "`Qvars'", "`absorb'", "`treatment'", "`touse'", "`wvar'", `tolerance', `maxiter', "`adj'")

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
      display as txt "warning: first-stage solver did not converge (max relative residual "
      display as txt "         " %8.1e scalar(__d2s_rmax) ", tolerance `tolerance')."
      display as txt "         This usually means the first-stage system is inconsistent, i.e. some"
      display as txt "         observations' fixed effects are not identified from untreated data."
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

    tempname V_final
    mata: st_matrix(st_local("V_final"), construct_V_final(V))

    matrix rownames `V_final' = `V_names'
    matrix colnames `V_final' = `V_names'

    ereturn clear
    ereturn post `b' `V_final', esample(`touse')
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
capture mata mata drop _d2s_flag_nosupp()
capture mata mata drop _d2s_Xv()
capture mata mata drop _d2s_XtWv()
capture mata mata drop _d2s_M0v()
capture mata mata drop _d2s_pcg()
capture mata mata drop _d2s_setup()
capture mata mata drop _d2s_hdfe_V()
mata: 

  struct _d2s_prob {
    real scalar     N, q, D, k1, kfe, rmax
    real colvector  y, w, w0, treat, first_u, dprec, ends, beta1
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

  /* Flag observations whose fixed-effect level has no untreated support in
     some dimension. Returns the count. */
  real scalar _d2s_flag_nosupp(string scalar fenames, string scalar treatname, string scalar tousename, string scalar wname, string scalar flagname)
  {
    string rowvector fes
    real colvector f, p, w, w0, treat, sw0g, supp, bad
    real matrix inf
    real scalar d, D, N, g, Lall

    fes = tokens(fenames)
    D = cols(fes)
    treat = st_data(., treatname, tousename)
    N = rows(treat)
    if (wname != "") {
      w = st_data(., wname, tousename)
    }
    else {
      w = J(N, 1, 1)
    }
    /* Support indicator: untreated AND positively weighted. Counting rows
       rather than summing weights avoids assuming weights are non-negative,
       which iweights do not guarantee. Identical to sum(w0) > 0 whenever
       they are. */
    w0 = (treat :== 0) :* (w :> 0)
    bad = J(N, 1, 0)

    for (d = 1; d <= D; d++) {
      f = st_data(., fes[d], tousename)
      p = order(f, 1)
      inf = panelsetup(f[p], 1)
      Lall = rows(inf)
      sw0g = _d2s_grpsum(w0, p, inf[., 2])
      supp = sw0g :> 0
      for (g = 1; g <= Lall; g++) {
        if (supp[g] == 0) {
          bad[p[|inf[g, 1] \ inf[g, 2]|]] = J(inf[g, 2] - inf[g, 1] + 1, 1, 1)
        }
      }
    }
    st_store(., flagname, tousename, bad)
    return(sum(bad))
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
     Starting from zero keeps the iterates in the range of the operator, so a
     consistent singular system converges to the minimum-norm solution. A
     divergence guard stops runaway iterations on inconsistent systems. */
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

  struct _d2s_prob scalar _d2s_setup(string scalar yname, string scalar qnames, string scalar fenames, string scalar treatname, string scalar tousename, string scalar wname, real scalar tol, real scalar maxiter, string scalar adjname)
  {
    struct _d2s_prob scalar s
    string rowvector fes
    real colvector f, p, s0, sw0g, supp, gid, lev, c1, beta, u, ends, gg, dq
    real matrix inf
    real scalar d, g, Lall, L

    s.y = st_data(., yname, tousename)
    s.N = rows(s.y)
    s.treat = st_data(., treatname, tousename)
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

    /* Support indicator, matching _d2s_flag_nosupp(). s.w0 remains the
       actual weight vector used in the normal equations. */
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

    /* Stage 1 */
    c1 = _d2s_XtWv(s, s.w0 :* s.y)
    beta = _d2s_pcg(s, c1, tol, maxiter)
    s.beta1 = beta

    u = s.y - _d2s_Xv(s, beta)
    st_store(., adjname, tousename, u)
    s.first_u = u :* (1 :- s.treat)

    st_numscalar("__d2s_rmax", s.rmax)
    return(s)
  }

  /* Corrected VCV: the construct_V() formula with (X10'WX10)^- applied
     iteratively.
       B = (X10'WX10)^- X1'W X2 ; Z = X1 B
       S = X2 .* (w u2) - Z .* (w u1)
       V = (X2'WX2)^-1 [sum_g S_g' S_g] (X2'WX2)^-1                       */
  real matrix _d2s_hdfe_V(struct _d2s_prob scalar s, string scalar x2names, string scalar u2name, string scalar clname, string scalar tousename, real scalar tol, real scalar maxiter)
  {
    real matrix X2, C, B, Z, S, A, info, meat, bread
    real colvector u2, clv, p
    real scalar k2, G, g

    X2  = st_data(., tokens(x2names), tousename)
    u2  = st_data(., u2name, tousename)
    clv = st_data(., clname, tousename)
    k2  = cols(X2)

    C = _d2s_XtWv(s, s.w :* X2)
    B = _d2s_pcg(s, C, tol, maxiter)
    Z = _d2s_Xv(s, B)

    S = (X2 :* (s.w :* u2)) - (Z :* (s.w :* s.first_u))

    p    = order(clv, 1)
    S    = S[p, .]
    info = panelsetup(clv[p], 1)
    G    = rows(info)
    A = J(G, k2, 0)
    for (g = 1; g <= G; g++) {
      A[g, .] = colsum(panelsubmatrix(S, g, info))
    }

    meat  = cross(A, A)
    bread = invsym(cross(X2, s.w, X2))

    st_numscalar("__d2s_rmax", s.rmax)
    return(bread * meat * bread)
  }

  /* ---- classic path, unchanged ---- */

  matrix construct_V(string scalar treatment_str, string scalar cluster_str, string scalar first_u_str, string scalar second_u_str, string scalar touse_str, string scalar vars_first_str, string scalar vars_second_str, string scalar weights_str, real scalar k2) {

    real colvector treat, cluster_var, first_u, second_u, weights, w0, p
    real matrix X1, X2, V, S, A, info, meat, cov, invX2X2
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

    /* Only calculate this part once */
    V = cross(X2, weights, X1) * invsym(cross(X1, w0, X1))

    /* Observation-level scores (N x k2):
         s_i = w_i X2_i u2_i - V w_i X1_i u1_i
       first_u is already zero on treated rows, so restricting the second
       term to untreated rows is automatic. Summing s_i within a cluster
       reproduces the W vector of the original loop exactly. */
    S = (X2 :* (weights :* second_u)) - (X1 :* (weights :* first_u)) * V'

    /* Sum scores within cluster: sort once, then one pass over the panels,
       instead of scanning all N rows once per cluster. */
    p    = order(cluster_var, 1)
    S    = S[p, .]
    info = panelsetup(cluster_var[p], 1)
    G    = rows(info)

    A = J(G, k2, 0)
    for(g = 1; g <= G; g++) {
      A[g, .] = colsum(panelsubmatrix(S, g, info))
    }

    meat = cross(A, A)

    invX2X2 = invsym(cross(X2, weights, X2))
    cov = invX2X2 * meat * invX2X2

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
