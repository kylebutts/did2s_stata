*! version 0.7  reghdfe first stage with optional absorbed fixed effects

* Main function
capture program drop did2s
program define did2s, eclass sortpreserve
  *-> Setup

    version 13
    syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw /], second_stage(varlist fv) TREATment(varname) cluster(varname) [ ABSorb(varlist numeric) first_stage(varlist fv) TOLerance(real 1e-10) MAXITer(integer 2000) ]
    if ("`absorb'" == "" & "`first_stage'" == "") {
      display as error "one of {bf:first_stage()} or {bf:absorb()} is required"
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
      display as error "cluster variable {bf:`cluster'} is not a numeric variable."
      exit(198)
    }

    capture noisily ms_get_version reghdfe, min_version("5.7.3")
    if (_rc != 0) {
      display as error "{bf:did2s} requires {bf:reghdfe} version 5.7.3 or newer"
      display as error "install or update it with {bf:ssc install reghdfe, replace}"
      exit(199)
    }

    * Materialize the weight once for reghdfe, regress, and Mata.
    tempvar wvar
    if ("`weight'" != "") {
      qui gen double `wvar' = `exp'
      local weightexp "[`weight'=`wvar']"
    }
    else {
      quietly gen double `wvar' = 1
      local weightexp ""
    }

    * Expand each factor-variable specification once. The generated variables
    * are used for both sample marking and the matrix calculations below.
    local Qvars ""
    if ("`first_stage'" != "") {
      fvrevar `first_stage' if `touse'
      local Qvars "`r(varlist)'"
    }
    fvrevar `second_stage' if `touse'
    local full_second_stage "`r(varlist)'"

    markout `touse' `varlist' `treatment' `cluster'
    if ("`absorb'" != "") markout `touse' `absorb'
    markout `touse' `wvar'
    if ("`Qvars'" != "") markout `touse' `Qvars'
    markout `touse' `full_second_stage'

    quietly count if `touse' & !inlist(`treatment', 0, 1)
    if (r(N) > 0) {
      display as error "option {bf:treatment()} must identify treatment with only 0 and 1"
      exit(459)
    }

    *-> Stage 1. Save each absorbed effect separately so it can be propagated
    *   to treated observations sharing a level observed in the untreated
    *   estimation sample.
    local absorb_vars "`absorb'"
    if ("`absorb_vars'" == "") {
      * reghdfe still supplies the intercept; the user-requested factor
      * variables remain explicit regressors, preserving the original syntax.
      tempvar intercept_fe
      quietly gen byte `intercept_fe' = 1 if `touse'
      local absorb_vars "`intercept_fe'"
    }

    local absorb_spec ""
    local fehatvars ""
    foreach fe of varlist `absorb_vars' {
      tempvar fehat
      local absorb_spec "`absorb_spec' `fehat'=`fe'"
      local fehatvars "`fehatvars' `fehat'"
    }

    tempvar first_u
    quietly reghdfe `varlist' `Qvars' `weightexp' ///
      if `touse' & `treatment' == 0, absorb(`absorb_spec') ///
      noconstant keepsingletons residuals(`first_u') verbose(-1)
    quietly replace `first_u' = 0 if `touse' & `treatment' == 1
    quietly replace `first_u' = 0 if `touse' & missing(`first_u') & `wvar' == 0
    quietly count if `touse' & `treatment' == 0 & missing(`first_u')
    if (r(N) > 0) {
      display as error "reghdfe excluded untreated observations from the first stage"
      exit(2000)
    }

    tempvar xb y0hat adj
    if ("`Qvars'" != "") {
      quietly predict double `xb' if `touse', xb
    }
    else {
      quietly gen double `xb' = 0 if `touse'
    }
    _did2s_recover_fe, touse(`touse') absorb(`absorb_vars') ///
      estimates(`fehatvars') xb(`xb') generate(`y0hat')

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

    capture mata: mata drop __D2S
    mata: __D2S = _d2s_setup("`Qvars'", "`absorb_vars'", "`treatment'", "`touse'", "`wvar'", "`first_u'")

    *-> Second Stage
    quietly reg `adj' `second_stage' `weightexp' if `touse', ///
      nocons

      tempname b_second omit
      matrix `b_second' = e(b)
      _ms_omit_info `b_second'
      matrix `omit' = r(omit)

      local n_expanded : word count `full_second_stage'
      if (`n_expanded' != colsof(`b_second')) {
        display as error "internal error matching expanded second-stage variables"
        exit(503)
      }

      local i = 1
      local newlist
      foreach var of local full_second_stage {
        if `omit'[1,`i'] == 0 {
          local newlist "`newlist' `var'"
        }
        local ++i
      }
      local vars_second "`newlist'"

    tempvar second_u
    quietly predict double `second_u' if `touse', residual
    tempname b
    matrix `b' = e(b)
    local V_names: colnames `b'
    local N = e(N)
    tempvar finalsample
    quietly gen byte `finalsample' = e(sample)

    *-> SE correction
    mata: V = construct_V(__D2S, "`vars_second'", "`second_u'", "`cluster'", "`touse'", `tolerance', `maxiter')

    if (scalar(__d2s_rmax) > `tolerance') {
      display as txt "warning: variance-correction solver did not converge (max relative residual "
      display as txt "         " %8.1e scalar(__d2s_rmax) ", tolerance `tolerance')."
      display as txt "         DO NOT report these standard errors. Try {bf:maxiter()}, a coarser"
      display as txt "         fixed-effect specification, or a cluster bootstrap."
    }
    capture scalar drop __d2s_rmax

    *-> Export
    tempname V_final
    mata: st_matrix(st_local("V_final"), construct_V_final(V))

    matrix rownames `V_final' = `V_names'
    matrix colnames `V_final' = `V_names'

    ereturn clear
    ereturn post `b' `V_final', esample(`finalsample') obs(`N') depname(`varlist')
    ereturn local cmdline `"`0'"'
    ereturn local cmd "did2s"
    ereturn local clustvar "`cluster'"
    ereturn local absorb "`absorb'"

    ereturn display

    capture mata: mata drop __D2S
end

* Carry separately saved fixed effects to observations with the same level.
* This is the part of prediction that reghdfe's predict cannot do outside
* its untreated estimation sample.
capture program drop _did2s_recover_fe
program define _did2s_recover_fe, sortpreserve
    version 13
    syntax, TOUSE(varname) ABSorb(varlist numeric) ///
      ESTimates(varlist numeric) XB(varname numeric) GENerate(name)

    local n_absorb : word count `absorb'
    local n_estimates : word count `estimates'
    if (`n_absorb' != `n_estimates') exit(503)

    quietly gen double `generate' = `xb' if `touse'
    forvalues i = 1/`n_absorb' {
      local fe : word `i' of `absorb'
      local fehat : word `i' of `estimates'
      gsort `fe' -`fehat'
      quietly by `fe': replace `fehat' = `fehat'[1] if missing(`fehat') & `touse'
      quietly replace `generate' = `generate' + `fehat' if `touse'
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
mata: 

  struct _d2s_prob {
    real scalar     N, q, D, k1, rmax
    real colvector  w, w0, treat, first_u, dprec, ends
    real matrix     Q, LEV, P2
    real rowvector  Ld, has0, eoff, elen
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

    s.k1 = s.q + sum(s.Ld)

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
  real matrix construct_V(struct _d2s_prob scalar s, string scalar x2names, string scalar u2name, string scalar clname, string scalar tousename, real scalar tol, real scalar maxiter)
  {
    real matrix X2, first_stage_rhs, first_stage_coef, first_stage_fitted
    real matrix influence, cluster_influence, info, x2tx2_inv, cov
    real colvector u2, clv, p

    X2  = st_data(., tokens(x2names), tousename)
    u2  = st_data(., u2name, tousename)
    clv = st_data(., clname, tousename)
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
    p = order(clv, 1)
    info = panelsetup(clv[p], 1)
    cluster_influence = _d2s_grpsum(influence, p, info[., 2])
    cov = cross(cluster_influence, cluster_influence)

    st_numscalar("__d2s_rmax", s.rmax)
    return(cov)
  }

  matrix construct_V_final(numeric matrix V_adj){
    matrix V_final, omit, idx
    omit = st_matrix(st_local("omit"))
    V_final = J(length(omit), length(omit), 0)

    /* index of non-omitted variables */
    idx = select(1..length(omit), omit :== 0)
    V_final[idx, idx] = V_adj
    
    return(V_final)
  }
end 
