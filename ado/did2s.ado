*! version 0.7  first stage solved by preconditioned conjugate gradients

* Main function
capture program drop did2s
program define did2s, eclass
  *-> Setup

    version 13
    syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw /], second_stage(varlist fv) TREATment(varname) cluster(varname) [ first_stage(varlist fv) ABSorb(varlist numeric) unit(varname) TOLerance(real 1e-12) MAXITer(integer 2000) ]

    if ("`first_stage'" == "" & "`absorb'" == "" & "`unit'" == "") {
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

    * confirm cluster is a numeric variable
    capture confirm numeric variable `cluster'
    if (_rc != 0) {
      display as error "cluster variable {bf:`cluster'} is not a numeric variable."
      exit(198)
    }

    * to use
    tempvar touse
    mark `touse' `if' `in'

    * Mata reads the weight with st_data(), which requires a variable name, so
    * the weight expression is materialized once and used everywhere.
    if ("`weight'" != "") {
      tempvar wvar
      qui gen double `wvar' = `exp' if `touse'
      local weightexp "`weight'=`wvar'"
    }
    else {
      local wvar ""
      local weightexp ""
    }

  *-> Fixed effects to absorb
    * Every first-stage fixed effect is absorbed rather than entered as dummies.
    * Plain factor main effects in first_stage() (i.x, ib#.x, ibn.x, ...) are
    * moved to the absorbed list, and unit(x) is the same as absorb(x). All
    * other first-stage terms -- continuous covariates, interactions, slopes,
    * level-restricted factors such as i(1/3).x -- stay as regressors.
    local fe_vars "`absorb'"
    local fe_show "`absorb'"
    local fs_rest ""
    if ("`first_stage'" != "") {
      fvunab fs_terms : `first_stage'
      local rest `"`fs_terms'"'
      while (`"`rest'"' != "") {
        gettoken term rest : rest, bind
        if (`"`term'"' == "") continue
        local route = 0
        local dot = strpos(`"`term'"', ".")
        if (`dot' > 1) {
          local pre = substr(`"`term'"', 1, `dot' - 1)
          local nm  = substr(`"`term'"', `dot' + 1, .)
          if ((`"`pre'"' == "i" | substr(`"`pre'"', 1, 2) == "ib") & strpos(`"`nm'"', ".") == 0 & strpos(`"`nm'"', "#") == 0) {
            capture confirm numeric variable `nm', exact
            if (_rc == 0) local route = 1
          }
        }
        if (`route') {
          local fe_vars "`fe_vars' `nm'"
          local fe_show "`fe_show' `nm'"
        }
        else {
          local fs_rest `"`fs_rest' `term'"'
        }
      }
      local fs_rest : list retokenize fs_rest
    }
    if ("`unit'" != "") {
      capture confirm numeric variable `unit'
      if (_rc != 0) {
        tempvar unitid
        qui egen long `unitid' = group(`unit') if `touse'
        local fe_vars "`fe_vars' `unitid'"
      }
      else {
        local fe_vars "`fe_vars' `unit'"
      }
      local fe_show "`fe_show' `unit'"
    }
    local fe_vars : list uniq fe_vars
    local fe_show : list uniq fe_show

    * With no fixed effect at all, absorb a constant, as the dummy first stage
    * of earlier versions always included one.
    if ("`fe_vars'" == "") {
      tempvar cons
      qui gen byte `cons' = 1
      local fe_vars "`cons'"
    }

  *-> Estimation sample
    markout `touse' `varlist' `treatment' `cluster' `fe_vars' `wvar'
    if (`"`fs_rest'"' != "") {
      fvrevar `fs_rest', list
      markout `touse' `r(varlist)'
    }
    fvrevar `second_stage', list
    markout `touse' `r(varlist)'

    capture assert inlist(`treatment', 0, 1) if `touse'
    if (_rc != 0) {
      display as error "treatment variable {bf:`treatment'} must be 0/1"
      exit(450)
    }

    qui count if `touse'
    if (r(N) == 0) error 2000

    local Qvars ""
    if (`"`fs_rest'"' != "") {
      fvrevar `fs_rest' if `touse'
      local Qvars "`r(varlist)'"
    }

  *-> First Stage
    * _d2s_setup() also screens the sample: observations in a fixed-effect
    * level with no untreated observations have no identified counterfactual,
    * and are flagged in `noid' and dropped. Dropping them leaves every other
    * level's untreated support unchanged, so one pass suffices.
    tempvar adj noid
    qui gen double `adj' = .
    qui gen byte `noid' = 0
    capture mata: mata drop __D2S
    mata: __D2S = _d2s_setup("`varlist'", "`Qvars'", "`fe_vars'", "`treatment'", "`touse'", "`wvar'", `tolerance', `maxiter', "`adj'", "`noid'")

    qui count if `touse' & `noid'
    local ndrop = r(N)
    if (`ndrop' > 0) {
      qui replace `touse' = 0 if `noid'
      display as txt "note: `ndrop' observation(s) dropped: they belong to a fixed-effect level with no"
      display as txt "      untreated observations, so their counterfactual Y(0) is not identified."
    }

    qui count if `touse'
    if (r(N) == 0) {
      capture mata: mata drop __D2S
      capture scalar drop __d2s_rmax
      display as error "no observations left: no fixed-effect level has untreated observations"
      exit(2000)
    }

  *-> Second Stage
    fvrevar `second_stage' if `touse'
    local full_second_stage `r(varlist)'

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
          if ("`var'" != "_cons") local newlist "`newlist' `var'"
        }
        local ++i
      }
      local vars_second "`newlist'"

    tempvar second_u
    quietly predict double `second_u' if `touse', residual

  *-> Standard Error Adjustment
    mata: V = _d2s_hdfe_V(__D2S, "`vars_second'", "`second_u'", "`cluster'", "`touse'", `tolerance', `maxiter')
    capture mata: mata drop __D2S

    if (scalar(__d2s_rmax) > `tolerance') {
      display as txt "warning: first-stage solver did not converge (max relative residual "
      display as txt "         " %8.1e scalar(__d2s_rmax) ", tolerance `tolerance')."
      display as txt "         This usually means some treated observation's counterfactual is not"
      display as txt "         identified from the untreated data, e.g. a first-stage covariate"
      display as txt "         that does not vary among untreated observations, or fixed effects"
      display as txt "         whose untreated observations fall into disconnected groups."
      display as txt "         DO NOT report these standard errors. Try {bf:maxiter()}, a coarser"
      display as txt "         first-stage specification, or a cluster bootstrap."
    }
    capture scalar drop __d2s_rmax

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
    ereturn local absorb "`fe_show'"
    ereturn scalar N = `N'
    ereturn scalar N_noid = `ndrop'

    ereturn display
end

* Mata Functions
version 13
capture mata mata drop construct_V_final()
capture mata mata drop _d2s_prob()
capture mata mata drop _d2s_grpsum()
capture mata mata drop _d2s_lev_vec()
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

  /* Level index for one FE: 1..L for levels with untreated support, 0
     otherwise. p2 receives a permutation sorted by that index (unsupported
     rows first), built by a stable partition of a single sort. */
  real colvector _d2s_lev_vec(real colvector f, real colvector s0, real colvector p2)
  {
    real colvector p, fp, newg, gidrow, ends, cs, sw, supp, gid, lev, z
    real scalar n, Lall, L

    n = rows(f)
    p = order(f, 1)
    fp = f[p]
    if (n > 1) newg = 1 \ (fp[|2 \ n|] :!= fp[|1 \ n - 1|])
    else       newg = 1
    gidrow = runningsum(newg)
    Lall = gidrow[n]
    if (n > 1) ends = selectindex(newg[|2 \ n|] \ 1)
    else       ends = 1
    cs = quadrunningsum(s0[p])
    if (Lall == 1) sw = cs[ends]
    else           sw = cs[ends] - (0 \ cs[ends[|1 \ Lall - 1|]])
    supp = sw :> 0
    L = sum(supp)
    gid = J(Lall, 1, 0)
    if (L > 0) gid[selectindex(supp)] = (1::L)
    lev = J(n, 1, 0)
    lev[p] = gid[gidrow, 1]
    z = (lev[p] :== 0)
    p2 = select(p, z) \ select(p, !z)
    return(lev)
  }

  struct _d2s_prob scalar _d2s_setup(string scalar yname, string scalar qnames, string scalar fenames, string scalar treatname, string scalar tousename, string scalar wname, real scalar tol, real scalar maxiter, string scalar adjname, string scalar flagname)
  {
    struct _d2s_prob scalar s
    string rowvector fes
    real colvector treat, w, s0, keep, kidx, map, lev, p, c1, beta, u, ends, gg, dq, adjout, y
    real matrix LV, PP, inf, Q
    real scalar d, n0

    treat = st_data(., treatname, tousename)
    n0 = rows(treat)
    if (wname != "") {
      w = st_data(., wname, tousename)
    }
    else {
      w = J(n0, 1, 1)
    }
    /* Support indicator: untreated AND positively weighted. Counting rows
       rather than summing weights avoids assuming non-negative weights. */
    s0 = (treat :== 0) :* (w :> 0)

    /* Level indexes, one sort per FE, on the full estimation sample */
    fes = tokens(fenames)
    s.D = cols(fes)
    LV = J(n0, s.D, .)
    PP = J(n0, s.D, .)
    for (d = 1; d <= s.D; d++) {
      lev = _d2s_lev_vec(st_data(., fes[d], tousename), s0, p)
      if (max(lev) == 0) {
        _error(2000, "did2s: a fixed effect has no levels with untreated observations")
      }
      LV[., d] = lev
      PP[., d] = p
    }

    /* Identification screen: a row is identified iff its level has untreated
       support in every FE. Untreated rows with positive weight always are, so
       every supported level keeps its rows and its index. */
    keep = (rowmin(LV) :> 0)
    st_store(., flagname, tousename, !keep)
    kidx = selectindex(keep)
    if (rows(kidx) == 0) {
      _error(2000, "did2s: no observations have identified first-stage predictions")
    }
    if (rows(kidx) < n0) {
      map = runningsum(keep)
      LV = LV[kidx, .]
      for (d = 1; d <= s.D; d++) {
        p = PP[., d]
        p = select(p, keep[p])
        PP[|1, d \ rows(kidx), d|] = map[p]
      }
      PP = PP[|1, 1 \ rows(kidx), s.D|]
    }

    y = st_data(., yname, tousename)
    s.y = y[kidx]
    s.N = rows(kidx)
    s.treat = treat[kidx]
    s.w = w[kidx]
    s.w0 = s.w :* (1 :- s.treat)

    if (qnames != "") {
      Q = st_data(., tokens(qnames), tousename)
      s.Q = Q[kidx, .]
      s.q = cols(s.Q)
    }
    else {
      s.Q = J(s.N, 0, 0)
      s.q = 0
    }

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
      lev = LV[., d]
      p = PP[., d]
      s.Ld[d] = max(lev)
      s.LEV[., d] = lev + (lev :== 0) :* (s.Ld[d] + 1)
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
    adjout = J(n0, 1, .)
    adjout[kidx] = u
    st_store(., adjname, tousename, adjout)
    s.first_u = u :* (1 :- s.treat)

    st_numscalar("__d2s_rmax", s.rmax)
    return(s)
  }

  /* Corrected VCV: Gardner's variance formula with (X10'WX10)^- applied
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
