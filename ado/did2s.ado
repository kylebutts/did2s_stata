*! version 0.6

* Main function
capture program drop did2s
program define did2s, eclass
  *-> Setup

    version 13
    syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw /], first_stage(varlist fv) second_stage(varlist fv) TREATment(varname) cluster(varname) [ unit(varname) ] 

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
mata: 
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
