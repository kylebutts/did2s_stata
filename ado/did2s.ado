*! version 0.1

capture program drop did2s
program define did2s, eclass
    *-> Setup

        version 13
        syntax varlist(min=1 max=1 numeric) [if] [in] [aw fw iw pw], first_stage(varlist fv) treat_formula(varlist fv) treat_var(varname) cluster(varname)

        * to use
        tempvar touse
        mark `touse' `if' `in'

        * confirm cluster is a numeric variable
        capture confirm numeric variable `cluster'

        if _rc != 0 { 
            display as error "cluster variable {bf:`clustervar'} is not a numeric variable."
            exit(198)
        }
        

    *-> First Stage 

        fvrevar `first_stage'
        local full_first_stage `r(varlist)'

        * First stage regression (with clustering and weights)
        qui reg `varlist' `full_first_stage' [`weight'`exp'] if `touse' & `treat_var' == 0, nocons vce(cluster `cluster')

        * Residualize outcome variable
        tempvar adj
        predict double `adj', residual

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
        gen `first_u' = `adj' * (1 - `treat_var')

    *-> Second Stage

        fvrevar `treat_formula'
        local full_second_stage `r(varlist)'

        * Second stage regression
        qui reg `adj' `full_second_stage' [`weight'`exp'] if `touse', nocons vce(cluster `cluster')

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

        **-> Create first_u, with 0s in row where D_it = 1
        tempvar second_u
        predict double `second_u', residual 
            

    *-> Standard Error Adjustment
        
        * Create initialized matrix
        mata: V = construct_V("`treat_var'", "`cluster'", "`first_u'", "`second_u'", "`touse'", "`vars_first'", "`vars_second'", `n_non_omit_second')

    *-> Export
        tempname b V_final

        * Second stage regression (with pretty display)
        qui reg `adj' `treat_formula' [`weight'`exp'] if `touse', nocons robust depname(`varlist')
        matrix `b' = e(b)
        local V_names: rownames e(V)
        local N = e(N)
        * local r2 = e(r2)
        * local r2_a = e(r2_a)
        * local F = e(F)

        * Fill in V for omitted variables
        mata: st_matrix(st_local("V_final"), construct_V_final(V)) 

        matrix rownames `V_final' = `V_names'
        matrix colnames `V_final' = `V_names'


        ereturn clear
        ereturn post `b' `V_final', esample(`touse')
        ereturn local cmdline `"`0'"'
        ereturn local cmd "did2s"
        ereturn local  vce      "`vce'"
        ereturn local  vcetype  "`vcetype'"
        ereturn local  clustvar "`cluster'"
        ereturn scalar N = `N'
        * ereturn scalar r2 = `r2'
        * ereturn scalar r2_a = `r2_a'
        * ereturn scalar F = `F'

        ereturn display
end


version 13

capture mata mata drop construct_V()
capture mata mata drop construct_V_final()
mata: 
    matrix construct_V(string scalar treat_var_str, string scalar cluster_str, string scalar first_u_str, string scalar second_u_str, string scalar touse_str, string scalar vars_first_str, string scalar vars_second_str, real scalar n2) {
        real colvector treat, cluster_var, first_u, second_u, cl, idx
        real matrix X1, X2, X10, V, meat, W, cov

        st_view(treat = ., ., treat_var_str, touse_str)
        st_view(cluster_var = ., ., cluster_str, touse_str)
        st_view(first_u = ., ., first_u_str, touse_str)
        st_view(second_u = ., ., second_u_str, touse_str)
        
        st_view(X1 = ., ., vars_first_str, touse_str)
        st_view(X2 = ., ., vars_second_str, touse_str)

        /* For Testing
        st_view(treat = ., ., "`treat_var'", "`touse'")
        st_view(cluster_var = ., ., "`cluster'", "`touse'")
        st_view(first_u = 0, ., "`first_u'", "`touse'")
        st_view(second_u = 0, ., "`second_u'", "`touse'")
        st_view(X1 = ., ., "`vars_first'", "`touse'")
        st_view(X2 = ., ., "`vars_second'", "`touse'")
        n2 = `n_non_omit_second'
        */
        
        /* Create X10 */
        X10 = X1
        for(i=1; i <= rows(X1); i++) {
            if(treat[i] == 1) {
                X10[i,] = X10[i,] :* 0
            }
        }

        /* Only calculate this part once */
        V = X2' * X1 * invsym(X10' * X10)

        cl = uniqrows(cluster_var)

        /* Initialize meat */
        meat = J(n2, n2, 0)

        /* real colvector temp_first_u, temp_second_u */
        for(i=1; i <= length(cl); i++) {
            idx = cluster_var :== cl[i]
        
            W = select(X2, idx)' * select(second_u , idx) - V * select(X10, idx)' * select(first_u, idx)

            meat = meat + W * W'
        }

        cov = invsym(X2'*X2) * meat * invsym(X2'*X2)

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
