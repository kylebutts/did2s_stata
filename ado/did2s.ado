*! version 0.1


* local varlist dep_var
* local first_stage i.state i.year
* local treat_formula i.treat
* local treat_var treat


capture program drop did2s
program define did2s, eclass
    *-> Setup

        version 13
        syntax varlist(min=1 max=1 numeric) [if] [in], [first_stage(varlist fv) treat_formula(varlist fv) treat_var(varname)]

        * to use
        tempvar touse
        mark `touse' `if' `in'

    *-> First Stage 

        fvrevar `first_stage'
        local full_first_stage `r(varlist)'

        * First stage regression
        qui reg `varlist' `full_first_stage' if `touse' & `treat_var' == 0, nocons robust

        * Store reg results
        tempname b_first V_first noomit omit
        matrix `b_first' = e(b)
        matrix `V_first' = e(V)

        * Residualize outcome variable
        tempvar adj
        predict double `adj', residual

        **-> Keeping only vcov of non-omitted
        * https://www.stata.com/support/faqs/programming/factor-variable-support/
            
            _ms_omit_info `b_first'
            local cols_first = colsof(`b_first')
            matrix `noomit' =  J(1,`cols_first',1) - r(omit)
            matrix `omit' = r(omit)
            
            * Store V_first in mata
            mata: V_first = select(st_matrix(st_local("V_first")),(st_matrix(st_local("noomit"))))
            mata: V_first = select(V_first, (st_matrix(st_local("noomit")))')

        **-> Get names of non-omitted variables
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

    *-> Second Stage

        fvrevar `treat_formula'
        local full_second_stage `r(varlist)'

        * Second stage regression
        qui reg `adj' `full_second_stage' if `touse', nocons robust

        * Store reg results
        tempname b_second V_second noomit omit
        matrix `b_second' = e(b)
        matrix `V_second' = e(V)

        *-> Keeping only vcov of non-omitted
        * https://www.stata.com/support/faqs/programming/factor-variable-support/
            
            _ms_omit_info `b_second'
            local cols_second = colsof(`b_second')
            matrix `noomit' =  J(1,`cols_second',1) - r(omit)
            matrix `omit' = r(omit)
            
            * Store V_second in mata
            mata: V_second = select(st_matrix(st_local("V_second")), (st_matrix(st_local("noomit"))))
            mata: V_second = select(V_second, (st_matrix(st_local("noomit")))')

            * Get non-omitted second stage variables
            mata: b_second = select(st_matrix(st_local("b_second")),(st_matrix(st_local("noomit"))))
            mata: st_matrix(st_local("b_second"), b_second)

        * Get names of non-omitted variables
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
            * correct for loop
            local n_non_omit_second: word count `vars_second'
            * disp "`vars_second'"

    *-> Standard Error Adjustment
        
        * Create initialized matrix
        mata: coef = J(0, `n_non_omit_second', .)

        * loop through first_stage variables
        foreach v in `vars_first' {
            qui reg `v' `vars_second', nocons
            mata: coef = coef \ st_matrix("e(b)")
        }

        * V_second + gmm adjustment
        mata: V_adj = V_second + coef' * V_first * coef


    *-> Export
        tempname b V_final

        * Second stage regression (with pretty display)
        qui reg `adj' `treat_formula' if `touse', nocons robust depname(`varlist')
        matrix `b' = e(b)
        local V_names: rownames e(V)
        local N = e(N)
        * local r2 = e(r2)
        * local r2_a = e(r2_a)
        * local F = e(F)

        * Fill in V for omitted variables
        mata: st_matrix(st_local("V_final"), construct_V_final(V_adj)) 

        matrix rownames `V_final' = `V_names'
        matrix colnames `V_final' = `V_names'


        ereturn clear
        ereturn post `b' `V_final', esample(`touse')
        ereturn local cmdline `"`0'"'
        ereturn local cmd "did2s"
        ereturn scalar N = `N'
        * ereturn scalar r2 = `r2'
        * ereturn scalar r2_a = `r2_a'
        * ereturn scalar F = `F'

        ereturn display
end

capture program drop Display
program Display
        version 9
        local version : di "version " string(_caller()) ":"
        syntax [, Level(cilevel) noHEader * ]
        _get_diopts diopts options , `options'
        if "`e(prefix)'" != "" {
                _prefix_display, level(`level') `header' `diopts' `options'
        }
        else {
                `version' _regress, level(`level') `header' `diopts' `options'
                _prefix_footnote
        }
end


version 13

capture mata mata drop construct_V_final()
mata:
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
