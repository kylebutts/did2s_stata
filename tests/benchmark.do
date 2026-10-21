*! benchmark.do -- rev 2 (2026-10-21): warm-up call before timing
*------------------------------------------------------------------------------
* did2s: v0.7 vs v0.6 timing on the same first_stage(i.unit i.year) command.
*
* Simulated staggered panel, 20 years, half the units never treated, static
* and event-study second stages, clustered on unit. v0.6 is skipped above
* 4,000 units, where its k1 x k1 dummy matrix gets slow and large.
*
* SETUP: run from the repository root or tests/. A few minutes.
*------------------------------------------------------------------------------

clear all
set more off

capture confirm file "ado/did2s.ado"
if (_rc == 0) global ROOT "."
else global ROOT ".."
run "$ROOT/ado/did2s.ado"
run "$ROOT/tests/reference/did2s_ref.ado"

capture program drop mkpanel
program define mkpanel
    args nunits
    clear
    set seed 12345
    qui set obs `nunits'
    gen long unit = _n
    gen double ufe = rnormal()
    gen int g = cond(runiform() < 0.5, 0, 2005 + floor(runiform() * 12))
    qui expand 20
    bysort unit: gen int year = 1999 + _n
    gen byte treat = (g > 0) & (year >= g)
    gen int rel = cond(g > 0, year - g + 20, 100)
    qui replace rel = 10 if rel < 10
    qui replace rel = 30 if rel > 30 & rel < 100
    gen double y = ufe + 0.05 * (year - 2000) + 2 * treat + rnormal()
end

capture program drop timeit
program define timeit, rclass
    timer clear 9
    timer on 9
    qui `0'
    timer off 9
    qui timer list 9
    return scalar t = r(t9)
end

* Warm-up: the first call in a session pays one-off costs (loading and
* compiling), which would otherwise land on the 500-unit row.
mkpanel 200
qui did2s     y, first_stage(i.unit i.year) second_stage(ib100.rel) treatment(treat) cluster(unit)
qui did2s_ref y, first_stage(i.unit i.year) second_stage(ib100.rel) treatment(treat) cluster(unit)

display as txt _n "{hline 76}"
display as txt "  units        N   second stage        v0.6 (s)    v0.7 (s)    speedup"
display as txt "{hline 76}"

foreach n in 500 1000 2000 4000 8000 16000 {
    mkpanel `n'
    foreach ss in "i.treat" "ib100.rel" {
        local cmd "y, first_stage(i.unit i.year) second_stage(`ss') treatment(treat) cluster(unit)"

        timeit did2s `cmd'
        local t7 = r(t)

        if (`n' <= 4000) {
            timeit did2s_ref `cmd'
            local t6 = r(t)
            local sp = `t6' / `t7'
            display as txt %7.0fc `n' %10.0fc _N "   " %-16s "`ss'" ///
                as res %10.2f `t6' %12.2f `t7' %10.1f `sp' "x"
        }
        else {
            display as txt %7.0fc `n' %10.0fc _N "   " %-16s "`ss'" ///
                as txt %10s "--" as res %12.2f `t7'
        }
    }
}
display as txt "{hline 76}"
