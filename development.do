use https://github.com/kylebutts/did2s_stata/raw/main/data/df_hom.dta, clear
gen rel_year_shift = rel_year + 20
replace rel_year_shift = 100 if rel_year_shift == .
gen popwt = runiform()

** Weights
local varlist dep_var
local first_stage i.year i.unit
local second_stage treat
local treatment treat
local cluster state
local weight = "pw"
local exp = "= popwt"
local vce "`r(vce)'"
tempvar touse
gen `touse' = 1 



/*
local varlist dep_var
local first_stage i.state i.year
local second_stage ib(last).rel_year_shift
local treatment treat
local cluster state
local weight = ""
local exp = ""
local vce "`r(vce)'"
tempvar touse
gen `touse' = 1
*/





/* use "https://github.com/scunning1975/mixtape/raw/master/castle.dta", clear

did2s l_homicide [pw=popwt], first_stage(i.year i.sid) second_stage(post) treatment(post) cluster(sid)
*/