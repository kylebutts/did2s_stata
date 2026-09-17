{smcl}
{it:v. 0.7}
{viewerjumpto "Anticipation" "did2s##anticipation"}{...}
{viewerjumpto "Absorbing fixed effects" "did2s##absorb"}{...}

{title:Two-Stage Difference-in-Differences}

{pstd}
{bf:did2s} - Estimates a TWFE model using the two-stage difference-in-differences approach from {browse "https://jrgcmu.github.io/2sdd_current.pdf":Gardner (2021)}

{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd:did2s} {depvar} {ifin} [{it:{help regress##weight:weight}}]{cmd:,} {cmdab:absorb(}{help varlist}{cmd:)} {cmdab:second_stage(}{help varlist}{cmd:)} {cmdab:treatment(}{help varname}{cmd:)} {cmdab:cluster(}{help varname}{cmd:)} [{cmdab:first_stage(}{help varlist}{cmd:)} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth absorb(varlist)}}Categorical variables whose fixed effects are absorbed with {help reghdfe}. Unit and time effects normally belong here.{p_end}
{synopt :{opth first_stage(varlist)}}Optional first-stage covariates other than the fixed effects listed in {opt absorb()}. Factor-variable notation is allowed.{p_end}
{synopt :{opth second_stage(varlist)}}List of treatment variables. This could be, for example a 0/1 treatment dummy, a set of event-study leads/lags, or a continuous treatment variable.{p_end}
{synopt :{opth treatment(varname)}}This must be a 0/1 dummy for when treatment is occuring; D_it = 1 if unit i is treated during period t. See {help did2s##anticipation:Anticipation} for details on how to deal with anticipation (shift D_it the number of periods of anticipation).{p_end}
{synopt :{opth cluster(varname)}}What variable to cluster on (use unit id if you don't want to cluster).{p_end}
{synopt :{opt tol:erance(#)}}Relative residual tolerance for the matrix-free variance correction. Default {cmd:tolerance(1e-10)}.{p_end}
{synopt :{opt maxit:er(#)}}Maximum iterations for each variance-correction solve. Default {cmd:maxiter(2000)}.{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{bf:did2s} implements Two-Staged Difference-in-Differences by Gardner (2021). A TWFE model for outcomes is given by unit/group fixed effects, time fixed effects, treatment variable (or variables in the case of event study), and potentially covariates. To avoid the problems with OLS estimation of difference-in-differences/event-studies in the presence of staggered treatment adoption, this method proceeds in two stages:{p_end}

{phang2}
1. This program estimates the unit/group fixed effects, time fixed effects, and potentially covariates using only untreated/not-yet-treated observations. This is used to predict counterfactual outcomes in all periods and residualize the observed outcome.{p_end}

{phang2}
2. Then regress the residualized outcome on the treatment variable(s) to estimate the treatment effects.{p_end}

{pstd}
More details can be found in {browse "https://jrgcmu.github.io/2sdd_current.pdf":Gardner (2021)} or more informally in {browse "https://kylebutts.com/blog/posts/2021-05-24-two-stage-difference-in-differences/":this blog post}. 
{p_end}


{marker anticipation}
{title:Anticipation}

{pstd}
This procedure works so long as μ_i and μ_t are consistently estimated. The key is to use only untreated/not-yet-treated observations to estimate the fixed effects. For example, if you used observations with D_it=1, you would attribute treatment effects τ as “fixed characteristics” and would combine μ_i with the treatment effects.
{p_end}

{pstd}
The fixed effects could be biased/inconsistent if there are anticipation effects, i.e. units respond before treatment starts. The fix is fairly simple, simply “shift” treatment date earlier by as many years as you suspect anticipation to occur (e.g. 2 years before treatment starts) and estimate on the subsample where the shifted treatment equals zero. If you suspect there is anticipation, pass this modified variable to {bf: treat_var}.
{p_end}

{marker absorb}
{title:Absorbing fixed effects}

{pstd}
{opt absorb()} requires {help reghdfe}. Install it with {cmd:ssc install reghdfe}.
The first stage is estimated on untreated observations with {cmd:reghdfe}. Saved
fixed-effect estimates are then propagated to treated observations with the
same fixed-effect levels. This reproduces the role of absorbed fixed effects in
the R implementation without materializing a dummy for every level.
{p_end}

{pstd}
Put ordinary categorical fixed effects in {opt absorb()} and other controls in
{opt first_stage()}. For example, R's
{cmd:first_stage = ~ x1 + x2 | unit + year} corresponds to
{cmd:first_stage(x1 x2) absorb(unit year)}.
{p_end}

{pstd}
The corrected variance calculation applies the first-stage cross-product as a
linear operator and does not construct the full fixed-effect dummy matrix.
{opt tolerance()} and {opt maxiter()} control this iterative solve.
{p_end}

{marker untreated}
{title:No Never-treated Group}

{pstd}
With {opt absorb()}, observations are dropped when any absorbed fixed-effect
level is absent from the untreated estimation sample. Their counterfactual
outcomes are not identified. {cmd:e(N)} and {cmd:e(sample)} reflect the reduced
sample.
{p_end}


{marker warning}
{title:Important Warning}

When some absorbed fixed-effect level is never observed untreated, its
counterfactual cannot be imputed. {cmd:did2s} leaves that prediction missing and
drops the affected observations with a note.

{title:Examples:}

{pstd}Setup{p_end}
{phang2}{cmd:. use https://github.com/kylebutts/did2s_stata/raw/main/data/df_hom.dta, clear}{p_end}

{pstd}Static Model{p_end}
{phang2}{cmd:. did2s dep_var, absorb(state year) second_stage(i.treat) treatment(treat) cluster(state)}{p_end}

{pstd}Event Study Model{p_end}
{phang2}{cmd:. gen rel_year_shift = rel_year + 20}{p_end}
{phang2}{cmd:. replace rel_year_shift = 100 if rel_year_shift == .}{p_end}
{phang2}{cmd:. did2s dep_var, absorb(state year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)}{p_end}

{pstd}With Covariates{p_end}
{phang2}{cmd:. use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear}{p_end}
{phang2}{cmd:. global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], first_stage($demo) absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}




{marker results}{...}
{title:Stored results}

{pstd}
{cmd:did2s} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:did2s}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed-effect variables{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}



{title:Program Author}

{pstd}
Kyle Butts   {break}
University of Colorado, Boulder      {break}
buttskyle96@gmail.com     {break}


{title:References}

{phang}
Gardner (2021), {browse "https://jrgcmu.github.io/2sdd_current.pdf":Two-stage Difference-in-differences}
{p_end}



