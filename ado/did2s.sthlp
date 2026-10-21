{smcl}
{it:v. 0.7} 
{viewerjumpto "Anticipation" "did2s##anticipation"}{...}
{viewerjumpto "Fixed effects" "did2s##absorb"}{...}
{viewerjumpto "Changes in version 0.7" "did2s##changes"}{...}

{title:Two-Stage Difference-in-Differences}

{pstd}
{bf:did2s} - Estimates a TWFE model using the two-stage difference-in-differences approach from {browse "https://jrgcmu.github.io/2sdd_current.pdf":Gardner (2021)}

{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd:did2s} {depvar} {ifin} [{it:{help regress##weight:weight}}]{cmd:,} {cmdab:second_stage(}{help varlist}{cmd:)} {cmdab:treatment(}{help varname}{cmd:)} {cmdab:cluster(}{help varname}{cmd:)} [{cmdab:absorb(}{help varlist}{cmd:)} {cmdab:first_stage(}{help varlist}{cmd:)} {it:options}]

{pstd}
One of {opt first_stage()} or {opt absorb()} is required. They may be combined.
{cmd:first_stage(i.unit i.year)} and {cmd:absorb(unit year)} are the same model and
give identical results; see {help did2s##absorb:Fixed effects}.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth first_stage(varlist)}}Model for the counterfactual Y_it(0), estimated on untreated observations. Factor-variable notation is allowed. Fixed effects written as plain factor main effects ({cmd:i.unit}, {cmd:ib2010.year}, {cmd:ibn.firm}, ...) are absorbed automatically; everything else (continuous covariates, interactions such as {cmd:i.state#c.year}, level-restricted factors such as {cmd:i(1/3).x}) enters as regressors. See {help did2s##absorb:Fixed effects}.{p_end}
{synopt :{opth absorb(varlist)}}Numeric categorical variables whose fixed effects are absorbed. {cmd:absorb(unit year)} is the same as {cmd:first_stage(i.unit i.year)}; use whichever reads better.{p_end}
{synopt :{opth second_stage(varlist)}}List of treatment variables. This could be, for example a 0/1 treatment dummy, a set of event-study leads/lags, or a continuous treatment variable.{p_end}
{synopt :{opth treatment(varname)}}This must be a 0/1 dummy for when treatment is occurring; D_it = 1 if unit i is treated during period t. See {help did2s##anticipation:Anticipation} for details on how to deal with anticipation (shift D_it the number of periods of anticipation).{p_end}
{synopt :{opth cluster(varname)}}What variable to cluster on (use unit id if you don't want to cluster).{p_end}
{synopt :{opth unit(varname)}}Kept for compatibility with earlier versions: the same as {cmd:absorb(}{it:varname}{cmd:)}. Numeric or string. It is no longer needed for speed, since every fixed effect is now absorbed.{p_end}

{syntab:Solver}
{synopt :{opt tol:erance(#)}}Relative residual at which the first-stage solver stops. Default {cmd:tolerance(1e-12)}.{p_end}
{synopt :{opt maxit:er(#)}}Maximum solver iterations. Default {cmd:maxiter(2000)}.{p_end}
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
The fixed effects could be biased/inconsistent if there are anticipation effects, i.e. units respond before treatment starts. The fix is fairly simple, simply “shift” treatment date earlier by as many years as you suspect anticipation to occur (e.g. 2 years before treatment starts) and estimate on the subsample where the shifted treatment equals zero. If you suspect there is anticipation, pass this modified variable to {opt treatment()}.
{p_end}

{marker absorb}
{title:Fixed effects}

{pstd}
Every first-stage fixed effect is absorbed: the first stage is solved by
preconditioned conjugate gradients from level indices and group sums, so no
dummy variables are created and no matrix with one row per level is ever
formed. How a first-stage term is handled depends only on its form:
{p_end}

{synoptset 24}{...}
{synopthdr:first-stage term}
{synoptline}
{synopt:{cmd:i.}{it:var}}absorbed fixed effect; so are {cmd:ib#.}{it:var}, {cmd:ib(first).}{it:var}, {cmd:ibn.}{it:var} and {cmd:i.(}{it:var1 var2}{cmd:)}{p_end}
{synopt:{cmd:absorb(}{it:var}{cmd:)}, {cmd:unit(}{it:var}{cmd:)}}absorbed fixed effect{p_end}
{synopt:{it:var}}continuous covariate: regressor{p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#c.}{it:var2}}group-specific slope: regressor{p_end}
{synopt:{cmd:c.}{it:var2}{cmd:##i.}{it:var1}}{cmd:i.}{it:var1} is absorbed; {cmd:c.}{it:var2} and the interaction are regressors{p_end}
{synopt:{cmd:i.}{it:var1}{cmd:#i.}{it:var2}}regressors, one per cell; for speed, {cmd:egen group()} first and absorb the result{p_end}
{synopt:{cmd:i(1/3).}{it:var}}level-restricted factor: regressor{p_end}
{synopt:(no fixed effect)}a constant is absorbed{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
The same variable named twice (say {cmd:i.unit} and {cmd:unit(unit)}) is
absorbed once. {cmd:e(absorb)} lists the absorbed dimensions.
{p_end}

{dlgtab:Speed}

{pstd}
Version 0.6 built the dummy design matrix and inverted its cross-product, so
its run time grows with the cube of the number of fixed-effect levels. The
conjugate-gradient solver is linear in the number of observations per
iteration.
{p_end}

{pstd}
Running the same {cmd:first_stage(i.unit i.year)} command in both versions;
20-period staggered panel, half the units never treated, clustered on unit;
event study with 21 relative-time indicators. Seconds, Stata 18 MP, 4 cores:
{p_end}

{col 9}{it:units}{col 18}{it:obs}{col 29}{it:static}{col 49}{it:event study}
{col 29}{it:0.6}{col 37}{it:0.7}{col 49}{it:0.6}{col 57}{it:0.7}
{col 9}{hline 54}
{col 9}  500{col 16} 10,000{col 27}  0.59{col 35}  0.07{col 47}  0.60{col 55}  0.21
{col 9}1,000{col 16} 20,000{col 27}  3.18{col 35}  0.08{col 47}  3.21{col 55}  0.38
{col 9}2,000{col 16} 40,000{col 27} 20.61{col 35}  0.15{col 47} 20.96{col 55}  0.77
{col 9}4,000{col 16} 80,000{col 27}182.96{col 35}  0.27{col 47}182.75{col 55}  1.36
{col 9}8,000{col 16}160,000{col 29}{c -}{col 35}  0.49{col 49}{c -}{col 55}  2.74
{col 9}16,000{col 16}320,000{col 29}{c -}{col 35}  1.01{col 49}{c -}{col 55}  5.24
{col 9}{hline 54}

{pstd}
Version 0.6 was not run beyond 4,000 units. Its time rose five- to ninefold
with each doubling, while version 0.7's roughly doubled.
{p_end}

{pstd}
The variance correction performs one solve per second-stage regressor, so an
event study costs more than a static specification, about five times as much
in the table above. Each solve is linear in the sample size, so this
changes the constant and not the scaling. {cmd:tests/benchmark.do} reproduces
the table.
{p_end}

{dlgtab:Accuracy}

{pstd}
The solve is iterative and stops at a relative residual of
{opt tolerance()}. At the default, point estimates matched the explicit-dummy
first stage of version 0.6 to 1e-11 or better in the package's tests, and to
6e-11 on real data with a covariate. The variance matrix agreed to between
1e-10 and 5e-8; that gap does not shrink when {opt tolerance()} is
tightened, because it originates in version 0.6, which takes a generalized
inverse of a singular matrix with {cmd:invsym()} and its tolerance-based rank
determination. If the solver does not reach {opt tolerance()} within
{opt maxiter()} iterations, {cmd:did2s} warns that the standard errors should
not be reported and sets {cmd:e(converged)} to 0. See
{help did2s##warning:Important Warning} for the main reason this happens.
{p_end}

{dlgtab:Correspondence with the R package}

{pstd}
In R the first stage is a single {cmd:fixest} formula with the absorbed
dimensions after the pipe. In Stata either form works:
{p_end}

{phang2}{cmd:first_stage = ~ 0 | unit + year}  becomes  {cmd:first_stage(i.unit i.year)} or {cmd:absorb(unit year)}{p_end}
{phang2}{cmd:first_stage = ~ x1 + x2 | unit + year}  becomes  {cmd:first_stage(i.unit i.year x1 x2)} or {cmd:first_stage(x1 x2) absorb(unit year)}{p_end}


{marker untreated}
{title:No Never-treated Group}

{pstd}
When all units are treated in a year or when some units are always treated, then the time or unit fixed effects can not be estimated.
{p_end}

{pstd}
Since version 0.7 this is detected automatically. Observations in a
fixed-effect level that has no untreated observations, such as an
always-treated unit or a year in which every unit is treated, have no
identified counterfactual. They are dropped with a note reporting how many;
{cmd:e(N_noid)} holds the count, and {cmd:e(N)} and {cmd:e(sample)} reflect the
reduced sample. The estimates are identical to those of version 0.6 run on the remaining
observations, to solver tolerance. This costs no flexibility: the check uses
only the absorbed dimensions, and any other first-stage specification, such as
gender-specific year effects, is still accepted.
{p_end}


{marker warning}
{title:Important Warning}

{pstd}
Before version 0.7, a fixed-effect level with no untreated observations (an
always-treated unit, or a period in which every unit is treated) was predicted
as 0, which could produce very odd results. Version 0.7 drops those
observations instead and reports how many. Results on such data therefore
differ from earlier versions; see {help did2s##untreated:No Never-treated Group}.
{p_end}

{pstd}
Two other cases are not screened in this way: a covariate or interaction term
(such as {cmd:i.state#c.year}) that does not vary among untreated observations,
and fixed effects whose untreated observations fall into disconnected groups.
Either can leave a treated observation's counterfactual unidentified. That
matters only if the second-stage estimates depend on it. When they do, the
variance solve cannot converge, so {cmd:did2s} warns that the standard errors
should not be reported and sets {cmd:e(converged)} to 0. When they do not, the
estimates are identified despite the gap, and {cmd:did2s} reports them
normally.
{p_end}


{marker changes}
{title:Changes in version 0.7}

{phang2}o  The first stage is solved by conjugate gradients with absorbed fixed effects. Existing {cmd:first_stage()} commands run unchanged, much faster when there are many fixed-effect levels.{p_end}
{phang2}o  New option {opt absorb()}, equivalent to entering {cmd:i.}{it:var} in {opt first_stage()}.{p_end}
{phang2}o  Observations whose fixed-effect level has no untreated observations are dropped instead of being predicted with that effect set to 0. New result {cmd:e(N_noid)}.{p_end}
{phang2}o  Weight expressions such as {cmd:[aw = pop*2]} now work.{p_end}
{phang2}o  {opt unit()} is kept as an alias for {opt absorb()} and accepts string variables.{p_end}
{phang2}o  New options {opt tolerance()} and {opt maxiter()} control the solver. New result {cmd:e(converged)}.{p_end}
{phang2}o  {cmd:e(depvar)} is now stored.{p_end}
{phang2}o  The treatment variable must be 0/1.{p_end}

{title:Examples:}

{pstd}Setup{p_end}
{phang2}{cmd:. use https://github.com/kylebutts/did2s_stata/raw/main/data/df_hom.dta, clear}{p_end}

{pstd}Static Model{p_end}
{phang2}{cmd:. did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)}{p_end}

{pstd}Event Study Model{p_end}
{phang2}{cmd:. gen rel_year_shift = rel_year + 20}{p_end}
{phang2}{cmd:. replace rel_year_shift = 100 if rel_year_shift == .}{p_end}
{phang2}{cmd:. did2s dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)}{p_end}

{pstd}With Covariates{p_end}
{phang2}{cmd:. use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear}{p_end}
{phang2}{cmd:. global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $demo) second_stage(i.post) treatment(post) cluster(sid)}{p_end}

{pstd}The same models written with {opt absorb()}; results are identical{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], first_stage($demo) absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}

{pstd}Unit rather than state fixed effects (1,000 levels){p_end}
{phang2}{cmd:. use https://github.com/kylebutts/did2s_stata/raw/main/data/df_het.dta, clear}{p_end}
{phang2}{cmd:. gen rel_year_shift = rel_year + 20}{p_end}
{phang2}{cmd:. replace rel_year_shift = 100 if rel_year_shift == .}{p_end}
{phang2}{cmd:. did2s dep_var, first_stage(i.unit i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)}{p_end}

{pstd}Three absorbed dimensions (schematic: {cmd:firm}, {cmd:worker} and {cmd:year} stand for your own variables){p_end}
{phang2}{cmd:. did2s y, absorb(firm worker year) second_stage(i.treat) treatment(treat) cluster(firm)}{p_end}




{marker results}{...}
{title:Stored results}

{pstd}
{cmd:did2s} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_noid)}}observations dropped because a fixed-effect level had no untreated observations{p_end}
{synopt:{cmd:e(converged)}}1 if the solver reached {opt tolerance()}, 0 otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:did2s}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed-effect dimensions, from {opt absorb()}, {opt unit()} and {opt first_stage()}{p_end}

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





