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
{cmd:did2s} {depvar} {ifin} [{it:{help regress##weight:weight}}]{cmd:,} {cmdab:second_stage(}{help varlist}{cmd:)} {cmdab:treatment(}{help varname}{cmd:)} {cmdab:cluster(}{help varname}{cmd:)} [{cmdab:absorb(}{help varlist}{cmd:)} {cmdab:first_stage(}{help varlist}{cmd:)} {it:options}]

{pstd}
One of {opt first_stage()} or {opt absorb()} is required. They may be combined.

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth absorb(varlist)}}Numeric categorical variables whose fixed effects are absorbed rather than entered as dummies. Any number of dimensions. {bf:Unit and time effects belong here.} Far faster than entering them in {opt first_stage()} once there are many levels; see {help did2s##absorb:Absorbing fixed effects}. May not be combined with {opt unit()}.{p_end}
{synopt :{opth first_stage(varlist)}}Covariates {it:other than} the absorbed fixed effects, used along with them to estimate the counterfactual Y_it(0): continuous controls, slope interactions such as {cmd:i.state#c.year}, and anything else absorption cannot express. Factor-variable notation is allowed, so fixed effects may also be entered here as dummies, which is what earlier versions required. One of {opt absorb()} or {opt first_stage()} is required; they may be combined.{p_end}
{synopt :{opth second_stage(varlist)}}List of treatment variables. This could be, for example a 0/1 treatment dummy, a set of event-study leads/lags, or a continuous treatment variable.{p_end}
{synopt :{opth treatment(varname)}}This must be a 0/1 dummy for when treatment is occuring; D_it = 1 if unit i is treated during period t. See {help did2s##anticipation:Anticipation} for details on how to deal with anticipation (shift D_it the number of periods of anticipation).{p_end}
{synopt :{opth cluster(varname)}}What variable to cluster on (use unit id if you don't want to cluster).{p_end}
{synopt :{opth unit(varname)}}Unit variable. If you use this option, *do not include unit fixed effects* in the first_stage regression. This is a trick that speeds up numerical computation that will speed things up a bunch if you have a lot of unit fixed effects. Note that the numerical equivalency only works if clutser is at the unit level (or higher!! e.g. units = people, cluster at state level) In most cases {opt absorb()} is preferable; see {help did2s##absorb:Absorbing fixed effects}.{p_end}

{syntab:Absorption}
{synopt :{opt tol:erance(#)}}Relative residual at which the first-stage solver stops. Default {cmd:tolerance(1e-12)}.{p_end}
{synopt :{opt maxit:er(#)}}Maximum solver iterations. Default {cmd:maxiter(2000)}.{p_end}
{synopt :{opt keepnoid}}Keep observations whose fixed-effect level has no untreated observations. Diagnostic only: the standard errors are not usable. See {help did2s##untreated:No Never-treated Group}.{p_end}
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
Put plain categorical fixed effects in {opt absorb()}, put everything else in
{opt first_stage()}, and use both together when the first stage needs both.
{opt absorb()} is never meaningfully slower and is dramatically faster once the
fixed effects run past a few hundred levels.
{p_end}

{p2colset 9 50 52 2}{...}
{p2col :{it:if the first stage contains}}{it:put it in}{p_end}
{p2line}
{p2col :categorical effects, any number of dimensions}{cmd:absorb(unit year ...)}{p_end}
{p2col :continuous covariates}{cmd:first_stage(x1 x2)}{p_end}
{p2col :both of the above}{cmd:first_stage(x1 x2) absorb(unit year)}{p_end}
{p2col :slope interactions, e.g. {cmd:i.state#c.year}}{cmd:first_stage(i.state#c.year)}{p_end}
{p2col :interacted effects, e.g. state-by-year cells}{helpb egen} {cmd:group()} first, then {cmd:absorb()}{p_end}
{p2line}
{p2colreset}{...}

{pstd}
Anything {opt first_stage()} accepted before it still accepts, and behaves
exactly as before. {opt absorb()} takes numeric categorical variables only, and
may not be combined with {opt unit()} — put the unit effect in {opt absorb()}
instead.
{p_end}

{dlgtab:Why absorb() is faster}

{pstd}
{opt first_stage()} builds the design matrix explicitly, so {cmd:i.unit i.year}
creates one column per level, and the first stage forms and inverts a matrix
whose dimension is the number of first-stage parameters. That is cubic in the
number of levels in time and quadratic in memory. {opt absorb()} names the
dimensions instead and solves the same normal equations by preconditioned
conjugate gradients using only level indices and group sums, so the matrix is
never formed. Memory and per-iteration cost are linear in the number of
observations.
{p_end}

{pstd}
Measured on a 20-period panel, event-study second stage, clustered on unit:
{p_end}

{p2colset 9 32 34 2}{...}
{p2col :{it:units}}{it:speedup}{p_end}
{p2line}
{p2col :50 (a state-year panel)}about even{p_end}
{p2col :500}4x{p_end}
{p2col :1,000}7x{p_end}
{p2col :2,000}24x{p_end}
{p2col :4,000}107x{p_end}
{p2col :8,000}456x{p_end}
{p2col :16,000}not comparable{p_end}
{p2line}
{p2colreset}{...}

{pstd}
At 16,000 units the explicit-dummy path is not slow but infeasible: the
first-stage cross-product matrix alone is about 2 GB before it is inverted, and
the variance calculation builds a second matrix of the same dimension.
{p_end}

{pstd}
{opt absorb()} performs one solve per second-stage regressor, so an event study
with forty leads and lags costs roughly forty times a static specification.
Each solve is linear in the sample size, so this changes the constant and not
the scaling.
{p_end}

{dlgtab:When the two give different answers}

{pstd}
On any design where every fixed-effect level has at least one untreated
observation, the two are the same estimator on the same sample and agree to
solver tolerance. Point estimates match to about 1e-12 at the default; the
variance matrix agrees to between 4e-10 and 5e-8. That variance gap does not
shrink when {opt tolerance()} is tightened, because it originates in the
explicit-dummy path, which takes a generalized inverse of a singular matrix
using {cmd:invsym()} and its tolerance-based rank determination. Absorption
resolves the same singularity by a different route. Agreement is at its worst
when weights are widely dispersed, which is what degrades the inverse.
{p_end}

{pstd}
Where some level has no untreated support, the two differ substantively, not
numerically, because they use different samples. See
{help did2s##untreated:No Never-treated Group}. This is worth checking before
switching an existing specification from one to the other.
{p_end}

{dlgtab:What absorb() gives up}

{pstd}
The fixed-effect estimates are never materialized, so unlike the explicit-dummy
path there is no design matrix to recover them from afterwards. Neither option
reports first-stage coefficients. {cmd:predict} is not supported after
{opt absorb()}.
{p_end}

{dlgtab:Correspondence with the R package}

{pstd}
In R the first stage is a single {cmd:fixest} formula with the absorbed
dimensions after the pipe. Stata's option syntax has no equivalent of the pipe,
so the same distinction is carried by two options:
{p_end}

{phang2}{cmd:first_stage = ~ 0 | unit + year}  becomes  {cmd:absorb(unit year)}{p_end}
{phang2}{cmd:first_stage = ~ x1 + x2 | unit + year}  becomes  {cmd:first_stage(x1 x2) absorb(unit year)}{p_end}


{marker untreated}
{title:No Never-treated Group}

{pstd}
When all units are treated in a year or when some units are always treated, then the time or unit fixed effects can not be estimated. The software could be designed to be able to detect this and automatically drop these observations. However, I have chosen not to do this as it would require a major change in API that would limit the flexibility of the program. As an example, requiring an option for "year fixed effects" would prevent the user from specifying gender-specific year fixed effects. The flexibility of the method (i.e. specifying a well-designed model for untreated outcomes) in my opinion outweighs the cost of not being able to detect this case (I think practically, there are relatively few cases of no never-treated group).
{p_end}

{pstd}
{opt absorb()} names the fixed-effect dimensions directly, so this detection no longer costs any flexibility: {opt first_stage()} still accepts gender-specific year effects or anything else, and behaves exactly as before. Within {opt absorb()}, levels containing no untreated observations are detected and the affected observations are dropped, with a message reporting how many. {cmd:e(N)} and {cmd:e(sample)} reflect the reduced sample. Specify {opt keepnoid} to keep them. Note that this does not reproduce what {opt first_stage()} returns on the same design: the two resolve the unidentified level inside different parameterizations, so the point estimates differ materially. {opt keepnoid} is a diagnostic. The point-estimate solve remains well posed under it, because the flagged rows carry no untreated weight, but the variance calculation does not: it uses the full-sample weight, so the right-hand side leaves the range of the operator and there is nothing for the solver to converge to. Point estimates are normalization-dependent and the standard errors must not be reported. {cmd:did2s} warns accordingly.
{p_end}


{marker warning}
{title:Important Warning}

Due to the way Stata handles the `predict` command, users of this package should be careful when they have always-treated units in their data or if all units end up treated before the end of the panel. 
If any value of `i.unit` or `i.year` are not estimated, then predict will estimate them as 0s which can predict very weird results. 

{pstd}
{opt absorb()} detects this case and reports it. See {help did2s##untreated:No Never-treated Group}.
{p_end}

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

{pstd}Absorbing the fixed effects instead. Equivalent to the above, and much faster when there are many levels{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}
{phang2}{cmd:. did2s l_homicide [aweight=popwt], first_stage($demo) absorb(sid year) second_stage(i.post) treatment(post) cluster(sid)}{p_end}

{pstd}A panel with many units, where the explicit-dummy first stage is impractical{p_end}
{phang2}{cmd:. did2s dep_var, absorb(unit year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(unit)}{p_end}

{pstd}Three absorbed dimensions{p_end}
{phang2}{cmd:. did2s y, absorb(firm worker year) second_stage(i.treat) treatment(treat) cluster(firm)}{p_end}




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
{synopt:{cmd:e(absorb)}}absorbed dimensions, if {opt absorb()} was used{p_end}

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





