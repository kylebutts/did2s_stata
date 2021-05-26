{smcl}
{it:v. 0.1.0} 
{viewerjumpto "Anticipation" "did2s##anticipation"}{...}

{title:Two-Stage Difference-in-Differences}

{pstd}
{bf:did2s} - Estimates a TWFE model using the two-stage difference-in-differences approach from {browse "https://jrgcmu.github.io/2sdd_current.pdf":Gardner (2021)}

{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd:did2s} {depvar}{cmd:,} {cmdab:first_stage(}{help varlist}{cmd:)} {cmdab:treat_formula(}{help varlist}{cmd:)} {cmdab:treat_var(}{help varname}{cmd:)} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opth first_stage(varlist)}}Fixed effects and covariates that will be used to estimate counterfactual Y_it(0){p_end}
{synopt :{opth treat_formula(varlist)}}List of treatment dummies. This could be for example a 0/1 Treat dummy or a set of event-study leads/lags.{p_end}
{synopt :{opth treat_var(varname)}}This must be a 0/1 dummy for when treatment is occuring. See {help did2s##anticipation:Anticipation} for details on how to deal with anticipation.{p_end}
{syntab:SE/Robust}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {bf: robust} or {bf:cluster}.{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{bf:did2s} implements Two-Staged Difference-in-Differences by Garener (2021). A TWFE model for outcomes is given by unit/group fixed effects, time fixed effects, treatment variable (or variables in the case of event study), and potentially covariates. To avoid the problems with OLS estimation of difference-in-differences/event-studies in the presence of staggered treatment adoption, this method proceeds in two stages:{p_end}

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
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:Program Author}

{pstd}
Kyle Butts   {break}
University of Colorado, Boulder      {break}
buttskyle96@gmail.com     {break}


{title:References}

{phang}
Gardner (2021), {browse "https://jrgcmu.github.io/2sdd_current.pdf":Two-stage Difference-in-differences}
{p_end}





