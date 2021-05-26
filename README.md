
# did2s

<!-- badges: start -->
<!-- badges: end -->

The goal of did2s is to estimate TWFE models without running into the
problem of staggered treatment adoption.

## Installation

You can install did2s from github with:

    net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")
    * ssc install did2s


## Two-stage Difference-in-differences (Gardner 2021)

Researchers often want to estimate either a static TWFE model,

<img src="man/figures/twfe.png" width="400px" height="100%">

where *μ*<sub>*i*</sub> are unit fixed effects, *μ*<sub>*t*</sub> are
time fixed effects, and *D* <sub>it</sub> is an indicator for receiving
treatment, or an event-study TWFE model

<img src="man/figures/es.png" width="450px" height="100%">

where *D* <sub>it</sub><sup>k</sup> are lag/leads of treatment (k
periods from initial treatment date). Sometimes researches use variants
of this model where they bin or drop leads and lags.

However, running OLS to estimate either model has been shown to not
recover an average treatment effect and has the potential to be severely
misleading in cases of treatment effect heterogeneity (Borusyak,
Jaravel, and Spiess 2021; Callaway and Sant’Anna 2018; Chaisemartin and
D’Haultfoeuille 2019; Goodman-Bacon 2018; Sun and Abraham 2020).

One way of thinking about this problem is through the FWL theorem. When
estimating the unit and time fixed effects, you create a residualized
*ỹ* <sub>it</sub> which is commonly said to be “the outcome variable
after removing time shocks and fixed units characteristics”, but you
also create a residulaized *D̃* <sub>it</sub> or *D̃*
<sub>it</sub><sup>k</sup>. To simplify the literature, this residualized
treatment indicators is what creates the problem of interpreting *τ* or
*τ*<sup>*k*</sup>, especially when treatment effects are heterogeneous.

That’s where Gardner (2021) comes in. What Gardner does to fix the
problem is quite simple: estimate *μ*<sub>*i*</sub> and
*μ*<sub>*t*</sub> seperately so you don’t residualize the treatment
indicators. In the absence of treatment, the TWFE model gives you a
model for (potentially unobserved) untreated outcomes

<img src="man/figures/twfe_count.png" width="350px" height="100%">

Therefore, if you can ***consistently*** estimate *y* <sub>it</sub> (0),
you can impute the untreated outcome and remove that from the observed
outcome *y* <sub>it</sub>. The value of *y* <sub>it</sub> $ - $
<sub>it</sub> (0) should be close to zero for control units and should
be close to *τ* <sub>it</sub> for treated observations. Then, regressing
*y* <sub>it</sub> $ - $ <sub>it</sub> (0) on the treatment variables
should give unbiased estimates of treatment effects (either static or
dynamic/event-study). This is the same logic as the new paper Borusyak,
Jaravel, and Spiess (2021)

The steps of the two-step estimator are:

1.  First estimate *μ*<sub>*i*</sub> and *μ*<sub>*t*</sub> using
    untreated/not-yet-treated observations, i.e. the subsample with *D*
    <sub>it</sub>  = 0. Residualize outcomes:

<img src="man/figures/resid.png" width="350px" height="100%">

1.  Regress *ỹ* <sub>it</sub> on *D* <sub>it</sub> or *D*
    <sub>it</sub><sup>k</sup>’s to estimate the treatment effect *τ* or
    *τ*<sup>*k*</sup>’s.

Some notes:

### Standard Errors

First, the standard errors on *τ* or *τ*<sup>*k*</sup>’s will be
incorrect as the dependent variable is itself an estimate. This is
referred to the generated regressor problem in econometrics parlance.
Therefore, Gardner (2021) has developed a GMM estimator that will give
asymptotically correct standard errors. Details are left to the paper,
but are implemented in the R package

### Anticipation

Second, this procedure works so long as *μ*<sub>*i*</sub> and
*μ*<sub>*t*</sub> are ***consistently*** estimated. The key is to use
only untreated/not-yet-treated observations to estimate the fixed
effects. For example, if you used observations with *D* <sub>it</sub> $
= 1$, you would attribute treatment effects *τ* as “fixed
characteristics” and would combine *μ*<sub>*i*</sub> with the treatment
effects.

The fixed effects could be biased/inconsistent if there are anticipation
effects, i.e. units respond before treatment starts. The fix is fairly
simple, simply “shift” treatment date earlier by as many years as you
suspect anticipation to occur (e.g. 2 years before treatment starts) and
estimate on the subsample where the shifted treatment equals zero. The R
package allows you to specify the variable *D* <sub>it</sub>, if you
suspect anticipation, provide the shifted variable to this option.

### Covariates

This method works with pre-determined covariates as well. Augment the
above step 1. to include *X*<sub>*i*</sub> and remove that from *y*
<sub>it</sub> along with the fixed effects to get *ỹ* <sub>it</sub>.

## Stata Package

I have created an R package with the help of John Gardner to estimate
the two-stage procedure. To install the package, run the following:

    net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")
    * ssc install did2s

To view the documentation, type `help did2s` into the console.


    ********************************************************************************
    * Static
    ********************************************************************************

    use data/df_het.dta
        
    * Manually (note standard errors are off)
    qui reg dep_var i.state i.year if treat == 0
    predict adj, residuals
    reg adj i.treat, vce(cluster state)


    * With did2s correction 
    did2s dep_var, first_stage(i.state i.year) treat_formula(i.treat) treat_var(treat) vce(cluster state)

    Linear regression                               Number of obs     =     31,000
                                                    F(1, 39)          =    2803.01
                                                    Prob > F          =     0.0000
                                                    R-squared         =     0.2896
                                                    Root MSE          =     1.7505

                                     (Std. Err. adjusted for 40 clusters in state)
    ------------------------------------------------------------------------------
                 |               Robust
             adj |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
         1.treat |   2.380156   .0449566    52.94   0.000     2.289222    2.471089
           _cons |   4.05e-10          .        .       .            .           .
    ------------------------------------------------------------------------------

                                      (Std. Err. adjusted for clustering on state)
    ------------------------------------------------------------------------------
         dep_var |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
         1.treat |   2.380208   .0504193    47.21   0.000     2.281388    2.479028
    ------------------------------------------------------------------------------

You can also do event-study by changing the `treat_formula`

    use data/df_het.dta

    * can't have negatives in factor variable
    gen rel_year_shift = rel_year + 20
    did2s dep_var, first_stage(i.state i.year) treat_formula(i.rel_year_shift) treat_var(treat) vce(cluster state)
    (11,408 missing values generated)

                                        (Std. Err. adjusted for clustering on state)
    --------------------------------------------------------------------------------
           dep_var |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    ---------------+----------------------------------------------------------------
    rel_year_shift |
                1  |    .155387    .095312     1.63   0.103    -.0314211    .3421951
                2  |   .0433077   .1103479     0.39   0.695    -.1729703    .2595857
                3  |   .0801822   .1093592     0.73   0.463    -.1341579    .2945223
                4  |   .1027144   .1056284     0.97   0.331    -.1043135    .3097423
                5  |   .2168214   .1087754     1.99   0.046     .0036255    .4300174
                6  |   .1711757   .1062672     1.61   0.107    -.0371042    .3794555
                7  |   .0943497   .0975148     0.97   0.333    -.0967758    .2854752
                8  |   .0902389   .1107409     0.81   0.415    -.1268092     .307287
                9  |   .1980108   .1155329     1.71   0.087    -.0284295    .4244511
               10  |   .1079317   .0742746     1.45   0.146    -.0376439    .2535073
               11  |   .0512958   .0798466     0.64   0.521    -.1052006    .2077922
               12  |   .0877925    .059598     1.47   0.141    -.0290174    .2046024
               13  |   .1542725   .0689907     2.24   0.025     .0190531    .2894918
               14  |   .0221227   .0693362     0.32   0.750    -.1137739    .1580192
               15  |   .0351602   .0679596     0.52   0.605    -.0980382    .1683586
               16  |  -.0508045   .0838868    -0.61   0.545    -.2152196    .1136105
               17  |  -.0093709   .0707269    -0.13   0.895    -.1479932    .1292513
               18  |    .008913    .075345     0.12   0.906    -.1387605    .1565865
               19  |   .1179371   .0720275     1.64   0.102    -.0232343    .2591084
               20  |    1.72709   .0825903    20.91   0.000     1.565216    1.888964
               21  |   1.752237   .0787708    22.24   0.000     1.597849    1.906625
               22  |   1.871322   .0872593    21.45   0.000     1.700296    2.042347
               23  |   1.918404   .0798623    24.02   0.000     1.761876    2.074931
               24  |   1.939901   .0822252    23.59   0.000     1.778743     2.10106
               25  |   2.145896   .0836287    25.66   0.000     1.981986    2.309805
               26  |   2.180405   .0853684    25.54   0.000     2.013086    2.347724
               27  |   2.347653   .0803612    29.21   0.000     2.190148    2.505158
               28  |   2.413051   .0730497    33.03   0.000     2.269876    2.556226
               29  |   2.619696   .0968117    27.06   0.000     2.429948    2.809443
               30  |   2.681013   .0868222    30.88   0.000     2.510844    2.851181
               31  |   2.712357   .1260245    21.52   0.000     2.465354    2.959361
               32  |   2.671891   .1386156    19.28   0.000      2.40021    2.943573
               33  |    2.65582   .1318246    20.15   0.000     2.397448    2.914191
               34  |   2.754776   .1318009    20.90   0.000     2.496451    3.013101
               35  |   2.823113   .1179188    23.94   0.000     2.591996    3.054229
               36  |   2.693967   .1125634    23.93   0.000     2.473347    2.914587
               37  |   2.896505   .1193611    24.27   0.000     2.662561    3.130448
               38  |   3.130011   .1191483    26.27   0.000     2.896485    3.363538
               39  |    3.23059   .1120913    28.82   0.000     3.010895    3.450285
               40  |   3.307945   .1195181    27.68   0.000     3.073694    3.542196
    --------------------------------------------------------------------------------

This method works with pre-determined covariates as well!


    ********************************************************************************
    * Castle Doctrine
    ********************************************************************************

    use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear

    * Define Covariates
    global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44

    * No Covariates
    did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) treat_formula(i.post) treat_var(post) vce(cluster sid)

    * Covariates
    did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $demo) treat_formula(i.post) treat_var(post) vce(cluster sid)

                                        (Std. Err. adjusted for clustering on sid)
    ------------------------------------------------------------------------------
      l_homicide |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
          1.post |   .1454483   .0824557     1.76   0.078    -.0161618    .3070584
    ------------------------------------------------------------------------------

                                        (Std. Err. adjusted for clustering on sid)
    ------------------------------------------------------------------------------
      l_homicide |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
          1.post |   .0802279   .0375127     2.14   0.032     .0067043    .1537516
    ------------------------------------------------------------------------------

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Borusyak_Jaravel_Spiess_2021">

Borusyak, Kirill, Xavier Jaravel, and Jann Spiess. 2021. “Revisiting
Event Study Designs: Robust and Efficient Estimation,” 48.

</div>

<div id="ref-Callaway_SantAnna_2018">

Callaway, Brantly, and Pedro H. C. Sant’Anna. 2018.
“Difference-in-Differences with Multiple Time Periods and an Application
on the Minimum Wage and Employment.” *arXiv:1803.09015 \[Econ, Math,
Stat\]*, August. <http://arxiv.org/abs/1803.09015>.

</div>

<div id="ref-deChaisemartin_DHaultfoeuille_2019">

Chaisemartin, Clement de, and Xavier D’Haultfoeuille. 2019. *Two-Way
Fixed Effects Estimators with Heterogeneous Treatment Effects*. w25904.
National Bureau of Economic Research. <https://doi.org/10.3386/w25904>.

</div>

<div id="ref-Gardner_2021">

Gardner, John. 2021. “Two-Stage Difference-in-Differences.” Working
Paper. <https://jrgcmu.github.io/2sdd_current.pdf>.

</div>

<div id="ref-Goodman-Bacon_2018">

Goodman-Bacon, Andrew. 2018. *Difference-in-Differences with Variation
in Treatment Timing*. w25018. National Bureau of Economic Research.
<https://doi.org/10.3386/w25018>.

</div>

<div id="ref-Sun_Abraham_2020">

Sun, Liyang, and Sarah Abraham. 2020. “Estimating Dynamic Treatment
Effects in Event Studies with Heterogeneous Treatment Effects,” 53.

</div>

</div>
