
# did2s

<!-- badges: start -->
<!-- badges: end -->

The goal of did2s is to estimate TWFE models without running into the
problem of staggered treatment adoption. For details on the methodology,
view this
[vignette](http://kylebutts.com/did2s/articles/Two-Stage-Difference-in-Differences.html)

## Installation

You can install did2s from github with:

    net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")
    * ssc install did2s

## Two-stage Difference-in-differences (Gardner 2021)

I have created an Stata package with the help of John Gardner to
estimate the two-stage procedure. The command is `did2s` which estimates
the two-stage did procedure. This function requires the following syntax

`did2s depvar [if] [in] [weight], first_stage(varlist) treat_formula(varlist) treat_var(varname) cluster(varname)`

-   `first_stage`: formula for first stage, can include fixed effects
    and covariates, but do not include treatment variable(s)!
-   `treat_formula`: Second stage, these should be the treatment
    indicator(s) (e.g. treatment variable or es leads/lags), use i() for
    factor variables, following fixest::feols.
-   `treat_var`: This has to be the 0/1 treatment variable that marks
    when treatment turns on for a unit. If you suspect anticipation, see
    note above for accounting for this.
-   `cluster`: Which variable to cluster on.

To view the documentation, type `help did2s` into the console.

## Example Usage


    ********************************************************************************
    * Static
    ********************************************************************************

    use data/df_het.dta
        
    * Manually (note standard errors are off)
    qui reg dep_var i.state i.year if treat == 0, nocons
    predict adj, residuals
    reg adj i.treat, cluster(state) nocons


    * With did2s standard error correction  
    did2s dep_var, first_stage(i.state i.year) treat_formula(i.treat) treat_var(treat) cluster(state)

    Linear regression                               Number of obs     =     31,000
                                                    F(1, 39)          =    2787.70
                                                    Prob > F          =     0.0000
                                                    R-squared         =     0.3776
                                                    Root MSE          =     1.7506

                                     (Std. Err. adjusted for 40 clusters in state)
    ------------------------------------------------------------------------------
                 |               Robust
             adj |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
         1.treat |   2.380208   .0450809    52.80   0.000     2.289024    2.471393
    ------------------------------------------------------------------------------

                                      (Std. Err. adjusted for clustering on state)
    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
         1.treat |   2.380208   .0614314    38.75   0.000     2.259805    2.500612
    ------------------------------------------------------------------------------

You can also do event-study by changing the `treat_formula`

    use data/df_het.dta

    * can not have negatives in factor variable
    gen rel_year_shift = rel_year + 20
    replace rel_year_shift = 100 if rel_year_shift == .

    did2s dep_var, first_stage(i.state i.year) treat_formula(ib100.rel_year_shift) treat_var(treat) cluster(state)
    (11,408 missing values generated)

    (11,408 real changes made)

                                        (Std. Err. adjusted for clustering on state)
    --------------------------------------------------------------------------------
                   |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    ---------------+----------------------------------------------------------------
    rel_year_shift |
                0  |   .0746601   .0839355     0.89   0.374    -.0898505    .2391707
                1  |    .155387   .0793007     1.96   0.050    -.0000395    .3108135
                2  |   .0433077   .0861965     0.50   0.615    -.1256343    .2122497
                3  |   .0801822   .0804671     1.00   0.319    -.0775303    .2378948
                4  |   .1027144    .088289     1.16   0.245    -.0703289    .2757576
                5  |   .2168214   .0947375     2.29   0.022     .0311394    .4025035
                6  |   .1711757   .0839522     2.04   0.041     .0066325    .3357189
                7  |   .0943497   .0806924     1.17   0.242    -.0638044    .2525039
                8  |   .0902389   .0839479     1.07   0.282     -.074296    .2547738
                9  |   .1980108   .0799579     2.48   0.013     .0412963    .3547253
               10  |   .1079317   .0650773     1.66   0.097    -.0196175     .235481
               11  |   .0512958   .0586111     0.88   0.381    -.0635799    .1661715
               12  |   .0877925   .0403538     2.18   0.030     .0087006    .1668845
               13  |   .1542725   .0439659     3.51   0.000     .0681009    .2404441
               14  |   .0221227   .0509763     0.43   0.664     -.077789    .1220343
               15  |   .0351602   .0489502     0.72   0.473    -.0607804    .1311009
               16  |  -.0508045   .0504791    -1.01   0.314    -.1497417    .0481326
               17  |  -.0093709    .049563    -0.19   0.850    -.1065126    .0877707
               18  |    .008913   .0564742     0.16   0.875    -.1017744    .1196004
               19  |   .1179371   .0515514     2.29   0.022     .0168982    .2189759
               20  |    1.72709   .0827019    20.88   0.000     1.564997    1.889183
               21  |   1.752237   .0798446    21.95   0.000     1.595744    1.908729
               22  |   1.871322   .0929648    20.13   0.000     1.689114    2.053529
               23  |   1.918404   .0755407    25.40   0.000     1.770347    2.066461
               24  |   1.939901   .0841578    23.05   0.000     1.774955    2.104847
               25  |   2.145896   .0846846    25.34   0.000     1.979917    2.311874
               26  |   2.180405   .0920294    23.69   0.000     2.000031    2.360779
               27  |   2.347653   .0818133    28.70   0.000     2.187302    2.508004
               28  |   2.413051   .0764681    31.56   0.000     2.263176    2.562925
               29  |   2.619696     .10755    24.36   0.000     2.408901     2.83049
               30  |   2.681013   .0954122    28.10   0.000     2.494008    2.868017
               31  |   2.712357   .1203332    22.54   0.000     2.476509    2.948206
               32  |   2.671891   .1532795    17.43   0.000     2.371469    2.972314
               33  |    2.65582   .1224423    21.69   0.000     2.415837    2.895802
               34  |   2.754776   .1293034    21.30   0.000     2.501346    3.008206
               35  |   2.823113   .1341072    21.05   0.000     2.560267    3.085958
               36  |   2.693967   .1199888    22.45   0.000     2.458793    2.929141
               37  |   2.896505   .1265275    22.89   0.000     2.648515    3.144494
               38  |   3.130011   .1160092    26.98   0.000     2.902638    3.357385
               39  |    3.23059   .1235021    26.16   0.000      2.98853    3.472649
               40  |   3.307945   .1119849    29.54   0.000     3.088458    3.527431
    --------------------------------------------------------------------------------

This method works with pre-determined covariates as well!


    ********************************************************************************
    * Castle Doctrine
    ********************************************************************************

    use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear

    * Define Covariates
    global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44

    * No Covariates
    did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) treat_formula(i.post) treat_var(post) cluster(sid)

    * Covariates
    did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $demo) treat_formula(i.post) treat_var(post) cluster(sid)

                                        (Std. Err. adjusted for clustering on sid)
    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
          1.post |   .1454483   .1012669     1.44   0.151    -.0530313    .3439279
    ------------------------------------------------------------------------------

                                        (Std. Err. adjusted for clustering on sid)
    ------------------------------------------------------------------------------
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
          1.post |   .0802279   .0540177     1.49   0.137    -.0256447    .1861006
    ------------------------------------------------------------------------------

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Gardner_2021">

Gardner, John. 2021. “Two-Stage Difference-in-Differences.” Working
Paper. <https://jrgcmu.github.io/2sdd_current.pdf>.

</div>

</div>
