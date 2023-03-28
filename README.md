
# did2s

<!-- badges: start -->
<!-- badges: end -->

The goal of did2s is to estimate TWFE models without running into the
problem caused by staggered treatment adoption. For details on the
methodology, view this
[vignette](http://kylebutts.com/did2s/articles/Two-Stage-Difference-in-Differences.html)

## Installation

You can install did2s from github with:

``` stata
net install did2s, replace from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")
* ssc install did2s
```

## Two-stage Difference-in-differences (Gardner 2021)

I have created an Stata package with the help of John Gardner to
estimate the two-stage procedure. The command is `did2s` which estimates
the two-stage did procedure. This function requires the following syntax

`did2s depvar [if] [in] [weight], first_stage(varlist) second_stage(varlist) treatment(varname) cluster(varname)`

- `first_stage`: formula for first stage, can include fixed effects and
  covariates, but do not include treatment variable(s)!
- `second_stage`: List of treatment variables. This could be, for
  example a 0/1 treatment dummy, a set of event-study leads/lags, or a
  continuous treatment variable
- `treatment`: This has to be the 0/1 treatment variable that marks when
  treatment turns on for a unit. If you suspect anticipation, see note
  above for accounting for this.
- `cluster`: Which variable to cluster on.

To view the documentation, type `help did2s` into the console.

## Example Usage

``` stata

********************************************************************************
* Static
********************************************************************************

use data/df_het.dta
    
* Manually (note standard errors are off)
qui reg dep_var i.state i.year if treat == 0, nocons
predict adj, residuals
reg adj i.treat, cluster(state) nocons

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
```

``` stata


* With did2s standard error correction  
did2s dep_var, first_stage(i.state i.year) second_stage(i.treat) treatment(treat) cluster(state)

                                  (Std. Err. adjusted for clustering on state)
------------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
     1.treat |   2.380156   .0614383    38.74   0.000     2.259739    2.500573
------------------------------------------------------------------------------
```

You can also do event-study by changing the `second_stage`

``` stata
use data/df_het.dta

* can not have negatives in factor variable
gen rel_year_shift = rel_year + 20
replace rel_year_shift = 100 if rel_year_shift == .

did2s dep_var, first_stage(i.state i.year) second_stage(ib100.rel_year_shift) treatment(treat) cluster(state)
(11,408 missing values generated)

(11,408 real changes made)

                                    (Std. Err. adjusted for clustering on state)
--------------------------------------------------------------------------------
               |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
---------------+----------------------------------------------------------------
rel_year_shift |
            0  |    .049467   .0795074     0.62   0.534    -.1063647    .2052986
            1  |   .1550051   .0793407     1.95   0.051    -.0004999      .31051
            2  |   .0429258   .0861871     0.50   0.618    -.1259978    .2118494
            3  |   .0798003   .0804802     0.99   0.321    -.0779379    .2375386
            4  |   .1023325   .0882446     1.16   0.246    -.0706237    .2752886
            5  |   .2164395   .0947508     2.28   0.022     .0307313    .4021477
            6  |   .1707938    .083863     2.04   0.042     .0064253    .3351622
            7  |   .0939678   .0805816     1.17   0.244    -.0639692    .2519048
            8  |    .089857   .0839759     1.07   0.285    -.0747327    .2544467
            9  |   .1976289   .0799969     2.47   0.013     .0408378      .35442
           10  |   .0952583   .0619518     1.54   0.124    -.0261651    .2166817
           11  |   .0512636   .0586073     0.87   0.382    -.0636046    .1661318
           12  |   .0877603   .0403517     2.17   0.030     .0086725    .1668482
           13  |   .1542402   .0439602     3.51   0.000     .0680799    .2404006
           14  |   .0220904   .0509833     0.43   0.665    -.0778349    .1220158
           15  |    .035128   .0489484     0.72   0.473    -.0608091    .1310651
           16  |  -.0508368   .0504833    -1.01   0.314    -.1497822    .0481087
           17  |  -.0094032   .0495642    -0.19   0.850    -.1065472    .0877408
           18  |   .0088808   .0564702     0.16   0.875    -.1017987    .1195602
           19  |   .1179048   .0515519     2.29   0.022      .016865    .2189447
           20  |   1.726992   .0826976    20.88   0.000     1.564907    1.889076
           21  |   1.752138   .0798351    21.95   0.000     1.595664    1.908612
           22  |   1.871223   .0929743    20.13   0.000     1.688997    2.053449
           23  |   1.918305   .0755331    25.40   0.000     1.770263    2.066347
           24  |   1.939803   .0841477    23.05   0.000     1.774876    2.104729
           25  |   2.145797   .0846879    25.34   0.000     1.979812    2.311782
           26  |   2.180307   .0920339    23.69   0.000     1.999923     2.36069
           27  |   2.347555   .0818049    28.70   0.000      2.18722    2.507889
           28  |   2.412952   .0764437    31.57   0.000     2.263125    2.562779
           29  |   2.619597   .1075448    24.36   0.000     2.408813    2.830381
           30  |   2.680793   .0954052    28.10   0.000     2.493802    2.867784
           31  |   2.712427    .120355    22.54   0.000     2.476536    2.948319
           32  |   2.671961   .1533243    17.43   0.000     2.371451    2.972471
           33  |    2.65589   .1224654    21.69   0.000     2.415862    2.895917
           34  |   2.754846   .1293217    21.30   0.000      2.50138    3.008312
           35  |   2.823183   .1341382    21.05   0.000     2.560277    3.086089
           36  |   2.694037   .1199969    22.45   0.000     2.458847    2.929226
           37  |   2.896575   .1265512    22.89   0.000     2.648539     3.14461
           38  |   3.130081   .1160177    26.98   0.000     2.902691    3.357472
           39  |    3.23066   .1235224    26.15   0.000      2.98856    3.472759
           40  |   3.308015   .1120092    29.53   0.000     3.088481    3.527549
--------------------------------------------------------------------------------
```

This method works with pre-determined covariates as well!

``` stata

********************************************************************************
* Castle Doctrine
********************************************************************************

use https://github.com/scunning1975/mixtape/raw/master/castle.dta, clear

* Define Covariates
global demo blackm_15_24 whitem_15_24 blackm_25_44 whitem_25_44

* No Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year) second_stage(i.post) treatment(post) cluster(sid)

* Covariates
did2s l_homicide [aweight=popwt], first_stage(i.sid i.year $demo) second_stage(i.post) treatment(post) cluster(sid)

                                    (Std. Err. adjusted for clustering on sid)
------------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      1.post |   .0751416   .0353795     2.12   0.034     .0057991    .1444842
------------------------------------------------------------------------------

                                    (Std. Err. adjusted for clustering on sid)
------------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      1.post |   .0760161   .0324715     2.34   0.019     .0123731    .1396591
------------------------------------------------------------------------------
```

### Large Datasets or Many Fixed Effects

There are some situations where standard errors can not be calculate
analytically in memory. The reason for this is that the analytic
standard errors require the creation of the matrix containing all the
fixed effects used in estimation. When there are a lot of observations
and/or many fixed effects, this matrix can’t be stored in memory.

In this case, it’s possible to obtain standard errors via bootstrapping
a custom program. Here is an example for the example data. You could
spend time to make the command more programmable with args, but I find
it easier to just write the estimation out.

``` stata
use data/df_het.dta, clear

egen unique_id = group(state unit)

capture program drop did2s_est

program did2s_est, rclass
    version 13.0
    regress dep_var i.new_id i.year if treat == 0
    tempvar dep_var_resid
    predict `dep_var_resid', residuals
    regress `dep_var_resid' ib0.treat, nocons
end

xtset unique_id year
sort unique_id year
bootstrap, cluster(state) idcluster(new_id) group(unique_id) reps(100): did2s_est
       panel variable:  unique_id (strongly balanced)
        time variable:  year, 1990 to 2020
                delta:  1 unit


(running did2s_est on estimation sample)

Bootstrap replications (100)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
..................................................    50
..................................................   100

Linear regression                               Number of obs     =     31,000
                                                Replications      =        100
                                                Wald chi2(1)      =    1568.60
                                                Prob > chi2       =     0.0000
                                                R-squared         =     0.3776
                                                Adj R-squared     =     0.3776
                                                Root MSE          =     1.7505

                                  (Replications based on 40 clusters in state)
------------------------------------------------------------------------------
             |   Observed   Bootstrap                         Normal-based
    __000001 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
     1.treat |   2.380156   .0600965    39.61   0.000     2.262369    2.497943
------------------------------------------------------------------------------
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Gardner_2021" class="csl-entry">

Gardner, John. 2021. “<span class="nocase">Two-Stage
Difference-in-Differences</span>.” Working Paper.
<https://jrgcmu.github.io/2sdd_current.pdf>.

</div>

</div>
