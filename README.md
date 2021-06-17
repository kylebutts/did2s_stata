
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

-   `first_stage`: formula for first stage, can include fixed effects
    and covariates, but do not include treatment variable(s)!
-   `second_stage`: List of treatment variables. This could be, for
    example a 0/1 treatment dummy, a set of event-study leads/lags, or a
    continuous treatment variable
-   `treatment`: This has to be the 0/1 treatment variable that marks
    when treatment turns on for a unit. If you suspect anticipation, see
    note above for accounting for this.
-   `cluster`: Which variable to cluster on.

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
     1.treat |   2.380156   .0614477    38.73   0.000      2.25972    2.500591
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
            0  |    .049467   .0780278     0.63   0.526    -.1034647    .2023986
            1  |   .1550051   .0793378     1.95   0.051    -.0004941    .3105043
            2  |   .0429258   .0861951     0.50   0.618    -.1260136    .2118651
            3  |   .0798003   .0804722     0.99   0.321    -.0779223     .237523
            4  |   .1023325   .0882875     1.16   0.246    -.0707078    .2753728
            5  |   .2164395   .0947135     2.29   0.022     .0308044    .4020746
            6  |   .1707938   .0838902     2.04   0.042     .0063719    .3352156
            7  |   .0939678   .0806006     1.17   0.244    -.0640064     .251942
            8  |    .089857   .0839978     1.07   0.285    -.0747756    .2544897
            9  |   .1976289   .0799682     2.47   0.013     .0408941    .3543637
           10  |   .0952583   .0640277     1.49   0.137    -.0302337    .2207503
           11  |   .0512636   .0586116     0.87   0.382     -.063613    .1661402
           12  |   .0877603   .0403556     2.17   0.030     .0086649    .1668558
           13  |   .1542402   .0439649     3.51   0.000     .0680706    .2404099
           14  |   .0220904   .0509844     0.43   0.665    -.0778372    .1220181
           15  |    .035128   .0489496     0.72   0.473    -.0608115    .1310675
           16  |  -.0508368   .0504788    -1.01   0.314    -.1497735       .0481
           17  |  -.0094032    .049567    -0.19   0.850    -.1065527    .0877463
           18  |   .0088808   .0564685     0.16   0.875    -.1017955    .1195571
           19  |   .1179048   .0515501     2.29   0.022     .0168685    .2189412
           20  |   1.726992   .0827157    20.88   0.000     1.564872    1.889111
           21  |   1.752138   .0798513    21.94   0.000     1.595633    1.908644
           22  |   1.871223   .0929891    20.12   0.000     1.688968    2.053478
           23  |   1.918305   .0755381    25.40   0.000     1.770253    2.066357
           24  |   1.939803   .0841515    23.05   0.000     1.774869    2.104737
           25  |   2.145797    .084697    25.33   0.000     1.979794      2.3118
           26  |   2.180307   .0920414    23.69   0.000     1.999909    2.360704
           27  |   2.347555   .0818194    28.69   0.000     2.187192    2.507918
           28  |   2.412952   .0764653    31.56   0.000     2.263083    2.562822
           29  |   2.619597   .1075586    24.36   0.000     2.408786    2.830408
           30  |   2.680793   .0954201    28.09   0.000     2.493773    2.867813
           31  |   2.712427   .1203408    22.54   0.000     2.476564    2.948291
           32  |   2.671961   .1533124    17.43   0.000     2.371475    2.972448
           33  |    2.65589   .1224591    21.69   0.000     2.415874    2.895905
           34  |   2.754846    .129317    21.30   0.000      2.50139    3.008303
           35  |   2.823183   .1341359    21.05   0.000     2.560281    3.086084
           36  |   2.694037   .1199964    22.45   0.000     2.458848    2.929225
           37  |   2.896575    .126548    22.89   0.000     2.648545    3.144604
           38  |   3.130081   .1160191    26.98   0.000     2.902688    3.357475
           39  |    3.23066   .1235199    26.15   0.000     2.988565    3.472754
           40  |   3.308015   .1120013    29.54   0.000     3.088496    3.527533
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
      1.post |   .0751416   .0572312     1.31   0.189    -.0370295    .1873128
------------------------------------------------------------------------------

                                    (Std. Err. adjusted for clustering on sid)
------------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      1.post |   .0760161   .0509904     1.49   0.136    -.0239232    .1759554
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
