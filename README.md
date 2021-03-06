MCMC_legacy
======================

Markov Chain Monte Carlo code for analyzing legacy data using tapered covariance multilevel model.
Details on the theoretical justification and MCMC procedure can be found in multilevelmodel_writeup.pdf.

Prerequisites
======================
1. Place Horizon.csv and Profile.csv in the code folder.

Execution:
======================
1. main_code.R: By changing the soil_property variable in the code, this code prepares data for the following MCMC procedure (Y, X, d.site);
2. Depending on the statistics model that you want to estimate, you have the following options:
    a. multilevel continuous normal model: MCMC_geostat_continuous_list.R
    b. multilevel truncated normal model: MCMC_geostat_truncnorm_list.R
    c. multilevel beta model: MCMC_geostat_beta_list_v2.R (still debugging)
    d. multilevel logistic model: MCMC_geostat_cat_list.R (still debugging)

After running step 2, a RData file will be saved in the current directory, which is needed for future prediction procedures.


Author: Jiehua Chen <jc3288@columbia.edu>
