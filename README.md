MCMC_legacy
===========

MCMC code for analyzing legacy data using taped covariance multilevel model.
The details for the theoretical justificaiton and MCMC procedure are in multilevelmodel_writeup.pdf.

prerequisite:

1. you need to have Horizon.csv, and Profile.csv saved in the git folder

steps to run:

1. step 1: main_code.R: by changing the soil_property variable in the code, this code prepares data for the following MCMC procedure (Y, X, d.site);
2. step 2: after running step 1, depending on the statistics model that you want to estimate, you have the following options:
    a. multilevel continuous normal model: MCMC_geostat_continuous_list.R
    b. multilevel truncated normal model: MCMC_geostat_truncnorm_list.R
    c. multilevel beta model: MCMC_geostat_beta_list_v2.R (still debugging)
    d. multilevel logistic model: MCMC_geostat_cat_list.R (still debugging)

After running step 2, a RData file will be saved in the current directory, which is needed for future prediction procedure.


Author: Jiehua Chen <jc3288@columbia.edu>
