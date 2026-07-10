# before you start
- always familiarise yourself with 
    - `../MSFR.pdf`, the main technical resource for this package and all functions concerning multi-study factor regression
    - `../FA.pdf`, the resource for all functions concerning factor analysis
    - `../FR.pdf`, the resource for all functions concerning factor regression
    - confirm you have read them

# file structure
- `R/ecm_msfr.R` is the main function for the multi-study factor regression algorithm
- all helpers which are shared across the three main modules `R/ecm_msfr.R`, `R/ecm_fr.R` and `R/ecm_fa.R` are stored in `R/helpers.R`
- helpers that are specific to only one of the main modules are stored in module's `.R` file.
- `R/start_fafr.R` is a function handling the initialization of parameters for the factor analysis and factor regression modules
- `R/start_msfr.R` handles parameter initialisation of the multi-study factor regression module. 
- `R/heat_plot.R` and `R/vcov_msfr.R` are both post-processing modules
- you must only modify files in the `R/`directory of this package project, unless specificly stated otherwise by me

# project requirements
- functions and modules must be optimized for long EM algorithm calculations
    - standard `niter = 50000`
- number of studies `S` is assumed to be small 
- number of observations per study `n_s` are assumed to be large
    - with small `n_s = 500` as a reference for a lower bound 
- number of function calls must to be kept at a minimum
- redundancy must be avoided
- direct BLAS calls are fastest for runtime
