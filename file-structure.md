# file structure
- `R/ecm_msfr.R` is the main function for the multi-study factor regression algorithm
- all helpers which are shared across the three main modules `R/ecm_msfr.R`, `R/ecm_fr.R` and `R/ecm_fa.R` are stored in `R/helpers.R`
- helpers that are specific to only one of the main modules are stored in module's `.R` file.
- `R/start_fafr.R` is a function handling the initialization of parameters for the factor analysis and factor regression modules
- `R/start_msfr.R` handles parameter initialisation of the multi-study factor regression module.
- `R/exp_values_msfr.R` is a module that handles the estimation of expected values for the `R/ecm_msfr.R` module
- `R/exp_values_fr.R` is a module that handles the estimation of expected values for the `R/ecm_fr.R` module
- `R/heat_plot.R` and `R/vcov_msfr.R` are both post-processing modules
- you must only modify files in the `R/`directory of this package project, unless specificly stated otherwise by me