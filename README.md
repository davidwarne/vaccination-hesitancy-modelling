# Bayesian analysis of vaccination hesitancy patterns
Modelling and Bayesian Analysis of stochastic epidemiological model including vaccine hesitancy behaviours. 

This code is provided as supplementary material to accompany the paper DJ Warne, A Varghese, AP Browning, MM Krell, C Drovandi, W Hu, A Mira, K Mengersen, AL Jenner (2022). ``Bayesian uncertainty quantification to identify population level vaccine hesitancy behaviours''. medRxiv.org (TBA)

The work also builds upon methods presented in the paper DJ Warne, A Ebert, D Drovandi, W Hu, A Mira, K Mengersen (2020) ``Hindsight is 2020 vision: a characterisation of the global response to the COVID-19 pandemic'' BMC Public Health 20:1868 (https://doi.org/10.1186/s12889-020-09972-z)

## Developers

 The following developers implemented these tools:

 - David J. Warne [1,2] (david.warne@qut.edu.au)
 - Abhishek Varghese [2] 
 - Christopher Drovandi [1,2]

 Affilations:

   1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
   2. Centre for Data Science, Queensland University of Technology


## Model summary
The model is a discrete state continuous-time Markov Process model that is based on an SEIR model with important extensions. These extension account for case identification processes and changes in population behaviours that affect virus transmission and vaccine uptake.

## Contents

The following files are provided:

* `run_smc_vax_scenarios.m` runs analysis pipeline (computes posteriors and samples prior predictive) given a vaccine hesitancy scenario.

* `simuldata_reg_fA_vax_h.m` forwards stochastic simulation of model given a vector of parameters.

* `smry.m` summary statistic function for usage in the ABC method.

* `smc_abc_rw.m` Adaptive sequential Monte Carlo sampler for ABC. This should not need to be modified even if the model completely changes.

* `TauLeapingMethod.m` routine for approximate stochastic simulation of discrete-state continuous-time Markov process (used in model simulation).

## Usage

1. Start MATLAB from the repository root directory and run `>> init` to add all relevant file paths.
2. To reproduce figures run `>> fig2`, `>> fig3`, `>> fig5`, `>> fig6`, and `>> figS1toS8`
3. To repeat the analysis, open `run_smc_vax_scenarios.m` and edit the entries in the parameter vector `theta` as appropriate for the the desired scenario. Then run `>> run_smc_vax_scenarios`

## Acknowledgement
The development of this model and software was supported by School of Mathematical Sciences and the Centre for Data Science at QUT.

