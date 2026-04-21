# SparseHybrids_TCstage1
Simulations for stage 1 testcrosses: GCA and SCA

Testcross_Scripts:
* 4Loop_testerchoice_dom.R: Script over different dominance levels and number of testers. Runs the main simulation + analysis, looping over dominance levels and number of testers, then summarises accuracy (e.g., GCA correlations/rankings) and writes outputs to the results folder. This is the updated version with graphs and results presented in the poster session of 20th April 2026 at Roslin :)

* Loop_dom_test_it.R: Script over different dominance levels and number of testers. Runs the main simulation + analysis, looping over dominance levels and number of testers, then summarises accuracy (e.g., GCA correlations/rankings) and writes outputs to the results folder.

* Step1_testcross_dominance.Rmd: Additive + dominance, but still simplified (no full looping; uses a simple dominance setting). Kept as an intermediate/reference.

* Step1_testcross.Rmd: First, simplest version: additive-only trait. Kept as a baseline/reference.
