# SparseHybrids_TCstage1
Simulations for stage 1 testcrosses: GCA and SCA

Testcross_Scripts:
* Loop_dom_test_it.R: Script over different dominance levels and number of testers. Runs the main simulation + analysis, looping over dominance levels and number of testers, then summarises accuracy (e.g., GCA correlations/rankings) and writes outputs to the results folder.

* Step1_testcross_dominance.Rmd: Additive + dominance, but still simplified (no full looping; uses a simple dominance setting). Kept as an intermediate/reference.

* Step1_testcross.Rmd: First, simplest version: additive-only trait. Kept as a baseline/reference.
