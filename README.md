# Random Functions for Population Genetics (probably in R)
Somewhere to store random, useful, functions/applications for popgen

## List of functions
* *vcfR2AF()* - Function takes a vcfR object and calculates allele frequencies for all populations in a provided popmap. Popmap should be a two column data.frame; first column = individual IDs, second column = population IDs. Function can parallelise with `mclapply` and returns either the "REF" (default) or "ALT" allele frequency as defined by the VCF. Output is a matrix of allele frequencies, with *m* columns (for *m* populations) and *n* rows (for *n* snps).
