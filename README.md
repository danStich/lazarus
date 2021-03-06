# lazarus

<br>
 
Mark-recapture survival models with state uncertainty for radio-tagged animals.

This repository contains `R` and `BUGS` code used for generalized estimation of organismal survival when observations of mortality contain false-positive identification of fish (e.g., in closed systems or for animals with unreliable mortality-sensing radio tags). This study came about because fish were frequently relocated long distances from previous locations, and were previously assumed to be dead. This occurred because they remained fairly still for long periods of time and frequently sent radio-tags into a burst rate that indicated they were dead.

<br>
 
## Files

<br>
 
`Lazarus_sims.R` contains all code necessary to simulate and analyze data for the validation study conducted in the citation provided below. [JAGS](http://sourceforge.net/projects/mcmc-jags/files/l) must be installed to run the code, and the package `R2jags` must be installed in R.

The code makes extensive use of R code from [Kéry and Schaub (2010) _Bayesian Population Analysis Using WinBUGS_](http:/www.vogelwarte.ch/de/projekte/publikationen/bpa/).

<br>
 
## Supporting publications

<br>
 
The citation for the study for which this code was written is:

**Stich, D. S.**, Y. Jiao, and B. R. Murphy. 2015. Life, death, and resurrection: accounting for state uncertainty in survival estimation from tagged grass carp. North American Journal of Fisheries Management 35:321-330.
