Empirical harvest control rules for ICES categories 4-6
================

This repository
([GA_MSE_cat456](https://github.com/shfischer/GA_MSE_cat456)) is a
mirror of [GA_MSE_HR](https://github.com/shfischer/GA_MSE_HR),
[GA_MSE](https://github.com/shfischer/GA_MSE), and
[GA_MSE_PA](https://github.com/shfischer/GA_MSE_PA) with the `cat456`
branch displayed as default branch.

## Introduction

This repository contains the code for an exploration of the ICES
approach for categories 4-6. It includes

- the typical ICES harvest control rule (constant catch with PA buffer,
  `MP = "constant_catch"`)

- a new harvest control rule that adjust the catch based on the slope of
  a length-based catch curve (`fcc`)

- several exploratory (but abandoned) harvest control rules (adjusting
  the constant catch based on mean catch length `CC_f`, guessing the
  stock trend from catch and mean catch length time series `CL`,
  guessing the stock trend from catch and length-based catch curves
  `CL_cc`)

- it also includes code for the ICES category 3 methods (hr and rfb
  rules, `rfb`, `chr`)

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and the Assessment for All (a4a)
standard MSE framework ([`FLR/mse`](github.com/FLR/mse)) developed
during the Workshop on development of MSE algorithms with R/FLR/a4a
([Jardim et al.,
2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

The simulation framework is based on the
([GA_MSE_HR](https://github.com/shfischer/GA_MSE_HR)) repository which
explores the use of harvest rates and contains the code for the
publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2022). Exploring a relative harvest rate strategy for moderately
> data-limited fisheries management. ICES Journal of Marine Science. 79:
> 1730-1741. <https://doi.org/10.1093/icesjms/fsac103>.

The operating models provided as an input are those from the repository
[shfischer/wklifeVII](https://github.com/shfischer/wklifeVII) as
described in:

> Fischer, S. H., De Oliveira, J. A. A., and Laurence T. Kell (2020).
> Linking the performance of a data-limited empirical catch rule to
> life-history traits. ICES Journal of Marine Science, 77: 1914-1926.
> <https://doi.org/10.1093/icesjms/fsaa054>.

but the operating models were updated
([`OM_lh.R`](https://github.com/shfischer/GA_MSE_cat456/blob/cat456/OM_lh.R)).

## Repository structure

The root folder contains the following R scripts:

- `OM_lh.R`: This script creates the (updated) operating models (OMs)
- `OM.R`: prepares the OMs for the MSE
- `funs_OM.R` contains functions and methods used for the creation of
  the operating models
- `funs.R` contains functions for running the MSE such as the OM and MP
  modules,
- `funs_GA.R` contains the function used in the optimisation procedure,
- `run_ms_cat456.R` is an R script for running MSE projections and is
  called from a job submission script
- `run*.pbs` are job submission scripts which are used on a high
  performance computing cluster and call `run_ms_cat456.R`
- \`run_ms_hr.R\`\` legacy script for running the chr rule
- `analysis_cat456.R` is for analysing the results

The following input files are provided:

- `input/stocks.csv` contains the stock definitions and life-history
  parameters
- `input/brps.rds` contains the FLBRP objects which are the basis for
  the OMs (old, not used anymore)
- `input/brps_new.rds` updated FLBRP objects created in `OM_lh.R`
- `input/brps_dome_new.rds` FLBRP objects with dome-shaped selectivity
  (based on `brps_new.R`

## R, R packages and version info

This is branch `cat456` in which R and R packages have been updated to R
4.4:

``` r
> sessionInfo()
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)
```

The package versions and their dependencies are recorded with the R
package [renv](https://rstudio.github.io/renv/) and stored in the file
[renv.lock](https://github.com/shfischer/GA_MSE_cat456/blob/cat456/renv.lock).
The exact package version can be restored by cloning this repository,
checking out the current branch (`cat456`), navigating into this folder
in R (or setting up a project), installing the renv package

``` r
install.packages("renv")
```

and calling

``` r
renv::restore()
```

See [renv](https://rstudio.github.io/renv/) and the package
documentation for details.

The framework is based on the Fisheries Library in R (FLR) framework and
uses the [FLR packages](https://flr-project.org/) `FLCore`, `FLasher`,
`FLBRP`, `ggplotFL`, `mse`. See
[renv.lock](https://github.com/shfischer/GA_MSE_cat456/blob/cat456/renv.lock)
for version details and sources.

The FLR package versions can also be installed manually with `remotes`
(requires suitable tools to compile R packages):

``` r
remotes::install_github(repo = "flr/FLCore", ref = "5dac55024c83fc6ee198d780d9cd810819c574a1")
remotes::install_github(repo = "flr/FLasher", ref = "84268b3bc21941bfde20dc9cbb0cfc3367f570a7")
remotes::install_github(repo = "flr/FLBRP", ref = "9fad7869462eb71456cd9e50a46e788a3b46f7f4")
remotes::install_github(repo = "shfischer/mse", ref = "fbea3bc9351a9efd831baa6be11104e0e26bb569")
```

For using MPI parallelisation, an MPI backend such as OpenMPI and the R
packages `Rmpi` and `doMPI` are required. This is not required and most
simulations also run with standard parallelisation using `doParallel`.
