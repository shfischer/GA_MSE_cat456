Exploring a relative harvest rate strategy for moderately data-limited
fisheries management
================

This repository ([GA_MSE_HR](https://github.com/shfischer/GA_MSE_HR)) is
a mirror of [GA_MSE](https://github.com/shfischer/GA_MSE) and
[GA_MSE_PA](https://github.com/shfischer/GA_MSE_PA) with the
`harvest_rate` branch displayed as default branch.

## Introduction

This repository contains the code for testing a harvest rate-based
management procedure and optimising it with a genetic algorithm. The
simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and the Assessment for All (a4a)
standard MSE framework ([`FLR/mse`](github.com/FLR/mse)) developed
during the Workshop on development of MSE algorithms with R/FLR/a4a
([Jardim et al.,
2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

This is the **`harvest_rate`** branch which explores the use of harvest
rates and contains the code for the publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2022). Exploring a relative harvest rate strategy for moderately
> data-limited fisheries management. ICES Journal of Marine Science. 79:
> 1730-1741. <https://doi.org/10.1093/icesjms/fsac103>.

The `master` branch ([GA_MSE](https://github.com/shfischer/GA_MSE))
contains the code for the publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2021). Using a genetic algorithm to optimise a data-limited catch
> rule. ICES Journal of Marine Science. 78: 1311-1323.
> <https://doi.org/10.1093/icesjms/fsab018>.

The `PA` branch ([GA_MSE_PA](https://github.com/shfischer/GA_MSE_PA))
includes the optimisation with specific risk limits for the ICES
precautionary approach (PA) and contains the code for the publication:

> Fischer, S. H., De Oliveira, J. A. A., Mumford, J. D., and Kell, L. T.
> (2021). Application of explicit precautionary principles in
> data-limited fisheries management. ICES Journal of Marine Science. 78:
> 2931-2942. <https://doi.org/10.1093/icesjms/fsab169>.

The operating models provided as an input are those from the repository
[shfischer/wklifeVII](https://github.com/shfischer/wklifeVII) as
described in:

> Fischer, S. H., De Oliveira, J. A. A., and Laurence T. Kell (2020).
> Linking the performance of a data-limited empirical catch rule to
> life-history traits. ICES Journal of Marine Science, 77: 1914-1926.
> <https://doi.org/10.1093/icesjms/fsaa054>.

## Repository structure

The root folder contains the following R scripts:

- `OM_hr.R`: This script creates the operating models (OMs),
- `funs.R` contains functions and methods used for the creation of the
  operating models and for running the MSE,
- `funs_GA.R` contains the function used in the optimisation procedure,
- `run_ms_hr.R` is an R script for running MSE projections and is called
  from a job submission script
- `run*.pbs` are job submission scripts which are used on a high
  performance computing cluster and call `run_ms.R`
- `analysis_hr.R` is for analysing the results

The following input files are provided:

- `input/stocks.csv` contains the stock definitions and life-history
  parameters
- `input/brps.rds` contains the FLBRP objects which are the basis for
  the OMs
- `input/catch_rates.rds` examples catch rates

Summarised outputs are provided in:

- `output/`

## R, R packages and version info

The MSE simulations were run on a high-performance computing cluster:

``` r
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)
```

This is branch `harvest_rate_R4.2` in which R and R packages have been
updated to R 4.2 (4.3 should also work):

``` r
> sessionInfo()
R version 4.2.3 (2023-03-15 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)
```

The package versions and their dependencies are recorded with the R
package [renv](https://rstudio.github.io/renv/) and stored in the file
[renv.lock](https://github.com/shfischer/GA_MSE_HR/blob/harvest_rate_R4.2/renv.lock).
The exact package version can be restored by cloning this repository,
checking out the current branch (`harvest_rate_4.2`), navigating into
this folder in R (or setting up a project), installing the renv package

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
uses the [FLR packages](https://flr-project.org/) `FLCore`, `FLash`,
`FLBRP`, `ggplotFL`, `mse`. See
[renv.lock](https://github.com/shfischer/GA_MSE_HR/blob/harvest_rate_R4.2/renv.lock)
for version details and sources.

The FLR package versions can also be installed manually with `remotes`
(requires suitable tools to compile R packages):

``` r
remotes::install_github(repo = "flr/FLCore", ref = "031f7742bca4b128f80ee3cc386c21d9c00ddcf7")
remotes::install_github(repo = "flr/FLash", ref = "1008c040b760f68df5101f8aafba28d9e61bba5e", INSTALL_opts = "--no-multiarch")
# INSTALL_opts = "--no-multiarch" to avoid issues in Windows
remotes::install_github(repo = "flr/FLBRP", ref = "5fa2b741d81f3d5c3e06b40cc221c98418d3a3c7", INSTALL_opts = "--no-multiarch")
remotes::install_github(repo = "shfischer/mse", ref = "18ce8c8d51362d5ced27b6c5978912f6da0ab5dc", INSTALL_opts = "--no-multiarch")
# and the GA package for running the genetic algorithm
remotes::install_github(repo = "shfischer/GA", ref = "48c1092437629b86a5310fa2873621621ff0b0e0")
```

For using MPI parallelisation, an MPI backend such as OpenMPI and the R
packages `Rmpi` and `doMPI` are required.
