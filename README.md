# Inference with context-dependent models of nucleotide substitution

## Description of the method

Words and math describing the method are in the subdirectory `writeups/`.

## R packages

The method is implemented in an efficient way in pure R through the use of sparse matrices,
making use of the fact that once certain structures are set up, 
the values of the relevant matrices can be computed using fast linear algebra.
R functions implementing various aspects of the method are in three packages:

- [`contextual`](contextual/) : the basic machinery for building generator matrices, transition matrices, and dealing with Tmer counts
- [`simcontext`](simcontext/) : simulation from context-dependent models
- [`contextutils`](contextutils/) : miscellaneous other utility functions used in this project

The code structure and key functions are described in [this file](model-desc.md).

To install these, run
```r
library(devtools)
install_github("petrelharp/context/contextual")
install_github("petrelharp/context/simcontext")
install_github("petrelharp/context/contextutils")
```
Simulation requires some Bioconductor packages;
see below for how to install those.

## Command-line scripts

The general strategy for inference and visualization is:

1. Configuration in `json` files
2. Computation with R scripts: run `Rscript (scriptname) --help` for options
3. Visualization using templated Rmarkdown files

These are in the [`scripts/`](scripts/) directory.
The most useful ones are:

Computation:

* [`sim-seq.R`](scripts/sim-seq.R) - simulate from the model
* [`count-seq.R`](scripts/count-seq.R) - count Tmer occurrences from simulated data
* [`make-genmat.R`](scripts/make-genmat.R) - precompute generator matrices (for a large model this can be time-consuming, but the results can be reused in subsequent steps)
* [`fit-model.R`](scripts/fit-model.R) and [`fit-tree-model.R`](scripts/fit-tree-model.R) - fit a model using maximum likelihood
* [`compute-resids.R`](scripts/compute-resids.R) - compute the difference in number of Tmers seen in data to those predicted under a fitted model
* [`mcmc-model.R`](scripts/mcmc-model.R) - find the posterior distribution of parameters using MCMC (note: can continue from previous MCMC runs)

Visualization:

* [`templated-Rmd.sh`](scripts/templated-Rmd.sh) - compile a Rmarkdown file into html using a particular data set
* [`simulation.Rmd`](scripts/simulation.Rmd) - visualize Tmer counts in a data file *and* model fit
* [`check-sim.Rmd`](scripts/check-sim.Rmd) - check simulation results match expected counts

Compilation and comparison of different models:

* [`gather-results.R`](scripts/gather-results.R) - compile information across many model fitting procedures (e.g., different window lengths)
* [`collect-params-results.R`](scripts/collect-params-results.R) - pull parameters from many json files into a table

## Example models

Full analysis pipelines, from simulation to inference and visualization,
are implemented for several example models (see writeup for descriptions).
The first few are motivated by statstical physics, not DNA,
but serve as good examples.
In each directory are shell scripts (usually `workflow.sh`) that demonstrate the workflow.

* [`tasep`](tasep/) : "totally asymmetric simple exclusion process" - two states, which are 'empty' and 'occupied' slots for moving particles
* [`ising`](ising/) : a "spin glass" model of up/down magnetic particles
* [`cpg`](cpg/) : a simple model of nucleotide evolution with an increased rate of `CG -> TG` mutations
* [`tree-cpg`](tree-cpg/) : the cpg model, but rather than inference using before-after observations, observations evolve in two sister taxa
* [`shape`](shape/) : a model of purifying selection on DNA shape


# Prerequisites:

To install the prerequisites separately:
```r
install.packages(c("expm", "mcmc", "stringdist", "optparse", "jsonlite", "ape", "rmarkdown", "ggplot2"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "IRanges", "S4Vectors"))
```
If the BioConductor installer does not work due to permissions, the following may:
```r
install.packages(c("Biostrings", "IRanges", "S4Vectors"), repos=biocinstallRepos())
```
