# Selezione Trattamento 
## Pacchetto treatppmx

This repo is useful for install the package and to reproduce all the simulation studies.

## Simulated patients
Simulation for patient level prognostic and predictive feature

Reference: *Bayesian personalized treatment selection strategies that integrate predictive with prognostic determinants*, **Biom.J.**, *Junsheng Ma, Francesco C. Stingo and Brian P. Hobbs*

This code is all written by Junsheng MA. I am just reorganizing it for my purposes.

The R scrips of "/data/data-raw/simu152pats.R" was used to implement the simulation of patient-level prognostic and predictive features (152 simulated patients each with 92 features). 

* The input file is: `/data/data-raw/bionmfBrunet.txt`
* The script file is: `/data/data-raw/simu152pats.R`
* The output file is: `data/Simu152pats.txt`. This data set contains 152 rows and 92 columns,
and served as the basis for different simulation scenarios.

Output is generated following the same reference and the script is `/data/data-raw/simout.R`. 

* The input file is: `data/simu152pats.txt`
* The 100 replicated dataset, along with all data, are stored in `data/SimuOutSce2.rda`

## NGG

NGG is now implemented. Need to perform a simulation study to choose parameters sigma and alpha. 
Run some comparison agains Ma and DP-like cohesion function.

### check

* Warning: devo documentare dati: fai dopo aver sistemato funzione x meccanismi generatori
* Note: sub-directories of 1Mb or more:
      libs   5.7Mb
      
      https://www.r-bloggers.com/2016/04/cran-check-note-sub-directories-of-1mb-or-more-libs/

## Installation

You can install the development version from [GitHub](https://CRAN.R-project.org), using the **devtools** package with:

``` r
devtools::install_github("mattpedone/treatppmx")
```

**NOTE** that this package depends on [vweights](https://github.com/mattpedone/vweights). 

It has only been tested on a PC running Ubuntu 20.04.2 LTS.


# to-do-list

- [ ] [References in the documentation](https://cran.r-project.org/web/packages/Rdpack/vignettes/Inserting_bibtex_references.pdf)

