# Selezione Trattamento 
## Pacchetto treatppmx

## NGG

NGG is now implemented. Need to perform a simulation study to choose parameters sigma and alpha. 
Run some comparison agains Ma and DP-like cohesion function.

### check

* Warning: devo documentare dati: fai dopo aver sistemato funzione x meccanismi generatori
* Note: sub-directories of 1Mb or more:
      libs   5.7Mb
      
      \url{https://www.r-bloggers.com/2016/04/cran-check-note-sub-directories-of-1mb-or-more-libs/}

## Installation

You can install the development version from [GitHub](https://CRAN.R-project.org), using the **devtools** package with:

``` r
devtools::install_github("mattpedone/treatppmx")
```

**NOTE** that this package depends on [vweights](https://github.com/mattpedone/vweights). 

It has only been tested on a PC running Ubuntu 20.04.2 LTS.


# to-do-list

- [ ] [References in the documentation](https://cran.r-project.org/web/packages/Rdpack/vignettes/Inserting_bibtex_references.pdf)

- [ ] sistema per bene tutta la documentazione (references)
- [ ] descrivi nel readme 

Queste cose falle nell'altro pacchetto:

- [ ] script per ottenere scenari
- [ ] salva scenario in `.RData`
- [ ] scipt con parallel per analisi scenario/sensitività (fissa seed)
- [ ] salva risultati in `.RData`

