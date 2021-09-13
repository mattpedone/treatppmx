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

- [ ] il parametro $\alpha$ della funzione R e Rcpp deve diventare un vettore di dimensioni 2. Contiene (seguendo la notazione di notex) alpha e sigma

-[ ] devo aggiungere un argomento alla funzione che sarà cohesion = {1, 2}

-[ ] calcola matrice triangolare V 1 colonna (integrazione numerica)

-[ ] calcola matrice triangolare V restante (formula ricorsiva)

-[ ] aggiungi pesi nuova coesione
