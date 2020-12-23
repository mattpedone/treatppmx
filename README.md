# Selezione Trattamento 
## Pacchetto treatppmx

Questa è una repo per condividere il codice e il *work in progress* del progetto FMR.

Ho iniziato costruendo il modello PPMx che è presentato nel paper *Page Quintana (2018) Stat Comp*. Ho guardato molto al codice di Page nel pacchetto, spero non sia un problema. 

Al momento il modello implementato è quello nell'equazione (15) di Page Quintana (2018), quindi la $y$ è ancora univariata. Le opzioni già implementate prevedono: 

  * **cohesion function**:

    - Dirichlet process 
    - Uniform cohesion 
  
  * **similarity function**:
  
    - Auxiliary similarity
    - Double dipper similarity
  
  * **calibration**:
  
    - no calibration
    - standardize similarity value for each covariate
    - coarsening 
  
  * **covariate**:
  
    - continue (N–N e N-NIG)
    - categoriche (DM).

I files ancora non sono *impacchettati* in maniera rigorosa ma lo script `demo/myscript.R` dovrebbe essere eseguibile senza grossi problemi (`devtools::check()` da un paio di WARNING e una NOTE, ma riferiti alla documentazione e al file `demo/00Index` che devo aggiungere). 
Il codice **sembra** funzionare abbastanza bene e i tempi credo siano accettabili. 

### Roadmap
Per promemoria (più che altro per me) segno i prossimi steps:

  - [ ] Finisci debug + risistema il codice
  - [ ] script su dati bear fatto meglio
  - [ ] testa il codice su scenari Page Quintana (2018)
  - [ ] rinomina y con $\eta$ per evitare di fare confusione andando avanti (segui notazione file `notex`)
  - [ ] fai diventare $\eta$ multivariata. vd appunti skype prof. (MVN-Inv Wishart)
  - [ ] studia e implementa Reuse algorithm Favaro e Teh
  - [ ] Inserisci il modello nel sampling scheme MD 
  - [ ] aggiungi le $\boldsymbol{Z}$
  - [ ] rendi $\boldsymbol{\beta}$ adattivi cfr paper Griffin e Walker
  
### Questioni:
Per trasferire oggetti dalla funzione Rcpp a R ho usato una Lista. forse è brutale come soluzione, ma "sicura". Sicuramente si potrebbe fare diversamente usando i puntatori, ma da quello che ho capito è sconsigliato usarli con Rcpp. In ogni caso il guadagno computazionale in questo caso sarebbe minimo, visto che servirebbero solo per restituire gli output.

nel file `notex` $\alpha$ e $M$ sono la stessa cosa? CHIEDI!!