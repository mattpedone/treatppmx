# Selezione Trattamento 
## Pacchetto treatppmx

Questa è una repo per condividere il codice e il *work in progress* del progetto FMR.

Ho iniziato costruendo il modello PPMx che è presentato nel paper *Page Quintana (2018) Stat Comp*. Ho guardato molto al codice di Page nel pacchetto, spero non sia un problema. 

Al momento il modello implementato è quello nell'equazione (15) di *Page Quintana (2018)*, quindi la y è ancora univariata. Le opzioni già implementate prevedono: 

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
  - [X] script su dati bear fatto meglio
  - [ ] testa il codice su scenari Page Quintana (2018)
  - [ ] rinomina y con &eta; per evitare di fare confusione andando avanti (segui notazione file `notex`)
  - [ ] fai diventare &eta; multivariata. vd appunti skype prof. (MVN-Inv Wishart)
  - [ ] studia e implementa Reuse algorithm Favaro e Teh
  - [ ] Inserisci il modello nel sampling scheme MD 
  - [ ] aggiungi le **Z**
  - [ ] rendi &beta; adattivi cfr paper Griffin e Walker
  
### Questioni:
Per trasferire oggetti dalla funzione `Rcpp` a `R` ho usato una Lista. forse è brutale come soluzione, ma "sicura". Sicuramente si potrebbe fare diversamente usando i puntatori, ma da quello che ho capito è sconsigliato usarli con Rcpp. In ogni caso il guadagno computazionale in questo caso sarebbe minimo, visto che servirebbero solo per restituire gli output.

nel file `notex` &alpha; e M sono la stessa cosa? CHIEDI!!

nel file `docs/adr1.md` ci stanno alcuni appunti per spiegare come sono fatti i vettori delle label e delle cardinalità dei singoli clusters.

Studiando il codice di Page e vedendo gli appunti nell'**Appendice A** di *Page Quintana (2018) Stat Comp* sembra che venga considerato un solo *empty cluster*. Dalla discussione di *Neal (2000) JCGS* mi sembra di capire che è una valida opzione, soprattutto dal punto di vista computazionale. In questo senso non avrebbe senso usare il *Reuse algorithm* di *Favaro and Teh (2013) Stat Sci* (visto che usando un solo *auxiliary parameter*) non ci sono *discarded draws*. D'altra parte *Neal (2000) JCGS* fa vedere come aumentare il numero di *auxiliary parameters* riduca l'autocorrelazione per il numero di clusters ed i relativi *cluster specific parameters*. Forse per evitare i costi computazionali dei *discarded draws* *Page Quintana (2018) Stat Comp* avevano scelto di usare un solo *auxiliary parameter*, considerando il *tradeoff* con l'autocorrelazione sostenibile. Usare *Reuse algorithm* potrebbe essere vantaggioso in questo senso.

Per adattare codice:

- [ ] estendere il numero di *auxiliary parameters* considerati da 1 a m
  - vedi appunti sul codice
  - studia bene ll. 841 ss.
  - implementa in nuovo branch
- [ ] introdurre *Reuse option*