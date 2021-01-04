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

Le *''estensioni''* apportate al modello di *Page Quintana (2018)* finora sono:

  * *auxiliary parameters*: è possibile scegliere un numero m &#62; 1 (unica opzione nel codice di Page).
  * *Reuse algorithm* di *Favaro & Teh (2013) Stat Sci*

I files ancora non sono *impacchettati* in maniera rigorosa ma lo script `demo/myscript.R` dovrebbe essere eseguibile senza grossi problemi (`devtools::check()` da un WARNING e una NOTE, ma riferiti alla documentazione e alla cartella `docs` che però è utile per documentare il codice). 
~~Il codice **sembra** funzionare abbastanza bene e i tempi credo siano accettabili.~~ sto scrivendo uno script per fare uno studio e verificare il corretto funzionamento dei metodi.

### Roadmap
Per promemoria (più che altro per me) segno i prossimi steps:

  - [ ] Finisci debug + risistema il codice
  - [X] script su dati bear fatto meglio
  - [ ] testa il codice su scenari Page Quintana (2018)
  - [ ] rinomina y con &eta; per evitare di fare confusione andando avanti (segui notazione file `notex`)
  - [ ] fai diventare &eta; multivariata. vd appunti skype prof. (MVN-Inv Wishart)
  - [X] studia e implementa Reuse algorithm Favaro e Teh
  - [ ] Inserisci il modello nel sampling scheme MD 
  - [ ] aggiungi le **Z**
  - [ ] rendi &beta; adattivi cfr paper Griffin e Walker
  
### Questioni:
Per trasferire oggetti dalla funzione `Rcpp` a `R` ho usato una Lista. forse è brutale come soluzione, ma "sicura". Sicuramente si potrebbe fare diversamente usando i puntatori, ma da quello che ho capito è sconsigliato usarli con Rcpp. In ogni caso il guadagno computazionale in questo caso sarebbe minimo, visto che servirebbero solo per restituire gli output.

nel file `notex` &alpha; e M sono la stessa cosa? CHIEDI!!

nel file `docs/adr1.md` ci stanno alcuni appunti per spiegare come sono fatti i vettori delle label e delle cardinalità dei singoli clusters.

Studiando il codice di Page e vedendo gli appunti nell'**Appendice A** di *Page Quintana (2018) Stat Comp* sembra che venga considerato un solo *empty cluster*. Dalla discussione di *Neal (2000) JCGS* mi sembra di capire che è una valida opzione, soprattutto dal punto di vista computazionale. In questo senso non avrebbe senso usare il *Reuse algorithm* di *Favaro and Teh (2013) Stat Sci* (visto che usando un solo *auxiliary parameter*) non ci sono *discarded draws*. D'altra parte *Neal (2000) JCGS* fa vedere come aumentare il numero di *auxiliary parameters* riduca l'autocorrelazione per il numero di clusters ed i relativi *cluster specific parameters*. **Forse** per evitare i costi computazionali dei *discarded draws* *Page Quintana (2018) Stat Comp* avevano scelto di usare un solo *auxiliary parameter*, considerando il *tradeoff* con l'autocorrelazione sostenibile. Usare *Reuse algorithm* potrebbe essere vantaggioso in questo senso.

Dai **primissimi risultati** in `demo/testing1.R`sembra che l'uso di m > 1 sia meglio in termini di autocorrelazione, come ci si aspetta (*Neal (2000)*) ed è ovviamente più costoso computazionalmente. L'uso del Reuse algorithm migliora i tempi computazionali per m>1, ma ha maggiore autocorrelazione... strano? Non mi torna/non capisco perché. Anche nel paper di *Favaro e Teh* Alg 8 (nel caso cs) ha ESS minore quando usa Reuse. 

Per adattare codice:

- [X] estendere il numero di *auxiliary parameters* considerati da 1 a m
  - vedi appunti sul codice
  - studia bene ll. 841 ss.
  - implementa in nuovo branch
- [X] **controlla** 
  - non sono sicuro che sia ok. sullo script sui dati bear all'aumentare di m aumentano i clusters individuati. è dovuto al fatto che diminuisce l'autocorrelazione? ha senso?
- [X] introdurre *Reuse option*
  - segui appunti su `myppmx.cpp`
- [ ] errore strano post predictive 
  * ricontrolla ciclo su p (pp)
  * adatta m auxiliary + reuse in predictive
  * se reuse = true, m deve essere >1 altrimenti errore
- [ ] controlla posterior predictive
- [ ] **confronto** *con e senza reuse* su dati bear e su scenari *Page Quintana (2018)*
    - m=1
    - m=3
    - m=10
    - m=30 
    sulla base 
    - tempo
    - autocorrelazione numero clusters
    - autocorrelazione parametri *cluster specific*
    - *effective sample size*
    - classificazione corretta (posterior predictive)
 
 