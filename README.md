# Selezione Trattamento 
## Pacchetto treatppmx

Questa è una repo per condividere il codice e il *work in progress* del progetto FMR.

Ho iniziato costruendo il modello PPMx che è presentato nel paper *Page Quintana (2018) Stat Comp*. Ho guardato molto al codice di Page nel pacchetto, spero non sia un problema. Ho aggiunto il Reuse Algorithm di *Favaro & Teh (2013) Stat Sci* e esteso il modello al caso multivariato.

Le opzioni già implementate prevedono: 

  * **cohesion function**:

    - Dirichlet process
  
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
    
  * **auxiliary variables**
  
    - parametro per scegliere numero di auxiliary variables
    
  * **Reuse**
  
    - parametro *''flag''* per scegliere se usare il reuse algorithm

Le *''estensioni''* (per punti) apportate al modello di *Page Quintana (2018)* finora sono:

  * *auxiliary parameters*: è possibile scegliere un numero m &#62; 1 (unica opzione nel codice di Page).
  * *Reuse algorithm* di *Favaro & Teh (2013) Stat Sci*
  * il modello è multivariato (Multivariate Normal - Inverse Wishart)

I files sono *impacchettati* ma `devtools::check()` da un po' WARNING e qualche NOTE, ma riferiti alla documentazione e alla funzione di *post-processing*.

### Roadmap
Per promemoria (più che altro per me) segno i prossimi steps:

  - [X] Finisci debug + risistema il codice
  - [X] script su dati bear fatto meglio
  - [ ] testa il codice su scenari Page Quintana (2018)
  - [ ] rinomina y con &eta; per evitare di fare confusione andando avanti (segui notazione file `notex`)
  - [X] fai diventare &eta; multivariata. vd appunti skype prof. (MVN-Inv Wishart)
  - [X] studia e implementa Reuse algorithm Favaro e Teh
  - [ ] Inserisci il modello nel sampling scheme MD 
  - [ ] aggiungi le **Z**
  - [ ] rendi &beta; adattivi cfr paper Griffin e Walker
  
### Questioni:

nel file `docs/adr1.md` ci stanno alcuni appunti per spiegare come sono fatti i vettori delle label e delle cardinalità dei singoli clusters.

