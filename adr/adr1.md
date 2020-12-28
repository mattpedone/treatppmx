## ADR 1

### 28/12/2020

#### How to initialize cluster counts

<!--- qui devo metterci immagine -->

##### Pseudocode

- n &in; N *number of observations*
- i, ii &in; {1, ..., n} *indices*
- nh, Si &in; R<sup>n</sup> *vectors*

> **nh** is the vector of clusters' cardinality, while **Si** is the vector of labels (indicates for each subject to which clusters it belongs to)

**for** i **in** n  
&nbsp;**for** ii **in** n  
&nbsp;&nbsp;**if** Si(i) == ii+1  
&nbsp;&nbsp;&nbsp; nh(ii) += 1

##### NB

Il vettore delle label `Si` indica i cluster non come `C++` (potrebbe creare confusione), quindi le etichette saranno inizieranno da 1. 
Per verificare se un cluster è ''popolato'' o è un singoletto userò il vettore `Si` ed il suo indice come indice:

> `if(nh(Si(i)-1) > 1)`  
> se questa condizione è vera l'elemento i-esimo appartiene ad un non-singleton