# da fare

* modifica calculate_gamma
* scrivi funzione update beta (come update_eta)
* prende in input beta come vettore $Q\times J$
* prende in input X come matrice $n\times Q$
* crea vettore beta_tmep dim $Q$ (lo userò per una categoria alla volta)

<!---
for(jj = 0 ; jj < n_cats ; jj++){

      // fill in beta_temp and update the proposal variance
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta_temp[kk] = beta[hh];
      }

      update_beta_jj(XX, JJ, loggamma, beta_temp, inclusion_indicator,
                     prop_per_beta, mu_be, sig_be, aa_hp, bb_hp, jj);

      // update beta with beta_temp and write to file
      for(kk = 0 ; kk < n_vars ; kk++){
        hh = kk + jj * n_vars;
        beta[hh] = beta_temp[kk];
        if((gg >= burn) & (gg % thin == 0)){
          accepted_beta[hh] = accepted_beta[hh] + accepted_beta_flag[hh];
          accepted_beta_flag[hh] = 0;
          fprintf(fout_beta, "%e ", beta[hh]); // space delimited output files
        }
      }
    }
-->
