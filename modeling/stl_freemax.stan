// Freemax
// Hierarchical model
// STL adapted for nmax as free parameter

functions {
  real partial_log_lik(array[] int slice_subj_idx,
                     int start,
                     int end,
                     array[,] int outcome,
                     array[,] int npumps,
                     array[,] int opportunity,
                     array[,] int balloon_color,
                     array[,,] int d,
                     vector vwin,
                     vector vloss,
                     vector omegaone, 
                     vector beta,
                     vector blue_max,
                     vector orange_max,
                     vector yellow_max) {
    
    real ll = 0;
    int ntrial = dims(outcome)[2];
    
    for (i in 1:size(slice_subj_idx)) {
      int subj = slice_subj_idx[i]; 
      
      // track omega dynamically by color
      vector[ntrial] omega;
      
      for (k in 1:ntrial) {
        int color = balloon_color[subj, k];
        real nmax_k = (balloon_color[subj, k] == 1) ? blue_max[subj] :
              (balloon_color[subj, k] == 2) ? orange_max[subj] :
                                             yellow_max[subj];
                                             
        if (k == 1) {  // First trial initialization
            omega[k] = nmax_k * omegaone[subj];
        }
        else {
            // Use color-specific nmax for updating omega based on previous trial
            int prev_color = balloon_color[subj, k-1];
            real prev_nmax = (prev_color == 1) ? blue_max[subj] :
                            (prev_color == 2) ? orange_max[subj] :
                                               yellow_max[subj];
            
            if (outcome[subj, k-1] == 1) {
              omega[k] = omega[k-1] *
                (1 - vloss[subj] * (1 - npumps[subj, k-1] * 1.0 / prev_nmax));
            } else {
              omega[k] = omega[k-1] *
                (1 + vwin[subj]  * (     npumps[subj, k-1] * 1.0 / prev_nmax));
            }
        }
      
        real beta_k = beta[subj];
        real omega_k = omega[k];
        int opp = opportunity[subj, k];
        
        // More robust safety check to avoid indexing errors
        if (opp > 0 && opp <= dims(d)[3]) {
          row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
          ll += bernoulli_logit_lpmf(d[subj, k, 1:opp] | -beta_k * (n_idx - omega_k));
        }
      }
    }

    return ll;
  }
}

data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  int<lower=1> maxpump;  // overall maximum pump opportunity (e.g., 128)
  array[nsub, ntrial] int<lower=1, upper=3> balloon_color;
    // blue is 1
    // orange is 2
    // yellow is 3
  array[nsub, ntrial, maxpump] int d;
}

parameters {
  // group-level parameters by color
  real<lower=0> mu_vwin;
  real<lower=0> mu_vloss;
  real<lower=0> mu_beta;
  real<lower=0,upper=1> mu_omegaone;
  real<lower=1, upper=100> mu_blue_max;
  real<lower=1,upper=100> mu_orange_max;
  real<lower=1,upper=100> mu_yellow_max;

  // subject level parameters by color
  vector<lower=0,upper=1>[nsub] vwin;
  vector<lower=0,upper=1>[nsub] vloss;
  vector<lower=0,upper=1>[nsub] omegaone;
  vector<lower=0,upper=3>[nsub] beta;
  vector<lower=1,upper=100>[nsub] blue_max;
  vector<lower=1,upper=100>[nsub] orange_max;
  vector<lower=1,upper=100>[nsub] yellow_max;

  // variance 
  vector<lower=0>[7] sigma;
}

model {
  // group priors
  mu_vwin       ~ normal(0.5, 0.25);
  mu_vloss      ~ normal(0.5, 0.25);
  mu_omegaone   ~ normal(0.5, 0.25);
  mu_beta       ~ normal(1, 0.5);
  mu_blue_max   ~ normal(20, 5);
  mu_orange_max ~ normal(20, 5);
  mu_yellow_max ~ normal(20, 5);

  // need different variance term for each parameter
  // sigma is a vector of those variance terms
  // (1:vwin, 2:vloss, 3:omegaone, 4:beta, 5:blue_nmax, 6:orange_nmax, 7:yellow_nmax)
  sigma[1:4] ~ inv_gamma(1, 1);
  sigma[5:7] ~ exponential(1);

  // subject level
  for (i in 1:nsub) {
    vwin[i]        ~ normal(mu_vwin, sigma[1]);
    vloss[i]       ~ normal(mu_vloss, sigma[2]);
    omegaone[i]    ~ normal(mu_omegaone, sigma[3]);
    beta[i]        ~ normal(mu_beta, sigma[4]);
    blue_max[i]    ~ normal(mu_blue_max, sigma[5]);
    orange_max[i]  ~ normal(mu_orange_max, sigma[6]);
    yellow_max[i]  ~ normal(mu_yellow_max, sigma[7]);
  }

  // likelihood
  {
    // Create an array of indices 1:nsub for reduce_sum
    array[nsub] int subj_idx;
    for (i in 1:nsub) {
      subj_idx[i] = i;
    }
    
    target += reduce_sum(
      partial_log_lik, subj_idx, 1,
      outcome, npumps, opportunity, balloon_color, d,  
      vwin, vloss, omegaone, beta,
      blue_max, orange_max, yellow_max  
    );
  }
}

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] real log_lik;
  array[nsub, ntrial] int npumps_pred;  // changed from npumps_saved and now used for predictions

  for (i in 1:nsub) {
    // We need to calculate the omega values again here
    vector[ntrial] omega;
    
    // Initialize all values
    for (k in 1:ntrial) {
      omega[k] = 0.0;
      log_lik[i,k] = 0.0;
      npumps_pred[i, k] = 0;  // initialize with 0
    }
    
    // Calculate all omega values
    for (k in 1:ntrial) {
      int color = balloon_color[i, k];
      real nmax_k = (color == 1) ? blue_max[i] :
                    (color == 2) ? orange_max[i] :
                                   yellow_max[i];
        
      if (k == 1) {  // First trial
        omega[k] = omegaone[i] * nmax_k;
      }
      else {
        // Use color-specific nmax for updating omega based on previous trial
        int prev_color = balloon_color[i, k-1];
        real prev_nmax = (prev_color == 1) ? blue_max[i] :
                        (prev_color == 2) ? orange_max[i] :
                                           yellow_max[i];
        
        if (outcome[i, k-1] == 1) {
          omega[k] = omega[k-1] *
            (1 - vloss[i] * (1 - npumps[i, k-1] * 1.0 / prev_nmax));
        } else {
          omega[k] = omega[k-1] *
            (1 + vwin[i]  * (     npumps[i, k-1] * 1.0 / prev_nmax));
        }
      }
    }
      
    // Now calculate log_lik, set omega_out, and simulate predictions
    for (k in 1:ntrial) {
      omega_out[i,k] = omega[k];
      int opp = opportunity[i,k];

      // Simulate predicted number of pumps
      int pred_pumps = 0;
      for (pump in 1:opp) {
        real pump_prob = inv_logit(-beta[i] * (pump - omega[k]));
        int decision = bernoulli_rng(pump_prob);
        if (decision == 1) {
          pred_pumps += 1;
        } else {
          break;
        }
      }
      npumps_pred[i, k] = pred_pumps;

      // Safety check to avoid array indexing issues
      if (opp > 0 && opp <= dims(d)[3]) {
        row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
        log_lik[i,k] = bernoulli_logit_lpmf(
          d[i,k,1:opp] |
          -beta[i] * ( n_idx - omega[k] )
        );
      } else {
        log_lik[i,k] = 0.0;  // Set to 0 if invalid opportunity
      }
    }
  }
}
