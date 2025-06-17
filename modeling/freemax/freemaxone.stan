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
                     array[,,] int d,
                     vector vwin,
                     vector vloss,
                     vector omegaone, 
                     vector beta,
                     vector nmaxone,
                     vector eta
                     ) {
    
    real ll = 0;
    int ntrial = dims(outcome)[2];
    
    for (i in 1:size(slice_subj_idx)) {
      int subj = slice_subj_idx[i]; 
      
      // track omega dynamically by color
      vector[ntrial] omega;
      vector[ntrial] nmax;
      for (k in 1:ntrial) {
        if (k == 1) {
          nmax[k] = nmaxone[subj];
          omega[k] = nmaxone[subj] * omegaone[subj];
        }
        else {
            if (outcome[subj, k-1] == 1) {
              nmax[k] = nmax[k-1] * (1 - (exp(-1.0 * eta[subj] * nmax[k-1])));
              omega[k] = omega[k-1] *
                (1 - vloss[subj] * (1 - npumps[subj, k-1] * 1.0 / nmax[k]));
            } else {
              nmax[k] = nmax[k-1] * (1 + (exp(-1.0 * eta[subj] * nmax[k-1])));
              omega[k] = omega[k-1] *
                (1 + vwin[subj]  * (     npumps[subj, k-1] * 1.0 / nmax[k]));
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
  array[nsub, ntrial, maxpump] int d;
}

parameters {
  // group-level parameters by color
  real<lower=0> mu_vwin;
  real<lower=0> mu_vloss;
  real<lower=0> mu_beta;
  real<lower=0, upper=1> mu_omegaone;
  real<lower=1, upper=100> mu_nmaxone;
  real<lower=0, upper=1> mu_eta;

  // subject level parameters by color
  vector<lower=0,upper=1>[nsub] vwin;
  vector<lower=0,upper=1>[nsub] vloss;
  vector<lower=0,upper=1>[nsub] omegaone;
  vector<lower=0,upper=3>[nsub] beta;
  vector<lower=1,upper=100>[nsub] nmaxone;
  vector<lower=0,upper=1>[nsub] eta;

  // variance 
  vector<lower=0>[7] sigma;
}

model {
  // group priors
  mu_vwin       ~ normal(0.5, 0.25);
  mu_vloss      ~ normal(0.5, 0.25);
  mu_omegaone   ~ normal(0.5, 0.25);
  mu_beta       ~ normal(1, 0.5);
  mu_nmaxone    ~ normal(50, 25);
  mu_eta        ~ normal(0.5, 0.25);

  // need different variance term for each parameter
  // sigma is a vector of those variance terms
  // (1:vwin, 2:vloss, 3:omegaone, 4:beta, 5: nmaxone, 6:eta)
  sigma[1:4] ~ inv_gamma(1, 1);
  sigma[5:6] ~ exponential(1);

  // subject level
  for (i in 1:nsub) {
    vwin[i]        ~ normal(mu_vwin, sigma[1]);
    vloss[i]       ~ normal(mu_vloss, sigma[2]);
    omegaone[i]    ~ normal(mu_omegaone, sigma[3]);
    beta[i]        ~ normal(mu_beta, sigma[4]);
    nmaxone[i]     ~ normal(mu_nmaxone, sigma[5]);
    eta[i]         ~ normal(mu_eta, sigma[6]);
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
      outcome, npumps, opportunity, d,  
      vwin, vloss, omegaone, beta,
      nmaxone, eta
    );
  }
}

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] real nmax_out; 
  array[nsub, ntrial] real log_lik;
  array[nsub, ntrial] int npumps_pred;  // changed from npumps_saved and now used for predictions

  for (i in 1:nsub) {
    // We need to calculate the omega values again here
    vector[ntrial] omega;
    vector[ntrial] nmax;

    // Initialize all values
    for (k in 1:ntrial) {
      omega[k] = 0.0;
      log_lik[i,k] = 0.0;
      npumps_pred[i, k] = 0;  // initialize with 0
    }
    
    // Calculate all omega values
    for (k in 1:ntrial) {
        
      if (k == 1) {
        nmax[k] = nmaxone[i];
        omega[k] = omegaone[i] * nmaxone[i];
      }
      else {
        if (outcome[i, k-1] == 1) {
          nmax[k] = nmax[k-1] * (1 - (exp(-1.0 * eta[nsub] * nmax[k-1])));
          omega[k] = omega[k-1] *
            (1 - vloss[i] * (1 - npumps[i, k-1] * 1.0 / nmax[k]));
        } else {
          nmax[k] = nmax[k-1] * (1 + (exp(-1.0 * eta[nsub] * nmax[k-1])));
          omega[k] = omega[k-1] *
            (1 + vwin[i]  * (     npumps[i, k-1] * 1.0 / nmax[k]));
        }
      }
    }
      
    // Now calculate log_lik, set omega_out, and simulate predictions
    for (k in 1:ntrial) {
      omega_out[i,k] = omega[k];
      nmax_out[i, k] = nmax[k]; 
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
