// The Scaled Target Learning Model
// Hierarchical model
// For 3 colors and a reversal

functions {
  real partial_log_lik(array[] int slice_subj_idx,
                     int start,
                     int end,
                     array[,] int outcome,
                     array[,] int npumps,
                     array[,] int opportunity,
                     array[,] int balloon_color,
                     array[,,] int d,
                     vector b_vwin_pre,
                     vector b_vwin_post,
                     vector b_vloss_pre,
                     vector b_vloss_post,
                     vector b_beta_pre,
                     vector b_beta_post,
                     vector b_omegaone,
                     vector o_vwin_pre,
                     vector o_vwin_post,
                     vector o_vloss_pre,
                     vector o_vloss_post,
                     vector o_beta_pre,
                     vector o_beta_post,
                     vector o_omegaone,
                     vector y_vwin,
                     vector y_vloss,
                     vector y_omegaone,
                     vector y_beta,
                     vector blue_max,
                     vector orange_max,
                     vector yellow_max) {
    
    real ll = 0;
    int ntrial = dims(outcome)[2];
    
    for (i in 1:size(slice_subj_idx)) {
      int subj = slice_subj_idx[i]; 
      
      // track omega dynamically by color
      vector[ntrial] omega_blue;
      vector[ntrial] omega_orange;
      vector[ntrial] omega_yellow;

      int last_blue = 0;
      int last_orange = 0;
      int last_yellow = 0;

      for (k in 1:ntrial) {
        omega_blue[k] = 0.0;
        omega_orange[k] = 0.0;
        omega_yellow[k] = 0.0;
      }

      for (k in 1:ntrial) {
        int color = balloon_color[subj, k];
        real nmax_k = (balloon_color[subj, k] == 1) ? blue_max[subj] :
              (balloon_color[subj, k] == 2) ? orange_max[subj] :
                                             yellow_max[subj];

        // blue = 1
        if (color == 1) {
          if (last_blue == 0)
            omega_blue[k] = nmax_k * b_omegaone[subj];
          else {
            real vloss = (k < 91 ? b_vloss_pre[subj] : b_vloss_post[subj]);
            real vwin  = (k < 91 ? b_vwin_pre[subj]  : b_vwin_post[subj]);
            if (outcome[subj, last_blue] == 1) {
              omega_blue[k] = nmax_k *
                (omega_blue[last_blue] / nmax_k) *
                (1 - vloss * (1 - npumps[subj, last_blue] * 1.0 / nmax_k));
            } else {
              omega_blue[k] = nmax_k *
                (omega_blue[last_blue] / nmax_k) *
                (1 + vwin  * (     npumps[subj, last_blue] * 1.0 / nmax_k));
            }
          }
          last_blue = k;
        }

        // orange = 2
        else if (color == 2) {
          if (last_orange == 0)
            omega_orange[k] = nmax_k * o_omegaone[subj];
          else {
            real vloss = (k < 91 ? o_vloss_pre[subj] : o_vloss_post[subj]);
            real vwin  = (k < 91 ? o_vwin_pre[subj]  : o_vwin_post[subj]);
            if (outcome[subj, last_orange] == 1) {
              omega_orange[k] = nmax_k *
                (omega_orange[last_orange] / nmax_k) *
                (1 - vloss * (1 - npumps[subj, last_orange] * 1.0 / nmax_k));
            } else {
              omega_orange[k] = nmax_k *
                (omega_orange[last_orange] / nmax_k) *
                (1 + vwin  * (     npumps[subj, last_orange] * 1.0 / nmax_k));
            }
          }
          last_orange = k;
        }

        // yellow = 3
        else if (color == 3) {
          if (last_yellow == 0)
            omega_yellow[k] = nmax_k * y_omegaone[subj];
          else {
            real vloss = y_vloss[subj];
            real vwin  = y_vwin[subj];
            if (outcome[subj, last_yellow] == 1) {
              omega_yellow[k] = nmax_k *
                (omega_yellow[last_yellow] / nmax_k) *
                (1 - vloss * (1 - npumps[subj, last_yellow] * 1.0 / nmax_k));
            } else {
              omega_yellow[k] = nmax_k *
                (omega_yellow[last_yellow] / nmax_k) *
                (1 + vwin  * (     npumps[subj, last_yellow] * 1.0 / nmax_k));
            }
          }
          last_yellow = k;
        }
      }

      // log likelihood
      for (k in 1:ntrial) {
        int color = balloon_color[subj, k];
        real beta_k = (k > 90) ? (
          color == 1 ? b_beta_post[subj] :
          color == 2 ? o_beta_post[subj] :
                       y_beta[subj]
        ) : (
          color == 1 ? b_beta_pre[subj] :
          color == 2 ? o_beta_pre[subj] :
                       y_beta[subj]
        );

        real omega_k = color == 1 ? omega_blue[k] :
                       color == 2 ? omega_orange[k] :
                                    omega_yellow[k];

        int opp = opportunity[subj, k];
        
        // safety check to avoid indexing errors
        if (opp > 0) {
          row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
          
          if (dims(d)[3] >= opp) {
            ll += bernoulli_logit_lpmf(d[subj, k, 1:opp] | -beta_k * (n_idx - omega_k));
          }
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

// add in a parameter that modulates ability to differentiate color contingencies 

parameters {
  // group-level parameters by color
  real<lower=0> b_mu_vwin_pre;
  real<lower=0> b_mu_vloss_pre;
  real<lower=0> b_mu_beta_pre;
  real<lower=0> b_mu_vwin_post;
  real<lower=0> b_mu_vloss_post;
  real<lower=0> b_mu_beta_post;
  real<lower=0,upper=1> b_mu_omegaone;

  
  real<lower=0> o_mu_vwin_pre;
  real<lower=0> o_mu_vloss_pre;
  real<lower=0> o_mu_beta_pre;
  real<lower=0> o_mu_vwin_post;
  real<lower=0> o_mu_vloss_post;
  real<lower=0> o_mu_beta_post;
  real<lower=0,upper=1> o_mu_omegaone;

  
  real<lower=0> y_mu_vwin;
  real<lower=0> y_mu_vloss;
  real<lower=0> y_mu_beta;
  real<lower=0,upper=1> y_mu_omegaone;

  real<lower=0.0001, upper=100> mu_blue_max;
  real<lower=0.0001,upper=100> mu_orange_max;
  real<lower=0.0001,upper=100> mu_yellow_max;

  // need different variance term for each parameter
  // sigma is a vector of those variance terms
  // 1-7 Blue (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 8-14 Orange (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 15-18 Yellow (vwin, vloss, omegaone, beta )
  vector<lower=0>[21] sigma;
  

  // subject level parameters by color
  vector<lower=0,upper=1>[nsub] b_vwin_pre;
  vector<lower=0,upper=1>[nsub] b_vwin_post;
  vector<lower=0,upper=1>[nsub] b_vloss_pre;
  vector<lower=0,upper=1>[nsub] b_vloss_post;
  vector<lower=0,upper=3>[nsub] b_beta_pre;
  vector<lower=0,upper=3>[nsub] b_beta_post;
  vector<lower=0,upper=1>[nsub] b_omegaone;

  vector<lower=0,upper=1>[nsub] o_vwin_pre;
  vector<lower=0,upper=1>[nsub] o_vwin_post;
  vector<lower=0,upper=1>[nsub] o_vloss_pre;
  vector<lower=0,upper=1>[nsub] o_vloss_post;
  vector<lower=0,upper=3>[nsub] o_beta_pre;
  vector<lower=0,upper=3>[nsub] o_beta_post;
  vector<lower=0,upper=1>[nsub] o_omegaone;

  vector<lower=0,upper=1>[nsub] y_vwin;
  vector<lower=0,upper=1>[nsub] y_vloss;
  vector<lower=0,upper=1>[nsub] y_omegaone;
  vector<lower=0,upper=3>[nsub] y_beta;
  
  vector<lower=0.0001,upper=100>[nsub] blue_max;
  vector<lower=0.0001,upper=100>[nsub] orange_max;
  vector<lower=0.0001,upper=100>[nsub] yellow_max;

}

model {
  // group priors
  b_mu_vwin_pre   ~ normal(0, 1);
  b_mu_vloss_pre  ~ normal(0, 1);
  b_mu_beta_pre   ~ normal(0, 1);
  b_mu_vwin_post  ~ normal(0, 1);
  b_mu_vloss_post ~ normal(0, 1);
  b_mu_beta_post  ~ normal(0, 1);
  b_mu_omegaone   ~ normal(0, 1);

  o_mu_vwin_pre   ~ normal(0, 1);
  o_mu_vloss_pre  ~ normal(0, 1);
  o_mu_beta_pre   ~ normal(0, 1);
  o_mu_vwin_post  ~ normal(0, 1);
  o_mu_vloss_post ~ normal(0, 1);
  o_mu_beta_post  ~ normal(0, 1);
  o_mu_omegaone   ~ normal(0, 1);

  y_mu_vwin       ~ normal(0, 1);
  y_mu_vloss      ~ normal(0, 1);
  y_mu_omegaone   ~ normal(0, 1);
  y_mu_beta       ~ normal(0, 1);
  
  blue_max        ~ normal(0, 1);
  orange_max      ~ normal(0, 1);
  yellow_max      ~ normal(0, 1);

  sigma ~ inv_gamma(1, 1);

  // subject level
  for (i in 1:nsub) {
    b_vwin_pre[i]    ~ normal(b_mu_vwin_pre, sigma[1]);
    b_vloss_pre[i]   ~ normal(b_mu_vloss_pre, sigma[2]);
    b_beta_pre[i]    ~ normal(b_mu_beta_pre, sigma[3]);
    b_vwin_post[i]   ~ normal(b_mu_vwin_post, sigma[4]);
    b_vloss_post[i]  ~ normal(b_mu_vloss_post, sigma[5]);
    b_beta_post[i]   ~ normal(b_mu_beta_post, sigma[6]);
    b_omegaone[i]    ~ normal(b_mu_omegaone, sigma[7]);

    o_vwin_pre[i]    ~ normal(o_mu_vwin_pre, sigma[8]);
    o_vloss_pre[i]   ~ normal(o_mu_vloss_pre, sigma[9]);
    o_beta_pre[i]    ~ normal(o_mu_beta_pre, sigma[10]);
    o_vwin_post[i]   ~ normal(o_mu_vwin_post, sigma[11]);
    o_vloss_post[i]  ~ normal(o_mu_vloss_post, sigma[12]);
    o_beta_post[i]   ~ normal(o_mu_beta_post, sigma[13]);
    o_omegaone[i]    ~ normal(o_mu_omegaone, sigma[14]);

    y_vwin[i]        ~ normal(y_mu_vwin, sigma[15]);
    y_vloss[i]       ~ normal(y_mu_vloss, sigma[16]);
    y_omegaone[i]    ~ normal(y_mu_omegaone, sigma[17]);
    y_beta[i]        ~ normal(y_mu_beta, sigma[18]);
    
    blue_max[i]      ~ normal(mu_blue_max, sigma[19]);
    orange_max[i]      ~ normal(mu_orange_max, sigma[20]);
    yellow_max[i]      ~ normal(mu_yellow_max, sigma[21]);

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
      b_vwin_pre, b_vwin_post, b_vloss_pre, b_vloss_post,  
      b_beta_pre, b_beta_post, b_omegaone,
      o_vwin_pre, o_vwin_post, o_vloss_pre, o_vloss_post,
      o_beta_pre, o_beta_post, o_omegaone,
      y_vwin, y_vloss, y_omegaone, y_beta,
      blue_max, orange_max, yellow_max  
    );
  }
}


// not saving y_pred as it is way too large. will simulate PPCs post model fitting
// saving log_lik for LOOIC
// saving omega_out for model output inflations per participant

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] real log_lik;
  array[nsub, ntrial] int npumps_saved;  
  
  for (i in 1:nsub) {
    // We need to calculate the omega values again here
    vector[ntrial] omega_blue;
    vector[ntrial] omega_orange;
    vector[ntrial] omega_yellow;
    
    int last_blue = 0;
    int last_orange = 0;
    int last_yellow = 0;
    
    // Initialize all omega values
    for (k in 1:ntrial) {
      omega_blue[k] = 0.0;
      omega_orange[k] = 0.0;
      omega_yellow[k] = 0.0;
      log_lik[i,k] = 0.0;  // Initialize log_lik
    }
    
    // First calculate all omega values
    for (k in 1:ntrial) {
      int color = balloon_color[i, k];
      real nmax_k = (balloon_color[i, k] == 1) ? blue_max[i] :
              (balloon_color[i, k] == 2) ? orange_max[i] :
                                             yellow_max[i];
                                             
      npumps_saved[i, k] = npumps[i, k];
        
      // blue = 1
      if (color == 1) {
        if (last_blue == 0)
          omega_blue[k] = nmax_k * b_omegaone[i];
        else {
          real vloss = (k < 91 ? b_vloss_pre[i] : b_vloss_post[i]);
          real vwin  = (k < 91 ? b_vwin_pre[i]  : b_vwin_post[i]);
          if (outcome[i, last_blue] == 1) {
            omega_blue[k] = nmax_k *
              (omega_blue[last_blue] / nmax_k) *
              (1 - vloss * (1 - npumps[i, last_blue] * 1.0 / nmax_k));
          } else {
            omega_blue[k] = nmax_k *
              (omega_blue[last_blue] / nmax_k) *
              (1 + vwin  * (     npumps[i, last_blue] * 1.0 / nmax_k));
          }
        }
        last_blue = k;
      }
      
      // orange = 2
      else if (color == 2) {
        if (last_orange == 0)
          omega_orange[k] = nmax_k * o_omegaone[i];
        else {
          real vloss = (k < 91 ? o_vloss_pre[i] : o_vloss_post[i]);
          real vwin  = (k < 91 ? o_vwin_pre[i]  : o_vwin_post[i]);
          if (outcome[i, last_orange] == 1) {
            omega_orange[k] = nmax_k *
              (omega_orange[last_orange] / nmax_k) *
              (1 - vloss * (1 - npumps[i, last_orange] * 1.0 / nmax_k));
          } else {
            omega_orange[k] = nmax_k *
              (omega_orange[last_orange] / nmax_k) *
              (1 + vwin  * (     npumps[i, last_orange] * 1.0 / nmax_k));
          }
        }
        last_orange = k;
      }
      
      // yellow = 3
      else if (color == 3) {
        if (last_yellow == 0)
          omega_yellow[k] = nmax_k * y_omegaone[i];
        else {
          real vloss = y_vloss[i];
          real vwin  = y_vwin[i];
          if (outcome[i, last_yellow] == 1) {
            omega_yellow[k] = nmax_k *
              (omega_yellow[last_yellow] / nmax_k) *
              (1 - vloss * (1 - npumps[i, last_yellow] * 1.0 / nmax_k));
          } else {
            omega_yellow[k] = nmax_k *
              (omega_yellow[last_yellow] / nmax_k) *
              (1 + vwin  * (     npumps[i, last_yellow] * 1.0 / nmax_k));
          }
        }
        last_yellow = k;
      }
    }

    // Now calculate log_lik and set omega_out
    for (k in 1:ntrial) {
      real beta_i;
      real omega_k;

      beta_i = (k > 90) ? (
        (balloon_color[i, k] == 1) ? b_beta_post[i] :
        (balloon_color[i, k] == 2) ? o_beta_post[i] :
                                     y_beta[i]
      ) : (
        (balloon_color[i, k] == 1) ? b_beta_pre[i] :
        (balloon_color[i, k] == 2) ? o_beta_pre[i] :
                                     y_beta[i]
      );

      omega_k = (balloon_color[i,k] == 1) ? omega_blue[k] :
                (balloon_color[i,k] == 2) ? omega_orange[k] :
                                           omega_yellow[k];

      omega_out[i,k] = omega_k;

      int opp = opportunity[i,k];
      
      // Safety check to avoid array indexing issues
      if (opp > 0 && dims(d)[3] >= opp) {
        row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
        log_lik[i,k] = bernoulli_logit_lpmf(
          d[i,k,1:opp] |
          -beta_i * ( n_idx - omega_k )
        );
      }
    }
  }
}
