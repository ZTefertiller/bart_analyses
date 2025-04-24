// The Scaled Target Learning Model
// Hierarchical model


data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  array[nsub, ntrial] int nmax;  // trial-specific maximum values
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
  real<lower=0> b_mu_vwin;
  real<lower=0> b_mu_vloss;
  real<lower=0,upper=1> b_mu_omegaone;
  real<lower=0> b_mu_beta;
  
  real<lower=0> y_mu_vwin;
  real<lower=0> y_mu_vloss;
  real<lower=0,upper=1> y_mu_omegaone;
  real<lower=0> y_mu_beta;
  
  real<lower=0> o_mu_vwin;
  real<lower=0> o_mu_vloss;
  real<lower=0,upper=1> o_mu_omegaone;
  real<lower=0> o_mu_beta;

  vector<lower=0>[4] sigma;
  

  // subject level parameters by color
  vector<lower=0,upper=1>[nsub] b_vwin;
  vector<lower=0,upper=1>[nsub] b_vloss;
  vector<lower=0,upper=1>[nsub] b_omegaone;
  vector<lower=0,upper=3>[nsub] b_beta;
  
  vector<lower=0,upper=1>[nsub] o_vwin;
  vector<lower=0,upper=1>[nsub] o_vloss;
  vector<lower=0,upper=1>[nsub] o_omegaone;
  vector<lower=0,upper=3>[nsub] o_beta;
  
  vector<lower=0,upper=1>[nsub] y_vwin;
  vector<lower=0,upper=1>[nsub] y_vloss;
  vector<lower=0,upper=1>[nsub] y_omegaone;
  vector<lower=0,upper=3>[nsub] y_beta;
}

model {
  //priors
  b_mu_vwin  ~ normal(0, 1);
  b_mu_vloss  ~ normal(0, 1);
  b_mu_omegaone  ~ normal(0, 1);
  b_mu_beta  ~ normal(0, 1);
  
  o_mu_vwin  ~ normal(0, 1);
  o_mu_vloss  ~ normal(0, 1);
  o_mu_omegaone  ~ normal(0, 1);
  o_mu_beta  ~ normal(0, 1);

  y_mu_vwin  ~ normal(0, 1);
  y_mu_vloss  ~ normal(0, 1);
  y_mu_omegaone  ~ normal(0, 1);
  y_mu_beta  ~ normal(0, 1);
  
  sigma ~ inv_gamma(1, 1);
  
  //likelihood
  for (i in 1:nsub) {
    vector[ntrial] omega;

    // Subject-level parameters 
    b_vwin[i]    ~ normal(b_mu_vwin, sigma[1]);
    b_vloss[i]   ~ normal(b_mu_vloss, sigma[2]);
    b_omegaone[i] ~ normal(b_mu_omegaone, sigma[3]); 
    b_beta[i]    ~ normal(b_mu_beta, sigma[4]);
    o_vwin[i]    ~ normal(o_mu_vwin, sigma[1]);
    o_vloss[i]   ~ normal(o_mu_vloss, sigma[2]);
    o_omegaone[i] ~ normal(o_mu_omegaone, sigma[3]); 
    o_beta[i]    ~ normal(o_mu_beta, sigma[4]);
    y_vwin[i]    ~ normal(y_mu_vwin, sigma[1]);
    y_vloss[i]   ~ normal(y_mu_vloss, sigma[2]);
    y_omegaone[i] ~ normal(y_mu_omegaone, sigma[3]); 
    y_beta[i]    ~ normal(y_mu_beta, sigma[4]);


    for (k in 1:ntrial) {
      
        // placeholders so that it can select by color and index conditionally
        real vwin_i;
        real vloss_i;
        real omegaone_i;
        real beta_i;
        
        // conditional color selection
        if (color_id[i, k] == 1) {
          vwin_i = b_vwin[i];
          vloss_i = b_vloss[i];
          omegaone_i = b_omegaone[i];
          beta_i = b_beta[i];
        } else if (color_id[i, k] == 2) {
          vwin_i = o_vwin[i];
          vloss_i = o_vloss[i];
          omegaone_i = o_omegaone[i];
          beta_i = o_beta[i];
        } else {
          vwin_i = y_vwin[i];
          vloss_i = y_vloss[i];
          omegaone_i = y_omegaone[i];
          beta_i = y_beta[i];
        }

      if (k < 2){
          omega[k] = nmax[i, k] * omegaone_i;
        }
  
        else{
          if (outcome[i, k-1] == 1) {
            omega[k] = nmax[i, k] * ((omega[k-1] / nmax[i, k]) * 
                      (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
          } else {
            omega[k] = nmax[i, k] * ((omega[k-1] / nmax[i, k]) * 
                      (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
          }
        }
  
        for (n in 1:opportunity[i,k]) {
          d[i,k,n] ~ bernoulli_logit(-beta_i * (n - omega[k]));
        }
      }
    }
  }










// generated quantities {
//   // Log-likelihood for model fit
//   array[nsub] real log_lik;
// 
//   // For posterior predictive check: allocate with fixed dimension maxpump
//   array[nsub, ntrial, maxpump] real y_pred;
// 
//   // Initialize y_pred with a placeholder
//   for (i in 1:nsub)
//     for (k in 1:ntrial)
//       for (l in 1:maxpump)
//         y_pred[i, k, l] = -1;
// 
//   { // Local section
//     for (i in 1:nsub) {
//       vector[ntrial] omega;
//       log_lik[i] = 0;
//       for (k in 1:ntrial) {
//         if (k < 2) {
//           omega[k] = nmax[i, k] * omegaone[i];
//         } else {
//           if (outcome[i, k-1] == 1) {
//             omega[k] = nmax[i, k] * ((omega[k-1] / nmax[i, k]) *
//                         (1 - (vloss[i] * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
//           } else {
//             omega[k] = nmax[i, k] * ((omega[k-1] / nmax[i, k]) *
//                         (1 + (vwin[i] * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
//           }
// 
//         }
//         for (n in 1:opportunity[i, k]) {
//           log_lik[i] += bernoulli_logit_lpmf(d[i, k, n] | -beta[i] * (n - omega[k]));
//         }
//         for (l in 1:nmax[i, k]) {
//           y_pred[i, k, l] = bernoulli_logit_rng(-beta[i] * (l - omega[k]));
//         }
//       }
//     }
//   }
// }
