// The Scaled Target Learning Model
// Hierarchical model


data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  array[nsub, ntrial] int color_max;  // trial-specific maximum values
  int<lower=1> maxpump;  // overall maximum pump opportunity (e.g., 128)
  vector[nsub] spq_z;  // centered and scaled SPQ scores
  vector[nsub] caps_z;  // centered and scaled SPQ scores
  vector[nsub] pdi_z;  // centered and scaled SPQ scores
  array[nsub, ntrial, maxpump] int d;
  // questionnaire predictors
}



parameters {
  // Group-level parameters
  real beta_spq;         // regression coefficient for SPQ
  real<lower=0> mu_vwin;
  real<lower=0> mu_vloss;
  real<lower=0,upper=1> mu_omegaone;
  real<lower=0> mu_beta;
  vector<lower=0>[4] sigma;
  
  real beta_vwin_spq;    // spq on win learning (vwin)
  real beta_vloss_spq;   // spq on loss learning (vloss)
  real beta_beta_spq;    // spq on behavioural consistency (beta)
  real beta_omegaone_spq; // spq on optimised pumps
  real beta_vwin_caps;    // caps
  real beta_vloss_caps;   
  real beta_beta_caps;    
  real beta_omegaone_caps; 
  real beta_vwin_pdi;    // pdi
  real beta_vloss_pdi;   
  real beta_beta_pdi;    
  real beta_omegaone_pdi; 

  vector<lower=0,upper=1>[nsub] vwin;
  vector<lower=0,upper=1>[nsub] vloss;
  vector<lower=0,upper=1>[nsub] omegaone;
  vector<lower=0,upper=3>[nsub] beta;
}

model {
  //priors
  mu_vwin  ~ normal(0, 1);
  mu_vloss  ~ normal(0, 1);
  mu_omegaone  ~ normal(0, 1);
  mu_beta  ~ normal(0, 1);
  sigma ~ inv_gamma(1, 1);
  beta_vwin_spq  ~ normal(0, 1);
  beta_vloss_spq ~ normal(0, 1);
  beta_beta_spq  ~ normal(0, 1);
  beta_omegaone_spq ~ normal(0, 1);
  beta_vwin_caps  ~ normal(0, 1);
  beta_vloss_caps ~ normal(0, 1);
  beta_beta_caps  ~ normal(0, 1);
  beta_omegaone_caps ~ normal(0, 1);
  beta_vwin_pdi  ~ normal(0, 1);
  beta_vloss_pdi ~ normal(0, 1);
  beta_beta_pdi  ~ normal(0, 1);
  beta_omegaone_pdi ~ normal(0, 1);
  
  //likelihood
  for (i in 1:nsub) {
    vector[ntrial] omega;

    // Subject-level parameters predicted by group intercepts and SPQ
    vwin[i]    ~ normal(mu_vwin  + beta_vwin_spq  * spq_z[i], sigma[1]);
    vloss[i]   ~ normal(mu_vloss + beta_vloss_spq * spq_z[i], sigma[2]);
    omegaone[i] ~ normal(mu_omegaone + beta_omegaone_spq * spq_z[i], sigma[3]); 
    beta[i]    ~ normal(mu_beta  + beta_beta_spq  * spq_z[i], sigma[4]);
    
    // CAPS
    vwin[i]    ~ normal(mu_vwin  + beta_vwin_caps  * caps_z[i], sigma[1]);
    vloss[i]   ~ normal(mu_vloss + beta_vloss_caps * caps_z[i], sigma[2]);
    omegaone[i] ~ normal(mu_omegaone + beta_omegaone_caps * caps_z[i], sigma[3]); 
    beta[i]    ~ normal(mu_beta  + beta_beta_caps  * caps_z[i], sigma[4]);
    
    // PDI
    vwin[i]    ~ normal(mu_vwin  + beta_vwin_pdi  * pdi_z[i], sigma[1]);
    vloss[i]   ~ normal(mu_vloss + beta_vloss_pdi * pdi_z[i], sigma[2]);
    omegaone[i] ~ normal(mu_omegaone + beta_omegaone_pdi * pdi_z[i], sigma[3]); 
    beta[i]    ~ normal(mu_beta  + beta_beta_pdi  * pdi_z[i], sigma[4]);

    for (k in 1:ntrial) {

      if (k < 2){
        omega[k] = color_max[i, k] * omegaone[i];
      }

      else{
        if (outcome[i, k-1] == 1) {
          omega[k] = color_max[i, k] * ((omega[k-1] / color_max[i, k]) * 
                    (1 - (vloss[i] * (1 - (npumps[i, k-1] * 1.0 / color_max[i, k])))));
        } else {
          omega[k] = color_max[i, k] * ((omega[k-1] / color_max[i, k]) * 
                    (1 + (vwin[i] * (npumps[i, k-1] * 1.0 / color_max[i, k]))));
        }
      }

      for (n in 1:opportunity[i,k]) {
        d[i,k,n] ~ bernoulli_logit(-beta[i] * (n - omega[k]));
      }
    }
  }
}
// 
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
//           omega[k] = color_max[i, k] * omegaone[i];
//         } else {
//           if (outcome[i, k-1] == 1) {
//             omega[k] = color_max[i, k] * ((omega[k-1] / color_max[i, k]) *
//                         (1 - (vloss[i] * (1 - (npumps[i, k-1] * 1.0 / color_max[i, k])))));
//           } else {
//             omega[k] = color_max[i, k] * ((omega[k-1] / color_max[i, k]) *
//                         (1 + (vwin[i] * (npumps[i, k-1] * 1.0 / color_max[i, k]))));
//           }
// 
//         }
//         for (n in 1:opportunity[i, k]) {
//           log_lik[i] += bernoulli_logit_lpmf(d[i, k, n] | -beta[i] * (n - omega[k]));
//         }
//         for (l in 1:color_max[i, k]) {
//           y_pred[i, k, l] = bernoulli_logit_rng(-beta[i] * (l - omega[k]));
//         }
//       }
//     }
//   }
// }
