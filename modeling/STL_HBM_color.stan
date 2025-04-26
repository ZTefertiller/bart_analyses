// The Scaled Target Learning Model
// Hierarchical model
// For 3 colors and a reversal

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
  real<lower=0> b_mu_vwin_pre;
  real<lower=0> b_mu_vloss_pre;
  real<lower=0,upper=1> b_mu_omegaone;
  real<lower=0> b_mu_beta;
  real<lower=0> b_mu_vwin_post;
  real<lower=0> b_mu_vloss_post;
  
  real<lower=0> o_mu_vwin_pre;
  real<lower=0> o_mu_vloss_pre;
  real<lower=0,upper=1> o_mu_omegaone;
  real<lower=0> o_mu_beta;
  real<lower=0> o_mu_vwin_post;
  real<lower=0> o_mu_vloss_post;
  
  real<lower=0> y_mu_vwin;
  real<lower=0> y_mu_vloss;
  real<lower=0,upper=1> y_mu_omegaone;
  real<lower=0> y_mu_beta;
  
  
  // need different variance term for each parameter
  // sigma is a vector of those variance terms
  // 1-6 Blue (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 7-12 Orange (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 13-16 Yellow (vwin, vloss, omegaone, beta )
  vector<lower=0>[16] sigma;
  

  // subject level parameters by color
  vector<lower=0,upper=1>[nsub] b_vwin_pre;
  vector<lower=0,upper=1>[nsub] b_vwin_post;
  vector<lower=0,upper=1>[nsub] b_vloss_pre;
  vector<lower=0,upper=1>[nsub] b_vloss_post;
  vector<lower=0,upper=1>[nsub] b_omegaone;
  vector<lower=0,upper=3>[nsub] b_beta;
  
  vector<lower=0,upper=1>[nsub] o_vwin_pre;
  vector<lower=0,upper=1>[nsub] o_vwin_post;
  vector<lower=0,upper=1>[nsub] o_vloss_pre;
  vector<lower=0,upper=1>[nsub] o_vloss_post;
  vector<lower=0,upper=1>[nsub] o_omegaone;
  vector<lower=0,upper=3>[nsub] o_beta;
  
  vector<lower=0,upper=1>[nsub] y_vwin;
  vector<lower=0,upper=1>[nsub] y_vloss;
  vector<lower=0,upper=1>[nsub] y_omegaone;
  vector<lower=0,upper=3>[nsub] y_beta;
}

model {
  //priors
  b_mu_vwin_pre  ~ normal(0, 1);
  b_mu_vloss_pre  ~ normal(0, 1);
  b_mu_omegaone  ~ normal(0, 1);
  b_mu_beta  ~ normal(0, 1);
  b_mu_vwin_post  ~ normal(0, 1);
  b_mu_vloss_post  ~ normal(0, 1);
  
  o_mu_vwin_pre  ~ normal(0, 1);
  o_mu_vloss_pre  ~ normal(0, 1);
  o_mu_omegaone  ~ normal(0, 1);
  o_mu_beta  ~ normal(0, 1);
  o_mu_vwin_post  ~ normal(0, 1);
  o_mu_vloss_post  ~ normal(0, 1);

  y_mu_vwin  ~ normal(0, 1);
  y_mu_vloss  ~ normal(0, 1);
  y_mu_omegaone  ~ normal(0, 1);
  y_mu_beta  ~ normal(0, 1);
  
  sigma ~ inv_gamma(1, 1);
  
  // subject-level parameters 
  for (i in 1:nsub) {
    
    // separate omegas for colors
    vector[ntrial] omega_blue;
    vector[ntrial] omega_orange;
    vector[ntrial] omega_yellow;

    // vwin vloss pre and post reversal
    b_vwin_pre[i]    ~ normal(b_mu_vwin_pre, sigma[1]);
    b_vloss_pre[i]   ~ normal(b_mu_vloss_pre, sigma[2]);
    b_omegaone[i] ~ normal(b_mu_omegaone, sigma[3]); 
    b_beta[i]    ~ normal(b_mu_beta, sigma[4]);
    b_vwin_post[i]    ~ normal(b_mu_vwin_post, sigma[5]);
    b_vloss_post[i]   ~ normal(b_mu_vloss_post, sigma[6]);
    
    o_vwin_pre[i]    ~ normal(o_mu_vwin_pre, sigma[7]);
    o_vloss_pre[i]   ~ normal(o_mu_vloss_pre, sigma[8]);
    o_omegaone[i] ~ normal(o_mu_omegaone, sigma[9]); 
    o_beta[i]    ~ normal(o_mu_beta, sigma[10]);
    o_vwin_post[i]    ~ normal(o_mu_vwin_post, sigma[11]);
    o_vloss_post[i]   ~ normal(o_mu_vloss_post, sigma[12]);
    
    y_vwin[i]    ~ normal(y_mu_vwin, sigma[13]);
    y_vloss[i]   ~ normal(y_mu_vloss, sigma[14]);
    y_omegaone[i] ~ normal(y_mu_omegaone, sigma[15]); 
    y_beta[i]    ~ normal(y_mu_beta, sigma[16]);


    for (k in 1:ntrial) {
      
      // placeholders so that it can select by color and index conditionally
      real vwin_i;
      real vloss_i;
      real omegaone_i;
      real beta_i;
      real omega_k;
        
      // conditional color selection "k" in balloon_color
          // blue is 1
          // orange is 2
          // yellow is 3
          // orange and blue are switching average burst point at trial 91
            
            
      // blue
      if (balloon_color[i, k] == 1) {  
    
        // learning rates pre/post reversal 
        if (k < 91) {
          vwin_i = b_vwin_pre[i];
          vloss_i = b_vloss_pre[i];
        } else {
          vwin_i = b_vwin_post[i];
          vloss_i = b_vloss_post[i];
        }
        
        // parameters for trial number 1 and beta values
        omegaone_i = b_omegaone[i];
        beta_i = b_beta[i];

        if (k == 1) { // first trial
        omega_blue[k] = nmax[i, k] * omegaone_i;
        
        } else if (balloon_color[i, k-1] == 1) {
          if (outcome[i, k-1] == 1) { // popped trial k-1
            omega_blue[k] = nmax[i, k] * ((omega_blue[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
          } else { // collected trial k-1
            omega_blue[k] = nmax[i, k] * ((omega_blue[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
          }
          
        } else {
          omega_blue[k] = nmax[i, k] * omegaone_i; // reset if color switched from last trial
        }
        omega_k = omega_blue[k];
      }



      // orange 
      else if (balloon_color[i, k] == 2) {  
        
        // learning rates pre/post reversal
        if (k < 91) {
          vwin_i = o_vwin_pre[i];
          vloss_i = o_vloss_pre[i];
        } else {
          vwin_i = o_vwin_post[i];
          vloss_i = o_vloss_post[i];
        }
        
        omegaone_i = o_omegaone[i];
        beta_i = o_beta[i];
        
        if (k == 1) { // first trial
          omega_orange[k] = nmax[i, k] * omegaone_i;
          
        } else if (balloon_color[i, k-1] == 2) {
          if (outcome[i, k-1] == 1) { // popped
            omega_orange[k] = nmax[i, k] * ((omega_orange[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
          } else { // collected
            omega_orange[k] = nmax[i, k] * ((omega_orange[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
          }
          
        } else {
          omega_orange[k] = nmax[i, k] * omegaone_i;
        }
        omega_k = omega_orange[k];
      }



      // yellow
      else if (balloon_color[i, k] == 3) {  
        
        // single learning rate for yellow as contingencies don't change
        vwin_i = y_vwin[i];
        vloss_i = y_vloss[i];
        omegaone_i = y_omegaone[i];
        beta_i = y_beta[i];
      
      // first trial
      if (k == 1) {
        omega_yellow[k] = nmax[i, k] * omegaone_i;
        
      } else if (balloon_color[i, k-1] == 3) {
        if (outcome[i, k-1] == 1) { // popped
          omega_yellow[k] = nmax[i, k] * ((omega_yellow[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
        } else { // collected
          omega_yellow[k] = nmax[i, k] * ((omega_yellow[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
        }
        
      } else {
        omega_yellow[k] = nmax[i, k] * omegaone_i;
      }
      omega_k = omega_yellow[k];
    }
    
    
    // finally, step 2, our bernoulli decision likelihood for the trial
      for (n in 1:opportunity[i, k]) {
        // remember at top, d is [participant, trial number, pump number]
        // omega_k was selected from the previous entry in the relevant omega color index
        d[i, k, n] ~ bernoulli_logit(-beta_i * (n - omega_k));
      }
    }
  }
}


// not saving y_pred as it is way too large. will simulate PPCs post model fitting
// saving log_lik for LOOIC
// saving omega_out for model output inflations per participant

generated quantities {
  
  array[nsub, ntrial] real omega_out;
  array[nsub] real log_lik;

  // same exact logic as in model
  for (i in 1:nsub) {
    
    log_lik[i] = 0; // except initializing log likelihood 
    
    vector[ntrial] omega_blue;
    vector[ntrial] omega_orange;
    vector[ntrial] omega_yellow;
    for (k in 1:ntrial) {
      real vwin_i;
      real vloss_i;
      real omegaone_i;
      real beta_i;
      real omega_k;
            
      // blue
      if (balloon_color[i, k] == 1) {  
        if (k < 91) {
          vwin_i = b_vwin_pre[i];
          vloss_i = b_vloss_pre[i];
        } else {
          vwin_i = b_vwin_post[i];
          vloss_i = b_vloss_post[i];
        }
        omegaone_i = b_omegaone[i];
        beta_i = b_beta[i];

        if (k == 1) { // first trial
        omega_blue[k] = nmax[i, k] * omegaone_i;
        
        } else if (balloon_color[i, k-1] == 1) {
          if (outcome[i, k-1] == 1) { // popped trial k-1
            omega_blue[k] = nmax[i, k] * ((omega_blue[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
          } else { // collected trial k-1
            omega_blue[k] = nmax[i, k] * ((omega_blue[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
          }
          
        } else {
          omega_blue[k] = nmax[i, k] * omegaone_i;
        }
        omega_k = omega_blue[k];
      }

      // orange 
      else if (balloon_color[i, k] == 2) {  
        if (k < 91) {
          vwin_i = o_vwin_pre[i];
          vloss_i = o_vloss_pre[i];
        } else {
          vwin_i = o_vwin_post[i];
          vloss_i = o_vloss_post[i];
        }
        
        omegaone_i = o_omegaone[i];
        beta_i = o_beta[i];
        
        if (k == 1) { // first trial
          omega_orange[k] = nmax[i, k] * omegaone_i;
          
        } else if (balloon_color[i, k-1] == 2) {
          if (outcome[i, k-1] == 1) { // popped
            omega_orange[k] = nmax[i, k] * ((omega_orange[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
          } else { // collected
            omega_orange[k] = nmax[i, k] * ((omega_orange[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
          }
          
        } else {
          omega_orange[k] = nmax[i, k] * omegaone_i;
        }
        omega_k = omega_orange[k];
      }

      // yellow
      else if (balloon_color[i, k] == 3) {  
        vwin_i = y_vwin[i];
        vloss_i = y_vloss[i];
        omegaone_i = y_omegaone[i];
        beta_i = y_beta[i];
      
      // first trial
      if (k == 1) {
        omega_yellow[k] = nmax[i, k] * omegaone_i;
        
      } else if (balloon_color[i, k-1] == 3) {
        if (outcome[i, k-1] == 1) { // popped
          omega_yellow[k] = nmax[i, k] * ((omega_yellow[k-1] / nmax[i, k]) * (1 - (vloss_i * (1 - (npumps[i, k-1] * 1.0 / nmax[i, k])))));
        } else { // collected
          omega_yellow[k] = nmax[i, k] * ((omega_yellow[k-1] / nmax[i, k]) * (1 + (vwin_i * (npumps[i, k-1] * 1.0 / nmax[i, k]))));
        }
        
      } else {
        omega_yellow[k] = nmax[i, k] * omegaone_i;
      }
      omega_k = omega_yellow[k];
    }
    
    
    // predicted inflations for participant i on trial k
    omega_out[i, k] = omega_k;
    
    
      for (n in 1:opportunity[i, k]) {
        log_lik[i] += bernoulli_logit_lpmf(d[i, k, n] | -beta_i * (n - omega_k));
      }
    }
  }
}

