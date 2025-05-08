// The Scaled Target Learning Model
// Hierarchical model
// For 3 colors and a reversal

data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  array[nsub, ntrial] int nmax;  // color-specific maximum values
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

  
  
  // need different variance term for each parameter
  // sigma is a vector of those variance terms
  // 1-7 Blue (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 8-14 Orange (vwin_pre, vloss_pre, omegaone, beta, vwin_post, vloss_post)
  // 15-18 Yellow (vwin, vloss, omegaone, beta )
  vector<lower=0>[18] sigma;
  

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
}
    
    
    // need conditional color selection "k" in balloon_color in the model
    // declaring omegas (optimised pumps that evolve over time deterministically) 
    // in transformed parameters for faster computing

transformed parameters {
  array[nsub] vector[ntrial] omega_blue;
  array[nsub] vector[ntrial] omega_orange;
  array[nsub] vector[ntrial] omega_yellow;

  for (i in 1:nsub) {
    // indices of the last trial of each color
    int last_blue   = 0;
    int last_orange = 0;
    int last_yellow = 0;

    for (k in 1:ntrial) {
      // --- BLUE ---
      if (balloon_color[i,k] == 1) {
        if (last_blue == 0) {
          // first-ever blue
          omega_blue[i,k] = nmax[i,k] * b_omegaone[i];
        } else {
          // update outcome and therefore learning rate/omega update from the last blue trial
          real vloss = (k < 91 ? b_vloss_pre[i]  : b_vloss_post[i]);
          real vwin  = (k < 91 ? b_vwin_pre[i]   : b_vwin_post[i]);
          if (outcome[i, last_blue] == 1) {
            omega_blue[i,k] = nmax[i,k] *
              (omega_blue[i,last_blue] / nmax[i,k]) *
              (1 - vloss * (1 - npumps[i,last_blue] / nmax[i,k]));
          } else {
            omega_blue[i,k] = nmax[i,k] *
              (omega_blue[i,last_blue] / nmax[i,k]) *
              (1 + vwin  * (     npumps[i,last_blue] / nmax[i,k]));
          }
        }
        last_blue = k;
      } else {
        // not blue → carry forward
        if (k == 1)
          omega_blue[i,k] = nmax[i,k] * b_omegaone[i];
        else
          omega_blue[i,k] = omega_blue[i,k-1];
      }

      // --- ORANGE ---
      if (balloon_color[i,k] == 2) {
        if (last_orange == 0) {
          omega_orange[i,k] = nmax[i,k] * o_omegaone[i];
        } else {
          real vloss = (k < 91 ? o_vloss_pre[i]  : o_vloss_post[i]);
          real vwin  = (k < 91 ? o_vwin_pre[i]   : o_vwin_post[i]);
          if (outcome[i, last_orange] == 1) {
            omega_orange[i,k] = nmax[i,k] *
              (omega_orange[i,last_orange] / nmax[i,k]) *
              (1 - vloss * (1 - npumps[i,last_orange] / nmax[i,k]));
          } else {
            omega_orange[i,k] = nmax[i,k] *
              (omega_orange[i,last_orange] / nmax[i,k]) *
              (1 + vwin  * (     npumps[i,last_orange] / nmax[i,k]));
          }
        }
        last_orange = k;
      } else {
        if (k == 1)
          omega_orange[i,k] = nmax[i,k] * o_omegaone[i];
        else
          omega_orange[i,k] = omega_orange[i,k-1];
      }

      // --- YELLOW ---
      if (balloon_color[i,k] == 3) {
        if (last_yellow == 0) {
          omega_yellow[i,k] = nmax[i,k] * y_omegaone[i];
        } else {
          real vloss = y_vloss[i];
          real vwin  = y_vwin[i];
          if (outcome[i, last_yellow] == 1) {
            omega_yellow[i,k] = nmax[i,k] *
              (omega_yellow[i,last_yellow] / nmax[i,k]) *
              (1 - vloss * (1 - npumps[i,last_yellow] / nmax[i,k]));
          } else {
            omega_yellow[i,k] = nmax[i,k] *
              (omega_yellow[i,last_yellow] / nmax[i,k]) *
              (1 + vwin  * (     npumps[i,last_yellow] / nmax[i,k]));
          }
        }
        last_yellow = k;
      } else {
        if (k == 1)
          omega_yellow[i,k] = nmax[i,k] * y_omegaone[i];
        else
          omega_yellow[i,k] = omega_yellow[i,k-1];
      }
    }
  }
}




model {
  //priors
  b_mu_vwin_pre  ~ normal(0, 1);
  b_mu_vloss_pre  ~ normal(0, 1);
  b_mu_beta_pre  ~ normal(0, 1);
  b_mu_vwin_post  ~ normal(0, 1);
  b_mu_vloss_post  ~ normal(0, 1);
  b_mu_beta_post  ~ normal(0, 1);
  b_mu_omegaone  ~ normal(0, 1);

  
  o_mu_vwin_pre  ~ normal(0, 1);
  o_mu_vloss_pre  ~ normal(0, 1);
  o_mu_beta_pre  ~ normal(0, 1);
  o_mu_vwin_post  ~ normal(0, 1);
  o_mu_vloss_post  ~ normal(0, 1);
  o_mu_beta_post  ~ normal(0, 1);
  o_mu_omegaone  ~ normal(0, 1);


  y_mu_vwin  ~ normal(0, 1);
  y_mu_vloss  ~ normal(0, 1);
  y_mu_omegaone  ~ normal(0, 1);
  y_mu_beta  ~ normal(0, 1);
  
  sigma ~ inv_gamma(1, 1);
  
  // subject-level parameters 
  for (i in 1:nsub) {
    
    // vwin vloss pre and post reversal
    b_vwin_pre[i]    ~ normal(b_mu_vwin_pre, sigma[1]);
    b_vloss_pre[i]   ~ normal(b_mu_vloss_pre, sigma[2]);
    b_beta_pre[i]    ~ normal(b_mu_beta_pre, sigma[3]);
    b_vwin_post[i]    ~ normal(b_mu_vwin_post, sigma[4]);
    b_vloss_post[i]   ~ normal(b_mu_vloss_post, sigma[5]);
    b_beta_post[i]    ~ normal(b_mu_beta_post, sigma[6]);
    b_omegaone[i] ~ normal(b_mu_omegaone, sigma[7]); 

    
    o_vwin_pre[i]    ~ normal(o_mu_vwin_pre, sigma[8]);
    o_vloss_pre[i]   ~ normal(o_mu_vloss_pre, sigma[9]);
    o_beta_pre[i]    ~ normal(o_mu_beta_pre, sigma[10]);
    o_vwin_post[i]    ~ normal(o_mu_vwin_post, sigma[11]);
    o_vloss_post[i]   ~ normal(o_mu_vloss_post, sigma[12]);
    o_beta_post[i]    ~ normal(o_mu_beta_post, sigma[13]);
    o_omegaone[i] ~ normal(o_mu_omegaone, sigma[14]); 

    
    y_vwin[i]    ~ normal(y_mu_vwin, sigma[15]);
    y_vloss[i]   ~ normal(y_mu_vloss, sigma[16]);
    y_omegaone[i] ~ normal(y_mu_omegaone, sigma[17]); 
    y_beta[i]    ~ normal(y_mu_beta, sigma[18]);

    for (k in 1:ntrial) {
      
          real beta_i;
          real omega_k;
          
          // choose beta once per trial
          beta_i = (k > 90) ? (
            (balloon_color[i, k] == 1) ? b_beta_post[i] :
            (balloon_color[i, k] == 2) ? o_beta_post[i] :
                                         y_beta[i]
             ) : (
                (balloon_color[i, k] == 1) ? b_beta_pre[i] :
                (balloon_color[i, k] == 2) ? o_beta_pre[i] :
                                             y_beta[i]
             );

          
          // grab the ω computed in transformed parameters
          omega_k = (balloon_color[i,k] == 1) ? omega_blue[i][k] :
                    (balloon_color[i,k] == 2) ? omega_orange[i][k] :
                                                omega_yellow[i][k];
              
          {
          int opp = opportunity[i, k];               
          row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);
          d[i, k, 1:opp] ~ bernoulli_logit(
                 -beta_i * ( n_idx - omega_k ) );

      }
    }
  }
}


// not saving y_pred as it is way too large. will simulate PPCs post model fitting
// saving log_lik for LOOIC
// saving omega_out for model output inflations per participant

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] real log_lik;

  for (i in 1:nsub) {
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

      omega_k = (balloon_color[i,k] == 1) ? omega_blue[i][k] :
                 (balloon_color[i,k] == 2) ? omega_orange[i][k] :
                                             omega_yellow[i][k];

      omega_out[i,k] = omega_k;

      int opp = opportunity[i,k];
      row_vector[opp] n_idx = linspaced_row_vector(opp, 1, opp);

      log_lik[i,k] = bernoulli_logit_lpmf(
        d[i,k,1:opp] |
        -beta_i * ( n_idx - omega_k )
      );
    }
  }
}

