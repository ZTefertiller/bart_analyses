data {
  int<lower=1> nsub;             // Number of subjects
  int<lower=1> ntrial;             // Maximum number of trials
  int<lower=0> Tsubj[nsub];      // Number of trials for each subject
  int<lower=2> P;             // Number of max pump + 1 ** CAUTION **
  int<lower=0> pumps[nsub, ntrial];   // Number of pump
  // int<lower=0> reward[N, T];  // Amount of rewards
  int<lower=0,upper=1> explosion[nsub, ntrial];  // Whether the balloon exploded (0 or 1)
}

transformed data {
  // Whether a subject pump the button or not (0 or 1)
  int d[nsub, ntrial, P];

  for (j in 1:nsub) {
    for (k in 1:Tsubj[j]) {
      for (l in 1:P) {
        if (l <= pumps[j, k])
          d[j, k, l] = 1;
        else
          d[j, k, l] = 0;
      }
    }
  }
}

parameters {
  // Group-level parameters
  vector[5] b_mu_pr;
  vector[5] o_mu_pr;
  vector[5] y_mu_pr;
  vector<lower=0>[15] sigma;
  // 5 sigma per color

  // Normally distributed error for Matt trick
  // blue
  vector[nsub] b_phi_pr;
  vector[nsub] b_beta_pr;
  vector[nsub] b_rho_pr;
  vector[nsub] b_tau_pr;
  vector[nsub] b_lambda_pr;
  // orange
  vector[nsub] o_phi_pr;
  vector[nsub] o_beta_pr;
  vector[nsub] o_rho_pr;
  vector[nsub] o_tau_pr;
  vector[nsub] o_lambda_pr;
  // yellow
  vector[nsub] y_phi_pr;
  vector[nsub] y_beta_pr;
  vector[nsub] y_rho_pr;
  vector[nsub] y_tau_pr;
  vector[nsub] y_lambda_pr;
  
}

transformed parameters {
  // Subject-level parameters with Matt trick
  
  //blue
  vector<lower=0,upper=1>[nsub] b_phi;
  vector<lower=0>[nsub] eta;
  vector<lower=-0.5,upper=0.5>[nsub] b_rho;
  vector<lower=0>[nsub] b_tau;
  vector<lower=0>[nsub] b_lambda;

  b_phi = Phi_approx(b_mu_pr[1] + sigma[1] * b_phi_pr);
  b_eta = Phi_approx(b_mu_pr[2] + sigma[2] * b_eta_pr);
  b_rho = 0.5 - Phi_approx(b_mu_pr[3] + sigma[3] * b_rho_pr);
  b_tau = exp(b_mu_pr[4] + sigma[4] * b_tau_pr);
  b_lambda = exp(b_mu_pr[5] + sigma[5] * b_lambda_pr);
  
  //orange
  vector<lower=0,upper=1>[nsub] o_phi;
  vector<lower=0>[nsub] eta;
  vector<lower=-0.5,upper=0.5>[nsub] o_rho;
  vector<lower=0>[nsub] o_tau;
  vector<lower=0>[nsub] o_lambda;

  o_phi = Phi_approx(o_mu_pr[1] + sigma[6] * o_phi_pr);
  o_eta = Phi_approx(o_mu_pr[2] + sigma[7] * o_eta_pr);
  o_rho = 0.5 - Phi_approx(o_mu_pr[3] + sigma[8] * o_rho_pr);
  o_tau = exp(o_mu_pr[4] + sigma[9] * o_tau_pr);
  o_lambda = exp(o_mu_pr[5] + sigma[10] * o_lambda_pr);
  
  //yellow
  vector<lower=0,upper=1>[nsub] y_phi;
  vector<lower=0>[nsub] eta;
  vector<lower=-0.5,upper=0.5>[nsub] y_rho;
  vector<lower=0>[nsub] y_tau;
  vector<lower=0>[nsub] y_lambda;

  y_phi = Phi_approx(y_mu_pr[1] + sigma[11] * y_phi_pr);
  y_eta = Phi_approx(y_mu_pr[2] + sigma[12] * y_eta_pr);
  y_rho = 0.5 - Phi_approx(y_mu_pr[3] + sigma[13] * y_rho_pr);
  y_tau = exp(y_mu_pr[4] + sigma[14] * y_tau_pr);
  y_lambda = exp(y_mu_pr[5] + sigma[15] * y_lambda_pr);
  
}

model {
  // Prior
  //blue
  b_mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2); // cauchy(0, 5);
  b_phi_pr ~ normal(0, 1);
  b_eta_pr ~ normal(0, 1);
  b_rho_pr ~ normal(0, 1);
  b_tau_pr ~ normal(0, 1);
  b_lambda_pr ~ normal(0, 1);
  
  //orange
  o_mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2); // cauchy(0, 5);
  o_phi_pr ~ normal(0, 1);
  o_eta_pr ~ normal(0, 1);
  o_rho_pr ~ normal(0, 1);
  o_tau_pr ~ normal(0, 1);
  o_lambda_pr ~ normal(0, 1);
  
  //yellow
  y_mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2); // cauchy(0, 5);
  y_phi_pr ~ normal(0, 1);
  y_eta_pr ~ normal(0, 1);
  y_rho_pr ~ normal(0, 1);
  y_tau_pr ~ normal(0, 1);
  y_lambda_pr ~ normal(0, 1);

  // Likelihood
  for (j in 1:nsub) {
    // Initialize n_succ and n_pump for a subject
    int n_succ = 0;  // Number of successful pumps
    int n_pump = 0;  // Number of total pumps
    real p_burst = phi[j];

    for (k in 1:Tsubj[j]) {
      real u_gain = 1;
      real u_loss;
      real u_pump;
      real u_stop = 0;
      real delta_u;

      for (l in 1:(pumps[j, k] + 1 - explosion[j, k])) {
        u_loss = (l - 1);

        u_pump = (1 - p_burst) * u_gain - lambda[j] * p_burst * u_loss +
        rho[j] * p_burst * (1 - p_burst) * (u_gain + lambda[j] * u_loss)^2;
        // u_stop always equals 0.

        delta_u = u_pump - u_stop;

        // Calculate likelihood with bernoulli distribution
        d[j, k, l] ~ bernoulli_logit(tau[j] * delta_u);
      }

      // Update n_succ and n_pump after each trial ends
      n_succ += pumps[j, k] - explosion[j, k];
      n_pump += pumps[j, k];

      if(n_pump>0){
        p_burst = phi[j] + (1 - exp(-n_pump * eta[j])) * ((0.0 + n_pump - n_succ) / n_pump - phi[j]);
      }
    }
  }
}

generated quantities {
  // Actual group-level mean
  real<lower=0> mu_phi = Phi_approx(mu_pr[1]);
  real<lower=0> mu_eta = Phi_approx(mu_pr[2]);
  real<lower=-0.5,upper=0.5> mu_rho = 0.5 - Phi_approx(mu_pr[3]);
  real<lower=0> mu_tau = exp(mu_pr[4]);
  real<lower=0> mu_lambda = exp(mu_pr[5]);

  // Log-likelihood for model fit
  real log_lik[nsub];

  // For posterior predictive check
  real y_pred[nsub, ntrial, P];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (j in 1:nsub)
    for (k in 1:ntrial)
      for(l in 1:P)
        y_pred[j, k, l] = -1;

  { // Local section to save time and space
    for (j in 1:nsub) {
      // Initialize n_succ and n_pump for a subject
      int n_succ = 0;  // Number of successful pumps
      int n_pump = 0;  // Number of total pumps
      real p_burst = phi[j];

      log_lik[j] = 0;

      for (k in 1:Tsubj[j]) {
        real u_gain = 1;
        real u_loss;
        real u_pump;
        real u_stop = 0;
        real delta_u;

        for (l in 1:(pumps[j, k] + 1 - explosion[j, k])) {
          // u_gain always equals r ^ rho.
          u_loss = (l - 1);

          u_pump = (1 - p_burst) * u_gain - lambda[j] * p_burst * u_loss +
          rho[j] * p_burst * (1 - p_burst) * (u_gain + lambda[j] * u_loss)^2;
          // u_stop always equals 0.

          delta_u = u_pump - u_stop;

          log_lik[j] += bernoulli_logit_lpmf(d[j, k, l] | tau[j] * delta_u);
          y_pred[j, k, l] = bernoulli_logit_rng(tau[j] * delta_u);
        }

        // Update n_succ and n_pump after each trial ends
        n_succ += pumps[j, k] - explosion[j, k];
        n_pump += pumps[j, k];

        if(n_pump>0){
          p_burst = phi[j] + (1 - exp(-n_pump * eta[j])) * ((0.0 + n_pump - n_succ) / n_pump - phi[j]);
        }
      }
    }
  }
}

