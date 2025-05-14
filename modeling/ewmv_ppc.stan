data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  int<lower=1> maxpump;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial, maxpump] int<lower=-1> d;
}

parameters {
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  vector[nsub] phi_pr;
  vector[nsub] eta_pr;
  vector[nsub] rho_pr;
  vector[nsub] tau_pr;
  vector[nsub] lambda_pr;
}


generated quantities {
  // Reconstruct subject-level parameters from group-level
  vector<lower=0, upper=1>[nsub] phi;
  vector<lower=0>[nsub] eta;
  vector<lower=-0.5, upper=0.5>[nsub] rho;
  vector<lower=0>[nsub] tau;
  vector<lower=0>[nsub] lambda;

  array[nsub, ntrial] int total_pumps;

  for (j in 1:nsub) {
    phi[j]    = Phi_approx(mu_pr[1] + sigma[1] * phi_pr[j]);
    eta[j]    = Phi_approx(mu_pr[2] + sigma[2] * eta_pr[j]);
    rho[j]    = 0.5 - Phi_approx(mu_pr[3] + sigma[3] * rho_pr[j]);
    tau[j]    = exp(mu_pr[4] + sigma[4] * tau_pr[j]);
    lambda[j] = exp(mu_pr[5] + sigma[5] * lambda_pr[j]);
  }

  for (j in 1:nsub) {
    int n_succ = 0;
    int n_pump = 0;
    real p_burst = phi[j];

    for (k in 1:ntrial) {
      real u_gain = 1;
      real u_stop = 0;
      real u_loss;
      real u_pump;
      real delta_u;

      int pumps_this_trial = 0;
      int max_decision = npumps[j, k] + 1 - outcome[j, k];

      for (l in 1:max_decision) {
        if (d[j, k, l] != 75) {
          u_loss = l - 1;
          u_pump = (1 - p_burst) * u_gain
                   - lambda[j] * p_burst * u_loss
                   + rho[j] * p_burst * (1 - p_burst) * square(u_gain + lambda[j] * u_loss);
          delta_u = u_pump - u_stop;

          int choice = bernoulli_logit_rng(tau[j] * delta_u);
          if (choice == 1) {
            pumps_this_trial += 1;
          } else {
            break;
          }
        }
      }

      total_pumps[j, k] = pumps_this_trial;

      n_succ += npumps[j, k] - outcome[j, k];
      n_pump += npumps[j, k];

      if (n_pump > 0) {
        p_burst = phi[j] + (1 - exp(-n_pump * eta[j])) * ((n_pump - n_succ) * 1.0 / n_pump - phi[j]);
      }
    }
  }
}
