data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  array[nsub] int<lower=0> Tsubj;
  int<lower=1> maxpump;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=0,upper=1> outcome;
}

transformed data {
  array[nsub, ntrial, maxpump] int d;
  for (j in 1:nsub) {
    for (k in 1:Tsubj[j]) {
      for (l in 1:maxpump) {
        d[j, k, l] = (l <= npumps[j, k]) ? 1 : 0;
      }
    }
  }
}

parameters {
  vector[4] mu_pr;
  vector<lower=0>[4] sigma;

  vector[nsub] phi_pr;
  vector[nsub] eta_pr;
  vector[nsub] gam_pr;
  vector[nsub] tau_pr;
}

transformed parameters {
  vector<lower=0.01, upper=0.99>[nsub] phi;
  vector<lower=0.01>[nsub] eta;
  vector<lower=0.01>[nsub] gam;
  vector<lower=0.01>[nsub] tau;

  for (j in 1:nsub) {
    phi[j] = 0.01 + 0.98 * Phi_approx(mu_pr[1] + sigma[1] * phi_pr[j]);
    eta[j] = 0.01 + exp(mu_pr[2] + sigma[2] * eta_pr[j]);
    gam[j] = 0.01 + exp(mu_pr[3] + sigma[3] * gam_pr[j]);
    tau[j] = 0.01 + exp(mu_pr[4] + sigma[4] * tau_pr[j]);
  }
}

generated quantities {
  array[nsub, ntrial] int npumps_pred;

  for (j in 1:nsub) {
    int n_succ = 0;
    int n_pump = 0;

    for (k in 1:Tsubj[j]) {
      real belief_component;
      real p_burst;
      real omega;

      if (n_pump == 0) {
        belief_component = phi[j];
      } else {
        belief_component = (phi[j] + eta[j] * n_succ) / (1 + eta[j] * n_pump);
      }

      belief_component = fmin(0.999, fmax(0.001, belief_component));
      p_burst = 1 - belief_component;
      p_burst = fmin(0.999, fmax(0.001, p_burst));

      real denom = log1m(p_burst);
      denom = fmin(-0.001, denom);
      omega = -gam[j] / denom;
      omega = fmin(1000, fmax(-1000, omega));

      int pred = 0;
      for (l in 1:maxpump) {
        real raw_pred = tau[j] * (omega - l);
        real bounded_pred = fmin(100, fmax(-100, raw_pred));
        if (bernoulli_rng(inv_logit(bounded_pred)) == 1) {
          pred += 1;
        } else {
          break;
        }
      }
      npumps_pred[j, k] = pred;

      n_succ += npumps[j, k] - outcome[j, k];
      n_pump += npumps[j, k];
    }

    for (k in (Tsubj[j] + 1):ntrial) {
      npumps_pred[j, k] = 0;
    }
  }
}

