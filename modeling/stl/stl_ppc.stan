data {
  int<lower=1> nsub;
  int<lower=1> ntrial;
  int<lower=1> maxpump;
  array[nsub, ntrial] int<lower=0, upper=1> outcome;
  array[nsub, ntrial] int<lower=0> npumps;
  array[nsub, ntrial] int<lower=1> opportunity;
  array[nsub, ntrial] int nmax;
}

parameters {                               // must match the fit CSV
  real<lower=0>              mu_vwin;
  real<lower=0>              mu_vloss;
  real<lower=0>              mu_beta;
  real<lower=0, upper=1>     mu_omegaone;
  vector<lower=0>[4]         sigma;
  vector<lower=0, upper=1>[nsub] vwin;
  vector<lower=0, upper=1>[nsub] vloss;
  vector<lower=0, upper=3>[nsub] beta;
  vector<lower=0, upper=1>[nsub] omegaone;
}

model {}                     // empty; Stan ignores it in GQ runs

generated quantities {
  array[nsub, ntrial] real omega_out;
  array[nsub, ntrial] int  npumps_pred;

  for (i in 1:nsub) {
    vector[ntrial] omega;
    int last = 0;
    for (k in 1:ntrial) {
      if (last == 0)
        omega[k] = nmax[i,k] * omegaone[i];
      else if (outcome[i,last] == 1)
        omega[k] = nmax[i,k] * (omega[last]/nmax[i,last]) *
                   (1 - vloss[i]*(1 - npumps[i,last]/nmax[i,last]));
      else
        omega[k] = nmax[i,k] * (omega[last]/nmax[i,last]) *
                   (1 + vwin[i]*(npumps[i,last]/nmax[i,last]));
      last = k;
    }

    for (k in 1:ntrial) {
      omega_out[i,k] = omega[k];
      int opp   = opportunity[i,k];
      int pred  = 0;
      for (n in 1:opp) {
        if (bernoulli_rng(inv_logit(-beta[i]*(n - omega[k])))) pred += 1;
        else break;
      }
      npumps_pred[i,k] = pred;
    }
  }
}
