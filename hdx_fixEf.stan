  data {
    int<lower=1> n;
    int<lower=0> Exposure[n];
    int<lower=0> State[n];
    real Mass[n];
    }
  parameters {
    vector[3] beta;
    //vector[S] s;
    real<lower=0> sigma_e;
    //real<lower=0> sigma_s;
  }
  model {
    real mu;
  //priors
   // s ~ normal(0,sigma_s);
  //likelihood
    for (i in 1:n){
	mu=beta[1] + beta[2]*State[i] + beta[3]*Exposure[i];
        Mass[i] ~ normal(mu,sigma_e); 

 }

  }
