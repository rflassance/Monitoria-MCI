data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real phi;
  real<lower=0> tau2;
}

model {
  for (n in 2:N)
    y[n] ~ normal(phi * y[n-1], tau2);
}

