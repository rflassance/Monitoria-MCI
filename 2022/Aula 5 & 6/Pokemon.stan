// N (escalar) - tamanho da amostra
// Y (vetor) - variavel resposta
// mu0 (escalar) - media a priori do vetor mu
data {
  int<lower=0> N;
  vector[N] Y;
  real mu0;
}
// mu (vetor) - media de cada grupo
// theta (escalar) - probabilidade de pertencer a cada grupo
// sigma (vetor) - variancia de caad grupo
parameters {
  ordered[2] mu;
  real<lower=0, upper=1> theta;
  vector<lower=0>[2] sigma;
}
// Parte 1: Calculo da log-verossimilhanca
// Parte 2: Elicitacao das prioris
model {
  for (n in 1:N){
    target += log_sum_exp(log(1 - theta) + normal_lpdf(Y[n] | mu[1], sqrt(sigma[1])),
                        log(theta) + normal_lpdf(Y[n] | mu[2], sqrt(sigma[2])));
  }
  theta ~ uniform(0,1);
  sigma ~ gamma(1,0.1);
  mu ~ normal(mu0, 10);
}

