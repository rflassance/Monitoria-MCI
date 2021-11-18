################################################################################
#################_____________________________________________##################
#################                                             ##################
#################  Lista de Fixação 5: Monte Carlo Sequencial ##################
#################_____________________________________________##################
#################                                             ##################
################################################################################



library(dplyr)
library(numbers)
#Exemplo do funcionamento da funcao modulo
mod(0:5,6)
mod(6,6)
mod(-1,6)

#Atualizacao do theta para cada etapa
theta.update <- function(theta, modulo, prob.reduz, prob.igual, prob.aumenta){
  
  #Ideia 2: sample
  mod(theta + sample(-1:1, 1, prob = c(prob.reduz, prob.igual, prob.aumenta)), modulo)
}

#Geracao dos dados
HMM <- function(N = 100, d = 2, sigma_x = 0.75, #Informacoes dos dados
                parm.space = 0:5, prob.reduz, prob.igual, prob.aumenta, #Informacoes do theta
                semente = NULL #Reproducibilidade
                ){
  set.seed(semente)
  modulo <- max(parm.space) + 1
  theta <- sample(parm.space, 1)
  X <- matrix(NA, nrow = N, ncol = d)
  X[1, ] <- rnorm(d, theta, sqrt(sigma_x))
  
  for(n in 2:N){
    theta[n] <- theta.update(theta = theta[n-1], modulo, prob.reduz, prob.igual, prob.aumenta)
    X[n, ] <- rnorm(d, theta[n], sqrt(sigma_x))
  }
  return(data.frame(X = X, theta = theta))
}

dados <- HMM(N = 100, d = 2, sigma_x = 0.75, parm.space = 0:5,
             prob.reduz = 0.4, prob.igual = 0.2, prob.aumenta = 0.4,
             semente = 42)

plot(dados$theta, type = "l", lty = 2, lwd = 2, ylim = c(min(dados), max(dados)),
     ylab = "Dados", xlab = "Tempo")
lines(dados$X.1, type = "l", col = "red")
lines(dados$X.2, type = "l", col = "blue")





################################################################################

SIR <- function(rprop, peso, d, B = 10^4){
  particles <- matrix(NA, nrow = B, ncol = d)
  for(jj in 1:d)
  {
    pesos <- rep(NA, B)
    for(ii in 1:B)
    {
      passado <- particles[ii, 1:(jj-1)]
      if(jj == 1) { passado <- c() }
      particles[ii,jj] <- rprop(passado)
      pesos[ii] <- peso(particles[ii,1:jj])
    }
    linhas = sample(1:B, B, replace = TRUE, prob = pesos)
    particles <- particles[linhas,]
  }
  particles
}

#Proposta: Uniforme ou condicional, dependendo do passado
rprop <- function(theta){
  modulo = 6
  prob.reduz = 0.4
  prob.igual = 0.2
  prob.aumenta = 0.4
  d = length(theta)
  if(d == 0) return(sample(0:(modulo-1), 1))
  return(theta.update(theta[d], modulo, prob.reduz, prob.igual, prob.aumenta))
}

#Pesos: Verossimilhanca
peso <- function(theta){
  sigma_x = 0.75
  d = length(theta)
  dnorm(X[d,1], theta[d], sqrt(sigma_x))*dnorm(X[d,2], theta[d], sqrt(sigma_x))
}

X <- dados[,1:2]
particulas <- SIR(rprop, peso, d = dim(dados)[1], B = 10^4)

colMeans(particulas[,99:100])
dados$theta[99:100]


cov(particulas[,99:100])
