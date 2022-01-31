#Exemplo: Maiores companhias do Brasil
#Fonte: https://www.kaggle.com/yamqwe/brazil-largest-companiese
#Todos os valores estao em bilhoes de dolares
library(GGally)
library(dplyr)

#Leitura dos dados
dados <- read.csv("Brazil_largest_companies.csv")
head(dados)

names(dados) <- c("RankGlobal", "Companhia", "Vendas", "Lucros", "Ativos", "ValorMercado")
ggpairs(dados[,-2]) #Removendo o nome da companhia do grafico de pares

#Tomando o logaritmo das covariaveia
dados <- dados %>%
  mutate(logVendas = log(Vendas), logAtivos = log(Ativos),
         logValor = log(ValorMercado)) %>%
  select(Companhia, RankGlobal, Lucros, logVendas, logAtivos, logValor)
head(dados)

ggpairs(dados[,-1])

#Regressao via metodos frequentistas
mod <- lm(Lucros ~ RankGlobal + logVendas + logAtivos + logValor, data = dados)
summary(mod)

################################################################################
################### Parte 1: MCMC com uma taxa de aceitação ####################
################################################################################
#Que valores iniciais escolher para o MCMC?
#Como controlar a taxa de convergencia da cadeia?

metropolis <- function(Y, X, mod, lvero, lpriori, ldpasso, rpasso, n.sims, semente = NULL){
  #Valores iniciais para os parametros
  beta <- matrix(NA, nrow = n.sims, ncol = dim(X)[2])
  beta[1,] <- as.numeric(mod$coefficients)
  tau <- as.numeric()
  tau[1] <- 1/var(mod$residuals)
  aceita = 0 #Iniciando a taxa de aceitacao
  set.seed(semente)
  for(i in 2:n.sims){
    prop <- rpasso(beta[i-1,], tau[i-1])
    beta.prop <- prop[-length(prop)]
    tau.prop <- tail(prop, 1)
    log.razao <-
      lvero(Y, X, beta.prop, tau.prop) - lvero(Y, X, beta[i-1,], tau[i-1]) + #Verossimilhanca
      lpriori(beta.prop, tau.prop) - lpriori(beta[i-1,], tau[i-1]) + #Priori
      ldpasso(beta[i-1,], tau[i-1], beta.prop, tau.prop) - ldpasso(beta.prop, tau.prop, beta[i-1,], tau[i-1]) #Passo Metropolis
    u <- runif(1)
    ifelse(log.razao > log(u),
           c(beta[i,] <- beta.prop, tau[i] <- tau.prop, aceita <- aceita + 1), #TRUE
           c(beta[i,] <- beta[i-1,], tau[i] <- tau[i-1])) #FALSE
  }
  return(list(beta = beta, tau = tau, aceita = aceita/n.sims))
}

lvero <- function(Y, X, beta, tau){
  sum(dnorm(Y, X%*%beta, sqrt(1/tau), log = T))
}

lpriori <- function(beta, tau){
  dgamma(tau, 1, 0.1, log = T) + sum(dnorm(beta, 0, sqrt(1000/tau), log = T))
}

ldpasso <- function(beta, tau, media, precisao, desvio = 0.000001){
  sum(dnorm(beta, media, sqrt(desvio), log = T)) +
    dgamma(tau, precisao^2/desvio, precisao/desvio, log=T)
}

rpasso <- function(beta, tau, desvio = 0.000001){
  beta.prop <- rnorm(length(beta), beta, sqrt(desvio))
  tau.prop <- rgamma(1, tau^2/desvio, tau/desvio)
  return(c(beta.prop,tau.prop))
}

Y <- dados$Lucros
X <- as.matrix(cbind(1,dados[,c("RankGlobal", "logVendas", "logAtivos", "logValor")]))

cadeias <- metropolis(Y, X, mod, lvero, lpriori, ldpasso, rpasso, n.sims = 30000, semente = 42)
cadeias$aceita #Taxa de aceitacao

#O que dizer da convergencia das cadeias?
par(mfrow = c(dim(cadeias$beta)[2]+1,1), mar = c(2,2,2,2))
for(i in 1:dim(cadeias$beta)[2]) plot(cadeias$beta[,i], type = 'l')
plot(cadeias$tau, type = 'l')


################################################################################
###################### Parte 2: MCMC com Gibbs sampling ########################
################################################################################

Gibbs <- function(Y, X, nu, alfa, gama, n.sims, semente = NULL){
  n <- dim(X)[1] #Tamanho da amostra
  p <- dim(X)[2] #Numero de covariaveis
  #Operacoes iniciais
  XX <- t(X)%*%X
  YY <- sum(Y*Y)
  V_inv <- XX + diag(1/nu, p)
  V <- chol2inv(chol(V_inv))
  chol_V <- t(chol(V))
  m <- V%*%t(X)%*%Y
  mVm <- t(m)%*%V_inv%*%m
  alfa_post <- (n+p)/2 + alfa
  #Definindo a cadeia
  beta <- matrix(NA, nrow = n.sims, ncol = p)
  beta[1,] <- 0
  tau <- as.numeric()
  tau[1] <- 1
  #Iteracoes
  set.seed(semente)
  for(i in 2:n.sims){
    gama_post <- gama + (YY - mVm + t(beta[i-1,] - m)%*%V_inv%*%(beta[i-1,] - m))
    tau[i] <- rgamma(1, alfa_post, gama_post)
    beta[i,] <- m + chol_V%*%rnorm(p)/sqrt(tau[i])
  }
  return(list(beta = beta, tau = tau))
}

cadeia_gibbs <- Gibbs(Y, X, nu = 1000, alfa = 1, gama = 0.1, n.sims = 30000, semente = 42)

par(mfrow = c(dim(cadeia_gibbs$beta)[2]+1,1), mar = c(2,2,2,2))
for(i in 1:dim(cadeia_gibbs$beta)[2]) plot(cadeia_gibbs$beta[,i], type = 'l')
plot(cadeia_gibbs$tau, type = 'l')

apply(cadeia_gibbs$beta, 2, quantile, probs = c(0.025, 0.975))
apply(cadeias$beta, 2, quantile, probs = c(0.025,0.975))
mod$coefficients
