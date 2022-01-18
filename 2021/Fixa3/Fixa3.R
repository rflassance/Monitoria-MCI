################################################################################
#################_____________________________________________##################
#################                                             ##################
#################        Lista de Fixação 3: MCMC e HMC       ##################
#################_____________________________________________##################
#################                                             ##################
################################################################################


############################# (a) Geração dos dados ############################

set.seed(42)

Y <- rnorm(1)
phi <- 0.5
tau2 = 1
n = 10^4

simulaAR <- function(n, Y, phi, tau2){
  ei <- rnorm(n-1,0,sqrt(tau2))
  for(i in 2:n){
    Y[i] <- phi*Y[i-1] + ei[i-1]
  }
  return(Y)
}

Y <- simulaAR(n, Y, phi, tau2)

plot(Y, type = 'l')
acf(Y) #Decaimento exponencial, indicando modelo AR



################################ (b) AR via MCMC ###############################


lveroAR <- function(Y, phi, tau2){
  lag.Y <- Y[1:(length(Y)-1)]
  log.vero <- dnorm(Y[1], log = T) + sum(dnorm(Y[-1], phi*lag.Y, sqrt(tau2), log=T))
}

metropolisAR <- function(M, Y, sdev = 1){
  #Inicializacao
  conta = 0
  phi <- 0
  tau2 <- .5
  for(i in 2:(M+1)){
    #Passo de phi
    phi.prop <- rnorm(1, phi[i-1], sdev)
    #Passo de tau2
    tau2.prop <- exp(rnorm(1, log(tau2[i-1]), sdev))
    #Log da razao
    lrazao <- lveroAR(Y, phi.prop, tau2.prop) - lveroAR(Y, phi[i-1], tau2[i-1]) +
      dnorm(phi.prop, log=T) - dnorm(phi[i-1], log=T) +
      dgamma(tau2.prop,1,1,log = T) - dgamma(tau2[i-1],1,1,log = T) +
      log(tau2.prop) - log(tau2[i-1])
    if(lrazao > log(runif(1))){
      phi[i] <- phi.prop
      tau2[i] <- tau2.prop
      conta = conta + 1
    } else{
      phi[i] <- phi[i-1]
      tau2[i] <- tau2[i-1]
    }
  }
  return(list(aceita = conta/M, cadeias = data.frame(phi,tau2)))
}

resultado <- metropolisAR(M = 6000, Y, sdev = .015)
resultado$aceita

par(mfrow = c(2,1), mar = c(2,5,2,2))
plot(resultado$cadeias[,1], type = 'l', main = names(resultado$cadeias)[1])
abline(h = phi, col= 'red', lty = 2)
plot(resultado$cadeias[,2], type = 'l', main = names(resultado$cadeias)[2])
abline(h = tau2, col= 'red', lty = 2)

#Burn-in
plot(resultado$cadeias[-(1:200),1], type = 'l', main = names(resultado$cadeias)[1])
abline(h = phi, col= 'red', lty = 2)
plot(resultado$cadeias[-(1:200),2], type = 'l', main = names(resultado$cadeias)[2])
abline(h = tau2, col= 'red', lty = 2)

#Autocorrelacao das cadeias
par(mar = c(2,4,4,2))
acf(resultado$cadeias[-(1:200),1]) #saltos de 15 parece bom
acf(resultado$cadeias[-(1:200),2])

#Removendo autocorrelacao
ind <- seq(from = 200, to = 6000, by = 15)
phi.MH <- resultado$cadeias[ind,1]
tau2.MH <- resultado$cadeias[ind,2]

acf(phi.MH) #Sem autocorrelacao
acf(tau2.MH)

mean(phi.MH)
quantile(phi.MH, probs = c(.025, .975))
mean(tau2.MH)
quantile(tau2.MH, probs = c(.025, .975))



############################# (c) AR via HMC (Stan) ############################
library(rstan)

aux = stan("StanAR.stan", data = list(N = n, y = Y))

plot(aux)

phi.resul <- extract(aux)$phi
mean(phi.resul)
quantile(phi.resul, probs = c(.025, .975))

tau2.resul <- extract(aux)$tau2
mean(tau2.resul)
quantile(tau2.resul, probs = c(.025, .975))



################################ (d) Comparação  ###############################

mean(phi.MH)
mean(phi.resul)
quantile(phi.MH, probs = c(.025, .975))
quantile(phi.resul, probs = c(.025, .975))

mean(tau2.MH)
mean(tau2.resul)
quantile(tau2.MH, probs = c(.025, .975))
quantile(tau2.resul, probs = c(.025, .975))
