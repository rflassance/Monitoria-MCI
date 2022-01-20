################################################################################
########################### Caso 1: Cadeia de Markov ###########################
################################################################################

#Matriz de transicao
DemRep <- matrix(c(.4,.6,.6,.4),ncol=2)
DemRep

passo.mat <- DemRep
prod.mat <- passo.mat%*%DemRep

while(abs(passo.mat[1,1] - prod.mat[1,1]) > 1e-10){
  passo.mat <-prod.mat
  prod.mat <- passo.mat%*%DemRep
}

#Limite da matriz de transicao
prod.mat

#Estacionaria?
c(0.5,0.5) %*% DemRep



################################################################################
########################## Caso 2: MCMC da dist Beta ###########################
################################################################################
#MCMC
metropolis <- function(X, log.vero, priori, dpasso, rpasso, n.sims, inicio,
                       semente = NULL){
  cadeia <- numeric()
  cadeia[1] <- inicio
  set.seed(semente)
  for(i in 2:n.sims){
    proposta <- rpasso(cadeia[i-1])
    log.razao <- log.vero(X, proposta) - log.vero(X, cadeia[i-1]) + #Verossimilhanca
      log(priori(proposta)) - log(priori(cadeia[i-1])) + #Priori
      log(dpasso(cadeia[i-1])) - log(dpasso(proposta)) #Passo Metropolis
    u <- runif(1)
    ifelse(log.razao > log(u), cadeia[i] <- proposta, cadeia[i] <- cadeia[i-1])
  }
  return(cadeia)
}
#Log-verossimilhanca
log.vero <- function(x, p) sum(dbinom(x, size = 1, prob = p, log=T))
#Dist priori
priori <- function(p) dunif(p, min = 0, max = 1)
#Passo de Metropolis
dpasso <- priori
rpasso <- function(p) runif(1, min = 0, max = 1) #Metropolis independente

dados <- c(rep(0,17), rep(1,19))
n.sims <- 50000
cadeia <- metropolis(dados, log.vero, priori, dpasso, rpasso, n.sims,
                     inicio = 0.5, semente = 42)

#Convergiu?
plot(cadeia, type = 'l')

#Dependencia da cadeia
acf(cadeia)
#Obtendo independencia
cadeia.ind <- cadeia[seq(from = 1, to = n.sims, by = 13)]

#Comparando
curve(dbeta(x, 20, 18), from = 0, to = 1, col = 'blue')
lines(density(cadeia.ind), type = 'l')
legend("topleft", c("Beta", "Metropolis"), lty = 1, col = c("blue", "black"),
       bty = "n")



################################################################################
####################### Caso 3: MCMC da dist Logit-normal ######################
################################################################################
#Priori logit-normal sem const de prop
priori <- function(p) exp( -(log(p/(1-p)))^2/(2*0.32^2) ) / (p*(1-p))

cadeia2 <- metropolis(dados, log.vero, priori, dpasso, rpasso, n.sims,
                      inicio = 0.5, semente = 42)
plot(cadeia2, type = 'l')
acf(cadeia2)
cadeia2.ind <- cadeia2[seq(from = 1, to = n.sims, by = 20)]

#Sensibilidade da priori
curve(dbeta(x, 20, 18), from = 0, to = 1, main = "Dist Posteriori",
      ylab = "Densidade", ylim = c(0,7))
lines(density(cadeia2.ind), type = 'l', col = 'red')
legend("topleft", c("Priori Uniforme", "Priori Logit-normal"), lty = 1,
       col = c("black", "red"), bty = "n")
