################################################################################
#################_____________________________________________##################
#################                                             ##################
#################  Lista de Fixação 6: Escore de Fisher e EM  ##################
#################_____________________________________________##################
#################                                             ##################
################################################################################

# Semestre: 2/2021
# Monitor: Rodrigo Ferrari Lucas Lassance
# Lista: https://www.overleaf.com/project/617c46fb34d1661acc0f77c3


###############################################################################
#######################   Questao 1 - Escore de Fisher  #######################
###############################################################################

##(a)
set.seed(42)
n <- 30
theta <- 0.5

X <- rexp(n, 1)
Y <- rpois(n, theta*X)

##(b)
## Escore: -SOMA(X) + SOMA(Y)/theta => EMV de theta = SOMA(Y)/SOMA(X)
## Matriz Informação de Fisher esperada: SOMA(X)/theta

##(c)
escore_fisher <- function(Y, X, iter, eps = 0.001, theta0){
  theta <- theta0
  n <- length(Y)
  S.X <- sum(X)
  S.Y <- sum(Y)
  
  for(i in 1:iter){
    theta <- theta + (-S.X + S.Y/theta)/(S.X/theta)
  }
  return(theta)
}

sum(Y)/sum(X)
escore_fisher(Y, X, iter = 100, theta0 = 1)

###############################################################################
#########################    Questao 2 - Lloyd e EM   #########################
###############################################################################
library(ggplot2)

##(a)
set.seed(24)
n <- 100
theta <- c(-.5, .75)
Z <- rbinom(n, 1, .5)
X <- rnorm(n, theta[Z+1])

ggplot(data = data.frame(X,Z), aes(x = X, fill = factor(Z))) + geom_density()
ggplot(data = data.frame(X,Z), aes(x = X)) + geom_density()

##(b)
lloyd <- function(X, theta0, iter){
  theta <- theta0
  for(i in 1:iter){
    Z <- abs(theta[2] - X) < abs(theta[1] - X)
    theta[1] <- sum((1-Z)*X)/sum(1-Z)
    theta[2] <- sum(Z*X)/sum(Z)
  }
  return(theta)
}
lloyd(X, theta0 = c(0, .1), iter = 10^4)

##(c)
#P(Zi = 1 | Xi, theta) = 1/((1-p)*exp(-( (Xi-theta[1])^2 - (Xi-theta[2])^2 )/(2*sigma^2))/p + 1)

##(d)
EM_norm <- function(X, theta0, p, sigma, iter){
  theta <- theta0
  
  for(i in 1:iter){
    p1 <- 1/((1-p)*exp(-( (X-theta[1])^2 - (X-theta[2])^2 )/(2*sigma^2))/p + 1)
    p0 <- 1-p1
    theta[1] <- sum(p0*X)/sum(p0)
    theta[2] <- sum(p1*X)/sum(p1)
  }
  return(theta)
}

EM_norm(X, theta0 = c(0, .1), p = .5, sigma = 1, iter = 100)



###############################################################################
########################    Questao 3 - Algoritmo EM   ########################
###############################################################################

##(a)
set.seed(420)

n <- 100
theta <- c(1, 3)

Z <- rbinom(n, 1, .5)
X <- rexp(n, theta[Z+1])


##(b)
#P(Z = 1 | X, theta) = 1 / ( (1-p)*theta[1]*exp(-(theta[1] - theta[2])*X)/(p*theta[2]) + 1 )
EM_exp <- function(X, theta0, p, iter){
  theta <- theta0
  for(i in 1:iter){
    p1 <- 1/( (1-p)*theta[1]*exp(-(theta[1] - theta[2])*X)/(p*theta[2]) + 1 )
    p0 <- 1-p1
    
    theta[1] <- sum(p0)/sum(p0*X)
    theta[2] <- sum(p1)/sum(p1*X)
  }
  return(theta)
}

EM_exp(X, theta0 = c(.1,.2), p = .5, iter = 100)
