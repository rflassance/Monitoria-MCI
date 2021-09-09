###############################################################################
###############################################################################
###########                                                         ###########
########### Lista Fixacao 2 - Metodos Computacionalmente Intensivos ###########
###########                                                         ###########
###############################################################################
###############################################################################

# Semestre: 2/2021
# Monitor: Rodrigo Ferrari Lucas Lassance
# Lista: https://www.overleaf.com/project/6123a4be8e1e62ce6f6cfcae


###############################################################################
#######################   Questao 1 - Gibbs c/ Normal   #######################
###############################################################################

##(a)
#theta2 | theta1 ~ N(theta1*rho, 1-rho^2)

##(b)
gibbs_normal <- function(M, rho, theta.ini = c(0,0)){
  
  #Criando os theta's
  theta1 <- theta2 <- numeric()
  
  #Adicionando o valor de inicialização
  theta1[1] <- theta.ini[1]
  theta2[1] <- theta.ini[2]
  
  #Gerando os novos valores via Gibbs
  for(i in 2:(M+1)){
    #theta1 -> theta2
    theta2[i] <- rnorm(1, theta1[i-1]*rho, sqrt(1 - rho^2))
    #theta2 -> theta1
    theta1[i] <- rnorm(1, theta2[i]*rho, sqrt(1 - rho^2)) #Sempre o valor mais recente
  }
  
  #Removendo os valores iniciais
  theta1 <- theta1[-1]
  theta2 <- theta2[-1]
  
  return(data.frame(theta1 = theta1, theta2 = theta2))
}

rhos <- c(.01, .5, .99) #Valores de rho de interesse
M <- 10000 #Numeros de simulacoes
cadeias <- list()
for(rho in rhos){
  Valor <- paste('rho = ', rho*100, '%', sep = '') #Nome para a lista
  cadeias[[Valor]] <- gibbs_normal(M, rho, theta.ini = c(0,0))
}

#Plot de cada configuracao
par(mfrow = c(3,2))
for(i in 1:3){
  #Convergencia de theta1
  plot(cadeias[[i]][,1], type = 'l', main = names(cadeias)[i],
       xlab = "Passos da cadeia", ylab = expression(theta[1]))
  #Convergencia de theta2
  plot(cadeias[[i]][,2], type = 'l', main = names(cadeias)[i],
       xlab = "Passos da cadeia", ylab = expression(theta[2]))
}
library(ggplot2)
library(dplyr)
#Distribuicao conjunta
conjunta <- bind_rows(cadeias)
conjunta$rho <- c(rep(.01,M), rep(.5,M), rep(.99,M))
ggplot(conjunta, aes(x=theta1, y=theta2) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlab(expression(theta[1])) +
  ylab(expression(theta[2])) +
  facet_wrap(~ rho)


##(c)
# O amostrador mais eficiente foi o terceiro, visto que
# esta mais concentrado. Isso faz sentido, ja que rho eh
# maior, portanto um parametro carrega mais informacao
# acerca do estado do outro




################################################################################
###################   Questao 2 - Gibbs hierarquico Normal   ###################
################################################################################

#(a)
##lambda | . ~ N[(mu1+mu2)/3, 1/3]
##mu1 | . ~ N[(n1*x1.barra + lambda)/(n1+1), 1/(n1+1)]
##mu2 | . ~ N[(n2*x2.barra + lambda)/(n2+1), 1/(n2+1)]


#(b)
gibbs_normal_hier <- function(M, x1.barra, x2.barra, n1, n2, parm.init = c(0,0,0)){
  lambda <- mu1 <- mu2 <- numeric()
  
  #Atribuindo os valores iniciais
  lambda[1] <- parm.init[1]
  mu1[1] <- parm.init[2]
  mu2[1] <- parm.init[3]
  
  #Atualizando os parametros
  for(i in 2:(M+1)){
    lambda[i] <- rnorm(1, (mu1[i-1]+mu2[i-1])/3, sqrt(1/3))
    mu1[i] <- rnorm(1, (n1*x1.barra + lambda[i])/(n1+1), sqrt(1/(n1+1)))
    mu2[i] <- rnorm(1, (n2*x2.barra + lambda[i])/(n2+1), sqrt(1/(n2+1)))
  }
  
  #Removendo valores iniciais
  lambda <- lambda[-1]
  mu1 <- mu1[-1]
  mu2 <- mu2[-1]
  
  return(data.frame(lambda = lambda, mu1 = mu1, mu2 = mu2))
}

x1.barra <- 10.2
x2.barra <- 5.7
n1 <- n2 <- 100
M <- 1000

cadeia <- gibbs_normal_hier(M, x1.barra, x2.barra, n1, n2)

par(mfrow = c(3,1))
##Convergencia de lambda
plot(cadeia[,1], type = 'l', main = names(cadeia)[1],
     xlab = "Passos da cadeia", ylab = expression(lambda))
##Convergencia de mu1
plot(cadeia[,2], type = 'l', main = names(cadeia)[2],
     xlab = "Passos da cadeia", ylab = expression(mu[1]))
##Convergencia de theta2
plot(cadeia[,3], type = 'l', main = names(cadeia)[3],
     xlab = "Passos da cadeia", ylab = expression(mu[2]))

##Distribuicao conjunta dos mu's
ggplot(cadeia, aes(x=mu1, y=mu2) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlab(expression(mu[1])) +
  ylab(expression(mu[2]))
cor(cadeia$mu1, cadeia$mu2) #Baixa ou nenhuma correlacao entre as amostras

#(c)
##lambda|x ~ N(m, v)
##d = 3*n1*n2 + 2*(n1+n2) + 1
##m = (n1*x1.barra*(n2+1) + n2*x2.barra*(n1+1))/d
##v = (n1+1)*(n2+1)/d

x1.barra <- 10.2
x2.barra <- 5.7
n1 <- n2 <- 100

d = 3*n1*n2 + 2*(n1+n2) + 1
m = (n1*x1.barra*(n2+1) + n2*x2.barra*(n1+1))/d
v = (n1+1)*(n2+1)/d

c(mean(cadeia$lambda[-(1:100)]), m)
c(var(cadeia$lambda[-(1:100)]), v)

par(mfrow = c(1,1))
plot(density(cadeia$lambda[-(1:100)]), main = "")
curve(dnorm(x, mean = m, sd = sqrt(v)), col='blue', add=T)
legend("topleft", c("Curva estimada", "Curva real"), lty = 1,
       col = c('black', 'blue'), bty = 'n')



###############################################################################
####################   Questao 3 - Gibbs c/ Ber/Beta/Ber   ####################
###############################################################################
#Passo 1: simular X4 e X6 de uma Ber(0.5)
#Se X4 = 0 e X6 = 0, theta0 ~ Beta(2, 4) e theta1 ~ Beta(3, 1)
#Se X4 = 0 e X6 = 1, theta0 ~ Beta(2, 3) e theta1 ~ Beta(3, 2)
#Se X4 = 1 e X6 = 0, theta0 ~ Beta(1, 4) e theta1 ~ Beta(4, 1)
#Se X4 = 1 e X6 = 1, theta0 ~ Beta(1, 3) e theta1 ~ Beta(4, 2)

#Para este caso, nao preciso de um valor inicial pois theta0 e theta1 sao
#condicionalmente independentes
Gibbs_BerBeta <- function(M){
  theta0 <- theta1 <- numeric()
  for(i in 1:M){
    X <- rbinom(2,1,.5)
    if(X[1] == 0 & X[2] == 0){
      theta0[i] <- rbeta(1,2,4)
      theta1[i] <- rbeta(1,3,1)
    } else if(X[1] == 0 & X[2] == 1){
      theta0[i] <- rbeta(1,2,3)
      theta1[i] <- rbeta(1,3,2)
    } else if(X[1] == 1 & X[2] == 0){
      theta0[i] <- rbeta(1,1,4)
      theta1[i] <- rbeta(1,4,1)
    } else{
      theta0[i] <- rbeta(1,1,3)
      theta1[i] <- rbeta(1,4,2)
    }
  }
  return(data.frame(theta0 = theta0, theta1 = theta1))
}

cadeia_theta <- Gibbs_BerBeta(M = 1000)

par(mfrow = c(2,1), mar = c(4,4,4,4))
##Convergencia de theta0
plot(cadeia_theta[,1], type = 'l', main = names(cadeia_theta)[1],
     xlab = "Passos da cadeia", ylab = expression(theta[0]))
##Convergencia de theta1
plot(cadeia_theta[,2], type = 'l', main = names(cadeia_theta)[2],
     xlab = "Passos da cadeia", ylab = expression(theta[1]))



###############################################################################
#######################   Questao 4 - Gibbs Beta/Gama   #######################
###############################################################################

MH_beta_gama <- function(M, x, parm.init = c(1,1), st.dev = c(.1,.1)){
  alfa1 <- alfa2 <- numeric()
  
  #Atribuindo os valores iniciais
  alfa1[1] <- parm.init[1]
  alfa2[1] <- parm.init[2]
  
  lvero <- function(x, a1, a2) sum(dbeta(x,a1,a2, log=T))
  
  conta <- c(0,0)
  
  #Iniciando MH
  for(i in 2:(M+1)){
    #Passo de alfa1
    alfa1.prop <- exp(rnorm(1, log(alfa1[i-1]), st.dev[1]))
    lrazao <- lvero(x, alfa1.prop, alfa2[i-1]) - lvero(x, alfa1[i-1], alfa2[i-1]) + #Verossimilhanca
      dgamma(alfa1.prop, 1, 1, log = T) - dgamma(alfa1[i-1], 1, 1, log = T) + #Priori
      log(alfa1.prop) - log(alfa1[i-1]) #Passo de MH
    ##Aceite
    ifelse(lrazao > log(runif(1)),
           c(alfa1[i] <- alfa1.prop, conta[1] <- conta[1] + 1),
           alfa1[i] <- alfa1[i-1])
    
    #Passo de alfa2
    alfa2.prop <- exp(rnorm(1, log(alfa2[i-1]), st.dev[2]))
    lrazao <- lvero(x, alfa1[i], alfa2.prop) - lvero(x, alfa1[i], alfa2[i-1]) + #Verossimilhanca
      dgamma(alfa2.prop, 1, 1, log = T) - dgamma(alfa2[i-1], 1, 1, log = T) + #Priori
      log(alfa2.prop) - log(alfa2[i-1]) #Passo de MH
    ##Aceite
    ifelse(lrazao > log(runif(1)),
           c(alfa2[i] <- alfa2.prop, conta[2] <- conta[2] + 1),
           alfa2[i] <- alfa2[i-1])
  }
  return(list(aceita = conta/M, cadeias = data.frame(alfa1, alfa2)))
}

x = c(0.21,0.35,0.23,0.27,0.31,0.29,0.24)
resultado <- MH_beta_gama(M = 10000, x, st.dev = c(.8,.8))

resultado$aceita #Checar se a taxa de aceitacao eh adequada

par(mfrow = c(2,1))
plot(resultado$cadeia$alfa1, ylab = expression(alfa[1]), xlab = "Passos de MH",
     type = 'l')
plot(resultado$cadeia$alfa2, ylab = expression(alfa[2]), xlab = "Passos de MH",
     type = 'l')



###############################################################################
#######################   Questao 5 - Gibbs p/ Imagem   #######################
###############################################################################
library(dplyr)

A <- 0:10

#Gera uma lista dos vizinhos de cada coordenada
vizinhos <- function(A){
  min.A <- min(A)
  max.A <- max(A)
  df <- data.frame(i=as.numeric(),
                   j=as.numeric(),
                   viz_i=as.numeric(),
                   viz_j=as.numeric())
  for(i in A){
    for(j in A){
      mexe.i <- c(i-1, i, i, i+1)
      mexe.j <- c(j, j-1, j+1, j)
      checa.viz <- which(mexe.i >= min.A & mexe.i <= max.A &
                         mexe.j >= min.A & mexe.j <= max.A)
      mexe.i <- mexe.i[checa.viz]
      mexe.j <- mexe.j[checa.viz]
      casos <- length(mexe.i)
      df <- bind_rows(df,
                      data.frame(i = rep(i, casos), j = rep(j, casos),
                                 viz_i = mexe.i, viz_j = mexe.j)
                      )
    }
  }
  return(df)
}

#Log da funcao de probabilidade sem constante de proporcionalidade
lprob_disc <- function(Z, lambda, i, j, viz, vals_Z){
  #Identificando os vizinhos
  viz_Z <- viz[which(viz[,1] == i & viz[,2] == j), 3:4]
  #Numero de vizinhos
  n_viz <- dim(viz_Z)[2]
  #Acumulando valores dos vizinhos
  soma_viz <- 0
  for(k in 1:n_viz){
    soma_viz <- soma_viz + Z[viz_Z[k,1] + 1, viz_Z[k,2] + 1]
  }
  return(-(4*lambda*n_viz)*(vals_Z - soma_viz/n_viz)^2/2)
}

#Log da funcao de densidade sem constante de proporcionalidade
lprob_cont <- function(Z, lambda, i, j, viz){
  #Identificando os vizinhos
  viz_Z <- viz[which(viz[,1] == i & viz[,2] == j), 3:4]
  #Numero de vizinhos
  n_viz <- dim(viz_Z)[1]
  #Acumulando valores dos vizinhos
  soma_viz <- 0
  for(k in 1:n_viz){
    soma_viz <- soma_viz + Z[viz_Z[k,1] + 1, viz_Z[k,2] + 1]
  }
  return(-(4*lambda*n_viz)*(Z[i,j] - soma_viz/n_viz)^2/2)
}


##(a)
Gibbs_discreto <- function(M, lambda, A, vals_Z = 0:1){
  #Aspectos de A
  tam.A <- length(A)
  viz <- vizinhos(A)
  #Configuracoes iniciais
  Z <- array(NA, c(tam.A, tam.A, M+1))
  Z.novos <- matrix(sample(vals_Z,tam.A^2,replace = T), nrow = tam.A, ncol = tam.A)
  Z[,,1] <- Z.novos
  #Inicializando a cadeia
  for(k in 2:(M+1)){ #Nivel da iteracao
    for(i in 1:tam.A){ #Linha da imagem
      for(j in 1:tam.A){ #Coluna da imagem
        lvero <- lprob_disc(Z = Z.novos, lambda, i-1, j-1, viz, vals_Z) #Coleta sua vero
        vero <- exp(lvero)
        probs = vero/sum(vero)
        Z.novos[i,j] <- sample(vals_Z, 1, prob = probs)
      }
    }
    Z[,,k] <- Z.novos
  }
  #Removendo o prmeiro caso
  Z <- Z[,,-1]
  return(Z)
}

m <- 5
cadeia_binaria <- Gibbs_discreto(M = m^2, lambda = 1, A = 0:10, vals_Z = 0:1)

par(mar=c(.5, .5, .5, .5), mfrow = c(m,m))
for(i in 1:dim(cadeia_binaria)[3]){
  image(cadeia_binaria[,,i], useRaster=TRUE, axes=FALSE)
}


##(b)
m <- 5
cadeia_gray <- Gibbs_discreto(M = m^2, lambda = 1, A = 0:10, vals_Z = 0:255)

par(mar=c(.5, .5, .5, .5), mfrow = c(m,m))
for(i in 1:dim(cadeia_binaria)[3]){
  image(cadeia_gray[,,i], useRaster=TRUE, axes=FALSE)
}

##(c)
Gibbs_continuo <- function(M, lambda, A){
  #Aspectos de A
  tam.A <- length(A)
  viz <- vizinhos(A)
  #Configuracoes iniciais
  Z <- array(NA, c(tam.A, tam.A, M+1))
  Z.novos <- matrix(0, nrow = tam.A, ncol = tam.A)
  Z[,,1] <- Z.novos
  #Inicializando a cadeia
  for(k in 2:(M+1)){ #Nivel da iteracao
    for(i in 1:tam.A){ #Linha da imagem
      for(j in 1:tam.A){ #Coluna da imagem
        #Identificando os vizinhos
        viz_Z <- viz[which(viz[,1] == i-1 & viz[,2] == j-1), 3:4]
        #Numero de vizinhos
        n_viz <- dim(viz_Z)[1]
        #Acumulando valores dos vizinhos
        soma_viz <- 0
        for(l in 1:n_viz){
          soma_viz <- soma_viz + Z.novos[viz_Z[l,1] + 1, viz_Z[l,2] + 1]
        }
        Z.novos[i,j] <- rnorm(1, soma_viz/n_viz, 1/sqrt(4*lambda*n_viz))
      }
    }
    Z[,,k] <- Z.novos
  }
  #Removendo o prmeiro caso
  Z <- Z[,,-1]
  return(Z)
}

n = 5
cadeia_cont <- Gibbs_continuo(M = 5^2, lambda = 1, A=0:10)

par(mar=c(.5, .5, .5, .5), mfrow = c(m,m))
for(i in 1:dim(cadeia_binaria)[3]){
  image(cadeia_cont[,,i], useRaster=TRUE, axes=FALSE)
}
