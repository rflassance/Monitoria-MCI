################################################################################
#################_____________________________________________##################
#################                                             ##################
#################   Lista de Fixação 4: MC Transdimensional   ##################
#################_____________________________________________##################
#################                                             ##################
################################################################################

#Regressao polinomial com dados normais
#Yi = soma(bj*Xi^j) + ei
#ei ~ N(0,1)
#beta ~ N(0, v2)

#Posteriori dado o modelo M - Normal multivariada
#V = (X'X + I/v2)^(-1)
#beta | Y, X, M ~ N(VX'Y, V)
#VX'Y = (X'X + I/v2)^(-1)X'Y => (X'X)^(-1)X'Y = EMQ

############################# (a) Geração dos dados ############################
gera_poly <- function(n, gerador.x = function(n) rnorm(n),
                      betas, semente){
  set.seed(semente)
  X <- gerador.x(n)
  
  X.acum = 0
  for(i in 1:length(betas)){
    X.acum = X.acum + betas[i]*X^(i-1)
  }
  
  Y <- rnorm(n, X.acum)
  
  return(data.frame(Y=Y, X=X))
}

dados <- gera_poly(100, betas = c(3,8,5), semente = 42)
lm(Y ~ X + I(X^2), dados)$coef #Testando se está gerando corretamente



############################ (b) MC transdimensional ###########################
poly_matrix <- function(X, grau){
  X_poly <- matrix(NA, nrow = length(X), ncol = grau+1)
  for(i in 0:grau){
    X_poly[,i+1] <- X^i
  }
  return(X_poly)
}

passo_Gibbs <- function(sims, Y, X_poly, v2){
  dim_X <- dim(X_poly)[2]
  chol_XX <- chol(
    t(X_poly)%*%X_poly + diag(1/v2, dim_X)
    )
  X_inv <- chol2inv(chol_XX)
  chol_inv <- chol(X_inv)
  beta_Gibbs <- matrix(NA, nrow = sims, ncol = dim_X)
  for(i in 1:sims){
    beta_Gibbs[i,] <- t(chol_inv)%*%rnorm(dim_X) + X_inv%*%t(X_poly)%*%Y
  }
  return(beta_Gibbs)
}

#Aqui, usaremos uma proposta de Gibbs para os betas sendo usados pelo modelo
#e adicionaremos um ruído branco aos outros
RJmetropolis <- function(B, #Número de simulações
                       dados, #Observações geradas (Y e X num dataframe)
                       grau, #Grau máximo do polinômio
                       v2, #Variância da priori dos betas
                       p, #Probabilidade de fazer uma proposta transdimensional
                       semente #Semente, para reproducibilidade
                       ){
  set.seed(semente)
  Y <- dados$Y
  X_poly <- poly_matrix(dados$X,grau)
  #Inicializando os parâmetros (modelo, intercepto e grau do polinômio)
  parm_RJ <- matrix(NA, nrow = B+1, ncol = grau + 2)
  parm_RJ[1,] <- c(0, lm(Y ~ -1 + X_poly)$coef)
  
  prop <- parms <- parm_RJ[1,]
  
  for(i in 2:(B+1)){
    
    #Etapa 1: Proposta
    if(runif(1) < p){
      prop[1] <- sample((0:grau)[-parms[1]-1],1) #Transdimensional
    } else{ #Interdimensional
      #Parte do modelo
      prop[2:(prop[1]+2)] <- passo_Gibbs(sims = 1, Y = Y,
                                         X_poly = X_poly[,1:(prop[1]+1), drop=F], v2 = v2)
      #Parte fora do modelo (se houver)
      if(prop[1] < 3){
        restantes <- (prop[1]+3):length(prop)
        prop[restantes] <- prop[restantes] + rnorm(length(restantes))
      }
    }
    
    #Etapa 2: Log da razão
    ##Obs1: A proposta de Gibbs não precisa entrar no cálculo
    ##Obs2: A proposta do ruído branco também pode ser desconsiderada no cálculo
    ##Obs3: Pelo fato da escolha do modelo ser uniforme, as probs são iguais
    ##Portanto, o aceite é automático, a não ser que seja um caso transdimensional
    if(prop[1] != parms[1]){
      lrazao <- log(p) - log(1-p) + 
        sum(
          dnorm(Y, X_poly[,1:(prop[1]+1), drop=F]%*%prop[2:(prop[1]+2)], log = T) -
            dnorm(Y, X_poly[,1:(parms[1]+1), drop=F]%*%parms[2:(parms[1]+2)], log = T)
        )
      if(lrazao > log(runif(1))){
        parms <- parm_RJ[i,] <- prop
      } else prop <- parm_RJ[i,] <- parms
    } else{
      parms <- parm_RJ[i,] <- prop
    }
  }
  
  return(parm_RJ)
}


colMeans(
  passo_Gibbs(sims=10000, Y = dados$Y, X_poly = poly_matrix(dados$X,2), v2=1000)
  )

cadeia <- RJmetropolis(B = 10000, dados = dados, grau = 3, v2 = 1000, p = 0.2,
                       semente = 24)
table(cadeia[,1])/10000



#################### (c) Estimativas do modelo mais provável ###################
modelo <- which.max(table(cadeia[,1]))
betas <- cadeia[which(cadeia[,1] == modelo),2:(modelo+2)]
for(i in 0:(dim(betas)[2]-1)){
  media <- round(mean(betas[,i+1]), digits = 3)
  quantis <- round(quantile(betas[,i+1], probs = c(.025, .975)), digits = 3)
  print(
    paste("beta", i, ": ", media,
          " (", quantis[1], ", ", quantis[2], ").", sep='')
    )
}



########################## (d) Valor esperado de Yn+1 ##########################

EY <- function(cadeia, x){
  #Etapa 1: Probabilidade de cada modelo
  prob_mod <- table(cadeia[,1])/dim(cadeia)[1]
  
  #Inicializando valores
  prev_Y <- 0
  mods <- as.numeric(names(prob_mod))
  poly_X <- poly_matrix(x, max(mods))
  for(mod in mods){
    betas <- colMeans(cadeia[which(cadeia[,1] == mod),2:(mod+2),drop=F])
    Xbeta <- poly_X[,1:(mod+1),drop=F]%*%betas
    prev_Y <- prev_Y + Xbeta*prob_mod[[as.character(mod)]]
  }
  return(prev_Y)
}

x = 2.1

c(EY(cadeia, x), 3 + 8*x + 5*x^2)
