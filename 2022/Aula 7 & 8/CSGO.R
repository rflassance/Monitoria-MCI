library(mvtnorm) #Densidade da normal multivariada (dmvnorm)
library(MASS) #Inversao de matriz generalizada (ginv)
library(dplyr) #Manipulacao de tabelas tibble
library(readr) #Leitura de dados como tabela tibble (read_csv)
library(tidyr) #Funcoes adicionais para tabelas tibble (drop_na)
library(bayesplot) #IC de cada parametro (mcmc_intervals)

CSGO_dados <- read_csv("tb_lobby_stats_player.csv")

CSGO <- CSGO_dados %>%
  #Removendo colunas da tabela
  select(-c(idLobbyGame, idRoom, descMapName, dtCreatedAt, qtRoundsPlayed,
            qt2Kill, qt3Kill, qt4Kill, qt5Kill, qtPlusKill, vlLevel)) %>%
  #Removendo observacoes faltantes (NA)
  drop_na() %>%
  #Agrupando os dados a partir do id do jogador
  group_by(idPlayer) %>%
  #Tomando a media de cada variavel
  summarise(across(everything(), mean)) %>%
  #Removendo casos extremados (sempre ganhou, sempre perdeu)
  filter(flWinner > 0, flWinner < 1) %>%
  #Tomando loglog da media de vitorias
  mutate(loglogWin = log(-log(flWinner))) %>%
  #Removendo id do jogador e media de vitorias
  select(-c(idPlayer, flWinner))

plot(density(CSGO$loglogWin))

#Variavel resposta
Y <- CSGO$loglogWin
#Variaveis explicativas (padronizadas para evitar influencia da escala)
X <- cbind(intercept = 1, scale(as.matrix(select(CSGO, -loglogWin))))

#MC transdimensional
rjmcmc <- function(n.sims, Y, X, p, desvio = 0.1, seed = NULL){
  #Informacoes gerais
  n <- dim(X)[1] #Tamanho da amostra
  d <- dim(X)[2] #Numero de covariaveis
  XX <- t(X)%*%X
  XX_inv <- ginv(t(X)%*%X)
  XX_inv_chol <- t(chol(XX_inv))
  beta_hat <- XX_inv%*%t(X)%*%Y
  sig_hat <- t(Y - X%*%beta_hat)%*%(Y - X%*%beta_hat)/n
  #Inicializacao dos objetos
  var_regis <- matrix(NA, nrow = n.sims, ncol = d)
  var_sel <- rep(1, d)
  var_regis[1,] <- var_sel
  sigma <- as.numeric()
  sigma[1] <- sig_hat
  beta <- matrix(NA, nrow = n.sims, ncol = d)
  colnames(beta) <- colnames(X)
  beta[1,] <- beta_hat
  aceita <- 0
  
  #Inicio das iteracoes
  set.seed(seed)
  for(i in 2:n.sims){
    
    #Parte 1: Proposta
    transdim <- rbinom(1, 1, p) #prop eh transdimensional?
    var.prop <- var_sel
    if(transdim){
      ind <- sample(1:d, 1)
      #tira a var ou a coloca de volta no modelo
      var.prop[ind] <- abs(var.prop[ind] - 1)
      sig.prop <- sigma[i-1]
      beta.prop <- beta[i-1,]
    } else{
      sig.prop <- 1/rgamma(1, n/2, n/2*sig_hat) #Simular da gama inversa
      beta.prop <- as.numeric(beta[i-1,] + XX_inv_chol%*%rnorm(d, 0, sqrt(desvio)))
    }
    
    #Parte 2: Log-razao
    Yhat <- X%*%(beta[i-1,]*var_sel)
    Yhat.prop <- X%*%(beta.prop*var.prop)
    lrazao <- transdim*(log(p) - log(1-p)) + #Entra se houver troca de modelo
      sum(dnorm(Y, Yhat.prop, sqrt(sig.prop), log = T) -
            dnorm(Y, Yhat, sqrt(sigma[i-1]), log = T)) + #Verossimilhanca
      log(sigma[i-1]) - log(sig.prop) + #Priori
      (n/2 + 1)*(log(sig.prop) - log(sigma[i-1])) +
      n*sig_hat/2*(1/sig.prop - 1/sigma[i-1]) #Passo Metropolis (den da gama inversa)
    
    u <- runif(1)
    ifelse(lrazao > log(u),
           c(beta[i,] <- beta.prop, sigma[i] <- sig.prop,
             var_regis[i,] <- var_sel <- var.prop, aceita <- aceita + 1 - transdim),
           c(beta[i,] <- beta[i-1,], sigma[i] <- sigma[i-1],
             var_regis[i,] <- var_sel))
  }
  
  return(list(beta = beta*var_regis,
              sigma = sigma,
              aceita = aceita/((1-p)*n.sims)))
}

#Execucao do MC transdimensional
cadeias <- rjmcmc(n.sims = 50000, Y, X, p = 0.3, desvio = .01, seed = 42)
cadeias$aceita

#Comportamento das cadeias
par(mfrow = c(4, 7), mar = c(4, 2, 2, 2))
for(i in 1:dim(cadeias$beta)[2]){
  plot(cadeias$beta[,i], type = 'l', xlab = colnames(X)[i])
  abline(h = 0, col = 'red')
}
plot(cadeias$sigma, type = 'l', xlab = expression(sigma^2))

#Avaliacao dos fatores individuais
mcmc_intervals(cadeias$beta)
