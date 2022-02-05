#Aulas 5 & 6: Modelagem de atributos de pokemon
#Fonte: https://www.kaggle.com/thiagoazen/all-pokemon-with-stats
library(rstan)

#Leitura do banco
pokemon <- read.csv("PokemonDb.csv")

summary(pokemon)
head(pokemon)

#Comportamento bimodal identificado
plot(density(pokemon$Total))
plot(density(log(pokemon$Total)))

#Informacoes para o Stan
Y <- log(pokemon$Total)
n <- length(Y)
mu0 <- mean(Y)

#Execucao do codigo Stan
aux = stan("Pokemon.stan", data = list(N = n, Y = Y, mu0 = mu0),
           iter = 40000, chains = 4, cores = 4, seed = 42, warmup = 30000,
           control = list(adapt_delta = 0.9))
plot(aux)

#Recuperacao da variavel Z
Z <- matrix(NA, nrow = n, ncol = dim(aux)[1])

cadeias <- extract(aux)

set.seed(42)
for(j in 1:dim(aux)[1]){
  d0 <- (1-cadeias$theta[j])*dnorm(Y, cadeias$mu[j,1], sqrt(cadeias$sigma[j,1]))
  d1 <- cadeias$theta[j]*dnorm(Y, cadeias$mu[j,2], sqrt(cadeias$sigma[j,2]))
  prob <- d1/(d0+d1)
  Z[,j] <- rbinom(n, 1, prob)
}

rowMeans(Z)

plot(density(log(pokemon$Total)))
abline(v = colMeans(cadeias$mu)[1], lty = 2)
abline(v = colMeans(cadeias$mu)[2], lty = 2)

#Gerando previsoes pontuais para Z (identificacao dos dois grupos)
Z.pred <- apply(Z, 1, function(x) as.numeric(names(which.max(table(x)))))

#Visualizacao inicial dos grupos
cbind(head(pokemon[Z.pred == 0, 1:2], 10), head(pokemon[Z.pred == 1, 1:2], 10))

#Quais pokemon sao mais caracteristicos de cada grupo?
head(pokemon[order(rowMeans(Z)), 1:2])
head(pokemon[order(rowMeans(Z), decreasing = T), 1:2])

#Quais sao os com maior status em cada grupo
poke1 <- pokemon[Z.pred == 0,]
ord1 <- order(poke1$Total, decreasing = T)
poke2 <- pokemon[Z.pred == 1,]
ord2 <- order(poke2$Total, decreasing = T)
cbind(head(poke1[ord1, 1:2], 10), head(poke2[ord2, 1:2], 10))
