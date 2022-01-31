#Exemplo Gibbs: Salários na área de dados na Índia
#Fonte: https://www.kaggle.com/iamsouravbanerjee/analytics-industry-salaries-2022-india

library(stringr)

salarios <- read.csv("Salary Dataset.csv")
summary(salarios)

head(salarios$Salary)
salarios$Salary <- gsub("₹", "", salarios$Salary)
salarios$Salary <- gsub(",", "", salarios$Salary)

mat_salarios <- matrix(unlist(str_split(salarios$Salary, "/")), ncol=2, byrow=T)
ind <- which(mat_salarios[,2] == "yr")
salarios$Salary <- as.numeric(mat_salarios[,1])
salarios$Salary[ind] <- salarios$Salary[ind]/12
ind <- which(mat_salarios[,2] != "hr")
salarios <- salarios[ind,]
salarios <- na.omit(salarios)

#Salario esta em rupias, vamos converter em reais (1 rupia = 0.071)
salarios$Salary <- salarios$Salary*0.071

table(salarios$Job.Title)

salarios_analista <- salarios$Salary[salarios$Job.Title == "Data Analyst"]
salarios_datasci <- salarios$Salary[salarios$Job.Title == "Data Scientist"]

par(mfrow = c(2,2), mar = c(2,2,2,2))
plot(density(salarios_analista))
plot(density(salarios_datasci))
plot(density(log(salarios_analista)))
plot(density(log(salarios_datasci)))

Gibbs <- function(Y, mu0, nu, alfa, beta, n.sims, semente = NULL){
  n <- length(Y)
  soma_Y2 <- sum(Y^2)
  Y_barra <- mean(Y)
  alfa_post <- (n+1)/2 + alfa
  mu0_post <- (n*Y_barra + mu0/nu)/(n + 1/nu)
  beta_part <- soma_Y2/2 + mu0^2/(2*nu) + beta - mu0_post^2*(n+1/nu)/2
  #Definir os parametros para o MCMC
  mu <- tau <- as.numeric()
  mu[1] <- mu0
  tau[1] <- alfa_post/beta_part
  set.seed(semente)
  for(i in 2:n.sims){
    beta_post <- beta_part + (n+1/nu)*(mu[i-1]-mu0_post)^2/2
    tau[i] <- rgamma(1, alfa_post, beta_post)
    mu[i] <- rnorm(1, mu0_post, sqrt(1/(tau[i]*(n+1/nu))))
  }
  return(list(mu=mu, tau=tau))
}

Gibbs_analista <- Gibbs(log(salarios_analista), mu0 =log(5000) , nu = 1000, alfa = 1, beta = 0.1,
                       n.sims = 10000, semente = 42)
Gibbs_datasci <- Gibbs(log(salarios_datasci), mu0 =log(7000) , nu = 1000, alfa = 1, beta = 0.1,
                      n.sims = 10000, semente = 42)


par(mfrow = c(2,1), mar = c(2,2,2,2))
plot(exp(Gibbs_analista$mu), type = 'l', ylim = range(c(exp(Gibbs_analista$mu), exp(Gibbs_datasci$mu))))
lines(exp(Gibbs_datasci$mu), type = 'l', col = 'red')

plot(Gibbs_analista$tau, type = 'l', ylim = range(c(Gibbs_analista$tau, Gibbs_datasci$tau)))
lines(Gibbs_datasci$tau, type = 'l', col = 'red')
