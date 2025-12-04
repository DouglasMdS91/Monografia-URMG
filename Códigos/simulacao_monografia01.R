rm(list = ls())
source("funcoes_definidas.R")
v <- c();N <- 10^5 #number of simulations
desvio <- 0
for (i in 1:N) {
  v <- c(v,SimularAteParar(nmax=10^4,z0=0,lambda=.25,
                           L=2.898,sigma=1,shift=desvio))
}
cat("Simulação com Shift=",desvio,"\n")
summary(v)
hist(v,col="lightblue",nclass = 30,
     main = paste0("Simulação com shift=",desvio), 
     ylab = "Frequência")
m <- mean(v)
abline(v = m,col="red")
text(
  x = m + 0.02 * diff(par("usr")[1:2]), 
  y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
  labels = paste0("Média = ", round(m, 2)),
  col = "red",
  adj = 0
)