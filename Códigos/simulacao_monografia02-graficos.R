rm(list=ls())
source("funcoes_definidas.R")

desvio <- 0
zero <- TRUE
for (desvio in c(0,.5,1,2,3,4)) {
  for (zero in c(TRUE,FALSE)) {
    nome <- paste0("ARL_estimado_",
                   ifelse(zero,"Zero_State_","Steady_State_"),
                   ifelse(desvio==.5,"_5",desvio),
                   ".png")
    png(filename = nome, width = 1000, height = 500, units = "px",res = 150)
    EstimateARL(lambda = .25,
                L = 2.898,
                shift = desvio,
                mu=0,
                sigma=1,
                zero_state=zero,
                gshow = TRUE)
    dev.off()
  }
  
}