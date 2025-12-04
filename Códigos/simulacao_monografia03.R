rm(list=ls())
source("funcoes_definidas.R")

df <- data.frame(Shift=seq(0,4,length.out=1000))
df$ARL <- numeric(nrow(df))
for (i in 1:nrow(df)) {
  df$ARL[i] <- EstimateARL(lambda=.25,
                        L=2.898,
                        shift = df$Shift[i],
                        mu=0,
                        sigma = 1,
                        zero_state = FALSE,
                        gshow = FALSE)
}
par(mfrow=c(2,1))
plot(df$Shift,df$ARL,type="l")
plot(df$Shift,exp(-df$Shift^2)/df$ARL,type="l")