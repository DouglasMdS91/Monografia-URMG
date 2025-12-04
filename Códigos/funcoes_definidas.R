SimularAteParar <- function(nmax,z0,lambda,L,sigma,shift){
  #Estimativa com shift + limites estacionários
  z <- z0;n <- 1
  lcl <- z0 - L*sigma*sqrt( lambda/(2-lambda))
  ucl <- z0 + L*sigma*sqrt( lambda/(2-lambda))
  while (n<=nmax) {
    x <- rnorm(1,z0+shift,sd=sigma)
    z <- c(z,lambda * x + (1-lambda)*z[length(z)])
    if(z[n+1] <lcl | z[n+1]>ucl){
      return(n)
    }
    n <- n+1
  }
  return(nmax)
}

GetTransitionMatrix <- function(m,lambda,L,shift,mu,sigma){
  t <- 2*m+1; delta <- L*sigma*sqrt( lambda/(2-lambda) )/t
  S <- mu+2*c(-m:m)*delta; P <- matrix(nrow = t+1,ncol = t+1)
  #Part R of P
  for (j in 1:t) {
    for (k in 1:t) {
      P[j,k] <- pnorm( ( S[k]+delta- (1-lambda)*S[j]-lambda*(mu+shift) )/
                         (lambda*sigma) ) - 
        pnorm( ( S[k]-delta- (1-lambda)*S[j]-lambda*(mu+shift) )/
                 (lambda*sigma) )
    }
  }
  #Part (I-R)1 of P
  for (j in 1:t) {
    P[j,t+1] <- 1-sum(P[j,1:t])
  }
  #Part 0 of P
  for (k in 1:t) {
    P[t+1,k] <- 0
  }
  #Part 1 of P
  P[t+1,t+1] <- 1
  return(P)
}

CalculateARL_T <- function(m,lambda,L,shift,mu,sigma,zero_state=FALSE){
  P <- GetTransitionMatrix(m,lambda,L,shift,mu,sigma)
  t <- 2*m+1; R <- P[1:t,1:t]
  Q <- GetTransitionMatrix(m,lambda,L,shift=0,mu,sigma)
  #O estado m+1 é o estado ZERO
  Q[t+1,t+1] <- 0; Q[t+1,m+1] <- 1
  
  ev <- eigen(t(Q))
  tol <- 1e-10
  idx <- which(Mod(ev$values - 1) < tol)
  
  p <- ev$vectors[, idx, drop = FALSE]
  
  if(all(Im(p)==0)){
    p <- Re(p)
  }else{
    stop("ERRO: Parte imaginária em P? Revisar.")
  }
  
  if(Re(p[1])<0){
    p <- -p
  }
  if(!all(p>0)){
    stop("ERRO: valor negativo em p. Contradiz teoria.")
  }
  
  p <- p/sum(p)
  #Eliminar última coordenada
  pss <- p[-(t+1)]
  #Normalizar probabilidade
  pss <- pss / sum(pss)
  
  It <- matrix(0,nrow = t,ncol = t)
  diag(It) <- 1
  arl <- pss %*% (solve(It-R)%*%rep(1,t))
  
  #Zero State
  if(zero_state){
    p0 <- rep(0,t)
    p0[m+1] <- 1
    arl0 <- p0 %*% (solve(It-R)%*%rep(1,t))
    return(arl0)
  }
  
  return(arl)
}

EstimateARL <- function(lambda,L,shift,mu,sigma,zero_state,gshow=FALSE){
  df <- data.frame(m=c(25,29,33,37,41)) #t=51 59 67 75 83
  df$t <- 2*df$m+1
  df$ARL <- numeric(nrow(df))
  
  for (i in 1:nrow(df)) {#i <- 1
    m <- df$m[i]
    df$ARL[i] <- CalculateARL_T(m,lambda,L,shift,mu,sigma,zero_state)
  }
  
  df$x1 <- 1/df$t
  df$x2 <- df$x1^2
  
  modelo <- lm(ARL ~ x1 + x2,data=df)
  
  arl_final <- modelo$coefficients[1]
  
  if(gshow){
    y_lim <- c(min(modelo$fitted.values,df$ARL,arl_final),
               max(modelo$fitted.values,df$ARL,arl_final))
    y_lim <- y_lim+.05*(y_lim[2]-y_lim[1])*c(-1,1)
    if(zero_state){
      titulo <- "Estado Zero"
    }else{
      titulo <- "Estado Estacionário"
    }
    plot(df$t, df$ARL, pch = 19,ylim=y_lim,
         xlab="t",ylab="ARL estimado",
         main=paste0(titulo," - ","shift = ",shift))
    curve(coef(modelo)[1] + coef(modelo)[2]/x + coef(modelo)[3]/x^2,
          add = TRUE, col = "red", lwd = 2)
    abline(h = arl_final,col="blue")
    text(
      x = 65, 
      y = arl_final,
      labels = paste0("ARL estimado = ", round(arl_final, 2)),
      col = "blue",
      pos=ifelse(arl_final>=min(df$ARL),1,3)
    )
  }
  return(round(modelo$coefficients[1],2))
}

