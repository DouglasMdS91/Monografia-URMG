rm(list=ls())
source("funcoes_definidas.R")

EstimateARL(lambda = .25,
            L = 2.898,
            shift = 3,
            mu=0,
            sigma=1,
            zero_state=FALSE,
            gshow = TRUE)
