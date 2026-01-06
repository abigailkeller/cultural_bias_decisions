library(nimble)
library(MCMCvis)
library(parallel)
library(viridis)
library(tidyverse)

# read in data
C_V <- readRDS("data/model_data/C_V.rds")
C_T <- readRDS("data/model_data/C_T.rds")

# calculate proportion
LPS_prop <- C_V / C_T

# write model
LPS_model <- nimbleCode({
  
  alpha_D ~ dunif(0, 100)
  beta_D ~ dunif(0, 100)
  
  for (i in 1:nprop) {
    LPS_prop[i] ~ dbeta(alpha_D, beta_D)
  }
  
})

# initial values
inits <- function() {
  list(
    alpha_D = runif(1, 0.1, 10),
    beta_D = runif(1, 0.1, 10)
  )
}

# bundle up data and constants
constants <- list (
  nprop = length(LPS_prop)
)

data <- list (
  LPS_prop = LPS_prop
)

cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("LPS_model", "inits", "data", "constants"))

# Create parallel function 
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # build model
  myModel <- nimbleModel(code = LPS_model,
                         data = data,
                         constants = constants,
                         inits = inits())
  
  
  # build the MCMC
  mcmcConf_myModel <- configureMCMC(
    myModel,
    monitors = c("alpha_D", "beta_D"),
    useConjugacy = FALSE, enableWAIC = TRUE)
  
  # build MCMC
  myMCMC <- buildMCMC(mcmcConf_myModel)
  
  # compile the model and MCMC
  CmyModel <- compileNimble(myModel)
  
  # compile the MCMC
  cmodel_mcmc <- compileNimble(myMCMC, project = myModel)
  
  # run MCMC
  cmodel_mcmc$run(2000, thin = 1,
                  reset = FALSE)
  
  samples <- as.mcmc(as.matrix(cmodel_mcmc$mvSamples))
  
  return(samples)
  
})

stopCluster(cl)

lower <- 200
upper <- 2001
sequence <- seq(lower, upper, 1)
out_sub <- list(
  out[[1]][sequence, ], out[[2]][sequence, ],
  out[[3]][sequence, ], out[[4]][sequence, ]
)
saveRDS(out_sub, "posterior_samples/LPS_posterior.rds")
