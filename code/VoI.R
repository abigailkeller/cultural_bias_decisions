library(tidyverse)
library(patchwork)
library(ggnewscale)
library(ggbrace)
library(grid)

source("code/timeseries_functions.R")

# read in MSFR posterior
posterior <- do.call(rbind, 
                     readRDS("posterior_samples/MSFR_posterior_1day.rds"))

# randomly select 1000 samples
set.seed(123)
index <- sample(nrow(posterior), 1000)

params <- list (
  # MSFR
  b = c(NA, NA), # attack rate [lamprey, salmon]
  q = NA, # attack rate density dependence
  h = c(NA, NA), # handling time [lamprey, salmon]
  # growth
  alpha_S = 250 * 0.9, # combined spawner mortality and intrinsic growth rate
  alpha_L = 250 * 0.65, # combined spawner mortality and intrinsic growth rate
  K_L = 3000000, # lamprey carrying capacity
  K_S = 3000000, # salmon carrying capacity
  # ocean survival
  mO_ddS = NA, # salmon ocean survival (function of lamprey density)
  mO_ddL = 5e-5, # lamprey ocean survival (function of salmon density)
  mO_surv = 8 / (8 + 125), # shared density-ind ocean survival
  # pinnipeds
  P_mean = 10000 # mean number of pinnipeds
)

######################
# get utility matrix #
######################

# simulate time series
ntime <- 200

burnin <- 50

# get s x a utilty matrix
alpha_change <- seq(params$alpha_S * 0.05, params$alpha_S * 0.45, 0.5)
utility <- as.data.frame(matrix(NA, nrow = length(alpha_change),
                                ncol = 3000))


# lamprey parasitism
L_para <- c(9e7, 5e6, 3e6)

for (a in seq_along(alpha_change)) {
  
  
  params$alpha_L <- alpha_change[a]
  
  for (i in seq_along(index)) {
    
    params$b <- c(exp(posterior[index[i], "b_log[1]"]), 
                  exp(posterior[index[i], "b_log[2]"]))
    params$q <- posterior[index[i], "q"]
    params$h <- c(posterior[index[i], "h[1]"], posterior[index[i], "h[2]"])
    
    for (j in seq_along(L_para)) {
      
      params$mO_ddS <- L_para[j]
      
      # get column
      column <- (i - 1) * length(L_para) + j
      
      df <- as.data.frame(matrix(NA, nrow = ntime, ncol = 7))
      colnames(df) <- c("t", "L", "S", "L_consumed", "S_consumed", 
                        "L_ocean", "S_ocean")
      df[1, c("L", "S")] <- c(60000, 500000)
      df$t <- 1:ntime
      
      for (t in 2:ntime) {
        
        # iterate through step function
        df[t, 2:7] <- step(L = df[t - 1, "L"], S = df[t - 1, "S"], 
                           params, N_scale = 10000)
      }
      
      # utility = salmon equilibria
      utility[a, column] <- mean(df[(burnin + 1):ntime, "S"])
    }
  }
}

# save utility
saveRDS(utility, "VOI/utility_fullposterior.rds")