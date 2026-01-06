library(nimble)
library(MCMCvis)
library(parallel)
library(viridis)
library(tidyverse)

################
# read in data #
################

# read in gastrointestinal data
gastro_CSL <- readRDS("data/model_data/gastro_data.rds")

# read in lamprey count data from years (2017, 2019) with LPS and night-time count
total_lamprey <- readRDS("data/model_data/total_lamprey.rds")

# read in bonneville dam visual count data
bonneville <- readRDS("data/model_data/bonneville_visualcount.rds")

# read in lamprey visual count passage proportion posterior summary
LPS_passage_posterior <- MCMCsummary(
  readRDS("posterior_samples/LPS_posterior.rds")
)

# read in unique dates
dates_df <- readRDS("data/model_data/unique_dates.rds")

# add time window to unique dates
dates_df <- dates_df %>% 
  mutate(start = Date - 1,
         end = Date + 3)

###############################################################
# calculate passed fish, L^P_T and S^P_t based on time window #
###############################################################

# get consumed at each time point
consumed_day <- gastro_CSL %>% 
  group_by(d_unq) %>% 
  summarize(total_L = sum(Lamprey),
            total_S = sum(Chinook))

# add passage data to gastro data
for (i in 1:nrow(dates_df)) {
  
  # get salmon passed
  temp <- bonneville[
    bonneville$Date >= dates_df[i, "start"] & 
      bonneville$Date <= dates_df[i, "end"], ]
  dates_df[i, "Chinook_P"] <- sum(abs(temp$Chin), na.rm = TRUE)
  
  # get lamprey passed - years with LPS and night-time count recorded
  if (dates_df[i, "Year"] %in% c("2017", "2019")) {
    temp_L <- total_lamprey[total_lamprey$Date >= dates_df[i, "start"] &
                             total_lamprey$Date <= dates_df[i, "end"], ]
    
    dates_df[i, "Lamprey_P"] <- sum(abs(temp_L$total), na.rm = TRUE)
    
  } else {
    # get lamprey passed - years without LPS and night-time count recorded
    
    temp_L <- bonneville[
      bonneville$Date >= dates_df[i, "start"] & 
        bonneville$Date <= dates_df[i, "end"], ]
    dates_df[i, "Lamprey_P"] <- sum(abs(temp_L$Lmpry), na.rm = TRUE)
  }
}


##############
# MSFR model #
##############

MSFR_model <- nimbleCode({
  
  ###############################
  # Chinook salmon prey abundance
  
  for (t in 1:n_time) {
    
    # Equation 8
    # get expected number of salmon passed
    S[t] ~ dpois(S_P[t])
    
    # Equation 9
    # get expected number of salmon available
    S_P[t] ~ dbinom(size = round(S_A[t] - F_day[t, 2]), 
                    prob = salmon_PE)
  }
  
  ########################
  # Lamprey prey abundance
  
  # Equation 10
  # get expected number of available lamprey passed - exact years
  for (t in 1:n_exact) {
    L_exact[t] ~ dpois(L_P[t_exact[t]])
  }
  
  # get expected number of available lamprey passed - inexact years
  for (t in 1:n_inexact) {
    # Equation 11
    L_inexact[t] ~ dbinom(size = L_P[t_inexact[t]], 
                          prob = p_D[t])
    # Equation 12
    p_D[t] ~ dbeta(alpha_D, beta_D)
  }
  
  for (t in 1:n_time) {
    # Equation 14
    # get expected number of lamprey available 
    L_P[t] ~ dbinom(size = round(L_A[t] - F_day[t, 1]), 
                    prob = lamprey_PE)
  }
  
  ###############
  # Consumed prey
  
  # number of prey consumed
  for (j in 1:n_obs) {
    
    for (i in 1:n_species) {
      
      # Equations 15-16
      consumed[j, i] ~ dpois(F_exp[obs_ref[j], i]) 
      
    }
  }
  
  ###################################
  # Multi-species functional response
  
  for (t in 1:n_time) {
    
    # Equations 17-18
    # expected number of prey consumed - function of available fish
    F_exp[t, 1:n_species] <- get_MSFR(L_A[t] / N_scale, 
                                      S_A[t] / N_scale, 
                                      q, h[1:n_species], 
                                      b_log[1:n_species])
  }
  
  
  ##########
  # priors #
  ##########
  
  # MSFR params
  for (i in 1:n_species) {
    h[i] ~ dunif(0, 1000) # handling time
    b_log[i] ~ dunif(-100, 100) # attack rate
  }
  
  q ~ dunif(-2, 2)
  
  # total expected number of lamprey and salmon available
  for (t in 1:n_time) {
    L_A[t] ~ dunif(0, 1000000)
    S_A[t] ~ dunif(0, 1000000)
  }
  
})


# initial values
inits <- function() {
  list(
    p_D = rep(0.2, length(which(!dates_df$Year %in% c("2017", "2019")))),
    h = c(1, 1),
    b_log = c(log(1000), log(200)),
    L_P = round((dates_df[, "Lamprey_P"] + 1) / 0.2),
    S_P = round(dates_df[, "Chinook_P"]),
    L_A = round((dates_df[, "Lamprey_P"] + 1) / 0.2 / 0.5),
    S_A = round(dates_df[, "Chinook_P"] / 0.3),
    q = 0
  )
}

# bundle up data and constants
constants <- list (
  salmon_PE = 0.84, # Frick report 2008
  lamprey_PE = 0.55, # Moser et al. 2002
  alpha_D = LPS_passage_posterior["alpha_D", "mean"],
  beta_D = LPS_passage_posterior["beta_D", "mean"],
  n_species = 2,
  N_scale = 10000,
  n_time = nrow(dates_df),
  t_exact = which(dates_df$Year %in% c("2017", "2019")),
  t_inexact = which(!dates_df$Year %in% c("2017", "2019")),
  n_exact = length(which(dates_df$Year %in% c("2017", "2019"))),
  n_inexact = length(which(!dates_df$Year %in% c("2017", "2019"))),
  obs_ref = gastro_CSL$d_unq,
  n_obs = nrow(gastro_CSL)
)

data <- list (
  L_exact = dates_df[dates_df$Year %in% c("2017", "2019"), "Lamprey_P"],
  L_inexact = dates_df[!dates_df$Year %in% c("2017", "2019"), "Lamprey_P"],
  S = dates_df[, "Chinook_P"],
  consumed = gastro_CSL[, c("Lamprey", "Chinook")],
  F_day = consumed_day[, c("total_L", "total_S")]
  
)

cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("MSFR_model", "inits", "data", "constants", "gastro_CSL",
                    "dates_df"))

# Create parallel function 
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # MSFR
  get_MSFR <- nimbleFunction (
    # input and output types
    run = function(L = double(0), S = double(0), 
                   q = double(0), h = double(1),
                   b_log = double(1))
    {
      returnType(double(1))
      
      # combine abundance into one vector
      x <- c(L, S)
      
      # calculate the denominator
      denominator <- 1 + (exp(b_log) * h) %*% x ^ (1 + q)
      
      L_prey <- exp(b_log[1]) * L ^ (1 + q) / denominator
      
      S_prey <- exp(b_log[2]) * S ^ (1 + q) / denominator
      
      return(c(L_prey, S_prey))
    }
  )
  assign("get_MSFR", get_MSFR, envir = .GlobalEnv)
  
  # build model
  myModel <- nimbleModel(code = MSFR_model,
                         data = data,
                         constants = constants,
                         inits = inits())
  
  
  # build the MCMC
  mcmcConf_myModel <- configureMCMC(
    myModel,
    monitors = c("p_D", "F_exp", "q",
                 "h", "b_log", "S_A", "L_A",
                 "S_P", "L_P"),
    useConjugacy = FALSE, enableWAIC = TRUE)
  
  # add block sampler for attack rate params
  mcmcConf_myModel$removeSamplers("b_log")
  mcmcConf_myModel$removeSamplers("q")
  mcmcConf_myModel$addSampler(c("b_log[1]", "b_log[2]", "q"),
                              type='RW_block')
  
  # build MCMC
  myMCMC <- buildMCMC(mcmcConf_myModel)
  
  # compile the model and MCMC
  CmyModel <- compileNimble(myModel)
  
  # compile the MCMC
  cmodel_mcmc <- compileNimble(myMCMC, project = myModel)
  
  # run MCMC
  cmodel_mcmc$run(50000000, thin = 10000,
                  reset = FALSE)
  
  
  samples <- as.mcmc(as.matrix(cmodel_mcmc$mvSamples))
  
  return(samples)
  
})

stopCluster(cl)

# remove burn-in
lower <- 200
upper <- 5001
sequence <- seq(lower, upper, 1)
out_sub <- list(
  out[[1]][sequence, ], out[[2]][sequence, ],
  out[[3]][sequence, ], out[[4]][sequence, ]
)

# save samples
saveRDS(out_sub, "posterior_samples/MSFR_posterior_3day.rds")
