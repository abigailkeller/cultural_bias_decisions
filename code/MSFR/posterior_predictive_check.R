library(nimble)

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
         end = Date + 1)

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


###########################
# build and compile model #
###########################

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

# function for calculating MSFR
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
                       constants = constants)

# compile the model
CmyModel <- compileNimble(myModel)

#################################
# get theta and read in samples #
#################################

# get theta
theta <- c("h", "b_log", "q", "S_A", "L_A")
theta <- myModel$expandNodeNames(theta)

# read in samples 
samples_in <- readRDS("posterior_samples/MSFR_posterior.rds")

# calculate with a subset of samples
lower <- 1
upper <- dim(samples_in[[1]])[1]
sub <- seq(lower, upper, 5)

# top nodes for posterior predictive samples likelihood
samples <- rbind(samples_in[[1]][sub, theta],
                 samples_in[[2]][sub, theta],
                 samples_in[[3]][sub, theta],
                 samples_in[[4]][sub, theta])

#########################################
# generate posterior predictive samples #
#########################################

set.seed(123)

# constants needed for indexing
t_exact <- which(dates_df$Year %in% c("2017", "2019"))
t_inexact <- which(!dates_df$Year %in% c("2017", "2019"))
cons_ref <- c(gastro_CSL$d_unq, gastro_CSL$d_unq + 46)

# function to generate posterior predictive samples and calculate discrepancies
ppSamplerNF <- nimbleFunction(
  setup = function(model, samples, theta) {
    
    theta <- model$topologicallySortNodes(theta)
    
    # nodes to simulate
    simNodes <- model$getNodeNames()[!(model$getNodeNames() %in% theta)]
    
    # data nodes
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    dataNodes_cons <- dataNodes[1:214]
    dataNodes_S <- dataNodes[215:260]
    dataNodes_Lexact <- dataNodes[c(261:271, 284:290)]
    dataNodes_Linexact <- dataNodes[c(272:283, 291:306)]
    
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    vars <- colnames(samples[, theta])
    
    # subset posterior samples to just theta (top nodes in graph)
    samples_sub <- samples[, theta]
    
  },
  run = function(samples = double(2), t_exact = double(1), t_inexact = double(1), cons_ref = double(1)) {
    
    nSamp <- dim(samples)[1]
    out <- matrix(nrow = nSamp, ncol = 8)
    
    for(i in 1:nSamp) {
      
      # add theta to model
      values(model, vars) <<- samples_sub[i, ]
      
      # update nodes based on theta
      model$simulate(simNodes, includeData = TRUE)
      
      # calculate deviance
      out[i, 1] <- 2 * model$calculate(dataNodes_cons) # consumed
      out[i, 2] <- 2 * model$calculate(dataNodes_S) # salmon count
      out[i, 3] <- 2 * model$calculate(dataNodes_Lexact) # lamprey count ("exact" measurement)
      out[i, 4] <- 2 * model$calculate(dataNodes_Linexact) # lamprey count ("inexact" measurement)

      # calculate chi-squared
      cons_exp <- values(model, "F_exp")[cons_ref]
      out[i, 5] <- sum(
        (values(model, dataNodes_cons) - cons_exp) ^ 2 / cons_exp
      )

      # salmon count
      out[i, 6] <- sum(
        (values(model, dataNodes_S) - values(model, "S_P")) ^ 2 / values(model, "S_P")
      )

      # lamprey count ("exact" measurement)
      exp_exact <- values(model, "L_P")[t_exact]
      out[i, 7] <- sum(
        (values(model, dataNodes_Lexact) - exp_exact) ^ 2 / (exp_exact + 1e-6)
      )
      
      # lamprey count ("inexact" measurement)
      exp_Linexact <- values(model, "L_P")[t_inexact] * values(model, "p_D")
      out[i, 8] <- sum(
        (values(model, dataNodes_Linexact) - exp_Linexact) ^ 2 / (exp_Linexact + 1e-6)
      )
    }
    returnType(double(2))
    return(out)
  })


# generate posterior predictive sampler
ppSampler <- ppSamplerNF(
  myModel, # uncompiled model
  samples, # posterior samples
  theta
)

# compile posterior predictive samples
cppSampler <- compileNimble(
  ppSampler, 
  project = CmyModel # compiled model
)

# run compiled function
post_predictive <- cppSampler$run(samples, t_exact, t_inexact, cons_ref)


##############################################
# calculate discrepancies with observed data #
##############################################

calc_disc_data <- nimbleFunction(
  setup = function(model, samples, theta) {
    
    # calculate model graph dependencies of theta to update
    deps <- model$getDependencies(theta, self = TRUE)
    
    # data nodes
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    dataNodes_cons <- dataNodes[1:214]
    dataNodes_S <- dataNodes[215:260]
    dataNodes_Lexact <- dataNodes[c(261:271, 284:290)]
    dataNodes_Linexact <- dataNodes[c(272:283, 291:306)]

    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    vars <- colnames(samples[, theta])
    
    # subset posterior samples to just theta (top nodes in graph)
    samples_sub <- samples[, theta]
    
  },
  run = function(samples = double(2), t_exact = double(1), t_inexact = double(1), cons_ref = double(1)) {
    
    nSamp <- dim(samples)[1]
    out <- matrix(nrow = nSamp, ncol = 8)
    
    for(i in 1:nSamp) {
      
      # add theta
      values(model, vars) <<- samples_sub[i, ]
      
      # update dependencies
      model$calculate(deps)
      
      # calculate deviance
      out[i, 1] <- 2 * model$calculate(dataNodes_cons) # consumed
      out[i, 2] <- 2 * model$calculate(dataNodes_S) # salmon count
      out[i, 3] <- 2 * model$calculate(dataNodes_Lexact) # lamprey count ("exact" measurement)
      out[i, 4] <- 2 * model$calculate(dataNodes_Linexact) # lamprey count ("inexact" measurement)

      # calculate chi-squared
      cons_exp <- values(model, "F_exp")[cons_ref]
      out[i, 5] <- sum(
        (values(model, dataNodes_cons) - cons_exp) ^ 2 / cons_exp
      )

      # salmon count
      out[i, 6] <- sum(
        (values(model, dataNodes_S) - values(model, "S_P")) ^ 2 / values(model, "S_P")
      )

      # lamprey count ("exact" measurement)
      exp_exact <- values(model, "L_P")[t_exact]
      out[i, 7] <- sum(
        (values(model, dataNodes_Lexact) - exp_exact) ^ 2 / (exp_exact + 1e-6)
      )
      
      # lamprey count ("inexact" measurement)
      exp_Linexact <- values(model, "L_P")[t_inexact] * values(model, "p_D")
      out[i, 8] <- sum(
        (values(model, dataNodes_Linexact) - exp_Linexact) ^ 2 / (exp_Linexact + 1e-6)
      )
    }
    returnType(double(2))
    return(out)
  })

# generate instance
disc_data <- calc_disc_data(
  myModel, # uncompiled model, contains data
  samples, # posterior samples
  theta
)

# compile
Cdisc_data <- compileNimble(
  disc_data, 
  project = CmyModel 
)

# run compiled function
data_discrepancies <- Cdisc_data$run(samples, t_exact, t_inexact, cons_ref)

######################
# calculate p values #
######################

pvalues <- data.frame(
  disc_type = c(rep("deviance", 4), rep("chisq", 4)),
  data = c("consumed", "salmon_A", "lamp_exact", "lamp_inexact", "consumed", "salmon_A", "lamp_exact", "lamp_inexact"),
  pvalue = rep(NA, 8)
)
for (i in seq_len(nrow(pvalues))) {
  pvalues[i, "pvalue"] <- mean(post_predictive[, i] > data_discrepancies[, i])
}

# save p values
saveRDS(pvalues, "data/ppp.rds")
