# multi-species functional response
get_MSFR <- function(L, S, params, N_scale) {
  
  # add small number of L and S and divide by N_scale
  L <- (L + 1) / N_scale
  S <- (S + 1) / N_scale
  
  x <- c(L, S)
  
  # calculate the denominator
  denominator <- 1 + (params$b * params$h) %*% x ^ (1 + params$q)
  
  L_prey <- params$b[1] * L ^ (1 + params$q) / denominator
  
  S_prey <- params$b[2] * S ^ (1 + params$q) / denominator
  
  return(c(L_prey, S_prey))
}

# growth - smolt production
# https://docs.salmonmse.com/articles/equations.html
growth <- function(params, L, S) {
  
  beta_L <- params$alpha_L / params$K_L
  beta_S <- params$alpha_S / params$K_S
  
  Lplus <- L * params$alpha_L / (1 + beta_L * L) 
  
  Splus <- S * params$alpha_S / (1 + beta_S * S)
  
  return(c(max(0, Lplus), max(0, Splus)))
}

# ocean survival
ocean <- function(L, S, params) {
  
  surv <- params$mO_surv
  
  Lplus <- L * (1 - exp(-params$mO_ddL * (S + 1e-6))) * surv
  Splus <- S * (1 - exp(-params$mO_ddS / (L + 1e-6))) * surv
  
  return(c(Lplus, Splus, (L - Lplus) / L, (S - Splus) / S))
  
}



# iterate one time step - 1d action
step <- function(L, S, params, N_scale) {
  
  # density dependent reproduction
  smolts <- growth(params, L, S)
  
  # ocean mortality
  ocean_out <- ocean(smolts[1], smolts[2], params)
  returns <- ocean_out[1:2]
  
  # get number of pinnipeds in estuary
  P_estuary <- params$P_mean
  
  # get MSFR
  consumed <- get_MSFR(L, S, params, N_scale) * P_estuary
  
  # get final numbers of spawners to system
  spawners_L <- max(returns[1] - consumed[1], 0)
  spawners_S <- max(returns[2] - consumed[2], 0)
  
  return(c(spawners_L, spawners_S, consumed[1], consumed[2], 
           ocean_out[3], ocean_out[4]))
  
}

# plot simulation
plot_sim <- function(df, burnin) {
  
  # convert to long
  df_N_long <- df[-c(1:burnin), c("L", "S", "t")] %>% 
    pivot_longer(cols = !t,
                 names_to = "species",
                 values_to = "N")
  df_C_long <- df[-c(1:burnin), c("L_consumed", "S_consumed", "t")] %>% 
    pivot_longer(cols = !t,
                 names_to = "species",
                 values_to = "N")
  
  df_O_long <- df[-c(1:burnin), c("L_ocean", "S_ocean", "t")] %>% 
    pivot_longer(cols = !t,
                 names_to = "species",
                 values_to = "N")
  
  # abundance plot
  N_plot <- ggplot(data = df_N_long) +
    geom_line(aes(x = t, y = N, color = species)) +
    ggtitle("Spawner abundance") +
    theme_minimal()
  
  # consumed plot
  C_plot <- ggplot(data = df_C_long) +
    geom_line(aes(x = t, y = N, color = species)) +
    ggtitle("Consumed") +
    theme_minimal()
  
  # ocean mortality plot
  O_plot <- ggplot(data = df_O_long) +
    geom_line(aes(x = t, y = N, color = species)) +
    ggtitle("Ocean mortality") +
    theme_minimal()
  
  final_plot <- N_plot + C_plot + O_plot + plot_layout(nrow = 1)
  
  return(final_plot)
}

