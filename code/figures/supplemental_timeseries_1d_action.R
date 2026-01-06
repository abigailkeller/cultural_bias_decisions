library(tidyverse)
library(patchwork)

source("code/timeseries_functions.R")

# read in MSFR posterior
samples <- do.call(rbind, 
                   readRDS("posterior_samples/MSFR_posterior_1day.rds"))
set.seed(123)
index <- sample(nrow(samples), 1000)
sub <- readRDS("posterior_samples/index_sub.rds")
posterior <- samples[index[sub], ]

#################
# deterministic #
#################

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

# simulate time series
ntime <- 200

########################
# get equilibria plots #
########################

burnin <- 50

df_equil <- as.data.frame(matrix(NA, nrow = 1, ncol = 5))
colnames(df_equil) <- c("alpha_L", "set", "L_eq", "S_eq", "mO_ddS")

# lamprey production action
alpha_change <- seq(params$alpha_S * 0.05, params$alpha_S * 0.4, 0.5)

# lamprey parasitism
# L_para <- c(9e7, 7e6, 5e6)
L_para <- c(9e7, 5e6, 3e6)

for (a in seq_along(alpha_change)) {
  
  params$alpha_L <- alpha_change[a]
  
  for (k in seq_along(L_para)) {
    
    params$mO_ddS <- L_para[k]
    
    for (j in seq_along(1:nrow(posterior))) {
      params$b <- c(exp(posterior[j, "b_log[1]"]), 
                    exp(posterior[j, "b_log[2]"]))
      params$q <- posterior[j, "q"]
      params$h <- c(posterior[j, "h[1]"], posterior[j, "h[2]"])
      
      # get index
      index <- (a - 1) * length(L_para) * nrow(posterior) + 
        (k - 1) * nrow(posterior) + 
        j
      
      params$alpha_L <- df_equil[index, "alpha_L"] <- alpha_change[a]
      
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
      
      # save equilibria
      df_equil[index, "L_eq"] <- mean(df[(burnin + 1):ntime, "L"])
      df_equil[index, "S_eq"] <- mean(df[(burnin + 1):ntime, "S"])
      
      # save params
      df_equil[index, "set"] <- j
      df_equil[index, "mO_ddS"] <- L_para[k]
    }
  }
}


para_labels <- c(
  "9e+07" = "D[S]==9e7",
  "7e+06" = "D[S]==7e6",
  "5e+06" = "D[S]==5e6"
)

# get optimal actions
optimal_actions <- df_equil %>% 
  group_by(set, mO_ddS) %>% 
  summarise(S_eq = max(S_eq))
optimal_actions <- left_join(optimal_actions, df_equil,
                             by = c("set", "mO_ddS", "S_eq"))


# plot
plot <- ggplot() +
  geom_line(data = df_equil,
            aes(x = alpha_L, y = S_eq, color = as.factor(set),
                linetype = as.factor(mO_ddS))) +
  geom_point(data = optimal_actions,
             aes(x = alpha_L, y = S_eq, color = as.factor(set))) +
  scale_color_manual(values = c("cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("1", "2", "3", "4")) +
  facet_wrap(~ mO_ddS, ncol = 1) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        labels = c("5e6", "7e6", "9e7")) +
  scale_x_continuous(breaks = c(250 * 0.9 * 0.1,
                                250 * 0.9 * 0.25,
                                250 * 0.9 * 0.4),
                     labels = c(
                       expression(alpha[L] == 0.1 * alpha[S]),
                       expression(alpha[L] == 0.25 * alpha[S]),
                       expression(alpha[L] == 0.4 * alpha[S])
                     )) +
  scale_y_continuous(breaks = c(100000, 120000, 140000)) +
  labs(x = "action (lamprey production)", 
       y = "utility (abundance of salmon returns)",
       color = "posterior\nsample", 
       linetype = "parasitism\nstrength") +
  theme_minimal(base_family = "Arial") + 
  theme(strip.text = element_blank(), 
        legend.title = element_text(hjust = 0.5)) 
ggsave("figures/supplemental/timeseries_equilibria.png",
       plot, dpi = 400, height = 6, width = 4.5, bg = "white")

