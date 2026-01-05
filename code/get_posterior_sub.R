library(tidyverse)
library(viridis)

source("final_code/code/timeseries_functions.R")

# read in utility
utility <- readRDS("final_code/VOI/utility_fullposterior.rds")

# get posterior
posterior <- do.call(rbind, 
                     readRDS("final_code/posterior_samples/MSFR_posterior_1day.rds"))
set.seed(123)
index <- sample(nrow(posterior), 1000)
posterior <- posterior[index, ]

# get 0.5 MSFR values
L <- seq(0, 150000, 1000)
S <- seq(150000, 0, -1000)
prey_L <- as.data.frame(matrix(NA, ncol = length(S), nrow = length(index)))
prey_S <- as.data.frame(matrix(NA, ncol = length(S), nrow = length(index)))

for (i in 1:length(index)) {
  
  params <- list(q = posterior[i, "q"],
                 h = c(posterior[i, "h[1]"],
                       posterior[i, "h[2]"]),
                 b = c(exp(posterior[i, "b_log[1]"]), 
                       exp(posterior[i, "b_log[2]"])))
  
  for (j in 1:length(L)) {
    
    prey_out <- get_MSFR(L = L[j], S = S[j], params = params,
                         N_scale = 10000)
    
    prey_L[i, j] <- prey_out[1] / (prey_out[1] + prey_out[2])
    prey_S[i, j] <- prey_out[2] / (prey_out[1] + prey_out[2])
    
  }
}

colnames(prey_L) <- L / (L + S)
colnames(prey_S) <- S / (L + S)

prey_L <- prey_L %>% 
  mutate(id = 1:1000) %>% 
  pivot_longer(cols = -id, 
               names_to = "x",
               values_to = "y")
prey_S <- prey_S %>% 
  mutate(id = 1:1000) %>% 
  pivot_longer(cols = -id, 
               names_to = "x",
               values_to = "y")

prey_L_0.5 <- prey_L[prey_L$x == 0.5, ]

# examine utility
utility_110000 <- utility[, seq(1, 3000, 3)]
index_110000 <- which(round(as.numeric(utility_110000[1, ]), -4) == 110000)

hist(prey_L_0.5[prey_L_0.5$id %in% index_110000, ]$y)
prey_L_0.5_sub <- prey_L_0.5[prey_L_0.5$id %in% index_110000, ]

utility_long <- utility_110000[, index_110000]
colnames(utility_long) <- index_110000
utility_long <- utility_long %>% 
  mutate(alpha_change = seq(250 * 0.9 * 0.05, 250 * 0.9 * 0.45, 0.5)) %>% 
  pivot_longer(cols = -alpha_change,
               names_to = "id",
               values_to = "utility")
utility_long$id <- as.numeric(utility_long$id)
utility_long <- left_join(utility_long, prey_L_0.5, by = "id")

ggplot(data = utility_long) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

ggplot(data = utility_long[utility_long$id %in% 
                             prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.7, ]$id, ]) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

ggplot(data = utility_long[utility_long$id %in% 
                             prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.6, ]$id, ]) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

ggplot(data = utility_long[utility_long$id %in% 
                             prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.5, ]$id, ]) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

ggplot(data = utility_long[utility_long$id %in% 
                             prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.41, ]$id, ]) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

#sub <- c(760, 445, 347, 294)
#sub <- c(760, 294, 445, 347)
# sub <- c(668, 815, 125, 299)
# sub <- c(prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.7, ]$id[2], 
#          prey_L_0.5_sub[round(prey_L_0.5_sub$y, 2) == 0.6, ]$id[12], 
#          125, 299)
# sub <- c(668, 996, 299, 125) # in order from highest y to lowest y
sub <- c(668, 125, 996, 299)

#which(index %in% sub)

ggplot(data = utility_long[utility_long$id %in% sub, ]) +
  geom_line(aes(x = alpha_change, y = utility, color = y, group = id)) +
  scale_color_viridis()

saveRDS(sub, "final_code/posterior_samples/index_sub.rds")
