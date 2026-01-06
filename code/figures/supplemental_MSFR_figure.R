library(tidyverse)
library(patchwork)
library(ggnewscale)
library(MCMCvis)
library(viridis)

source("code/timeseries_functions.R")

#######################
# prepare gastro data #
#######################

# read in gastro data
data_gastro <- read.csv("data/MSFR/BON_gastro_data.csv") %>% 
  select(-Removal_location)

# convert date columns
data_gastro$Date <- as.Date(data_gastro$Date, format = "%m/%d/%Y")

# replace NA with 0
cols <- setdiff(colnames(data_gastro), c("Date", "Pinniped_species"))
data_gastro[cols] <- lapply(data_gastro[cols], function(x) {
  x[is.na(x)] <- 0
  return(x)
})

# subset to just chinook and lamprey (and assume unidentified salmon = Chinook)
data_gastro_sub <- data_gastro %>% 
  mutate(Chinook = Salmonid_unidentified_adult + Chinook_adult) %>% 
  filter(Pinniped_species == "CSL")

# sort by date
data_gastro_sub <- data_gastro_sub[order(as.Date(data_gastro_sub$Date, 
                                                 format = "%d/%m/%Y")), ]

# get unique dates
dates_df <- data.frame(Date = as.Date(unique(data_gastro_sub$Date)),
                       d_unq = 1:length(unique(data_gastro_sub$Date)))

data_gastro_sub <- left_join(data_gastro_sub, dates_df, by = "Date")

# get posterior summaries
posterior <- do.call(rbind, 
                     readRDS("posterior_samples/MSFR_posterior_1day.rds"))
posterior_sub <- posterior[, c(1:46, 93:138)]
propN_L <- matrix(NA, nrow = nrow(posterior_sub), 
                  ncol = ncol(posterior_sub) / 2)
propN_S <- matrix(NA, nrow = nrow(posterior_sub), 
                  ncol = ncol(posterior_sub) / 2)
for (i in 1:(ncol(posterior_sub) / 2)) {
  L <- posterior_sub[, i]
  S <- posterior_sub[, i + ncol(posterior_sub) / 2]
  propN_L[, i] <- L / (L + S)
  propN_S[, i] <- S / (L + S)
}
colnames(propN_L) <- colnames(posterior)[1:46]
colnames(propN_S) <- colnames(posterior)[93:138]
propN_L_sum <- MCMCsummary(propN_L) %>% 
  mutate(d_unq = 1:46)
data_gastro_sub <- left_join(data_gastro_sub, propN_L_sum, by = "d_unq") %>% 
  rename(propNmean_L = mean,
         propNlow_L = `2.5%`,
         propNhigh_L = `97.5%`)
propN_S_sum <- MCMCsummary(propN_S) %>% 
  mutate(d_unq = 1:46)
data_gastro_sub <- left_join(data_gastro_sub, propN_S_sum, by = "d_unq") %>% 
  rename(propNmean_C = mean,
         propNlow_C = `2.5%`,
         propNhigh_C = `97.5%`)

# add expected available to consumed data
data <- data_gastro_sub %>% 
  mutate(propC_L = Lamprey / (Lamprey + Chinook),
         propC_C = Chinook / (Lamprey + Chinook)) %>% 
  filter(!is.na(propC_C)) %>% 
  select(propNmean_L, propNmean_C, propNlow_L, 
         propNlow_C, propNhigh_L, propNhigh_C, 
         propC_L, propC_C) %>% 
  mutate(id = 1:101) %>% 
  pivot_longer(cols = -id, 
               names_to = c("Type", "Fish_species"),
               values_to = "Prop",
               names_sep = "_") %>% 
  pivot_wider(id_cols = c(Fish_species, id),
              names_from = Type,
              values_from = Prop)

##################
# get model fits #
##################

############
# proportion L/S

L <- seq(0, 150000, 1000)
S <- seq(150000, 0, -1000)

set.seed(123)
index <- sample(nrow(posterior), 1000)

prey_L <- as.data.frame(matrix(NA, ncol = length(S), nrow = length(index)))
prey_S <- as.data.frame(matrix(NA, ncol = length(S), nrow = length(index)))

for (i in seq_along(index)) {
  
  params <- list(q = posterior[index[i], "q"],
                 h = c(posterior[index[i], "h[1]"],
                       posterior[index[i], "h[2]"]),
                 b = c(exp(posterior[index[i], "b_log[1]"]), 
                       exp(posterior[index[i], "b_log[2]"])))
  
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

# convert to numeric
prey_L$x <- as.numeric(prey_L$x)
prey_S$x <- as.numeric(prey_S$x)

# plot without MSFR fit, with abundance error bars
plot_nofit_errorbars_L <- ggplot() +
  geom_errorbarh(data = data[data$Fish_species == "L", ],
                 aes(xmin = propNlow, xmax = propNhigh,
                     y = propC), height = 0.01) +
  geom_point(data = data[data$Fish_species == "L", ],
             aes(x = propNmean, y = propC)) +
  ggtitle("Pacific lamprey") + 
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = NULL,
       y = "proportion of species in predator diet") +
  theme_minimal(base_family = "Arial")
plot_nofit_errorbars_S <- ggplot() +
  geom_errorbarh(data = data[data$Fish_species == "C", ],
                 aes(xmin = propNlow, xmax = propNhigh,
                     y = propC), height = 0.01) +
  geom_point(data = data[data$Fish_species == "C", ],
             aes(x = propNmean, y = propC)) +
  ggtitle("Chinook salmon") + 
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = NULL,
       y = "") +
  theme_minimal(base_family = "Arial")

final_nofit_errorbars <- 
  plot_nofit_errorbars_L + plot_nofit_errorbars_S + plot_layout(nrow = 1) +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0.5, family = "Arial",
                                  size = 11)
    ),
    caption = "estimated proportion of species in prey community"
  )
ggsave("figures/supplemental/MSFR_nofit_error.png", 
       final_nofit_errorbars, 
       dpi = 400, height = 4, width = 8)

# plot with MSFR fit
plot_L <- ggplot() +
  geom_line(data = prey_L,
            aes(x = as.numeric(x), y = y, group = id),
            alpha = 0.02, color = "black") +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_jitter(data = data[data$Fish_species == "L", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Pacific lamprey") + 
  labs(x = NULL,
       y = "proportion of species in predator diet") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "None")

plot_S <- ggplot() +
  geom_line(data = prey_S,
            aes(x = as.numeric(x), y = y, group = id),
            alpha = 0.02, color = "black") +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_jitter(data = data[data$Fish_species == "C", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Chinook salmon") +
  labs(x = NULL,
       y = NULL,
       color = "sample") +
  theme_minimal(base_family = "Arial")

final_plot <- 
  plot_L + plot_S + plot_layout(nrow = 1) +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0.5, family = "Arial",
                                  size = 11)
    ),
    caption = "estimated proportion of species in prey community"
  )
ggsave("figures/supplemental/MSFR_fit.png", final_plot,
       dpi = 400, height = 4, width = 8)

# plot MSFR fit, with full posterior colored
L_0.5 <- prey_L[prey_L$x == 0.5, c("id", "y")]
colnames(L_0.5) <- c("id", "color")
prey_L <- left_join(prey_L, L_0.5, by = "id")
prey_S <- left_join(prey_S, L_0.5, by = "id")

plot_L_color <- ggplot() +
  geom_line(data = prey_L,
            aes(x = as.numeric(x), y = y, group = id, color = color),
            alpha = 0.4) +
  scale_color_viridis() +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_jitter(data = data[data$Fish_species == "L", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Pacific lamprey") + 
  labs(x = NULL,
       y = "proportion of species in predator diet") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "None")

plot_S_color <- ggplot() +
  geom_line(data = prey_S,
            aes(x = as.numeric(x), y = y, group = id, color = color),
            alpha = 0.4) +
  scale_color_viridis() +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_jitter(data = data[data$Fish_species == "C", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Chinook salmon") + 
  labs(x = NULL,
       y = NULL) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "None")
final_plot_color <- 
  plot_L_color + plot_S_color + plot_layout(nrow = 1) +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0.5, family = "Arial",
                                  size = 11)
    ),
    caption = "estimated proportion of species in prey community"
  )
ggsave("figures/supplemental/MSFR_fit_fullpost.svg", 
       final_plot_color,
       dpi = 400, height = 4, width = 8)

