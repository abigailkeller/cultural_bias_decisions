library(tidyverse)
library(patchwork)
library(ggnewscale)
library(MCMCvis)
library(viridis)

source("final_code/code/timeseries_functions.R")

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
                     readRDS("final_code/posterior_samples/MSFR_posterior_1day.rds"))
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

# get data densities
data_density_L <- data[data$Fish_species == "L", ] %>% 
  mutate(group = case_when(propNmean < 0.05 ~ 0.025,
                           propNmean >= 0.05 & propNmean < 0.1 ~ 0.075,
                           propNmean >= 0.1 & propNmean < 0.15 ~ 0.125,
                           propNmean >= 0.15 & propNmean < 0.2 ~ 0.175,
                           propNmean >= 0.2 & propNmean < 0.25 ~ 0.225,
                           propNmean >= 0.25 & propNmean < 0.3 ~ 0.275,
                           )) %>% 
  group_by(group) %>% 
  summarise(count = n())
data_density_S <- data[data$Fish_species == "C", ] %>% 
  mutate(group = case_when(propNmean < 0.7 ~ 0.675,
                           propNmean >= 0.7 & propNmean < 0.75 ~ 0.725,
                           propNmean >= 0.75 & propNmean < 0.8 ~ 0.775,
                           propNmean >= 0.8 & propNmean < 0.85 ~ 0.825,
                           propNmean >= 0.85 & propNmean < 0.9 ~ 0.875,
                           propNmean >= 0.9 & propNmean < 0.95 ~ 0.925,
                           propNmean >= 0.95 & propNmean < 1 ~ 0.975,
  )) %>% 
  group_by(group) %>% 
  summarise(count = n())

datadens_L <- ggplot(data_density_L) +
  geom_tile(aes(x = group, y = 1, alpha = log(count + 1e-6)),
            fill = "cornflowerblue") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(legend.position = "None",
        axis.title = element_blank(), 
        axis.text = element_blank())

datadens_S <- ggplot(data_density_S) +
  geom_tile(aes(x = group, y = 1, alpha = log(count + 1e-6)),
            fill = "cornflowerblue") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_alpha_continuous(breaks = c(log(3), log(10), log(25), log(50)),
                         labels = c("3", "10", "25", "50")) +
  labs(alpha = "data\ndensity") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        legend.title = element_text(hjust = 0.5))

ggsave("final_code/figures/data_density.svg", datadens_L + datadens_S,
       dpi = 400, width = 6, height = 3)

# plot without MSFR fit
plot_nofit_L <- ggplot() +
  geom_vline(xintercept = seq(0, 1, length.out = 5), 
             color = "grey90", linewidth = 0.3) + 
  geom_vline(xintercept = seq(0, 1, length.out = 9)[-c(1, 9)], 
             color = "grey92", linewidth = 0.2) + 
  geom_hline(yintercept = seq(0, 1, length.out = 5), 
             color = "grey90", linewidth = 0.3) + 
  geom_hline(yintercept = seq(0, 1, length.out = 9)[-c(1, 9)], 
             color = "grey92", linewidth = 0.2) + 
  geom_jitter(data = data[data$Fish_species == "L", ],
             aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Pacific lamprey") + 
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(x = NULL,
       y = "proportion of species\nin predator diet") +
  theme_minimal(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12))
plot_nofit_S <- ggplot() +
  geom_vline(xintercept = seq(0, 1, length.out = 5), 
             color = "grey90", linewidth = 0.3) + 
  geom_vline(xintercept = seq(0, 1, length.out = 9)[-c(1, 9)], 
             color = "grey92", linewidth = 0.2) + 
  geom_hline(yintercept = seq(0, 1, length.out = 5), 
             color = "grey90", linewidth = 0.3) + 
  geom_hline(yintercept = seq(0, 1, length.out = 9)[-c(1, 9)], 
             color = "grey92", linewidth = 0.2) + 
  geom_jitter(data = data[data$Fish_species == "C", ],
             aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Chinook salmon") + 
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  labs(x = NULL,
       y = NULL) +
  theme_minimal(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12))

final_nofit <- 
  plot_nofit_L + plot_nofit_S + plot_layout(nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 1), "cm")) +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0.5, family = "Arial",
                                  size = 11)
    ),
    caption = "estimated proportion of species in prey community"
  )

# ggsave("final_code/figures/MSFR_nofit.png", final_nofit, 
#        dpi = 400, height = 4, width = 8)
ggsave("final_code/figures/MSFR_nofit.svg", final_nofit, 
       dpi = 400, height = 4, width = 8)

# plot MSFR fit with posterior subset highlighted
sub <- readRDS("final_code/posterior_samples/index_sub.rds")

prey_L_ordered <- prey_L %>%
  filter(id %in% sub) %>%
  mutate(id = factor(id, levels = sub))
prey_S_ordered <- prey_S %>%
  filter(id %in% sub) %>%
  mutate(id = factor(id, levels = sub))

plot_L_sub <- ggplot() +
  geom_vline(xintercept = seq(0, 1, length.out = 5),
             color = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = seq(0, 1, length.out = 9)[-c(1, 9)],
             color = "grey92", linewidth = 0.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 5),
             color = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = seq(0, 1, length.out = 9)[-c(1, 9)],
             color = "grey92", linewidth = 0.2) +
  geom_jitter(data = data[data$Fish_species == "L", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Pacific lamprey") + 
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  geom_line(data = prey_L,
            aes(x = as.numeric(x), y = y, group = id),
            alpha = 0.04, color = "black") +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_line(data = prey_L_ordered,
            aes(x = as.numeric(x), y = y, group = id, color = as.factor(id)),
            linewidth = 1) +
  scale_color_manual(values = c("cornflowerblue", "green4", 
                                "goldenrod", "coral3")) +
  labs(x = NULL,
       y = "proportion of species\nin predator diet",
       color = "sample", linetype = "") +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "None") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12))

plot_S_sub <- ggplot() +
  geom_vline(xintercept = seq(0, 1, length.out = 5),
             color = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = seq(0, 1, length.out = 9)[-c(1, 9)],
             color = "grey92", linewidth = 0.2) +
  geom_hline(yintercept = seq(0, 1, length.out = 5),
             color = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = seq(0, 1, length.out = 9)[-c(1, 9)],
             color = "grey92", linewidth = 0.2) +
  geom_jitter(data = data[data$Fish_species == "C", ],
              aes(x = propNmean, y = propC), alpha = 0.5, size = 3) +
  ggtitle("Chinook salmon") + 
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  geom_line(data = prey_S,
            aes(x = as.numeric(x), y = y, group = id),
            alpha = 0.04, color = "black") +
  geom_line(aes(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), linetype = "dashed") +
  geom_line(data = prey_S_ordered,
            aes(x = as.numeric(x), y = y, group = id, color = as.factor(id)),
            linewidth = 1) +
  scale_color_manual(values = c("cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("1", "2", "3", "4")) +
  labs(x = NULL,
       y = NULL,
       color = "posterior\nsample") +
  theme_minimal(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12))

final_plot_sub <- 
  plot_L_sub + plot_S_sub + plot_layout(nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 1), "cm")) +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0.5, family = "Arial",
                                  size = 11)
    ),
    caption = "estimated proportion of species in prey community"
  )
# ggsave("final_code/figures/MSFR_fit_sub.png", final_plot_sub,
#        dpi = 400, height = 4, width = 8)
ggsave("final_code/figures/MSFR_fit_sub.svg", final_plot_sub,
       dpi = 400, height = 4, width = 8)


