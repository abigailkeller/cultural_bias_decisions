library(tidyverse)
library(patchwork)

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
  mO_ddS = 2e7, # salmon ocean survival (function of lamprey density)
  mO_ddL = 5e-5, # lamprey ocean survival (function of salmon density)
  mO_surv = 8 / (8 + 125), # shared density-ind ocean survival
  # pinnipeds
  P_mean = 10000 # mean number of pinnipeds
)


# plot parasitism
dd <- c(9e7, 5e6, 3e6)
L <- seq(0, params$K_L, 10000)
surv <- params$mO_surv
para_data <- data.frame(
  L = L,
  dd1 = (1 - exp(-dd[1] / (L + 1e-6))) * surv,
  dd2 = (1 - exp(-dd[2] / (L + 1e-6))) * surv,
  dd3 = (1 - exp(-dd[3] / (L + 1e-6))) * surv
) %>% 
  pivot_longer(cols = -L,
               values_to = "survival",
               names_to = "dd_param")

para_plot <- ggplot(data = para_data) +
  geom_line(aes(x = L, y = survival, linetype = as.factor(dd_param)),
            size = 1) +
  labs(x = expression("lamprey juvenile abundance"), 
       y = expression("salmon ocean survival"), 
       linetype = expression("parasitism strength (" * D[S] * ")")) +
  scale_x_continuous(breaks = c(0, 1.5e6, 3e6),
                     labels = c("0", expression("0.5" * K[L]), 
                                expression(K[L]))) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"),
                        labels = c("9e7", "5e6", "3e6")) +
  theme_minimal(base_family = "Arial")
ggsave("final_code/figures/parasitism.svg",
       para_plot, dpi = 400, height = 3, width = 5)


