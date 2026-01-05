library(tidyverse)
library(wesanderson)

pal <- wes_palette("Royal2")[1]
pal <- wes_palette("Royal2", 5, type = "discrete")

prob_plot <- ggplot() +
  geom_density(aes(x = rnorm(100000, 0, 1)),
               fill = pal[3], color = NA, alpha = 0.8) +
  geom_density(aes(x = rnorm(100000, -0.5, 0.3)),
               fill = pal[5], color = NA, alpha = 0.8) +
  labs(x = "ecological quantity\n of interest",
       y = "probability") +
  scale_x_continuous(limits = c(-2.7, 2.7)) +
  theme_minimal() +
  theme(axis.text = element_blank())
ggsave("final_code/figures/uneven_prob.svg", prob_plot, dpi = 400, 
       width = 2.5, height = 2)


x <- seq(-1, 1, 0.01)

ggplot() +
  geom_line(aes(x = x, y = x)) +
  geom_line(aes(x = x, y = x ^ 2))
  
