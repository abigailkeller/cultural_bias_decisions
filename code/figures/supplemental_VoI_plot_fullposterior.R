library(tidyverse)
library(patchwork)
library(ggnewscale)
library(ggbrace)
library(grid)
library(ggtext)
library(viridis)

# read in utility
utility <- readRDS("final_code/VOI/utility_fullposterior.rds")

# actions
alpha_change <- seq(250 * 0.9 * 0.05, 250 * 0.9 * 0.45, 0.5)

# get utility of mean parasitism
utility_sub <- utility[, seq(2, 3000, 3)]

###############
# expectation #
###############

# term 2
bethedge_utility <- rowMeans(utility)[which.max(rowMeans(utility))]

# bet-hedging strategy
bethedge_action <- alpha_change[which.max(rowMeans(utility))]

# perfect information
evpi <- rep(NA, ncol(utility_sub))
evpi_policy <- rep(NA, ncol(utility_sub))
for (i in 1:length(evpi)) {
  evpi_policy[i] <- which.max(utility_sub[, i])
  evpi[i] <- max(utility_sub[, i])
}
evpi_action <- alpha_change[evpi_policy]

evpi_change <- as.data.frame(matrix(NA, nrow = length(alpha_change) * 
                                      ncol(utility) / 3,
                                     ncol = 3))
colnames(evpi_change) <- c("id", "alpha_change", "evpi")
for (i in 1:ncol(utility_sub)) {
  row_seq <- seq((i - 1) * length(alpha_change) + 1, 
                 i * length(alpha_change), 1)
  evpi_change[row_seq, "id"] <- rep(i, length(alpha_change))
  evpi_change[row_seq, "alpha_change"] <- alpha_change
  evpi_change[row_seq, "evpi"] <- utility_sub[, i]
}

index_color <- readRDS("final_code/VOI/index_color.rds")
evpi_change <- left_join(evpi_change, index_color, by = "id")

# optimal actions
actions <- data.frame(
  id = 1:1000,
  action = evpi_action,
  utility = evpi
)
actions <- left_join(actions, index_color, by = "id")

evpi_expect <- ggplot() +
  geom_line(aes(x = alpha_change, y = rowMeans(utility))) +
  geom_line(data = evpi_change, aes(x = alpha_change, y = evpi,
                                     group = id, color = color),
            alpha = 0.2) +
  geom_point(data = actions,
             aes(x = action, y = utility, color = color),
             show.legend = FALSE, size = 2, alpha = 0.5) +
  geom_line(aes(x = alpha_change, y = rowMeans(utility))) +
  scale_color_viridis() +
  ggtitle("Full posterior: expectation") +
  labs(x = "action", 
       y = "expected utility*") +
  geom_text(label = "bet-hedging strategy",
            aes(x = 60, y = 124000), angle = -10,
            size = 3, family = "Arial") +
  geom_point(aes(x = bethedge_action, y = bethedge_utility),
             show.legend = FALSE, size = 2) +
  scale_x_continuous(breaks = c(250 * 0.9 * 0.1, 
                                250 * 0.9 * 0.2, 
                                250 * 0.9 * 0.3),
                     limits = c(250 * 0.9 * 0.075,
                                250 * 0.9 * 0.35),
                     labels = c(
                       expression(alpha[L] == 0.1 * alpha[S]),
                       expression(alpha[L] == 0.2 * alpha[S]),
                       expression(alpha[L] == 0.3 * alpha[S])
                     )) +
  theme_minimal(base_family = "Arial") +
  theme(legend.position = "None",
        element_text(size = 11))


ggsave("final_code/figures/supplemental/evpi_fullposterior.svg", evpi_expect, 
       dpi = 400, height = 4, width = 4, bg = "white")

