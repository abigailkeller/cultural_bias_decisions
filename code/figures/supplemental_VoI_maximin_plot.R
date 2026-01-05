library(tidyverse)
library(patchwork)
library(ggnewscale)
library(ggbrace)
library(grid)
library(ggtext)

# read in utility
utility <- readRDS("final_code/VOI/utility_fullposterior.rds")

# actions
alpha_change <- seq(250 * 0.9 * 0.05, 250 * 0.9 * 0.45, 0.5)

# posterior subset
sub <- readRDS("final_code/posterior_samples/index_sub.rds")



###########
# maximin #
###########

# term 2
bethedge_utility_min <- apply(utility, 
                              1, min)[which.max(apply(utility, 1, min))]

# bet-hedging strategy
bethedge_action_min <- alpha_change[which.max(apply(utility, 1, min))]

# perfect information
evpxi_index <- sub * 3 - 1
evpi <- rep(NA, ncol(utility))
evpi_policy <- rep(NA, ncol(utility))
for (i in 1:length(evpi)) {
  evpi_policy[i] <- which.max(utility[, i])
  evpi[i] <- max(utility[, i])
}
evpi_utility <- evpi[evpxi_index]
evpi_action <- alpha_change[evpi_policy[evpxi_index]]

evpi_sub <- data.frame(
  sample = c(rep(1, length(alpha_change)), rep(2, length(alpha_change)), 
             rep(3, length(alpha_change)), rep(4, length(alpha_change))), 
  alpha_change = rep(alpha_change, 4), 
  evpi = c(utility[, evpxi_index[1]], utility[, evpxi_index[2]], 
           utility[, evpxi_index[3]], utility[, evpxi_index[4]])
)

actions <- data.frame(
  sample = c("bet", "1", "2", "3", "4"),
  action = c(bethedge_action,
             evpi_action),
  utility = c(bethedge_utility,
              evpi_utility)
)
actions$sample <- factor(actions$sample, 
                         levels = c("bet", "1", "2", "3", "4"))

# partial information
evpxi_index_min <- seq(1, ncol(utility) + 1, 3)
evpxi_min <- rep(NA, length(evpxi_index_min) - 1)
evpxi_policy_min <- rep(NA, length(evpxi_index_min) - 1)
for (i in 1:(length(evpxi_index_min) - 1)) {
  sic <- utility[, evpxi_index_min[i]:(evpxi_index_min[i + 1] - 1)]
  evpxi_policy_min[i] <- which.max(apply(sic, 1, min))
  evpxi_min[i] <- apply(sic, 1, min)[evpxi_policy_min[i]]
}
evpxi_utility_min <- evpxi_min[sub]
evpxi_action_min <- alpha_change[evpxi_policy_min[sub]]

evpxi_sub_min <- data.frame(
  sample = c(rep(1, length(alpha_change)), rep(2, length(alpha_change)), 
             rep(3, length(alpha_change)), rep(4, length(alpha_change))), 
  alpha_change = rep(alpha_change, 4), 
  evpxi_min = c(
    apply(utility[, evpxi_index_min[sub[1]]:(evpxi_index_min[sub[1] + 1] - 1)], 
          1, min),
  apply(utility[, evpxi_index_min[sub[2]]:(evpxi_index_min[sub[2] + 1] - 1)], 
        1, min),
  apply(utility[, evpxi_index_min[sub[3]]:(evpxi_index_min[sub[3] + 1] - 1)], 
        1, min),
  apply(utility[, evpxi_index_min[sub[4]]:(evpxi_index_min[sub[4] + 1] - 1)], 
        1, min))
)

# optimal actions
actions_min <- data.frame(
  sample = c("bet", "1", "2", "3", "4"),
  action = c(bethedge_action_min,
             evpxi_action_min),
  utility = c(bethedge_utility_min,
              evpxi_utility_min)
)

actions_min$sample <- factor(actions_min$sample,
                             levels = c("bet", "1", "2", "3", "4"))

evpxi_plot_min <- ggplot() +
  geom_line(aes(x = alpha_change, y = apply(utility, 1, min))) +
  geom_line(data = evpxi_sub_min, aes(x = alpha_change, y = evpxi_min,
                                      color = as.factor(sample))) +
  scale_color_manual(values = c("cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("1", "2", "3", "4")) +
  ggtitle("A.") +
  labs(x = "action", y = "min utility*",
       color = "posterior\nsample") +
  geom_text(label = "bet-hedging strategy",
            aes(x = 70, y = 81500, angle = -17), 
            size = 3, family = "Arial") +
  new_scale_color() +
  scale_y_continuous(limits = c(min(apply(utility, 1, min)), max(evpxi_min))) +
  geom_point(data = actions_min,
             aes(x = action, y = utility, color = sample),
             show.legend = FALSE, size = 3) +
  scale_color_manual(values = c("black", "cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("bet", "1", "2", "3", "4")) +
  scale_x_continuous(breaks = c(250 * 0.9 * 0.1, 
                                250 * 0.9 * 0.25, 
                                250 * 0.9 * 0.4),
                     labels = c(
                       expression(alpha[L] == 0.1 * alpha[S]),
                       expression(alpha[L] == 0.25 * alpha[S]),
                       expression(alpha[L] == 0.4 * alpha[S])
                     )) +
  theme_minimal(base_family = "Arial")

violin_min <- ggplot() +
  geom_violin(aes(x = rep(1, length(evpxi_min)), y = evpxi_min),
              alpha = 0.4, fill = "gray") +
  geom_point(aes(x = rep(1, length(actions_min$utility)),
                 y = actions_min$utility, 
                 color = as.factor(actions_min$sample)), size = 3) +
  scale_color_manual(values = c("black", "cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("bet", "1", "2", "3", "4")) +
  scale_y_continuous(limits = c(min(apply(utility, 1, min)), max(evpxi_min))) +
  geom_segment(aes(x = 0.5, xend = 1.5, y = mean(evpxi_min),
                   yend = mean(evpxi_min)), color = "purple3", 
               linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = 1.5, 
                   y = bethedge_utility_min,
                   yend = bethedge_utility_min), 
               color = "black", linetype = "dashed") +
  annotate("text",
           x = 0.82, y = bethedge_utility_min,
           label = "bet-hedging strategy", color = "black",
           vjust = -0.7, size = 2.5, family = "Arial") +
  annotate("text",
           x = 0.8, y = mean(evpxi_min),
           label = "expected maximin utility\nafter F uncertainty\nhas been resolved", 
           color = "purple3",
           vjust = -0.1, size = 2.5, family = "Arial") +
  annotate("text",
           x = 1.7, 
           y = mean(c(mean(evpxi_min), bethedge_utility_min)),
           label = "EVPXI", color = "purple3", size = 3, 
           family = "Arial") +
  ggtitle("B.") +
  stat_brace(
    aes(x = c(2, 2),
        y = c(mean(evpxi_min), bethedge_utility_min)),
    rotate = 90,
    outside = TRUE,
    width = 0.1,
    outerstart = 1.52,
    color = "purple3"
  ) + 
  labs(x = "", 
       y = "maximin utility*") +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "None", text = element_text(family = "Arial"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(hjust = 0.4))

final_VOI_min_plot <- evpxi_plot_min + violin_min

ggsave("final_code/figures/supplemental/evpxi_maximin.png", final_VOI_min_plot, 
       dpi = 400,
       height = 4.5, width = 11, bg = "white")

