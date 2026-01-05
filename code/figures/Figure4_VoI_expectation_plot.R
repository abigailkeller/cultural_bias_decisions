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

###############
# expectation #
###############

# term 2
bethedge_utility <- rowMeans(utility)[which.max(rowMeans(utility))]

# bet-hedging strategy
bethedge_action <- alpha_change[which.max(rowMeans(utility))]

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

evpi_expect <- ggplot() +
  geom_line(aes(x = alpha_change, y = rowMeans(utility))) +
  geom_line(data = evpi_sub, aes(x = alpha_change, y = evpi,
                                  color = as.factor(sample))) +
  scale_color_manual(values = c("cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("1", "2", "3", "4")) +
  ggtitle("A.") +
  labs(x = "action (lamprey production)", 
       y = "expected utility* (abundance of salmon returns)",
       color = "posterior\nsample") +
  geom_text(label = "bet-hedging strategy",
            aes(x = 55, y = 125000), angle = -9,
            size = 3, family = "Arial") +
  new_scale_color() +
  geom_point(data = actions,
             aes(x = action, y = utility, color = sample),
             show.legend = FALSE, size = 3) +
  scale_color_manual(values = c("black", "cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("bet", "1", "2", "3", "4")) +
  scale_y_continuous(limits = c(min(evpi), max(evpi))) +
  scale_x_continuous(breaks = c(250 * 0.9 * 0.15, 
                                250 * 0.9 * 0.2, 
                                250 * 0.9 * 0.25),
                     limits = c(250 * 0.9 * 0.125,
                                250 * 0.9 * 0.275),
                     labels = c(
                       expression(alpha[L] == 0.15 * alpha[S]),
                       expression(alpha[L] == 0.2 * alpha[S]),
                       expression(alpha[L] == 0.25 * alpha[S])
                     )) +
  theme_minimal(base_family = "Arial")


violin_expect <- ggplot() +
  geom_violin(aes(x = rep(1, length(evpi)), y = evpi),
              alpha = 0.4, fill = "gray") +
  geom_point(aes(x = rep(1, length(actions$utility)),
                 y = actions$utility, color = as.factor(actions$sample)),
             size = 3) +
  scale_color_manual(values = c("black", "cornflowerblue", "green4", 
                                "goldenrod", "coral3"),
                     labels = c("bet", "1", "2", "3", "4")) +
  scale_y_continuous(limits = c(min(evpi), max(evpi))) +
  geom_segment(aes(x = 0.5, xend = 1.5, y = mean(evpi),
                   yend = mean(evpi)), color = "#662d91", linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = 1.5, 
                   y = bethedge_utility,
                   yend = bethedge_utility), 
               color = "black", linetype = "dashed") +
  ggtitle("B.") +
  annotate("text",
           x = 0.78, y = mean(evpi),
           label = "expected utility\nafter uncertainty\nhas been resolved", 
           color = "#662d91",
           vjust = -0.2, size = 2.5, family = "Arial") +
  annotate("text",
           x = 0.84, y = bethedge_utility - 8500,
           label = "bet-hedging\nstrategy", color = "black",
           vjust = -0.7, size = 2.5, family = "Arial") +
  annotate("text",
           x = 1.6, 
           y = mean(c(mean(evpi), bethedge_utility)),
           label = "EVPI", color = "#662d91", size = 3, 
           family = "Arial") +
  geom_segment(
    aes(x = 1.52, y = bethedge_utility, xend = 1.52, yend = mean(evpi)),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "#662d91"
  ) +
  # stat_brace(
  #   aes(x = c(2, 2),
  #       y = c(mean(evpi), bethedge_utility)),
  #   rotate = 90,
  #   outside = TRUE,
  #   width = 0.1,
  #   outerstart = 1.52,
  #   color = "purple3"
  # ) + 
  labs(x = "", 
       y = "maximized expected utility*") +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "None", text = element_text(family = "Arial"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(hjust = 0.4))


final_VOI_plot <- evpi_expect + violin_expect

ggsave("final_code/figures/Figure4_evpi.svg", final_VOI_plot, dpi = 400,
       height = 4.5, width = 9, bg = "white")

