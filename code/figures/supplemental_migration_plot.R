library(tidyverse)
library(patchwork)

# read in bonneville dam data
file_list <- list.files(path = "data/Bonneville_dam_returns", 
                        pattern = "adultdaily_", full.names = TRUE)

# read all CSV files into a combined df
combined_data <- do.call(rbind, lapply(file_list, read.csv))

# convert date columns
combined_data$Date <- as.Date(combined_data$Date, 
                              format = "%Y-%m-%d")

# create a new column for the year and month extracted from the date
combined_data$Year <- format(combined_data$Date, "%Y")
combined_data$Month <- format(combined_data$Date, "%m")
combined_data$MonthDay <- format(combined_data$Date, "%m-%d")
combined_data$MonthDay <- as.Date(combined_data$MonthDay, 
                                  format = "%m-%d")

# subset to spring months (April and May)
combined_spring <- combined_data[c(combined_data$MonthDay > as.Date("03-31", 
                                                                    format = "%m-%d") &
                                     combined_data$MonthDay < as.Date("06-15", 
                                                                      format = "%m-%d")), ]

# remove years without lamprey data
total_lamprey_annual <- combined_spring %>% 
  group_by(Year) %>% 
  summarize(total_lamp = sum(Lmpry, na.rm = TRUE))
years <- total_lamprey_annual[which(total_lamprey_annual$total_lamp > 0), ]$Year
combined_spring <- combined_spring[combined_spring$Year %in% years, ]

plots <- list()
for(i in 1:length(years)) {
  data_sub <- combined_spring[!is.na(combined_spring$Year) & 
                                combined_spring$Year == years[i], ] %>% 
    select(Date, Chin,
           Lmpry) %>% 
    pivot_longer(cols = -Date,
                 names_to = 'Species',
                 values_to = 'N')
  
  plots[[i]] <- ggplot()+
    geom_line(data=data_sub,
              aes(x=Date, y = N,
                  color = Species))+
    scale_color_manual(labels = c('Chinook', 'Lamprey'),
                       values = c('dodgerblue', 'violetred'))+
    labs(x='Date', y='Count in Bonneville Dam')+
    ggtitle(years[i]) +
    theme_minimal()
}

combined_spring$Chin[is.na(combined_spring$Chin)] <- 0
combined_spring$Lmpry[is.na(combined_spring$Lmpry)] <- 0
combined_spring$Shad[is.na(combined_spring$Shad)] <- 0

# summarize early and later years
# early
data_early <- combined_spring[combined_spring$Year %in% 
                                years[as.numeric(years) < 1970], ] %>% 
  group_by(MonthDay) %>% 
  summarize(mean_S = mean(Chin),
            lower_S = max(0, mean_S - 0.5 * sd(Chin)),
            upper_S = mean_S + 0.5 * sd(Chin),
            mean_L = mean(Lmpry),
            lower_L = max(0, mean_L - 0.5 * sd(Lmpry)),
            upper_L = mean_L + 0.5 * sd(Lmpry),
            mean_Sh = mean(Shad),
            lower_Sh = max(0, mean_Sh - 0.5 * sd(Shad)),
            upper_Sh = mean_Sh + 0.5 * sd(Shad)) %>% 
  pivot_longer(
    cols = -MonthDay,
    names_to = c("stat", "species"),
    names_sep = "_",
    values_to = "value"
  ) %>% 
  pivot_wider(id_cols = c(MonthDay, species),
              values_from = value,
              names_from = stat)

early_all <- ggplot(data = data_early) +
  geom_ribbon(aes(x = MonthDay,
                  ymin = lower, ymax = upper, fill = species),
              alpha = 0.2) + 
  geom_line(aes(x = MonthDay, y = mean, color = species)) +
  scale_color_manual(values = c("green4", "dodgerblue", "violetred"),
                     labels = c("Lamprey", "Chinook", "Shad")) +
  scale_fill_manual(values = c("green4", "dodgerblue", "violetred"),
                    labels = c("Lamprey", "Chinook", "Shad")) +
  labs(x = "Date", y = "Count") +
  ggtitle("Bonneville dam: 1946 - 1969") +
  theme_minimal()
early_noshad <- ggplot(data = data_early[data_early$species != "Sh", ]) +
  geom_ribbon(aes(x = MonthDay,
                  ymin = lower, ymax = upper, fill = species),
              alpha = 0.2) + 
  geom_line(aes(x = MonthDay, y = mean, color = species)) +
  scale_color_manual(values = c("green4", "dodgerblue"),
                     labels = c("Lamprey", "Chinook")) +
  scale_fill_manual(values = c("green4", "dodgerblue"),
                    labels = c("Lamprey", "Chinook")) +
  labs(x = "Date", y = "Count") +
  ggtitle("Bonneville dam: 1946 - 1969") +
  theme_minimal()

# later
data_later <- combined_spring[combined_spring$Year %in% 
                                years[as.numeric(years) > 1970], ] %>% 
  group_by(MonthDay) %>% 
  summarize(mean_S = mean(Chin),
            lower_S = max(0, mean_S - 0.5 * sd(Chin)),
            upper_S = mean_S + 0.5 * sd(Chin),
            mean_L = mean(Lmpry),
            lower_L = max(0, mean_L - 0.5 * sd(Lmpry)),
            upper_L = mean_L + 0.5 * sd(Lmpry),
            mean_Sh = mean(Shad),
            lower_Sh = max(0, mean_Sh - 0.5 * sd(Shad)),
            upper_Sh = mean_Sh + 0.5 * sd(Shad)) %>% 
  pivot_longer(
    cols = -MonthDay,
    names_to = c("stat", "species"),
    names_sep = "_",
    values_to = "value"
  ) %>% 
  pivot_wider(id_cols = c(MonthDay, species),
              values_from = value,
              names_from = stat)

later_all <- ggplot(data = data_later) +
  geom_ribbon(aes(x = MonthDay,
                  ymin = lower, ymax = upper, fill = species),
              alpha = 0.2) + 
  geom_line(aes(x = MonthDay, y = mean, color = species)) +
  scale_color_manual(values = c("green4", "dodgerblue", "violetred"),
                     labels = c("Lamprey", "Chinook", "Shad")) +
  scale_fill_manual(values = c("green4", "dodgerblue", "violetred"),
                    labels = c("Lamprey", "Chinook", "Shad")) +
  labs(x = "Date", y = "Count") +
  ggtitle("Bonneville dam: 1999 - 2024") +
  theme_minimal()
later_noshad <- ggplot(data = data_later[data_later$species != "Sh", ]) +
  geom_ribbon(aes(x = MonthDay,
                  ymin = lower, ymax = upper, fill = species),
              alpha = 0.2) + 
  geom_line(aes(x = MonthDay, y = mean, color = species)) +
  scale_color_manual(values = c("green4", "dodgerblue"),
                     labels = c("Lamprey", "Chinook")) +
  scale_fill_manual(values = c("green4", "dodgerblue"),
                    labels = c("Lamprey", "Chinook")) +
  labs(x = "Date", y = "Count") +
  ggtitle("Bonneville dam: 1999 - 2024") +
  theme_minimal()


# combine
noshad_plot <- early_noshad + later_noshad + plot_layout(ncol = 1)
ggsave("figures/supplemental/migration.png", 
       noshad_plot, dpi = 400,
       width = 5, height = 6)
