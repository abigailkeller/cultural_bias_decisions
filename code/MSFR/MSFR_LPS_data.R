library(tidyverse)

############
# LPS data #
############

file_list_LPS <- c(
  "data/Bonneville_dam_LPS/WA_AWS_2017.csv",
  "data/Bonneville_dam_LPS/WA_AWS_2019.csv",
  "data/Bonneville_dam_LPS/BI_AWS_2017.csv",
  "data/Bonneville_dam_LPS/BI_AWS_2019.csv",
  "data/Bonneville_dam_LPS/CI_ENT_2017.csv",
  "data/Bonneville_dam_LPS/CI_ENT_2019.csv",
  "data/Bonneville_dam_LPS/BradfordIsland_2017.csv",
  "data/Bonneville_dam_LPS/BradfordIsland_2019.csv",
  "data/Bonneville_dam_LPS/WashingtonShore_2017.csv",
  "data/Bonneville_dam_LPS/WashingtonShore_2019.csv"
)

# read all CSV files into a combined df
combined_data_LPS <- do.call(rbind, lapply(file_list_LPS, read.csv))

# convert date columns
combined_data_LPS$CountDate <- as.Date(combined_data_LPS$CountDate, 
                                       format = "%m/%d/%Y")

# add day and month columns
combined_data_LPS$Day <- format(combined_data_LPS$CountDate, "%d")
combined_data_LPS$Day <- as.numeric(combined_data_LPS$Day)
combined_data_LPS$Month <- format(combined_data_LPS$CountDate, "%m")
combined_data_LPS$Month <- as.numeric(combined_data_LPS$Month)

# create 3-day interval column
int <- seq(0, 31, 3) + 1
combined_data_LPS <- combined_data_LPS %>% 
  mutate(interval = case_when(
           Day >= int[1] & Day < int[2] ~ 1, Day >= int[2] & Day < int[3] ~ 2,
           Day >= int[3] & Day < int[4] ~ 3, Day >= int[4] & Day < int[5] ~ 4,
           Day >= int[5] & Day < int[6] ~ 5, Day >= int[6] & Day < int[7] ~ 6,
           Day >= int[7] & Day < int[8] ~ 7, Day >= int[8] & Day < int[9] ~ 8,
           Day >= int[9] & Day < int[10] ~ 9, Day >= int[10] ~ 10)) 

data_LPS_2017 <- combined_data_LPS %>% 
  filter(year == 2017) %>% 
  select(Location, CountDate, LPSLamprey, interval, Month) %>% 
  pivot_wider(id_cols = c(CountDate, interval, Month),
              names_from = Location,
              values_from = LPSLamprey) %>% 
  rename(
    Date = "CountDate",
    wa_aws_lps = "WA AWS LPS Mech",
    bi_aws_lps = "BI AWS LPS Mech",
    ci_aws = "Cascades Island AWS Traps",
    BI = "Bonneville Bradford Island",
    WA = "Bonneville Washington Shore"
  ) %>% 
  mutate(total = wa_aws_lps + bi_aws_lps + ci_aws + BI + WA)

data_LPS_2019 <- combined_data_LPS %>% 
  filter(year == 2019) %>% 
  select(Location, CountDate, LPSLamprey, interval, Month) %>% 
  pivot_wider(id_cols = c(CountDate, interval, Month),
              names_from = Location,
              values_from = LPSLamprey) %>% 
  rename(
    Date = "CountDate",
    wa_aws_lps = "WA AWS LPS Mech",
    bi_aws_lps = "BI AWS LPS Mech",
    ci_aws = "Cascades Island AWS Traps",
    BI = "Bonneville Bradford Island",
    WA = "Bonneville Washington Shore"
  ) %>% 
  mutate(total = wa_aws_lps + bi_aws_lps + ci_aws + BI + WA)

###############
# Window data #
###############

# get adult passage data
file_list_adult_passage <- c(
  "data/Bonneville_dam_returns/adultdaily_2017.csv",
  "data/Bonneville_dam_returns/adultdaily_2019.csv"
)

# read all CSV files into a combined df
combined_data_passage <- do.call(rbind, 
                                 lapply(file_list_adult_passage, read.csv))

# convert date columns
combined_data_passage$Date <- as.Date(combined_data_passage$Date, 
                                      format = "%Y-%m-%d")

# remove rows with na
combined_data_passage <- combined_data_passage[
  !is.na(combined_data_passage$Date), ]

# create a new column for the year extracted from the date
combined_data_passage$Year <- format(combined_data_passage$Date, "%Y")


##########################
# combine window and LPS #
##########################

data_LPS_2017 <- left_join(data_LPS_2017, 
                           combined_data_passage[, c("Date", "Lmpry")], 
                           by = "Date") 

data_LPS_2019 <- left_join(data_LPS_2019, 
                           combined_data_passage[, c("Date", "Lmpry")], 
                           by = "Date") 

# get interval proportion
interval_prop_2017 <- data_LPS_2017 %>% 
  group_by(interval, Month) %>% 
  summarize(total = sum(total),
            Lmpry = sum(Lmpry),
            date = first(Date)) %>% 
  mutate(interval_prop = Lmpry/ total,
         interval_prop = case_when(
           total < 20 ~ NA,
           TRUE ~ interval_prop
         ))

interval_prop_2019 <- data_LPS_2019 %>% 
  group_by(interval, Month) %>% 
  summarize(total = sum(total),
            Lmpry = sum(Lmpry),
            date = first(Date)) %>% 
  mutate(interval_prop = Lmpry/ total,
         interval_prop = case_when(
           total < 20 ~ NA,
           TRUE ~ interval_prop
         ))

# combine 2017 and 2019
interval_prop <- rbind(interval_prop_2017, interval_prop_2019)

# remove NA
interval_prop <- interval_prop[!is.na(interval_prop$interval_prop), ]

# get 2017 - 2019 total lamprey
total_lamprey <- rbind(data_LPS_2017[, c("Date", "total")],
                       data_LPS_2019[, c("Date", "total")])

# save as model data
saveRDS(interval_prop$Lmpry,
        "data/model_data/C_V.rds")
saveRDS(interval_prop$total,
        "data/model_data/C_T.rds")
saveRDS(total_lamprey,
        "data/model_data/total_lamprey.rds")
