# plots inflow using VT and WVWA pressure transducers 
# CEB
# Oct 7 2024

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)

#load discharge data 
discharge <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/202/12/aae7888d68753b276d1623680f81d5de")

# separate DateTime column so you can grab just year of interest
discharge <- mutate(discharge, 
                    DateTime = ymd_hms(DateTime),
                    Year = year(DateTime),
                    Date = date(DateTime),
                    Hour = hour(DateTime))

# filter to just look at years we have sed trap data for
discharge <- discharge %>% 
  filter(Year>=2020)


# mean flow for every hour: first groups by the date
# so we get mean flow for everyday of 2020s and every hour of that day
discharge_mean <- discharge %>% 
  group_by(Date) %>% 
  summarise(mean_WVWA_cms = mean(coalesce(WVWA_Flow_cms, VT_Flow_cms), na.rm = TRUE)) %>% 
  mutate(Year = year(Date)) %>% 
  ungroup() %>% 
  mutate( # this creates dates that align with ISCO data
    Date_limited = case_when(
      Year == 2020 ~ if_else(Date >= as.Date("2020-05-11") & Date <= as.Date("2020-11-09"), Date, NA_Date_),
      Year == 2021 ~ if_else(Date >= as.Date("2021-04-16") & Date <= as.Date("2021-11-22"), Date, NA_Date_),
      Year == 2022 ~ if_else(Date >= as.Date("2022-04-19") & Date <= as.Date("2022-12-12"), Date, NA_Date_),
      Year == 2023 ~ if_else(Date >= as.Date("2023-05-08") & Date <= as.Date("2023-10-31"), Date, NA_Date_)))
# end date for 2023 is 12/4 but the ISCO data looks like it only plots until 10/31


## Add precip data
met_vars <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")
# changes to the proper format
met_vars$DateTime <- as.POSIXct(met_vars$DateTime)

met <- met_vars %>% 
  mutate(Year = year(DateTime)) %>% # mutate creates new Year variable from DateTime variable to be able to look at specific years
  mutate(Date = date(DateTime)) %>% # time doesnt matter so got rid of it
  filter(Year >= 2020) %>%
  select(Year, Date, Rain_Total_mm)

met_sum <- met %>% # groups by date and sums total rain for that day
  group_by(Date) %>%
  summarize(
    Rain_Total_mm = sum(Rain_Total_mm, na.rm = TRUE)) %>%
  ungroup()

met_discharge <- left_join(discharge_mean, met_sum, by="Date") %>% 
  mutate(Year= year(Date))

## this aligns with ISCO dates 
# Compute max, sum, and associated dates
met_discharge_summary_ISCO <- met_discharge %>%
  filter(!is.na(Date_limited)) %>%  # Keep only rows where Date_limited is not NA
  group_by(Year) %>%
  summarise(
    max_mean_WVWA_cms = max(mean_WVWA_cms, na.rm = TRUE),  # Find max mean_WVWA_cms
    date_max_mean_WVWA = Date[which.max(mean_WVWA_cms)],  # Get date of max mean_WVWA_cms
    sum_mean_WVWA_cms = sum(mean_WVWA_cms, na.rm = TRUE),  # Sum mean_WVWA_cms
    
    max_rain_total_mm = max(Rain_Total_mm, na.rm = TRUE),  # Find max Rain_Total_mm
    date_max_rain = Date[which.max(Rain_Total_mm)],  # Get date of max Rain_Total_mm
    sum_rain_total_mm = sum(Rain_Total_mm, na.rm = TRUE)   # Sum Rain_Total_mm
  )

write_csv(met_discharge_summary, 'Rain_Discharge_Summary_ISCO.csv')

## summary for entire year 
# Compute max, sum, and associated dates
met_discharge_summary <- met_discharge %>%
  group_by(Year) %>%
  summarise(
    max_mean_WVWA_cms = max(mean_WVWA_cms, na.rm = TRUE),  # Find max mean_WVWA_cms
    date_max_mean_WVWA = Date[which.max(mean_WVWA_cms)],  # Get date of max mean_WVWA_cms
    sum_mean_WVWA_cms = sum(mean_WVWA_cms, na.rm = TRUE),  # Sum mean_WVWA_cms
    
    max_rain_total_mm = max(Rain_Total_mm, na.rm = TRUE),  # Find max Rain_Total_mm
    date_max_rain = Date[which.max(Rain_Total_mm)],  # Get date of max Rain_Total_mm
    sum_rain_total_mm = sum(Rain_Total_mm, na.rm = TRUE)   # Sum Rain_Total_mm
  )

write_csv(met_discharge_summary, 'Rain_Discharge_Summary.csv')


##PLOT
ggplot(discharge_mean) +
  theme_bw() + theme(panel.grid = element_blank()) +
  facet_wrap(facets=vars(Year), nrow = 2, scales = 'free_x') +
  scale_y_continuous(limits=c(0,0.173)) +
  scale_x_date(date_labels = '%b', breaks = '1 month') +
  geom_point(aes(as.Date(Date_limited), mean_WVWA_cms), size = 0.75) +
  theme(plot.title = element_text(hjust = 0.5)) + xlab('Date') + ylab ('Discharge (cms)')

max(discharge_mean$mean_WVWA_cms, na.rm=T)


ggplot(met_discharge)+
  theme_bw()+ theme(panel.grid = element_blank())+
  facet_wrap(facets = vars(Year), nrow = 2, scales = 'free_x') +
  geom_bar(aes(Date_limited, Rain_Total_mm), stat = 'identity')+
  scale_y_reverse()+
  geom_point(aes(Date_limited, mean_WVWA_cms))

ggplot(met_discharge) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y.right = element_text(color = "blue")) +  # Make right axis label blue for clarity
  facet_wrap(vars(Year), nrow = 2, scales = 'free_x') +
  
  # Left axis: Discharge (mean_WVWA_cms) as points
  geom_point(aes(x = Date_limited, y = mean_WVWA_cms), color = "black", size = 0.75) + #change x = Date_limited for ISCO 
  # Right axis: Rain (Rain_Total_mm) as bars, scaled
  geom_bar(aes(x = Date_limited, y = Rain_Total_mm * 0.002149186), #change x = Date_limited for ISCO 
           stat = 'identity', fill = "blue", alpha = 0.5) +
  # Scale left axis normally, right axis transformed
  scale_y_continuous(
    name = "Discharge (cms)",  # Left y-axis label
    sec.axis = sec_axis(~ . / 0.002149186, name = "Rainfall (mm)"))+  # Right y-axis label
  # Optional: Improve axis readability
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = NULL)

max(met_discharge$mean_WVWA_cms, na.rm = T) / max(met_discharge$Rain_Total_mm, na.rm = T)
max(met_discharge$mean_WVWA_cms, na.rm=TRUE)
max(met_discharge$Rain_Total_mm, na.rm = T)


### ISCO Grey boxes 
## Data for full year
grey_boxes <- met_discharge %>%
  filter(!is.na(Date_limited)) %>%
  group_by(Year) %>%
  summarize(
    xmin = min(Date_limited),
    xmax = max(Date_limited)
  )

ggplot(met_discharge) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y.right = element_text(color = "blue")) +
  facet_wrap(vars(Year), nrow = 2, scales = 'free_x') +
  
  # Add shaded box per year
  geom_rect(data = grey_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE) +
  
  # Discharge points
  geom_point(aes(x = Date, y = mean_WVWA_cms), color = "black", size = 0.75) +
  
  # Rain bars
  geom_bar(aes(x = Date, y = Rain_Total_mm * 0.002149186),
           stat = 'identity', fill = "blue", alpha = 0.5) +
  
  # Axes
  scale_y_continuous(
    name = "Discharge (cms)",
    sec.axis = sec_axis(~ . / 0.002149186, name = "Rainfall (mm)")
  ) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = NULL)
