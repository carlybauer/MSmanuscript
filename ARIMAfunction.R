# function for processing metal, environmental, precip, and flow data for ARIMA
# CEB
# 19Feb2025
# does BVR and FCR for 2020 - 2023 - separates by year
# have to give it start and end dates of interest
# comment out flow info for BVR data if you want or ignore when running ARIMA

## FOR FCR: Depths = 1.6 or 8
##          CTD = 1.1 - 2.1 or 7.5 - 8.5
## FOR BVR: Depths = 3 or 9
##          CTD = 2.5 - 3.5 or 8.5 - 9.5
## Use the same precip data for both
###FCR
## 2020 dates: 01-31-2020, 11-06-2020
## 2023 dates: 04-04-2023 12-05-2023
###BVR
## 2020 dates: 03-30-2020, 11-16-2020
## 2023 dates: 03-21-2023 12-05-2023


process_reservoir_data <- function(reservoir, year, start_date, end_date) {
  # library(dplyr)
  # library(tidyr)
  # library(lubridate)
  # library(tsibble)
  # library(forecast)
  
# Define depth ranges for metals and environmental data
# options for hypolimnion 
  depth_ranges <- list(
    FCR = list(metal = c(8), env = list(c(7.5, 8.5))),
    BVR = list(metal = c(9), env = list(c(8.5, 9.5)))
  )
  
  metal_depths <- depth_ranges[[reservoir]]$metal
  env_depths <- depth_ranges[[reservoir]]$env
  
# Load Metals Data
  metal <- read_csv("https://raw.githubusercontent.com/carlybauer/Reservoirs/refs/heads/master/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLmetals/2024/Metals_2014_2023.csv")
  metal$DateTime <- as.POSIXct(metal$DateTime)
  
  # # this gets rid of any values that we flagged as abnormally high in the data
  # # not including this because it gets rid of a lot of BVR values
  # metalsflag <- metal %>% 
  #   filter(
  #     Flag_TAl_mgL != 8 & Flag_SAl_mgL != 8 & Flag_TAl_mgL != 68 & Flag_SAl_mgL != 68,
  #     Flag_TFe_mgL != 8 & Flag_SFe_mgL != 8 & Flag_TFe_mgL != 68 & Flag_SFe_mgL != 68,
  #     Flag_TMn_mgL != 8 & Flag_SMn_mgL != 8 & Flag_TMn_mgL != 68 & Flag_SMn_mgL != 68,
  #     Flag_TCu_mgL != 8 & Flag_SCu_mgL != 8 & Flag_TCu_mgL != 68 & Flag_SCu_mgL != 68,
  #     Flag_TSr_mgL != 8 & Flag_SSr_mgL != 8 & Flag_TSr_mgL != 68 & Flag_SSr_mgL != 68,
  #     Flag_TBa_mgL != 8 & Flag_SBa_mgL != 8 & Flag_TBa_mgL != 68 & Flag_SBa_mgL != 68)
  
  
  metals <- metal %>%
    mutate(Year = year(DateTime), 
           Date_real = date(DateTime)) %>%
    filter(Reservoir == reservoir, 
           Site == 50, 
           Year == year, 
           Depth_m %in% metal_depths) %>%
    select(Date_real, Year, Depth_m, TFe_mgL, SFe_mgL, TMn_mgL, SMn_mgL, SBa_mgL, TBa_mgL, 
           SCu_mgL, TCu_mgL, SSr_mgL, TSr_mgL, SAl_mgL, TAl_mgL) %>%
    drop_na()
  
 
  
  metals <- as_tsibble(metals) %>% #recognize data as time series 
    mutate(Week = week(Date_real))
  
# Create weekly date sequence
  weekly_dates <- data.frame(Date_fake = seq.Date(from = as.Date(start_date), to = as.Date(end_date), by = "week")) %>%
    mutate(Year = year(Date_fake), Week = week(Date_fake))
  
# Merge Metals with Weekly Dates
  metal_weekly <- weekly_dates %>% 
    left_join(metals, by = c("Year", "Week"))
  
# get rid of duplicates by taking the mean 
  metal_weekly <- metal_weekly %>% 
    group_by(Date_fake) %>% 
    summarise(across(where(is.numeric),  \(x) mean(x, na.rm = TRUE))) %>%
    ungroup()
    
  
# Interpolate missing values
  metal_interp <- metal_weekly %>% 
    as_tsibble(index = Date_fake) %>% 
    dplyr::mutate(dplyr::across(TFe_mgL:TAl_mgL, ~ forecast::na.interp(.)))
  

# Load Environmental Data
  
  # Handheld sensor data load and process 
  hand_env_vars <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/198/13/e50a50d062ee73f4d85e4f20b360ce4f")

  hand_env_vars$DateTime <- as.POSIXct(hand_env_vars$DateTime)
  
  hand_env <- hand_env_vars %>%
    mutate(Year = year(DateTime),
           Date_real = date(DateTime)) %>%
    filter(Reservoir == reservoir, Site == 50, Year == year, Depth_m %in% metal_depths) %>% 
    select(Year, Date_real, Depth_m, DO_mgL, Temp_C, pH) %>% 
    mutate(Week = week(Date_real))
  
  hand_env_weekly <- weekly_dates %>%
    left_join(hand_env, by = c("Year", "Week"))
  
  #CTD Data load and process
  env_vars <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a")
  
  env_vars$DateTime <- as.POSIXct(env_vars$DateTime)
  
  env <- env_vars %>% 
    mutate(Year = year(DateTime), 
           Date_real = date(DateTime)) %>% 
    filter(Reservoir == reservoir, Site == 50, Year == year, 
           Depth_m >= depth_ranges[[reservoir]]$env[[1]][1] & Depth_m <= depth_ranges[[reservoir]]$env[[1]][2]) %>%
    select(Year, Date_real, Depth_m, DO_mgL, Temp_C, Turbidity_NTU, pH)
  
  env_avg <- env %>% 
    group_by(Date_real) %>% 
    summarize(
      DO_mgL = mean(DO_mgL, na.rm = TRUE),
      Temp_C = mean(Temp_C, na.rm = TRUE),
      #DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
      Turbidity_NTU = mean(Turbidity_NTU, na.rm = TRUE),
      pH = mean(pH, na.rm = TRUE)) 
  
  env_avg <- env_avg %>% 
    mutate(Week = week(Date_real), 
           Year = year(Date_real))
  
  env_weekly_0 <- weekly_dates %>% 
    left_join(env_avg, by = c("Year", "Week"))
  
  env_weekly <- env_weekly_0 %>%
    left_join(hand_env_weekly %>% 
                select(Date_fake, DO_mgL, Temp_C, pH), by = "Date_fake", relationship = "many-to-many") %>%
    mutate(
      DO_mgL = coalesce(DO_mgL.x, DO_mgL.y),
      Temp_C = coalesce(Temp_C.x, Temp_C.y),
      pH = coalesce(pH.x, pH.y)) %>%
    select(-DO_mgL.x, -DO_mgL.y, -Temp_C.x, -Temp_C.y, -pH.x, -pH.y)
  
  env_weekly <- env_weekly %>% 
    group_by(Date_fake) %>% 
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup()
  
  print(head(env_weekly)) # Check for unexpected NA values
  print(colnames(env_weekly)) # Ensure DO_mgL exists
  
  
  env_interp <- env_weekly %>% 
    as_tsibble(index = Date_fake) %>% 
    dplyr::mutate(dplyr::across(Turbidity_NTU:pH, ~ forecast::na.interp(.)))

# Load and Process Precipitation Data
  met_vars <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")
  
  met_vars$DateTime <- as.POSIXct(met_vars$DateTime)
  
  met <- met_vars %>% 
    mutate(Year = year(DateTime), 
           Date_real = date(DateTime)) %>%
    filter(Year == year) %>%
    select(Year, Date_real, Rain_Total_mm) %>%
    group_by(Date_real) %>% 
    summarise(Rain_Total_mm = sum(Rain_Total_mm, na.rm = TRUE)) %>%  # sums rain for the day
    ungroup()
  
  met_sum <- met %>% 
    mutate(Week = week(Date_real), 
           Year = year(Date_real)) %>%
    group_by(Year, Week) %>% 
    summarise(Rain_Total_mm = sum(Rain_Total_mm, na.rm = TRUE), #sums rain for the week
              Date_real = min(Date_real)) %>%  # takes first date in week to save date
    mutate(Lag1_Rain_Total_mm = lag(Rain_Total_mm, n = 1)) %>%  #lags rain by 1 week
    ungroup()
  
  met_weekly <- weekly_dates %>% 
    left_join(met_sum, by = c("Year", "Week")) %>% 
    as_tsibble(index = Date_fake)

# Load and Process Flow Data
  q_vars <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/202/12/aae7888d68753b276d1623680f81d5de")
  
  q_vars$DateTime <- as.POSIXct(q_vars$DateTime)
  
  flow <- q_vars %>% 
    mutate(Year = year(DateTime), 
           Date_real = date(DateTime)) %>%
    filter(Year == year, Reservoir == "FCR", Site == 100) %>%
    select(Year, Date_real, WVWA_Flow_cms, VT_Flow_cms)
  
  q_mean <- flow %>% 
    group_by(Date_real) %>% 
    summarise(MeanDay_WVWA_cms = mean(coalesce(WVWA_Flow_cms, VT_Flow_cms), na.rm = TRUE), #averages q for the date using WVWA and VT transducers
              SumDay_WVWA_cms = sum(coalesce(WVWA_Flow_cms, VT_Flow_cms), na.rm = TRUE)) %>%  # sums q for date using VT and WVWA, finds first non-missing element and uses it
    ungroup()
  
  q_sum <- q_mean %>% 
    mutate(Week = week(Date_real), 
           Year = year(Date_real)) %>%
    group_by(Week, Year) %>% 
    summarise(Q_Mean_cms = sum(MeanDay_WVWA_cms, na.rm = TRUE), #sums q for the week using mean of day, this is original way of calculating 
              Q_Sum_cms = sum(SumDay_WVWA_cms, na.rm = TRUE), #sums q for week using sum of the day, new calculation, doesnt help 
              Date_real = min(Date_real)) %>% 
    ungroup() %>% 
    mutate(Lag1_Q_Mean_cms = lag(Q_Mean_cms, n = 1),
           Lag1_Q_Sum_cms = lag(Q_Sum_cms, n = 1)) #lags q by 1 week 
  
  q_weekly <- weekly_dates %>% 
    left_join(q_sum, by = c("Year", "Week")) %>% 
    as_tsibble(index = Date_fake)

# Merge all data
  merged_df <- env_interp %>% 
    left_join(metal_interp, by = "Date_fake") %>%
    left_join(met_weekly, by = "Date_fake") %>%
    left_join(q_weekly, by = "Date_fake") %>%
    select(-starts_with("Year"), -starts_with("Week"), -Depth_m, -Date_real.x, -Date_real.y)

#   return(merged_df)
# Standardize data
  merged_standard <- merged_df %>%
    mutate(across(where(is.numeric) & !c(Date_fake), scale)) %>%
    as.data.frame(lapply(merged_standard, as.vector)) %>%
    as_tsibble(index = Date_fake)

  return(merged_standard)
}
