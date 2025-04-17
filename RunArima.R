# runs ARIMA 
# CEB
#19Feb2025
# sources ARIMAfunction and runs ARIMA using Fable:Forecast 

## FOR FCR: Depths = 1.6 or 8
##          CTD = 1.1 - 2.1 or 7.5 - 8.5
## FOR BVR: Depths = 3 or 9
##          CTD = 2.5 - 3.5 or 8.5 - 9.5
## Use the same precip data for both
###FCR
## 2020 dates: "2020-01-31" "2020-11-06"
## 2021 dates: "2021-02-08" "2021-12-06"
## 2022 dates: "2022-01-31" "2022-12-12"
## 2023 dates: "2023-04-04" "2023-12-05"

###BVR
## 2020 dates: "2020-03-30" "2020-11-16"
## 2021 dates: "2021-03-08" "2021-12-06"
## 2022 dates: "2022-01-26" "2022-12-16"
## 2023 dates: "2023-03-21" "2023-12-05"


library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)
library(forecast)
library(fable)
library(tsibble)

# runs function - should only need to do this once if data is saved 
# source the ARIMA function from folder 
source("https://raw.githubusercontent.com/carlybauer/MSmanuscript/refs/heads/main/ARIMAfunction.R")


FCR2020 <- process_reservoir_data(reservoir = "FCR", year = 2020, 
                                  start_date = "2020-06-12", end_date = "2020-11-06")

FCR2021 <- process_reservoir_data(reservoir = "FCR", year = 2021, 
                                  start_date = "2021-06-12", end_date = "2021-11-06")

FCR2022 <- process_reservoir_data(reservoir = "FCR", year = 2022, 
                                  start_date = "2022-06-12", end_date = "2022-11-06")

FCR2023 <- process_reservoir_data(reservoir = "FCR", year = 2023, 
                                  start_date = "2023-06-12" , end_date = "2023-11-06")
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
# This will add discharge (Q) data for BVR - ignore these values in this dataset 
# since we don't actually collect Q at BVR

BVR2020 <- process_reservoir_data(reservoir = "BVR", year = 2020, 
                                  start_date = "2020-06-12" , end_date = "2020-11-06")

BVR2021 <- process_reservoir_data(reservoir = "BVR", year = 2021, 
                                  start_date = "2021-06-12" , end_date = "2021-11-06")

BVR2022 <- process_reservoir_data(reservoir = "BVR", year = 2022, 
                                  start_date = "2022-06-12" , end_date = "2022-11-06")

BVR2023 <- process_reservoir_data(reservoir = "BVR", year = 2023, 
                                  start_date = "2023-06-12" , end_date = "2023-11-06")
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
# ##### if running function to save merged data, NOT STANDARDIZED
# # standardize data that is processed above 
#   merged_standard <- FCR2023 %>%
#     mutate(across(where(is.numeric) & !c(Date_fake), scale)) %>%
#     as.data.frame(lapply(merged_standard, as.vector)) %>%
#     as_tsibble(index = Date_fake)
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
##### if running function to save standardized data then save as .RData and load 
save(BVR2023, file = "BVR2023_ARIMA_0612_1106.RData")
load("BVR2020_ARIMA_0612_1106.RData")
  
  
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
# Fit ARIMA model with metal as the response and other variables as regressors
fit <- BVR2020 %>%
  model(ARIMA(TBa_mgL ~ DO_mgL + Turbidity_NTU+ Lag1_Rain_Total_mm + SFe_mgL))

# View model summary, including coefficients
fit_report <- report(fit)

