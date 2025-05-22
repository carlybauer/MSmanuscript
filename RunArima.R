# runs ARIMA 
# CEB
#19Feb2025
# sources ARIMAfunction and runs ARIMA using Fable
# interpolates using Forecast

## FOR FCR: Depths = 1.6 or 8
##          CTD = 1.1 - 2.1 or 7.5 - 8.5
## FOR BVR: Depths = 3 or 9
##          CTD = 2.5 - 3.5 or 8.5 - 9.5
## Use the same precip data for both
## Do not use Q data for BVR


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
#source("https://raw.githubusercontent.com/carlybauer/MSmanuscript/refs/heads/main/ARIMAfunction.R")
source("ARIMAfunction.R")

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
save(BVR2023, file = "BVR2023_ARIMA_0612_1106_all.RData")
load("BVR2023_ARIMA_0612_1106_all.RData")
  
  
##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
# variables to include in model run
# FCR 2020: Tmetal_mgL ~ DO_mgL + Turbidity_NTU + Lag1_Rain_Total_mm + Lag1_Q_Sum_cms
# FCR 2021: Tmetal_mgL ~ DO_mgL + pH + Lag1_Rain_Total_mm + Lag1_Q_Sum_cms
# FCR 2022: Tmetal_mgL ~ DO_mgL + pH + Turbidity_NTU + Lag1_Rain_Total_mm + Lag1_Q_Sum_cms
# FCR 2023: Tmetal_mgL ~ DO_mgL + pH + Turbidity_NTU + Lag1_Rain_Total_mm + Lag1_Q_Sum_cms

# BVR 2020: Tmetal_mgL ~ DO_mgL + Turbidity_NTU + Lag1_Rain_Total_mm
# BVR 2021: nothing because too many interpolated metal and env variables
# BVR 2022: nothing because too many interpolated metal and env variables
# BVR 2023: Tmetal_mgL ~ DO_mgL + pH + Turbidity_NTU + Lag1_Rain_Total_mm
# Fit ARIMA model with metal as the response and other variables as regressors
fit <- BVR2023 %>%
  model(ARIMA(TSr_mgL ~ DO_mgL + pH + Turbidity_NTU + Lag1_Rain_Total_mm))

# View model summary, including coefficients
fit_report <- report(fit)
#tidy(fit)
#glance(fit)
#accuracy(fit)
##::::::::::::::::::::::::::::::::::::::::::::##
# if ARIMA (0,0,0) then can fit to lm 
lm_fit <- lm(TBa_mgL ~ DO_mgL  + Turbidity_NTU + Lag1_Rain_Total_mm, data = BVR2020)
summary(lm_fit)
# FCR 2020: 
  ##TAl p value = 7.509e-06 R2 = 0.8041
  ##TSr p-value = 8.037e-07 R2 = 0.8503
# FCR 2021: 
  ##TAl p-value = 0.003761 R2 = 0.5794
  ##TBa p-value = 5.678e-06 R2 = 0.8106
  ##TCu p-value = 5.067e-06 R2 = 0.8132
# FCR 2022: 
  ##TAl p-value = 6.388e-05 R2 = 0.7877
  ##TCu p-value = 0.0007022 R2 = 0.7086
  ##TSr p-value = 6.375e-06 R2 = 0.8425
# FCR 2023: 
  ##TBa p-value = 2.319e-10 R2 = 0.957
  ##TCu p-value = 8.042e-05 R2 = 0.7812
  ##TSr p-value = 3.781e-06 R2 = 0.8528
# BVR 2023: 
  ##TAl p-value = 0.4751  R2 = 0.1779
  ##TBa p-value = 0.000321  R2 = 0.6906
  ##TCu p-value = 0.5909  R2 = 0.1446
##::::::::::::::::::::::::::::::::::::::::::::##
#if ARIMA (2,0,0) or (1,0,0) then
# FCR 2020: 
  ##TBa: p-value =   R2 = 0.8405348 (1,0,0)
  ##TCu: p-value =   R2 = 0.8757426 (2,0,0)
# FCR 2021:
  ##TSr: p-value =   R2 = 0.8123439 (1,0,0)
# FCR 2022: 
  ##TBa: p-value =   R2 = 0.7894672 (2,0,0)
# FCR 2023: 
  ##TAl: p-value =   R2 = 0.4828279 (1,0,0)
# BVR 2020: 
  ##TCu: p-value =   R2 = 0.6143354 (2,0,0)
# BVR 2023: 
  ##TSr: p-value =   R2 = 0.7715246 (1,0,0)

# Add fitted values and compute R²
augmented <- augment(fit)

# Compute pseudo-R² manually
pseudo_r2 <- cor(augmented$.fitted, augmented$TSr_mgL)^2
print(pseudo_r2)
##::::::::::::::::::::::::::::::::::::::::::::##
# if ARIMA (1,1,0)
# BVR 2020: TAl, 

##::::::::::::::::::::::::::::::::::::::::::::##
# if ARIMA (0,1,0) 
#BVR: 2020: TBa, TSr
  ## TBa: p-value = 5.941e-05  R2 = 0.6999
  ## TSr: p-value = 0.0001833  R2 = 0.6589
# Create differenced versions of all variables
df_diff <- data.frame(
  dmetal = diff(BVR2020$TSr_mgL),
  dDO = diff(BVR2020$DO_mgL),
  dTurbidity = diff(BVR2020$Turbidity_NTU),
  dRain = diff(BVR2020$Lag1_Rain_Total_mm)
)

# Fit the linear model on differenced data (no intercept)
lm_diff <- lm(dmetal ~ dDO + dTurbidity + dRain - 1, data = df_diff)
summary(lm_diff)

