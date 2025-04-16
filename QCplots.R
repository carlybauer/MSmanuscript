## Q-C plot for discharge data 
## 16Aug2024
## Carly Bauer

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)
library(ggpubr)

# load packages
pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, cowplot, reshape2, scales
)

# read in metals data from EDI
met <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/455/9/9a072c4e4af39f96f60954fc4f7d8be5")

#getting the metals into the right format
metals <- met %>%
  mutate(Year = year(DateTime),
         Week = week(DateTime),
         Date = date(DateTime)) %>% 
  filter(Reservoir == 'FCR',
         Site == 50,
         Year >= 2020 & Year < 2024) %>% # can change to look at different years 
  mutate(Date=date(DateTime))

# removes random high TAl >0.6
metals$TAl_mgL[metals$Depth_m == 1.6 & metals$Date == ("2020-04-13")] <- NA



metals <- metals %>% 
  # this gets rid of any values that we flagged as abnormally high in the data
  filter(Flag_TCu_mgL != 8 & Flag_SCu_mgL != 8 & Flag_TCu_mgL != 68 & Flag_SCu_mgL != 68,
         Flag_TSr_mgL != 8 & Flag_SSr_mgL != 8 & Flag_TSr_mgL != 68 & Flag_SSr_mgL != 68,
         Flag_TBa_mgL != 8 & Flag_SBa_mgL != 8 & Flag_TBa_mgL != 68 & Flag_SBa_mgL != 68,
         Flag_TAl_mgL != 8 & Flag_SAl_mgL != 8 & Flag_TAl_mgL != 68 & Flag_SAl_mgL != 68)



#read in discharge data from EDI
dis <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/202/12/aae7888d68753b276d1623680f81d5de")

# getting discharge into right format
discharge <- dis %>% 
  mutate(DateTime = ymd_hms(DateTime),# separate DateTime column 
         Year = year(DateTime),
         Date = date(DateTime),
         Hour = hour(DateTime),
         Week = week(DateTime)) %>% 
  filter(Year >= 2020 & Year < 2024) %>% # change for different year
  group_by(Year, Week) %>% # groups by date then week and calculates the average inflow for the day or week, don't sum because its a rate
  summarize(meanQ = mean(coalesce(VT_Flow_cms, WVWA_Flow_cms), na.rm = TRUE)) %>% 
  ungroup()

#creates lags for discharge coming in at 100 to be seen in site 50
discharge <- discharge %>% 
 mutate(Lag1_meanQ = lag(meanQ, n = 1), # lags discharge by 1 week 
        Lag2_meanQ = lag(meanQ, n = 2)) # lags discharge by 2 weeks

# dataframe joining metal concentrations with mean discharge data by the Date
qc <- left_join(metals, discharge, by=c('Year','Week')) %>% # or Date
  mutate(Month = month(Date)) %>% 
  select(-(starts_with("Flag")))


# Q-C Plots without regression line 
ggplot(qc) + 
  theme_bw()+
  geom_point(aes(x = Lag1_meanQ, y = TAl_mgL, color = Month, shape = as.factor(Depth_m))) + # Change this based on metal you want to see
  scale_color_gradientn(
    colours = c("blue","green", "red","green", "blue"), 
    values = scales::rescale(c(1, 3, 6, 9, 12)),  # 1 = January (blue), 6 = June (green), 12 = December (blue)
    na.value = "gray") +
  scale_shape_manual(values = 1:10, name = "Depth_m") +  # Adjust as needed for different depths
  ggtitle("FCR TAl QC 2020 and 2023")+
  labs(x ="meanQ_lag1week_cms")


### Q-C with regression 
# Function to plot metal vs Lag1_meanQ with regression line and equation
plot_metal_lm <- function(df, metal_col, title_text) {
  # Fit linear model
  formula <- as.formula(paste(metal_col, "~ Lag1_meanQ"))
  lm_model <- lm(formula, data = df)
  
  # Create equation text
  eq <- substitute(italic(y) == a + b * italic(x)~~","~~italic(R)^2~"="~r2, 
                   list(a = format(as.numeric(coef(lm_model)[1]), digits = 3),
                        b = format(as.numeric(coef(lm_model)[2]), digits = 3),
                        r2 = format(summary(lm_model)$r.squared, digits = 3)))
  
  # Build ggplot
  p <- ggplot(df, aes_string(x = "Lag1_meanQ", y = metal_col)) + 
    theme_bw() +
    geom_point(aes(color = Month, shape = as.factor(Depth_m))) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    annotate("text", 
             x = max(df$Lag1_meanQ, na.rm = TRUE), 
             y = max(df[[metal_col]], na.rm = TRUE), 
             label = as.character(as.expression(eq)), 
             parse = TRUE, hjust = 1, vjust = 1) +
    scale_color_gradientn(
      colours = c("blue", "green", "red", "green", "blue"), 
      values = scales::rescale(c(1, 3, 6, 9, 12)),
      na.value = "gray"
    ) +
    scale_shape_manual(values = 1:10, name = "Depth_m") +
    ggtitle(title_text) +
    labs(x = "Mean Discharge (cms) 1-week Lag",
         y = paste0(gsub("_mgL", "", metal_col), " (mg/L)")) +
    theme(plot.title = element_text(size=10, hjust = 0),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10))
  
  return(p)
}

# make plots and arrange together
Al <- plot_metal_lm(qc, "TAl_mgL", "FCR TAl Q-C 2020–2023")
Ba <- plot_metal_lm(qc, "TBa_mgL", "FCR TBa Q-C 2020–2023")
Cu <- plot_metal_lm(qc, "TCu_mgL", "FCR TCu Q-C 2020–2023")
Sr <- plot_metal_lm(qc, "TSr_mgL", "FCR TSr Q-C 2020–2023")

ggarrange(Al, Ba, Cu, Sr, ncol = 2, nrow = 2, common.legend = TRUE, 
          legend = "right", labels = c("A", "B", "C", "D"))


