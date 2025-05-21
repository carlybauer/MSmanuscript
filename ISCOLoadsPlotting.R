# ISCO Loads Plotting
# 26June2024
# Carly Bauer
# uses csv from original workflow
# uses csv from stage ISCO data package 

rm(list=ls(all=TRUE))

library(tidyverse)
library(tidyr)
library(ggplot2)
library(paletteer)
library(lubridate)
library(dplyr)
library(readxl)
library(broom)

setwd("/Users/carlybauer/Documents/R/GitHub/Metals")
# # original dataset 
imetals <- read_csv("ISCO_metal_loads.csv")

imetals <- imetals %>% 
  mutate(End_time = as_date(End_time)) %>% 
  mutate(Start_time = as_date(Start_time))

imetals <- imetals %>%
  filter(variable != "Depth_m") %>% 
  mutate(Year = year(Start_time)) %>% 
  filter(!is.na(Start_time),
         Year >= 2020 & Year < 2024)

metals4 <- imetals %>% 
  filter(variable == "TBa_mgL" | variable == "TAl_mgL" | 
           variable == "TSr_mgL" | variable == "TCu_mgL")
##::::::::::::::::::::::::::::::::::::::::::::::::::##
# thesis plot 
ggplot(metals4) +
  ggtitle("Total Metal Loads")+
  theme_bw() + theme(panel.grid = element_blank()) + 
  facet_wrap(vars(Year), nrow = 2, scales = 'free_x' ) +
  scale_x_date(date_labels = '%b', breaks = '1 month') +
  geom_point(aes(x = Start_time, y = ISCO_load_r_kg, colour = variable, shape = "ISCO_load_r_kg"), size = 2) +
  scale_shape_manual(values = c("ISCO_load_r_kg" = 16))+
  scale_color_paletteer_d("ltc::fernande") +
  xlab(NULL)+
  ylab('Metal Load (kg)') +
  theme(plot.title = element_text(size=10),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  guides(shape = "none") # Adjust size for legend shapes













# # change this link after data is published
# imetals <- 	read_csv("https://pasta-s.lternet.edu/package/data/eml/edi/1308/5/a3faa5a1c31d6eb10ac698a8fdfb9dd2")
# 
# 
# imetals <- imetals %>%
#   mutate(Year = year(Collection_start_time)) %>% 
#   filter(!is.na(Collection_start_time),
#     Year >= 2020 & Year < 2024)
# 
# imetals <- imetals %>% 
#   mutate(Collection_end_time = as_date(Collection_end_time)) %>% 
#   mutate(Collection_start_time = as_date(Collection_start_time))
# 
# metals4 <- imetals %>% 
#   filter(Element_name == "TBa_mgL" | Element_name == "TAl_mgL" | 
#            Element_name == "TSr_mgL" | Element_name == "TCu_mgL")
# 
# 
# # thesis plot 
# ggplot(metals4) +
#   ggtitle("Total Metal Loads")+
#   theme_bw() + theme(panel.grid = element_blank()) + 
#   facet_wrap(vars(Year), nrow = 2, scales = 'free_x' ) +
#   scale_x_date(date_labels = '%b', breaks = '1 month') +
#   geom_point(aes(x = Collection_start_time, y = Load_kg, colour = Element_name), size = 2) +
#   scale_shape_manual(values = c("Load_kg" = 16))+
#   scale_color_paletteer_d("ltc::fernande") +
#   xlab(NULL)+
#   ylab('Metal Load (kg)') +
#   theme(plot.title = element_text(size=10),
#         axis.text=element_text(size=10), #change font size of axis text
#         axis.title=element_text(size=10), #change font size of axis titles
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10)) +
#   guides(shape = "none") # Adjust size for legend shapes




