# Statistics and plotting for sed traps 
# Author: CEB (data cleaning from CEW)
# 31 Jan 2025

## script calculates average sed flux for metals 
## calculates cumulative average sed flux for metals
## filters for specific years 
## tests for normality using shapiro test between hypo FCR and BVR sed fluxes
## tests for variance using var.test 
## runs t.test to test significance between FCR and BVR sed flux and cumul sed flux 

## also creates variables for bar plots to look at mean total sed flux for each year 
## with standard deviations calculated and plotted as error bars
## tests normality and uses appropriate two-sample t-test to test if mean are 
## statistically significantly different between FCR and BVR

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)
library(forcats)
library(ggrepel)

#read in from EDI 
#getting the data into the right format
sedtraps <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/1474/8/7cc135d14d92911eb7e8f6ceb096a29d")
sedtraps$Date <- as.Date(sedtraps$Date)

fluxes <- sedtraps %>% 
  select(Reservoir, Date, Depth_m, Duration_days, matches("^T.*g$"),  # Variables that start with 'T' and end with 'g'
         ends_with("_gm2d"), # Variables that end with '_gm2d'
         starts_with("ICPT")   # Variables that end with '_mgL'
  ) %>% 
  filter(Date >= "2020-01-01") %>% 
  group_by(Date, Reservoir, Depth_m) %>% 
  summarise(Duration_days = mean(Duration_days, na.rm = TRUE),
            # AvgSedMass_g = mean(SedMass_g), # not finding this anywhere 
            AvgTLi_g = mean(TLi_g, na.rm = TRUE), # calculates the mean concentration of each metal over the entire dataframe (2019-2023)
            AvgTNa_g = mean(TNa_g, na.rm = TRUE),
            AvgTMg_g = mean(TMg_g, na.rm = TRUE),
            AvgTAl_g = mean(TAl_g, na.rm = TRUE),
            AvgTK_g = mean(TK_g, na.rm = TRUE),
            AvgTCa_g = mean(TCa_g, na.rm = TRUE),
            AvgTFe_g = mean(TFe_g, na.rm = TRUE), #note that we are avging across reps here
            AvgTMn_g = mean(TMn_g, na.rm = TRUE), 
            AvgTCu_g = mean(TCu_g, na.rm = TRUE),
            AvgTSr_g = mean(TSr_g, na.rm = TRUE),
            AvgTBa_g = mean(TBa_g, na.rm = TRUE),
            AvgTLiFlux_gm2d = mean(TLiFlux_gm2d, na.rm = TRUE), #calculates mean flux of each metal over entire dataframe (2020-2023)
            AvgTNaFlux_gm2d = mean(TNaFlux_gm2d, na.rm = TRUE),
            AvgTMgFlux_gm2d = mean(TMgFlux_gm2d, na.rm = TRUE),
            AvgTAlFlux_gm2d = mean(TAlFlux_gm2d, na.rm = TRUE),
            AvgTKFlux_gm2d = mean(TKFlux_gm2d, na.rm = TRUE),
            AvgTCaFlux_gm2d = mean(TCaFlux_gm2d, na.rm = TRUE),
            AvgTFeFlux_gm2d = mean(TFeFlux_gm2d, na.rm = TRUE), 
            AvgTMnFlux_gm2d = mean(TMnFlux_gm2d, na.rm = TRUE),
            AvgTCuFlux_gm2d = mean(TCuFlux_gm2d, na.rm = TRUE),
            AvgTSrFlux_gm2d = mean(TSrFlux_gm2d, na.rm = TRUE),
            AvgTBaFlux_gm2d = mean(TBaFlux_gm2d, na.rm = TRUE),
            AvgICPTFe_mgL = mean(ICPTFe_mgL, na.rm = TRUE)) %>% 
  mutate(Layer = if_else(Depth_m <= 5, 'Epilimnion', 'Hypolimnion'),
         Year = year(Date),
         DOY = yday(Date)) %>% 
  ungroup() %>%
  mutate(across(c('AvgTLi_g', 'AvgTNa_g', 'AvgTMg_g','AvgTAl_g','AvgTK_g','AvgTCa_g','AvgTFe_g', 'AvgTMn_g', 
                  'AvgTCu_g', 'AvgTSr_g', 'AvgTBa_g', 'AvgTLiFlux_gm2d', 'AvgTNaFlux_gm2d', 'AvgTMgFlux_gm2d', 
                  'AvgTAlFlux_gm2d', 'AvgTKFlux_gm2d', 'AvgTCaFlux_gm2d', 'AvgTFeFlux_gm2d', 'AvgTMnFlux_gm2d', 
                  'AvgTCuFlux_gm2d','AvgTSrFlux_gm2d','AvgTBaFlux_gm2d'),
                ~ifelse(is.nan(.), NA, .)), #need to get rid of NaNs created by taking mean during summarise step
         Layer = as.factor(Layer)) %>% 
  arrange(Date)


fluxes2 <- sedtraps %>%
  select(Reservoir, Date, Depth_m, Duration_days, matches("^T.*g$"), ends_with("_gm2d"), starts_with("ICPT")) %>% 
  mutate(Layer = if_else(Depth_m <= 5, 'Epilimnion', 'Hypolimnion'),
         Year = year(Date),
         DOY = yday(Date)) %>% 
  filter(Year >= 2020)
#na.omit()



#now we use group_by to break into reservoir, year, and depth and calculate cumulative mass
fluxes <- fluxes %>% 
  group_by(Reservoir, Year, Layer) %>% 
  mutate(
    CumulLi_g = cumsum(AvgTLi_g), 
    CumulNa_g = cumsum(AvgTNa_g), 
    CumulMg_g = cumsum(AvgTMg_g), 
    CumulAl_g = cumsum(AvgTAl_g), 
    CumulK_g = cumsum(AvgTK_g), 
    CumulCa_g = cumsum(AvgTCa_g), 
    CumulFe_g = cumsum(AvgTFe_g), 
    CumulMn_g = cumsum(AvgTMn_g),
    CumulCu_g = cumsum(AvgTCu_g), 
    CumulSr_g = cumsum(AvgTSr_g), 
    CumulBa_g = cumsum(AvgTBa_g), 
    CumulFeMn = CumulFe_g/CumulMn_g, 
    
    #CumulSed_g = cumsum(AvgSedMass_g), #not sure where this variable comes from 
    
    CumulLiFlux_gm2 = cumsum(AvgTLiFlux_gm2d*Duration_days),
    CumulNaFlux_gm2 = cumsum(AvgTNaFlux_gm2d*Duration_days),
    CumulMgFlux_gm2 = cumsum(AvgTMgFlux_gm2d*Duration_days),
    CumulAlFlux_gm2 = cumsum(AvgTAlFlux_gm2d*Duration_days),
    CumulKFlux_gm2 = cumsum(AvgTKFlux_gm2d*Duration_days),
    CumulCaFlux_gm2 = cumsum(AvgTCaFlux_gm2d*Duration_days),
    CumulFeFlux_gm2 = cumsum(AvgTFeFlux_gm2d*Duration_days),
    CumulMnFlux_gm2 = cumsum(AvgTMnFlux_gm2d*Duration_days),
    CumulCuFlux_gm2 = cumsum(AvgTCuFlux_gm2d*Duration_days),
    CumulSrFlux_gm2 = cumsum(AvgTSrFlux_gm2d*Duration_days),
    CumulBaFlux_gm2 = cumsum(AvgTBaFlux_gm2d*Duration_days)) %>% 
  ungroup()

#let's add turnover dates

fluxesO2 <- fluxes %>% 
  mutate(TO = NA,
         #2020
         TO = ifelse(Year == 2020 & Reservoir == 'FCR' & Date < as.Date('2020-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'BVR' & Date < as.Date('2020-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'FCR' & Date >= as.Date('2020-11-02'), 'Post', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'BVR' & Date >= as.Date('2020-11-02'), 'Post', TO),
         #2021
         TO = ifelse(Year == 2021 & Reservoir == 'FCR' & Date < as.Date('2021-11-03'), 'Pre', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'BVR' & Date < as.Date('2021-11-05'), 'Pre', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'FCR' & Date >= as.Date('2021-11-03'), 'Post', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'BVR' & Date >= as.Date('2021-11-05'), 'Post', TO),
         #2022
         TO = ifelse(Year == 2022 & Reservoir == 'FCR' & Date < as.Date('2022-10-20'), 'Pre', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'BVR' & Date < as.Date('2022-10-19'), 'Pre', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'FCR' & Date >= as.Date('2022-10-20'), 'Post', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'BVR' & Date >= as.Date('2022-10-19'), 'Post', TO),
         #2023
         TO = ifelse(Year == 2023 & Reservoir == 'FCR' & Date < as.Date('2023-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'BVR' & Date < as.Date('2023-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'FCR' & Date >= as.Date('2023-11-02'), 'Post', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'BVR' & Date >= as.Date('2023-11-02'), 'Post', TO)) 

fluxes2 <- fluxes2 %>% 
  mutate(TO = NA,
         #2020
         TO = ifelse(Year == 2020 & Reservoir == 'FCR' & Date < as.Date('2020-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'BVR' & Date < as.Date('2020-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'FCR' & Date >= as.Date('2020-11-02'), 'Post', TO),
         TO = ifelse(Year == 2020 & Reservoir == 'BVR' & Date >= as.Date('2020-11-02'), 'Post', TO),
         #2021
         TO = ifelse(Year == 2021 & Reservoir == 'FCR' & Date < as.Date('2021-11-03'), 'Pre', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'BVR' & Date < as.Date('2021-11-05'), 'Pre', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'FCR' & Date >= as.Date('2021-11-03'), 'Post', TO),
         TO = ifelse(Year == 2021 & Reservoir == 'BVR' & Date >= as.Date('2021-11-05'), 'Post', TO),
         #2022
         TO = ifelse(Year == 2022 & Reservoir == 'FCR' & Date < as.Date('2022-10-20'), 'Pre', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'BVR' & Date < as.Date('2022-10-19'), 'Pre', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'FCR' & Date >= as.Date('2022-10-20'), 'Post', TO),
         TO = ifelse(Year == 2022 & Reservoir == 'BVR' & Date >= as.Date('2022-10-19'), 'Post', TO),
         #2023
         TO = ifelse(Year == 2023 & Reservoir == 'FCR' & Date < as.Date('2023-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'BVR' & Date < as.Date('2023-11-02'), 'Pre', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'FCR' & Date >= as.Date('2023-11-02'), 'Post', TO),
         TO = ifelse(Year == 2023 & Reservoir == 'BVR' & Date >= as.Date('2023-11-02'), 'Post', TO))

#
#
#
#
#

# filter for finding total cumulative flux in each reservoir by year 
# first make dataframe of just point pre turnover for FCR and BVR in the hypo
# removing dates where we have extra points sampled in FCR compared to BVR based on code below that is commented out
FluxHypPre <- fluxesO2 %>% 
  filter(Layer == 'Hypolimnion' & TO == 'Pre' & Date != as.Date("2021-06-07") & Date != as.Date("2022-10-17")) 

# calculate total number of days included in sampling year for each reservoir
# to confirm we are using around the same number of collection days to compare with
total_duration <- FluxHypPre %>%
  group_by(Year, Reservoir) %>%
  summarise(Total_Duration_Days = sum(Duration_days, na.rm = TRUE), .groups = "drop")

print(total_duration)# 2020 had 9 day difference, 2021 had 7 day difference, 2022 and 2023 same days

## FIND MAXES and save into dataframe 
# Define the metals of interest
metals <- c("CumulAlFlux_gm2", "CumulBaFlux_gm2", "CumulCuFlux_gm2", "CumulSrFlux_gm2")

# Compute max values for each year and reservoir
max_fluxes <- FluxHypPre %>%
  group_by(Year, Reservoir) %>%
  summarise(across(all_of(metals), ~ max(.x, na.rm = TRUE), .names = "max_{.col}"), .groups = "drop")

#write_csv(max_fluxes,'maxThesisFluxes_gm2.csv')

## test for normality on above dataset (max metal flux for each year and reservoir)
# Check normality using Shapiro-Wilk test for each metal
metals <- c("max_CumulAlFlux_gm2", "max_CumulBaFlux_gm2", 
            "max_CumulCuFlux_gm2", "max_CumulSrFlux_gm2")

normality_results <- lapply(metals, function(metal) {
  shapiro.test(max_fluxes[[metal]])
})

# Print results
names(normality_results) <- metals
normality_results
## Al not normal pvalue = 0.01861 Wilcoxon signed-rank test
## Ba normal pvalue = 0.6713 paired t test
## Cu normal pvalue = 0.5057 paired t test
## Sr normal pvalue = 0.05052 paired t test

# Run paired t-test for each metal 
t_test_results <- lapply(metals, function(metal) {
  t.test(max_fluxes[[metal]][max_fluxes$Reservoir == "BVR"], 
         max_fluxes[[metal]][max_fluxes$Reservoir == "FCR"], 
         paired = TRUE)
})

# Print results
names(t_test_results) <- metals
t_test_results
## valid only for Ba Cu Sr
## Ba p-value = 0.1001 no significant difference
## Cu p-value = 0.01096 significant difference
## Sr p-value = 0.2277 no significant difference

# Run Wilcoxon signed-rank test for each metal
wilcoxon_results <- lapply(metals, function(metal) {
  wilcox.test(max_fluxes[[metal]][max_fluxes$Reservoir == "BVR"], 
              max_fluxes[[metal]][max_fluxes$Reservoir == "FCR"], 
              paired = TRUE)
})

# Print results
names(wilcoxon_results) <- metals
wilcoxon_results
## valid only for Al
## Al p-value = 0.375 no significant difference

# Calculate mean and standard deviation of max flux for each metal per reservoir across all years
summary_max_fluxes <- max_fluxes %>%
  group_by(Reservoir) %>%
  summarise(
    across(starts_with("max_Cumul"), list(mean = ~ mean(.x, na.rm = TRUE), 
                                          sd = ~ sd(.x, na.rm = TRUE)), .names = "{.fn}_{.col}"),
    .groups = "drop"
  )


#write_csv(summary_max_fluxes,'summaryThesisFluxes_gm2.csv')
#####PLOT 
##::::::::::::::::::::::::::::::::::::::::::::::##
# create dot plot showing max cumulative removal for each year then star where mean is
library(tidyverse)
library(ggpubr)

# Reshape max_fluxes to long format
plot_data <- max_fluxes %>%
  pivot_longer(cols = starts_with("max_Cumul"),
               names_to = "Metal",
               values_to = "Flux") %>%
  mutate(Metal = str_remove(Metal, "max_Cumul"),
         Metal = str_remove(Metal, "Flux_gm2"),
         Metal = str_trim(Metal))

# Calculate mean max flux per Reservoir and Metal
mean_flux <- plot_data %>%
  group_by(Reservoir, Metal) %>%
  summarise(MeanFlux = mean(Flux), .groups = "drop")


library(tidyverse)
library(ggrepel)
library(ggpubr)

# Step 1: Jitter manually by converting Reservoir to numeric
plot_data <- plot_data %>%
  mutate(Reservoir_num = as.numeric(factor(Reservoir, levels = c("BVR", "FCR"))),
         jittered_x = jitter(Reservoir_num, amount = 0.1))

# Update mean_flux with numeric x for means
mean_flux <- mean_flux %>%
  mutate(Reservoir_num = as.numeric(factor(Reservoir, levels = c("BVR", "FCR"))))

# Define color palette
fill_colors <- c("FCR" = "#56B4E9", "BVR" = "#D55E00")

# Function to plot for a given metal
plot_metal <- function(metal_name, y_label, plot_title) {
  metal_data <- filter(plot_data, Metal == metal_name)
  mean_data <- filter(mean_flux, Metal == metal_name)
  
  metal_data$Reservoir <- factor(metal_data$Reservoir, levels = c("BVR", "FCR"))
  mean_data$Reservoir <- factor(mean_data$Reservoir, levels = c("BVR", "FCR"))
  
  ggplot(metal_data, aes(x = jittered_x, y = Flux, color = Reservoir)) +
    geom_point(size = 3, alpha = 1) +
    geom_text_repel(aes(label = Year), size = 3, show.legend = FALSE, max.overlaps = 50) +
    geom_point(data = mean_data,
               aes(x = Reservoir_num, y = MeanFlux),
               shape = 18, color = "black", size = 4, stroke = 1.2, alpha = 0.7) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = c("BVR", "FCR")
    ) +
    coord_cartesian(xlim = c(0.6, 2.4)) +  # Adjust limits to tighten spacing
    scale_color_manual(values = fill_colors) +
    labs(x = NULL, y = y_label, title = plot_title) +
    theme_bw() +
    guides(color = "none")
}

# Generate each plot
Al <- plot_metal("Al", "Cumulative Al Removal (g/m²)", "Aluminum")
Ba <- plot_metal("Ba", "Cumulative Ba Removal (g/m²)", "Barium")
Cu <- plot_metal("Cu", "Cumulative Cu Removal (g/m²)", "Copper")
Sr <- plot_metal("Sr", "Cumulative Sr Removal (g/m²)", "Strontium")

# Combine into one plot
ggarrange(Al, Ba, Cu, Sr,
          ncol = 2, nrow = 2,
          labels = c("A", "B", "C", "D"))



## this doesnt jitter nicely. overlaps year and points 
# # Define color palette
# fill_colors <- c("FCR" = "#56B4E9", "BVR" = "#D55E00")
# 
# # Plot for each metal 
# Al <- ggplot(filter(plot_data, Metal == "Al"), aes(x = Reservoir, y = Flux, color = Reservoir)) +
#   geom_jitter(width = 0.1, size = 3, alpha = 1) +
#   geom_text_repel(aes(label = Year), size = 3, show.legend = FALSE) +
#   geom_point(data = filter(mean_flux, Metal == "Al"), aes(y = MeanFlux),
#              shape = 18, color = "black", size = 4, stroke = 1.2, alpha = 0.7) +
#   scale_color_manual(values = fill_colors) +
#   labs(x = NULL, y = "Cumulative Al Flux (g/m²)", title = "Aluminum") +
#   theme_bw() + 
#   guides(color = "none")
# 
# 
# # Ba
# Ba <- ggplot(filter(plot_data, Metal == "Ba"), aes(x = Reservoir, y = Flux, color = Reservoir)) +
#   geom_jitter(width = 0.1, size = 3, alpha = 1) +
#   geom_text_repel(aes(label = Year), size = 3, show.legend = FALSE) +
#   geom_point(data = filter(mean_flux, Metal == "Ba"), aes(y = MeanFlux),
#              shape = 18, color = "black", size = 4, stroke = 1.2, alpha = 0.7) +
#   scale_color_manual(values = fill_colors) +
#   labs(x = NULL, y = "Cumulative Ba Flux (g/m²)", title = "Barium") +
#   theme_bw() + guides(color = "none")
# 
# # Cu
# Cu <- ggplot(filter(plot_data, Metal == "Cu"), aes(x = Reservoir, y = Flux, color = Reservoir)) +
#   geom_jitter(width = 0.1, size = 3, alpha = 1) +
#   geom_text_repel(aes(label = Year), size = 3, show.legend = FALSE) +
#   geom_point(data = filter(mean_flux, Metal == "Cu"), aes(y = MeanFlux),
#              shape = 18, color = "black", size = 4, stroke = 1.2, alpha = 0.7) +
#   scale_color_manual(values = fill_colors) +
#   labs(x = NULL, y = "Cumulative Cu Flux (g/m²)", title = "Copper") +
#   theme_bw() + guides(color = "none")
# 
# # Sr
# Sr <- ggplot(filter(plot_data, Metal == "Sr"), aes(x = Reservoir, y = Flux, color = Reservoir)) +
#   geom_jitter(width = 0.1, size = 3, alpha = 1) +
#   geom_text_repel(aes(label = Year), size = 3, show.legend = FALSE) +
#   geom_point(data = filter(mean_flux, Metal == "Sr"), aes(y = MeanFlux),
#              shape = 18, color = "black", size = 4, stroke = 1.2, alpha = 0.7) +
#   scale_color_manual(values = fill_colors) +
#   labs(x = NULL, y = "Cumulative Sr Flux (g/m²)", title = "Strontium") +
#   theme_bw() + guides(color = "none")
# 
# # Combine into one plot
# ggarrange(Al, Ba, Cu, Sr,
#           ncol = 2, nrow = 2,
#           labels = c("A", "B", "C", "D"))
