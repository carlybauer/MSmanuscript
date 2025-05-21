# D:T Ratio Heatmap plotting for multiple years and both reservoirs
# CEB
# 13 Feb 2025
# calculates S/T ratio of metals, first checking if T>S then swapping if true 
# for that particular metal (not the whole row)
# sources heatmap ratio function CEW added and interpolates only between 0 and 1
# shouldn't be above 1 because S/T = 1 means all of the sample is soluble 

setwd("/Users/carlybauer/Documents/R/GitHub/Metals")
# Read in Packages
pacman::p_load(tidyverse, gridExtra)
library(ggpubr)

# source the heatmap function from folder because edited 
source("Heatmap_function_multipleyears.R")

# Read in the data you want to use. You will have to change the path of this file because it doesn't exsist. 
#metals <- read_csv("./Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLmetals/Metals_2014_2023.csv")
#metals <- read_csv("metals_2014_2024.csv")

metals <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/455/9/9a072c4e4af39f96f60954fc4f7d8be5")


# filter the data to the time you want. Right now this is just for 2021.

met <- metals |>
  mutate(Year = year(DateTime)) %>% 
  filter(Year >= 2020, Year < 2024,
         Depth_m <= 9.0)

# this gets rid of any values that we flagged as abnormally high in the data
met <- met %>% 
  filter(Flag_TCu_mgL != 8 & Flag_SCu_mgL != 8 & Flag_TCu_mgL != 68 & Flag_SCu_mgL != 68,
         Flag_TSr_mgL != 8 & Flag_SSr_mgL != 8 & Flag_TSr_mgL != 68 & Flag_SSr_mgL != 68,
         Flag_TBa_mgL != 8 & Flag_SBa_mgL != 8 & Flag_TBa_mgL != 68 & Flag_SBa_mgL != 68,
         Flag_TAl_mgL != 8 & Flag_SAl_mgL != 8 & Flag_TAl_mgL != 68 & Flag_SAl_mgL != 68,
         Flag_TNa_mgL != 8 & Flag_SNa_mgL != 8 & Flag_TNa_mgL != 68 & Flag_SNa_mgL != 68,
         Flag_TCa_mgL != 8 & Flag_SCa_mgL != 8 & Flag_TCa_mgL != 68 & Flag_SCa_mgL != 68,
         Flag_TMg_mgL != 8 & Flag_SMg_mgL != 8 & Flag_TMg_mgL != 68 & Flag_SMg_mgL != 68)



### FUNCTION CALCULATES RATIOS

# Create a function to correct values and maintain full dataset - dont change this
correct_metal <- function(df, total_col, soluble_col) {
  df %>%
    mutate(
      # Calculate the difference
      Diff = !!sym(total_col) - !!sym(soluble_col),
      
      # Swap the values only if the total is less than the soluble value
      new_total = ifelse(Diff < 0, !!sym(soluble_col), !!sym(total_col)),
      new_soluble = ifelse(Diff < 0, !!sym(total_col), !!sym(soluble_col))
    ) %>%
    # Replace the old columns with the new swapped values
    mutate(
      !!sym(total_col) := new_total,
      !!sym(soluble_col) := new_soluble
    ) %>%
    select(-Diff, -new_total, -new_soluble)  # Remove intermediate columns
}

# Apply function to each metal while preserving the full dataset 
# this is where you change to input your metals
met_corrected <- met %>%
  correct_metal("TAl_mgL", "SAl_mgL") %>%
  correct_metal("TBa_mgL", "SBa_mgL") %>%
  correct_metal("TCu_mgL", "SCu_mgL") %>%
  correct_metal("TSr_mgL", "SSr_mgL") %>%
  correct_metal("TNa_mgL", "SNa_mgL") %>%
  correct_metal("TCa_mgL", "SCa_mgL") %>%
  correct_metal("TMg_mgL", "SMg_mgL") %>%
  mutate(
    STAl = SAl_mgL / TAl_mgL,
    STBa = SBa_mgL / TBa_mgL,
    STCu = SCu_mgL / TCu_mgL,
    STSr = SSr_mgL / TSr_mgL,
    STNa = SNa_mgL / TNa_mgL,
    STCa = SCa_mgL / TCa_mgL,
    STMg = SMg_mgL / TMg_mgL
  )

###END OF FUNCTION AND DATA CHANGE TO CALC RATIOS



# ### BELOW CALCULATES S/T RATIOS the long way - don't use this 
# #create diff column
# metdiff <- met %>% 
#   group_by(Reservoir, Site, DateTime, Depth_m) %>% 
#   mutate(AlDiff = TAl_mgL - SAl_mgL,
#          BaDiff = TBa_mgL - SBa_mgL,
#          CuDiff = TCu_mgL - SCu_mgL, 
#          SrDiff = TSr_mgL - SSr_mgL)
# 
# 
# #create a data frame with the correct values
# met_correct <- metdiff %>%
#   filter(AlDiff > 0,
#          BaDiff > 0, 
#          CuDiff > 0, 
#          SrDiff > 0)
# 
# #create a data frame with the incorrect values and switch them
# Alissues <- metdiff %>% 
#   select(Reservoir, Site, DateTime, Depth_m, TAl_mgL, SAl_mgL, AlDiff) %>% 
#   filter(AlDiff < 0) %>% 
#   rename('test1' = 'TAl_mgL',
#          'test2' = 'SAl_mgL') %>% 
#   rename('TAl_mgL' = 'test2',
#          'SAl_mgL' = 'test1') %>% 
#   mutate(AlDiff = TAl_mgL - SAl_mgL)
# 
# Baissues <- metdiff %>% 
#   select(Reservoir, Site, DateTime, Depth_m, TBa_mgL, SBa_mgL, BaDiff) %>% 
#   filter(BaDiff < 0) %>% 
#   rename('test1' = 'TBa_mgL',
#          'test2' = 'SBa_mgL') %>% 
#   rename('TBa_mgL' = 'test2',
#          'SBa_mgL' = 'test1') %>% 
#   mutate(BaDiff = TBa_mgL - SBa_mgL)
# 
# Cuissues <- metdiff %>% 
#   select(Reservoir, Site, DateTime, Depth_m, TCu_mgL, SCu_mgL, CuDiff) %>% 
#   filter(CuDiff < 0) %>% 
#   rename('test1' = 'TCu_mgL',
#          'test2' = 'SCu_mgL') %>% 
#   rename('TCu_mgL' = 'test2',
#          'SCu_mgL' = 'test1') %>% 
#   mutate(CuDiff = TCu_mgL - SCu_mgL)
# 
# Srissues <- metdiff %>% 
#   select(Reservoir, Site, DateTime, Depth_m, TSr_mgL, SSr_mgL, SrDiff) %>% 
#   filter(SrDiff < 0) %>% 
#   rename('test1' = 'TSr_mgL',
#          'test2' = 'SSr_mgL') %>% 
#   rename('TSr_mgL' = 'test2',
#          'SSr_mgL' = 'test1') %>% 
#   mutate(SrDiff = TSr_mgL - SSr_mgL)
# 
# #rbind to rejoin and calculate S/T ratio
# metals_correct <- rbind(met_correct, Alissues, Baissues, Cuissues, Srissues) %>% 
#   select(Reservoir, Site, DateTime, Depth_m, SAl_mgL, TAl_mgL, SBa_mgL, TBa_mgL, SCu_mgL, TCu_mgL, TSr_mgL, SSr_mgL) %>%
#   mutate(STAl = SAl_mgL/TAl_mgL,
#          STBa = SBa_mgL/TBa_mgL,
#          STCu = SCu_mgL/TCu_mgL,
#          STSr = SSr_mgL/TSr_mgL)
# 
# ### END OF CALCULATING RATIOS by long way

# sources heatmap_ratio function and makes basic plots 
#plots ratio
heatAl <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STAl")
heatBa <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STBa")
heatCu <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STCu")
heatSr <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STSr")
heatNa <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STNa")
heatCa <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STCa")
heatMg <- heatmap_ratio(data=met_corrected, reservoirs =c("FCR", "BVR"), site=50, z="STMg")


# look at the basic plot
print(heatAl)

print(heatBa)

print(heatCu)

print(heatSr)

print(heatNa)

print(heatCa)

print(heatMg)

# ### HOx on dates 
# # load csv to add lines to fancy plot below
# #create new frames to add TO lines
# HOxBoxes <- read_csv('HOxOnDates.csv') %>% 
#   mutate(Group = rownames(.),
#          Reservoir = 'FCR',
#          Year = year(HOxOn),
#          DOYon = yday(HOxOn),
#          DOYoff = yday(HOxOff))

### Turnover dates 
# load csv to add lines to fancy plot below
#create new frames to add TO lines
TO <- read_csv('TO_Dates.csv')

# ratio plots 
Ba <- heatBa + 
  # add a title
  ggtitle("Ba D:T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Ba)

Al <- heatAl + 
  # add a title
  ggtitle("Al D:T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Al)

Cu <- heatCu + 
  # add a title
  ggtitle("Cu D:T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Cu)

Sr <- heatSr + 
  # add a title
  ggtitle("Sr D:T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Sr)

Na <- heatNa + 
  # add a title
  ggtitle("Na S/T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = subset(HOxBoxes, Reservoir == "FCR"), aes(xintercept = as.numeric(HOxOn)), 
             linetype = "solid", color = 'black', linewidth = 0.5) +  
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Na)

Ca <- heatCa + 
  # add a title
  ggtitle("Ca S/T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = subset(HOxBoxes, Reservoir == "FCR"), aes(xintercept = as.numeric(HOxOn)), 
             linetype = "solid", color = 'black', linewidth = 0.5) +  
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Ca)

Mg <- heatMg + 
  # add a title
  ggtitle("Mg S/T ratio") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  # scale_fill_gradientn(colours = blue2green2red(60),
  #                      limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
  #                      na.value = "gray")+
  geom_vline(data = subset(HOxBoxes, Reservoir == "FCR"), aes(xintercept = as.numeric(HOxOn)), 
             linetype = "solid", color = 'black', linewidth = 0.5) +  
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text=element_text(size=12), #change font size of axis text
        axis.title=element_text(size=14), #change font size of axis titles
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=12)) #change font size of legend title   
print(Mg)
