# Heatmap plotting for multiple years and both reservoirs
# CEB
# 17 Sept 2024
# heatmap total and soluble metals based on year and reservoir 

setwd("/Users/carlybauer/Documents/R/GitHub/Metals")
# Read in Packages
pacman::p_load(tidyverse, gridExtra)
library(ggpubr)


# source the heatmap function from folder because edited 
# source("Heatmap_function_multipleyears.R")
source("https://raw.githubusercontent.com/carlybauer/MSmanuscript/refs/heads/main/Heatmap_function_multipleyears.R")

# Read in metals data from EDI package 
metals <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/455/9/9a072c4e4af39f96f60954fc4f7d8be5")

# filter the data to the time you want. Right now this is just for 2020-2023.

met <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Year >= 2020, Year < 2024,
         Depth_m <= 9.0)

# this gets rid of any values that we flagged as abnormally high in the data
met <- met %>% 
  filter(Flag_TCu_mgL != 8 & Flag_SCu_mgL != 8 & Flag_TCu_mgL != 68 & Flag_SCu_mgL != 68,
         Flag_TSr_mgL != 8 & Flag_SSr_mgL != 8 & Flag_TSr_mgL != 68 & Flag_SSr_mgL != 68,
         Flag_TBa_mgL != 8 & Flag_SBa_mgL != 8 & Flag_TBa_mgL != 68 & Flag_SBa_mgL != 68,
         Flag_TAl_mgL != 8 & Flag_SAl_mgL != 8 & Flag_TAl_mgL != 68 & Flag_SAl_mgL != 68)


# find ranges of min and max to make color bar consistant
# so we can merge the T and S into one plot and confirm gradients are same 
# code is with each metal below



### Make the baseheat map using the heatmap_EDI function. 
# You will have to run this for the different sites and elements.
# data = the name of your data frame with filtered data
# reservoir = the name of the reservoir which is used for filtering
# site = the number of the site you want which is used for filtering
# z = the column you want for your heatmap
# It will give you a warning because the points are all very similar. That's ok!

#plots normal T and S separate
heatTBa <- heatmap_EDI(data=met, reservoirs =c("FCR", "BVR"), site=50, z="TBa_mgL")

heatSBa <- heatmap_EDI(data=met, reservoirs = c("FCR", "BVR"), site = 50, z = "SBa_mgL")

heatTSr <- heatmap_EDI(data=met, reservoirs =c("FCR", "BVR"), site=50, z="TSr_mgL")

heatSSr <- heatmap_EDI(data=met, reservoirs = c("FCR", "BVR"), site = 50, z = "SSr_mgL")

heatTCu <- heatmap_EDI(data=met, reservoirs =c("FCR", "BVR"), site=50, z="TCu_mgL")

heatSCu <- heatmap_EDI(data=met, reservoirs = c("FCR", "BVR"), site = 50, z = "SCu_mgL")

heatTAl <- heatmap_EDI(data=met, reservoirs =c("FCR", "BVR"), site=50, z="TAl_mgL")

heatSAl <- heatmap_EDI(data=met, reservoirs = c("FCR", "BVR"), site = 50, z = "SAl_mgL")

# look at the basic plot
print(heatTAl)

print(heatBa)

print(heatCu)

print(heatSr)

### HOx on dates 
# load csv to add lines to fancy plot below
#create new frames to add TO lines
HOxBoxes <- read_csv('HOxOnDates.csv') %>% 
  mutate(Group = rownames(.),
         Reservoir = 'FCR',
         Year = year(HOxOn),
         DOYon = yday(HOxOn),
         DOYoff = yday(HOxOff))

### Turnover dates 
# load csv to add lines to fancy plot below
#create new frames to add TO lines
TO <- read_csv('TO_Dates.csv')



### Extra things you may want to change on the plot. 

TBa <- heatTBa + 
  # add a title
  ggtitle("Total Ba") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title   
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(TBa)

SBa <- heatSBa + 
  # add a title
  ggtitle("Dissolved Ba") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0030, 0.112300000), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  # adjust the size and the position of the title
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title 
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(SBa)

#library(ggpubr)
ggarrange(TBa, SBa, ncol = 1, common.legend = TRUE, 
          legend = "right")

## calculations for max and min and setting color gradient scale
# range(met$TBa_mgL, na.rm = T) #0.0091 0.1123
# range(met$SBa_mgL, na.rm = T) #0.0051 0.1098 
# range(heatTBa$data$z, na.rm = TRUE) #max = 0.112300000
# range(heatSBa$data$z, na.rm = TRUE) #min = 0.0030
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
TSr <- heatTSr + 
  # add a title
  ggtitle("Total Sr") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.005287647, 0.0267615), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10)) +#change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(TSr)

SSr <- heatSSr + 
  # add a title
  ggtitle("Dissolved Sr") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.005287647, 0.0267615), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  # adjust the size and the position of the title
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(SSr)

ggarrange(TSr, SSr, ncol = 1, common.legend = TRUE, 
          legend = "right")

## calculations for max and min and setting color gradient scale
# range(met$TSr_mgL, na.rm = T) #0.0057 0.0518
# range(met$SSr_mgL, na.rm = T) #0.0055 0.0512
# range(heatTSr$data$z, na.rm = TRUE) #max = 0.0267615
# range(heatSSr$data$z, na.rm = TRUE) #min = 0.005287647
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##
TCu <- heatTCu + 
  # add a title
  ggtitle("Total Cu") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0008807913, 0.0043056011), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10)) +#change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(TCu)

SCu <- heatSCu + 
  # add a title
  ggtitle("Dissolved Cu") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0008807913, 0.0043056011), # makes same gradient with min and max T and S
                       na.value = "gray")+
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  # adjust the size and the position of the title
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(SCu)

ggarrange(TCu, SCu, ncol = 1, common.legend = TRUE, 
          legend = "right")

## calculations for max and min and setting color gradient scale
# range(met$TCu_mgL, na.rm = T) #0.0010 0.0042
# range(met$SCu_mgL, na.rm = T) #0.0010 0.0043
# range(heatTCu$data$z, na.rm = TRUE) #max = 0.0043056011
# range(heatSCu$data$z, na.rm = TRUE) #min = 0.0008807913
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::##

TAl <- heatTAl + 
  # add a title
  ggtitle("Total Al") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0010, 0.6884), # makes same gradient with min and max T and S
                       na.value = "gray", 
                       values=rescale(c(0.0010,0.1,0.6884))) +
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+  
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(TAl)

SAl <- heatSAl + 
  # add a title
  ggtitle("Dissolved Al") +  
  theme(panel.grid = element_blank()) +
  # change the x axis labels depending on how many breaks you want and how you want it labeled
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  scale_fill_gradientn(colours = blue2green2red(60),
                       limits = c(0.0010, 0.6884), # makes same gradient with min and max T and S
                       na.value = "gray", 
                       values=rescale(c(0.0010,0.1,0.6884))) +
  geom_vline(data = TO, aes(xintercept = as.numeric(Date)), 
             linetype = "solid", color = 'black', linewidth = 0.5)+ 
  # adjust the size and the position of the title
  theme(plot.title = element_text(size=10, hjust = 0),
        axis.text=element_text(size=10), #change font size of axis text
        axis.title=element_text(size=10), #change font size of axis titles
        legend.text=element_text(size=10), #change font size of legend text
        legend.title=element_text(size=10))+ #change font size of legend title
  theme(strip.text = element_text(size = 8, margin = margin(0.8,0.8,0.8,0.8)),  # Reduce text size and margins
        strip.background = element_rect(size = 0.2))  # Reduce box thickness
print(SAl)

ggarrange(TAl, SAl, ncol = 1, common.legend = TRUE, 
          legend = "right")

## calculations for max and min and setting color gradient scale
# range(met$TAl_mgL, na.rm = T) #0.0012 0.6884
# range(met$SAl_mgL, na.rm = T) #0.0010 0.2066
# range(heatTAl$data$z, na.rm = TRUE) #max = 0.6884
# range(heatSAl$data$z, na.rm = TRUE) #min = 0.0010



# ### save the heatmap
# # change the name of the file you are saving
# # decide on where the file will be saved, file type, and the size of the image
# ggsave("heatCu2.png", width = 19, height = 13, units = "cm")
# 
# ### Putting multiple plots together into one file
# # You can add extra things here too. See https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
# #heat_all <- gridExtra::grid.arrange(heatTFe, heatTFeB, nrow = 2)
# 
# ### save plot
# # change the name of the file you are saving
# # decide on where the file will be saved, file type, and the size of the image
# ggsave("All_FCR.pdf", width = 20, height = 20, units = "cm")
# 
