# Title: heatmap_EDI function for multiple years and reservoirs
# Authors: Mary Lofton and Adrienne Breef-Pilz
# Created: 04 April 2024
# Edited: CEB 17Sept2024

# Heatmap function for inspection plots. 
# Took Mary's code from Fluoroprobe
# creates facet grid of reservoir and year. multiple inputs allowed for reservoirs
# function used for metal heatmaps and ratio metal heatmaps 

heatmap_EDI <- function(data, reservoirs, site, z){
  
  # load packages
  pacman::p_load(tidyverse, lubridate, akima, reshape2, 
                 gridExtra, grid, colorRamps, RColorBrewer, cowplot, reshape2, scales
                 )
  
  # Initialize final dataframe for both reservoirs
  df_interp_final <- data.frame()
  
  # Loop through each reservoir
  for (reservoir in reservoirs) {
    
    # Subset to relevant data for the current reservoir
    fp <- data %>%
      filter(Reservoir == reservoir & Site == site) %>%
      select(Reservoir, DateTime, Depth_m, {{z}}) %>%
      mutate(Date = as.Date(DateTime),
             Year = year(DateTime))
    
    # Select depth range based on the reservoir
    if (reservoir == "FCR"){
      depths = seq(0.1, 9.3, by = 0.3)
    } else if (reservoir == "BVR"){
      depths = seq(0.1, 10, by = 0.3)
    }
    
    # Initialize dataframe for interpolated data
    df_interp <- data.frame()
    
    # Loop through each depth and slice data
    for (i in 1:length(depths)){
      
      fp_layer <- fp %>%
        group_by(Date) %>%
        slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df_interp <- bind_rows(df_interp, fp_layer) %>%
        dplyr::distinct()
    }
    
    # Wrangle final dataframe for interpolation
    fp_new <- arrange(df_interp, DateTime) %>%
      drop_na({{z}})
    
    # Round each extracted depth to the nearest 10th
    fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5)
    
    # Interpolate data for the current reservoir
    interp <- interp(x = as.numeric(fp_new$Date), y = fp_new$Depth_m, z = unlist(fp_new[z]),
                     xo = seq(min(fp_new$Date), max(fp_new$Date), by = "day"),
                     yo = seq(min(fp_new$Depth_m), max(fp_new$Depth_m), by = 0.1),
                     extrap = TRUE, linear = TRUE, duplicate = "strip")
    
    interp <- interp2xyz(interp, data.frame = TRUE)
    
    # Add 'Reservoir' and 'Year' to the interpolated data
    interp <- interp %>%
      mutate(Date = as.Date(x, origin = "1970-01-01"),
             Year = year(Date),
             Reservoir = reservoir)  # Keep track of which reservoir
    
    # Bind the interpolated data for each reservoir
    df_interp_final <- bind_rows(df_interp_final, interp)
  }
  
  # Create title for plot
  fig_title <- paste(paste(reservoirs, collapse = " & "), "Site", site, z, sep = " ")

  # HOxBoxes <- read_csv('HOxOnDates.csv') %>% 
  #   mutate(Group = rownames(.),
  #          Reservoir = 'FCR',
  #          Year = year(HOxOn),
  #          DOYon = yday(HOxOn),
  #          DOYoff = yday(HOxOff))
  
  # Plot the heatmap with facet_grid by Reservoir and Year
  p1 <- ggplot(df_interp_final, aes(x = Date, y = y)) +
    geom_raster(aes(fill = z)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_date(expand = c(0, 0)) +
    scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray")+
    
    facet_grid(rows = vars(Reservoir), cols = vars(Year), scales = 'free_x', space = 'free_x') +
    labs(x = NULL, y = "Depth (m)", title = fig_title, fill = ("mg/L")) +
    theme_bw()
  
  return(p1)
}


heatmap_ratio <- function(data, reservoirs, site, z){
  
  # load packages
  pacman::p_load(tidyverse, lubridate, akima, reshape2, 
                 gridExtra, grid, colorRamps, RColorBrewer, cowplot, reshape2, scales
  )
  
  # Initialize final dataframe for both reservoirs
  df_interp_final <- data.frame()
  
  # Loop through each reservoir
  for (reservoir in reservoirs) {
    
    # Subset to relevant data for the current reservoir
    fp <- data %>%
      filter(Reservoir == reservoir & Site == site) %>%
      select(Reservoir, DateTime, Depth_m, {{z}}) %>%
      mutate(Date = as.Date(DateTime),
             Year = year(DateTime))
    
    # Select depth range based on the reservoir
    if (reservoir == "FCR"){
      depths = seq(0.1, 9.3, by = 0.3)
    } else if (reservoir == "BVR"){
      depths = seq(0.1, 10, by = 0.3)
    }
    
    # Initialize dataframe for interpolated data
    df_interp <- data.frame()
    
    # Loop through each depth and slice data
    for (i in 1:length(depths)){
      
      fp_layer <- fp %>%
        group_by(Date) %>%
        slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df_interp <- bind_rows(df_interp, fp_layer) %>%
        dplyr::distinct()
    }
    
    # Wrangle final dataframe for interpolation
    fp_new <- arrange(df_interp, DateTime) %>%
      drop_na({{z}})
    
    # Round each extracted depth to the nearest 10th
    fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5)
    
    # Interpolate data for the current reservoir
    interp <- interp(x = as.numeric(fp_new$Date), y = fp_new$Depth_m, z = unlist(fp_new[z]),
                     xo = seq(min(fp_new$Date), max(fp_new$Date), by = "day"),
                     yo = seq(min(fp_new$Depth_m), max(fp_new$Depth_m), by = 0.1),
                     extrap = TRUE, linear = TRUE, duplicate = "strip")
    
    # interp <- interp2xyz(interp, data.frame = TRUE)
    interp <- interp2xyz(interp, data.frame = TRUE) %>%
      mutate(z = pmax(0, pmin(z, 1)))  # Constrain values between 0 and 1
    
    
    # Add 'Reservoir' and 'Year' to the interpolated data
    interp <- interp %>%
      mutate(Date = as.Date(x, origin = "1970-01-01"),
             Year = year(Date),
             Reservoir = reservoir)  # Keep track of which reservoir
    
    # Bind the interpolated data for each reservoir
    df_interp_final <- bind_rows(df_interp_final, interp)
  }
  
  # Create title for plot
  fig_title <- paste(paste(reservoirs, collapse = " & "), "Site", site, z, sep = " ")
  
  # HOxBoxes <- read_csv('HOxOnDates.csv') %>% 
  #   mutate(Group = rownames(.),
  #          Reservoir = 'FCR',
  #          Year = year(HOxOn),
  #          DOYon = yday(HOxOn),
  #          DOYoff = yday(HOxOff))
  
  # Plot the heatmap with facet_grid by Reservoir and Year
  p1 <- ggplot(df_interp_final, aes(x = Date, y = y)) +
    geom_raster(aes(fill = z)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_date(expand = c(0, 0)) +
    scale_fill_gradientn(colours = blue2red(60), na.value = "gray") +
    facet_grid(rows = vars(Reservoir), cols = vars(Year), scales = 'free_x', space = 'free_x') +
    # geom_rect(data = HOxBoxes, 
    #           aes(xmin = HOxOn, xmax = HOxOff, ymin = -Inf, ymax = Inf, fill = NULL, color = NULL), 
    #           alpha = 0.2) +  
    #geom_vline(xintercept = FCRHOx$Date, linetype = 1, colour = '#D55E00') + #HOx lines!
    labs(x = NULL, y = "Depth (m)", title = fig_title, fill = ("Dissolved:Total")) +
    theme_bw() 
  
  return(p1)
}

