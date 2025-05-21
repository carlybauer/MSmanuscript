# NMDS Metals
# Authors: CCC and CEB
# April 14 2025

# NMDS for metals in appendix / SI 

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)
library(FactoMineR)
library(ggfortify)
library(ggpubr)

# load metal data from EDI using version through 2023 because not including 2024 in MS
metals <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/455/8/9c8c61b003923f4f03ebfe55cea8bbfd")

# this gets rid of any values that we flagged as abnormally high in the data
metals <- metals %>% 
  filter(Flag_TLi_mgL != 8 & Flag_SLi_mgL != 8 & Flag_TLi_mgL != 68 & Flag_SLi_mgL != 68,
         Flag_TNa_mgL != 8 & Flag_SNa_mgL != 8 & Flag_TNa_mgL != 68 & Flag_SNa_mgL != 68,
         Flag_TMg_mgL != 8 & Flag_SMg_mgL != 8 & Flag_TMg_mgL != 68 & Flag_SMg_mgL != 68,
         Flag_TAl_mgL != 8 & Flag_SAl_mgL != 8 & Flag_TAl_mgL != 68 & Flag_SAl_mgL != 68,
         Flag_TSi_mgL != 8 & Flag_SSi_mgL != 8 & Flag_TSi_mgL != 68 & Flag_SSi_mgL != 68,
         Flag_TK_mgL != 8 & Flag_SK_mgL != 8 & Flag_TK_mgL != 68 & Flag_SK_mgL != 68,
         Flag_TCa_mgL != 8 & Flag_SCa_mgL != 8 & Flag_TCa_mgL != 68 & Flag_SCa_mgL != 68,
         Flag_TFe_mgL != 8 & Flag_SFe_mgL != 8 & Flag_TFe_mgL != 68 & Flag_SFe_mgL != 68,
         Flag_TMn_mgL != 8 & Flag_SMn_mgL != 8 & Flag_TMn_mgL != 68 & Flag_SMn_mgL != 68,
         Flag_TCu_mgL != 8 & Flag_SCu_mgL != 8 & Flag_TCu_mgL != 68 & Flag_SCu_mgL != 68,
         Flag_TSr_mgL != 8 & Flag_SSr_mgL != 8 & Flag_TSr_mgL != 68 & Flag_SSr_mgL != 68,
         Flag_TBa_mgL != 8 & Flag_SBa_mgL != 8 & Flag_TBa_mgL != 68 & Flag_SBa_mgL != 68)


## format data so its just the 12 metals, nothing else in dataframe for FCR and BVR
#CCC's attempt at an NMDS, following the excellent tutorial here: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
install.packages("vegan")
library(vegan)

FCRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "FCR", Site == 50, Year >= 2020) %>%
  select(Reservoir, starts_with("T")) %>%
  na.omit()  # Remove rows with any NA values

# save new df that only includes BVR site 50 data and only T metal concentrations 
BVRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
  select(Reservoir, starts_with("T")) %>%
  na.omit()  # Remove rows with any NA values

data1 <- merge(BVRmetals1, FCRmetals1, all.x=T, all.y=T)  

original_rownames <- data1$Reservoir

data2 = data1 |> 
  select(starts_with("T"))

example_NMDS=metaMDS(data2, 
                     k=2, distance="bray")
stressplot(example_NMDS)
plot(example_NMDS)

#I JUST RAN OUT OF TIME, BUT THE NEXT STEP IS TO COLOR THE POINTS BY RESERVOIR-
# YOU CAN CLEARLY SEE CLUSTERING OF FCR VS BVR HERE. YOU CAN THEN ALSO COMPARE THE POLYGONS
#SIMILAR TO WHAT HEATHER DID STATISTICALLY.

# Extract site scores and add Reservoir info
nmds_scores <- as.data.frame(scores(example_NMDS, display = "sites"))
nmds_scores$Reservoir <- original_rownames  # your saved group info

## CALCULATE CENTROID FOR EACH RESERVOIR
#gives (x,y) coordinates for each group's center
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Reservoir, data = nmds_scores, FUN = mean) 
# For each reservoir, it calculates the mean NMDS1 and mean NMDS2 coordinate (centroid of all its points in NMDS space)
print(centroids) 


## CALCULATE EUCLIDEAN DISTANCE BETWEEN TWO CENTROIDS
# Use dist() if just two groups:
centroid_distance <- dist(centroids[, c("NMDS1", "NMDS2")])
print(centroid_distance)

# # Plot with ggplot2
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Reservoir)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(geom = "polygon", aes(fill = Reservoir), alpha = 0.2) +
  geom_point(data = centroids, aes(x = NMDS1, y = NMDS2, color = "black"),  # Add centroids
             size = 5, shape = 4, stroke = 2) +  # shape 4 = X
  theme_minimal() +
  labs(title = "NMDS of Metal Concentrations by Reservoir",
       x = "NMDS1", y = "NMDS2") +
  scale_colour_manual(values = c("FCR" = '#56B4E9', "BVR" = '#D55E00')) +
  scale_fill_manual(values = c("FCR" = '#56B4E9', "BVR" = '#D55E00')) +
  theme(legend.position = "right")


# Save grouping variable (Reservoir) for plotting
groups <- original_rownames  # should be same length as nrow(data2)

# Plot
groups <- as.character(data1$Reservoir)
ordiplot(example_NMDS, type = "n")
ordihull(example_NMDS, groups = groups, draw = "polygon", col = c("#D55E00", "#56B4E9"), label = TRUE)
points(example_NMDS, display = "sites", col = ifelse(groups == "FCR", "#56B4E9", "#D55E00"), pch = 16) # for sites
# Add centroids
points(centroids$NMDS1, centroids$NMDS2, 
       pch = 4, cex = 1.5, lwd = 1, 
       col = c("black"))  # Adjust to match order


### FOLLOWING SIMILAR TO HEATHER
## not sure what method is right
# Ensure 'Reservoir' is a character vector
data1$Reservoir <- as.character(data1$Reservoir)

var_results <- data.frame(
  fcr_disp = rep(NA, 500),
  bvr_disp = rep(NA, 500),
  centroid_dist = rep(NA, 500)
)

set.seed(42)

for (i in 1:500) {
  
  sub_data <- data1 %>%
    group_by(Reservoir) %>%
    slice_sample(n = 10) %>%
    ungroup()
  
  sub_metals <- sub_data %>% select(starts_with("T"))
  sub_metals_hel <- decostand(sub_metals, method = "hellinger")
  sub_euc <- vegdist(sub_metals_hel, method = "euclidean")
  
  disp <- betadisper(sub_euc, group = sub_data$Reservoir, type = "centroid") 
  
  # Save distances to centroid for each group
  dist_df <- data.frame(Reservoir = sub_data$Reservoir,
                        dist = disp$distances)
  
  var_results$fcr_disp[i] <- mean(dist_df$dist[dist_df$Reservoir == "FCR"])
  var_results$bvr_disp[i] <- mean(dist_df$dist[dist_df$Reservoir == "BVR"])
  
  # Distance between centroids
  var_results$centroid_dist[i] <- dist(disp$centroids)[1]
}

# Compare dispersion
wilcox.test(var_results$fcr_disp, var_results$bvr_disp)
## p-value < 2.2e-16 using euclidian, but this doesnt make sense since I'm distancing using bray in NMDS
## p-value = 0.6857 using bray, not significantly different 

#summarize/visualize
disp_long <- var_results %>%
  select(fcr_disp, bvr_disp) %>%
  pivot_longer(everything(), names_to = "Reservoir", values_to = "Dispersion")

ggboxplot(disp_long, x = "Reservoir", y = "Dispersion", fill = "Reservoir")

#Plot centroid distances
ggplot(var_results, aes(x = centroid_dist)) +
  geom_histogram(binwidth = 0.01, fill = "#5E435D") +
  labs(x = "Distance between centroids (FCR vs BVR)",
       y = "Frequency",
       title = "Bootstrapped centroid distances")

##### NOW DO IT FOR OBSERVING METAL PATTERNS COMBINED
# Combine both
data_all <- data1 |> 
  select(starts_with("T"))

# Transpose to analyze metal relationships
# This makes rows = metals, columns = samples
# Transpose and preserve metal names
data_transposed <- as.data.frame(t(data_all))
data_transposed$Metal <- rownames(data_transposed)  # Save rownames (metal names) as a column
#rownames(data_transposed) <- NULL  # Optional: clear rownames


# Compute distance matrix among metals
cor_dist <- vegdist(data_transposed |> select(-Metal), method = "bray")

# Run NMDS on transposed data
metal_NMDS <- metaMDS(cor_dist, k = 2, trymax = 20)

# Extract scores (now each point is a metal)
metal_scores <- as.data.frame(scores(metal_NMDS, display = "sites"))
metal_scores$Metal <- data_transposed$Metal

# Plot
ggplot(metal_scores, aes(x = NMDS1, y = NMDS2, label = Metal)) +
  geom_point(color = "darkblue", size = 4) +
  geom_text(vjust = -1.2, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(title = "NMDS of Metal Similarity",
       x = "NMDS1", y = "NMDS2")

### OBSERVING METAL PATTERNS SEPARATE IN RESERVOIRS
# repetition of above
FCRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "FCR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%  # Keep only metal columns
  na.omit()

BVRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  na.omit()

# Transpose so rows = metals, columns = samples
FCR_t <- as.data.frame(t(FCRmetals1))
FCR_t$Metal <- rownames(FCR_t)
rownames(FCR_t) <- NULL

BVR_t <- as.data.frame(t(BVRmetals1))
BVR_t$Metal <- rownames(BVR_t)
rownames(BVR_t) <- NULL

# FCR NMDS
# distancing matrix applied
FCR_dist <- vegdist(FCR_t |> select(-Metal), method = "bray")
FCR_NMDS <- metaMDS(FCR_dist, k = 2, trymax = 20)

FCR_scores <- as.data.frame(scores(FCR_NMDS, display = "sites"))
FCR_scores$Metal <- FCR_t$Metal

# BVR NMDS
# distancing matrix applied
BVR_dist <- vegdist(BVR_t |> select(-Metal), method = "bray")
BVR_NMDS <- metaMDS(BVR_dist, k = 2, trymax = 20)

BVR_scores <- as.data.frame(scores(BVR_NMDS, display = "sites"))
BVR_scores$Metal <- BVR_t$Metal

# Load ggrepel if you haven’t already
library(ggrepel)

# FCR Plot with ggrepel
ggplot(FCR_scores, aes(x = NMDS1, y = NMDS2, label = Metal)) +
  geom_point(color = "#56B4E9", size = 3) +
  geom_text_repel(size = 4, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "FCR",
       x = "NMDS1", y = "NMDS2")

# BVR Plot
ggplot(BVR_scores, aes(x = NMDS1, y = NMDS2, label = Metal)) +
  geom_point(color = "#D55E00", size = 3) +
  geom_text_repel(size = 4, max.overlaps = 100) +
  theme_minimal() +
  labs(title = "BVR",
       x = "NMDS1", y = "NMDS2")

# --- Assign metal categories to group ---
metal_groups <- data.frame(
  Metal = c("TLi_mgL", "TNa_mgL", "TMg_mgL", "TK_mgL", "TCa_mgL", 
            "TFe_mgL", "TMn_mgL", "TAl_mgL", "TCu_mgL", "TBa_mgL", 
            "TSr_mgL", "TSi_mgL"),
  Group = c("Physical", "Chemical", "Chemical", "Chemical", "Chemical", 
            "Redox", "Redox", "Physical", "Redox", "Redox", 
            "Chemical", "Chemical"))

# --- Plot with assigned groupings
FCR_scores <- FCR_scores %>% left_join(metal_groups, by = "Metal")
ggplot(FCR_scores, aes(x = NMDS1, y = NMDS2, color = Group, label = Metal)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE) +
  geom_point(size = 3) +
  geom_text(vjust = -1.2, hjust = 0.5, size = 3) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(title = "NMDS of Metal Similarity - FCR",
       subtitle = "Polygons group metals by category",
       x = "NMDS1", y = "NMDS2")


##PLOT WITH ASSIGNED GROUPINGS FOR BOTH RES COMBINED
# STEP 1: Clean and select metal data from both reservoirs
FCRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "FCR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  na.omit()

BVRmetals1 <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  na.omit()

# STEP 2: Combine both datasets (optional: log-transform if needed)
data_all <- bind_rows(FCRmetals1, BVRmetals1)

# STEP 3: Transpose so rows = metals, columns = samples
data_transposed <- as.data.frame(t(data_all))
data_transposed$Metal <- rownames(data_transposed)
rownames(data_transposed) <- NULL

# STEP 4: Manually define metal groups
metal_groups <- data.frame(
  Metal = c("TLi_mgL", "TNa_mgL", "TMg_mgL", "TK_mgL", "TCa_mgL", 
            "TFe_mgL", "TMn_mgL", "TAl_mgL", "TCu_mgL", "TBa_mgL", 
            "TSr_mgL", "TSi_mgL"),
  Group = c("Physical", "Chemical", "Chemical", "Chemical", "Chemical", 
            "Redox", "Redox", "Physical", "Redox", "Redox", 
            "Chemical", "Chemical")
)

# STEP 5: Join groups to transposed data
data_transposed <- left_join(data_transposed, metal_groups, by = "Metal")

# STEP 6: NMDS on metal rows (drop metadata columns first)
cor_dist <- vegdist(data_transposed |> select(-Metal, -Group), method = "bray")
metal_NMDS <- metaMDS(cor_dist, k = 2, trymax = 20)

# STEP 7: Extract NMDS scores + add back group info
metal_scores <- as.data.frame(scores(metal_NMDS, display = "sites"))
metal_scores$Metal <- data_transposed$Metal
metal_scores$Group <- data_transposed$Group

# STEP 8: Plot with polygons or ellipses by group
ggplot(metal_scores, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21, color = "black") +
  geom_text(aes(label = Metal), vjust = -1.2, hjust = 0.5, size = 3) +
  stat_ellipse(aes(group = Group), geom = "polygon", alpha = 0.2) +
  theme_minimal() +
  labs(title = "NMDS of Metal Similarity by Grouping",
       subtitle = "Points = Metals; Polygons = Manual Groupings",
       x = "NMDS1", y = "NMDS2")