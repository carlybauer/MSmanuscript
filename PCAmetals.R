# PCA metals
# 29 May 2024
# Carly Bauer

install.packages("factoextra")
install.packages("FactoMineR")
install.packages("ggfortify")

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)
library(stringr)
library(ggplot2)
library(readr)
library(factoextra) # for visualizing PCA 
library(stats) # prcomp to perform PCA
library(FactoMineR)
library(ggfortify)

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
# save new df that only includes FCR site 50 data and only T metal concentrations 
FCRmetals <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "FCR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  #select(-Site) %>% 
  mutate(across(starts_with("T"), log10)) %>% # log10 concentrations to help linearity, not assumption just helps increase strength
  na.omit()  # Remove rows with any NA values

# save new df that only includes BVR site 50 data and only T metal concentrations 
BVRmetals <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  #select(-Site) %>% 
  mutate(across(starts_with("T"), log10)) %>% # log10 concentrations to help linearity, not assumption just helps increase strength
  na.omit()  # Remove rows with any NA values


#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# do PCA
PCAmetals <- prcomp(BVRmetals, scale=TRUE) # rerun this line and below for other reservoir
summary(PCAmetals)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::
## scree plot to visualize percentages of variance explained by each principal component
 fviz_eig(PCAmetals, addlabels = TRUE)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# what variables contribute the most to PC1 and PC2 
fviz_cos2(PCAmetals, choice = "var", axes =1:2 )

 
 ## graph of variables - need to rerun line 58 for BVR results 
 #::::::::::::::::::::::::::::::::::::::::::::::::::::::
 # positive correlated variables point to same direction 
 # negative correlated variables point to opposite direction of other variables
 pca_var_plot<-fviz_pca_var(PCAmetals,
                            col.var = "contrib", # Color by contributions to the PC
                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # those with similar cos2 values will plot similar colors 
                            repel = TRUE     # Avoid text overlapping
 )
 library(ggplot2)
 pca_var_plot + ggtitle("PCA BVR Metal Concentrations 2020 - 2023")

##################
 #CCC's attempt at an NMDS, following the excellent tutorial here: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
 install.packages("vegan")
 library(vegan)
 
 FCRmetals1 <- metals %>%
   mutate(Year = year(DateTime)) %>% 
   filter(Reservoir == "FCR", Site == 50, Year >= 2020) %>%
   select(Reservoir, starts_with("T")) %>%
   #mutate(across(starts_with("T"), log10)) %>% # log10 concentrations to help linearity, not assumption just helps increase strength
   na.omit()  # Remove rows with any NA values
 
 # save new df that only includes BVR site 50 data and only T metal concentrations 
 BVRmetals1 <- metals %>%
   mutate(Year = year(DateTime)) %>% 
   filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
   select(Reservoir, starts_with("T")) %>%
   #mutate(across(starts_with("T"), log10)) %>% # log10 concentrations to help linearity, not assumption just helps increase strength
   na.omit()  # Remove rows with any NA values
 
 data1 <- merge(BVRmetals1, FCRmetals1, all.x=T, all.y=T)  

original_rownames <- data1$Reservoir

data2 = data1 |> 
  select(starts_with("T"))
 
 example_NMDS=metaMDS(data2, 
                      k=2, distance="bray")
 stressplot(example_NMDS)
 plot(example_NMDS)

#I JUST RAN OUT OF TIME, BUT CARLY THE NEXT STEP IS TO COLOR THE POINTS BY RESERVOIR-
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
  #sub_bray <- vegdist(sub_metals_hel, method = "bray")
  
  
  disp <- betadisper(sub_euc, group = sub_data$Reservoir, type = "centroid") #change for method 
  
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

# Optional: log-transform to help linearity (commented out if not needed)
# data_transposed <- log10(data_transposed + 1e-6)

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
 

## CEB stops here, stuff below was just play when first figuring out PCA
 #::::::::::::::::::::::::::::::::::::::::::::::::::::::
 #::::::::::::::::::::::::::::::::::::::::::::::::::::::
 # if wanting to looking at both T and S on same PCA, can then use identifier as
 # BVR FCR. Making new dataframe for that
 # allMetals <- metals %>% 
 #   filter(Reservoir =="FCR"| Reservoir =="BVR" & Site ==50) %>% 
 #   select(starts_with("S"), starts_with("T")) %>% 
 #   select(-Site) %>% 
 #   na.omit()  # Remove rows with any NA values
 
 # ResMetals <- metals %>% 
 #   filter(Reservoir =="FCR"| Reservoir =="BVR" & Site ==50) %>% 
 #   select("Reservoir", starts_with("S"), starts_with("T")) %>% 
 #   select(-Site) %>% 
 #   na.omit()  # Remove rows with any NA values
 
 
 # when doing soluble PCA, some metals have constant concentrations which doesn't
 # allow for the PCA to scale, so we can remove those constant columns with code
 # below. only use this for soluble
 # BVRmetals<- BVRmetals %>% 
 #     select_if(~n_distinct(.) > 1)
 # FCRmetals <- FCRmetals %>% 
 #     select_if(~n_distinct(.)>1)
 
## graph with both BVR and FCR as groupings 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# scores <- as.data.frame(PCAmetals$x)
# scores$Reservoir <- ResMetals$Reservoir
# 
# ggplot(scores, aes(x = PC1, y = PC2, color = Reservoir)) +
#   geom_point() +
#   theme_minimal() +
#   labs(title = "PCA of Metals", x = "PC1", y = "PC2")

#####
# library(ggfortify)
# autoplot(PCAmetals, data = ResMetals, colour = 'Reservoir',
#          loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3,
#          loadings.label.colour = "black")

## graph of individuals
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# individuals with a similar profile are grouped together
# doesn't work well because rows have no meaning 
# fviz_pca_ind(PCAmetals,
#              col.ind = "cos2",  # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE       # Avoid text overlapping
# )



## another way to plot results 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# loading library 
library(ggfortify) 
autoplot(PCAmetals, 
        data = FCRmetals) 

## contributions of variables plotted
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
fviz_contrib(PCAmetals, choice = "var", axes = 1, top = 10)
fviz_contrib(PCAmetals, choice = "var", axes = 2, top = 10)
fviz_contrib(PCAmetals, choice = "var", axes = 3, top = 10)
fviz_contrib(PCAmetals, choice = "var", axes = 4, top = 10)

## access to PCA results 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Eigenvalues
eig.val <- get_eigenvalue(PCAmetals)
eig.val

# Results for Variables
res.var <- get_pca_var(PCAmetals)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(PCAmetals)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

## PCA results for variables 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# helper function
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# compute coordinates
loadings <- PCAmetals$rotation
sdev <- PCAmetals$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])

#compute Cos2
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

# Compute contributions
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])

#::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::
# K means

Kmetals <- scale(FCRmetals)
head(Kmetals)
kmMetals<-kmeans(Kmetals, 4, iter.max = 10, nstart = 25)
print(kmMetals)

#load library
library("FactoMineR")

metals.pca <- PCA(FCRmetals, graph = FALSE)

kmMetals <- kmeans(Kmetals, centers = 3, nstart = 25)
grp <- as.factor(kmMetals$cluster)

fviz_pca_var(metals.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster FCR PCA")
             
             
             
             
             
             
