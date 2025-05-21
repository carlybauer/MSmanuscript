# PCA metals
# 29 May 2024
# Carly Bauer


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

# load metal data using version through 2023 because not including 2024 in MS
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

FCRmetals4 <- FCRmetals %>% 
  select(TAl_mgL, TBa_mgL, TCu_mgL, TSr_mgL)

# save new df that only includes BVR site 50 data and only T metal concentrations 
BVRmetals <- metals %>%
  mutate(Year = year(DateTime)) %>% 
  filter(Reservoir == "BVR", Site == 50, Year >= 2020) %>%
  select(starts_with("T")) %>%
  #select(-Site) %>% 
  mutate(across(starts_with("T"), log10)) %>% # log10 concentrations to help linearity, not assumption just helps increase strength
  na.omit()  # Remove rows with any NA values

BVRmetals4 <- BVRmetals %>% 
  select(TAl_mgL, TBa_mgL, TCu_mgL, TSr_mgL)
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
pca_var_plot + ggtitle("PCA BVR")




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
             
             
             
             
             
             