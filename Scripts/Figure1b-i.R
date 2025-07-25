#############  Figure 1b-i ###################
setwd("./Data/Figure1/")

# Load required libraries

library(tidyverse)      # Data manipulation & ggplot2
library(dplyr)          # Data manipulation
library(ggplot2)        # Plotting
library(FactoMineR)     # PCA
library(factoextra)     # PCA visualization
library(missMDA)        # Missing data imputation
library(vegan)          # PERMANOVA
library(lme4)           # Mixed effects models
library(emmeans)        # Estimated marginal means
library(gridExtra)      # Arranging plots
library(reshape2)       # Data reshaping


# Define sample list

samples_to_keep <- c(
  "MGr3I", "MGr3H", "MGr3G", "MGr2F", "MGr2E", "MGr2D", "MGr1C", "MGr1B", "MGr1A",
  "CRr3I", "CRr3H", "CRr3G", "CRr2F", "CRr2E", "CRr2D", "CRr1C", "CRr1B", "CRr1A",
  "MHr3I", "MHr3H", "MHr3G", "MHr2F", "MHr2E", "MHr2D", "MHr1C", "MHr1B", "MHr1A",
  "BOr3I", "BOr3H", "BOr3G", "BOr2F", "BOr2E", "BOr2D", "BOr1C", "BOr1B", "BOr1A",
  "BAr3I", "BAr3H", "BAr3G", "BAr2F", "BAr2E", "BAr2D", "BAr1C", "BAr1B", "BAr1A",
  "SEr3I", "SEr3H", "SEr3G", "SEr2F", "SEr2E", "SEr2D", "LASCr2F", "LASCr2E", "LASCr2D",
  "LABRr1C", "LABRr1B", "LABRr1A", "LASAr3I", "LASAr3G", "LASAr3H"
)

# Load data
# Oxygen data

O2 <- read.csv("oxygen_probe.csv") %>%
  filter(Sample %in% samples_to_keep) %>%
  mutate(Depth = factor(Depth, levels = c("0_5", "5_10", "40_45", "45_50"))) %>%
  spread(key = Depth, value = O2) %>%
  mutate(
    oxic = rowMeans(select(., "0_5", "5_10"), na.rm = TRUE),
    anoxic = rowMeans(select(., "40_45", "45_50"), na.rm = TRUE)
  )

# Moisture & bulk density

bulk <- read.csv("Bulk_Density.csv") %>%
  filter(Sample %in% samples_to_keep) %>%
  mutate(
    volume        = Volume_After - Volume_before,
    wet_weight    = weight_before - Crucible_weight,
    dry_weight    = weight_after - Crucible_weight,
    moisture      = ((wet_weight - dry_weight) / wet_weight) * 100,
    bulk_density  = dry_weight / volume,
    CUE_addition  = 0.334 * (moisture / 100) * 300,
    wet_bulk      = wet_weight / Volume_before
  )

bulk_mean <- bulk %>%
  group_by(Treatment) %>%
  summarise(moisture_mean = mean(moisture, na.rm = TRUE))


# Vegetation data

veg <- read.csv("vegetation_percentage.csv") %>%
  left_join(read.csv("Sample_ID.csv"), by = "Sample") %>%
  mutate(vascular = eric + gram)

TC <- read.csv("TC_TN.csv") %>%
  left_join(read.csv("Sample_ID.csv"), by = "Sample") %>%
  filter(Sample %in% samples_to_keep)

# pH

ph <- read.csv("Ph.csv") %>%
  filter(Sample %in% samples_to_keep & Depth == "0_10")

# Water table

wt_raw <- read.csv("Site_Data.csv")

wt_overall <- wt_raw %>%
  group_by(ID, Site, Region) %>%
  summarise(mean_wt = mean(WTD.mm., na.rm = TRUE) / 10)

wt_eco_index <- read.csv("wt_eco_index.csv")

wt_mean <- wt_eco_index %>%
  group_by(Treatment, Site) %>%
  summarise(mean_wt = mean(mean_wt, na.rm = TRUE))

dna_prod<- read.csv("CUE_calc_ashish.csv")

eco_index_mean <- dna_prod %>%
  group_by(Treatment, Site) %>%
  summarise(mean_eco = mean(degradation_index, na.rm = TRUE))

merged_data <- wt_mean %>%
  left_join(eco_index_mean, by = c("Treatment", "Site"))

#FTIR
ftir_pca1<- read.csv("FTIR_PCA1_export.csv")

# PCA of environmental variables
ftir_pca1<- ftir_pca1 %>%
  select(Sample, Dim.1, Dim.2)


bulk_pca<- bulk %>%
  select(Sample, Site, Treatment, moisture, bulk_density)

veg_pca<- veg %>%
  select(Sample, sphag, vascular)

o2_pca<- O2_combine %>%
  select(Sample, oxic)

TC_pca<- TC %>%
  select(Sample, X.C, X.N)

TC_pca<- TC %>%
  select(Sample, X.C, X.N, CN)

ph_pca<- ph %>%
  select(Sample, ph)

# Figure 1b
# Merge environmental datasets
env_pca<- bulk_pca %>%
  left_join(x=.,y=ph_pca, by=c("Sample")) %>%
  left_join(x=.,y=veg_pca, by=c("Sample"))   %>%
  left_join(x=.,y=o2_pca, by=c("Sample")) %>%
  left_join(x=.,y=TC_pca, by=c("Sample")) %>%
  left_join(x=.,y=ftir_pca1, by=c("Sample"))

# Prepare data for PCA
env_pca_m <- env_pca %>%
  select(-Sample, -Treatment, -Site) %>%
  rename(
    pH = ph,
    CN = CN,
    FTIR_axis_1 = Dim.1,
    O2 = oxic,
    Moss_cover = sphag,
    Moisture = moisture
  ) %>%
  mutate(across(where(is.numeric), scale))

# center and scale the data
for (i in 1:length(colnames(env_pca_m))){
  
  if (is.numeric(env_pca_m[, i])==TRUE)
    
    env_pca_m[, i] <- as.numeric(scale(env_pca_m[, i]))
  
  else
    
    env_pca_m[, i] <- env_pca_m[, i]
  
}

# Impute missing data and run PCA
nb <- estim_ncpPCA(env_pca_m, ncp.max = 5)
env_pca_m_imputed <- imputePCA(env_pca_m, ncp = nb$ncp)$completeObs

#extract PCA scores for individuals and loadings for variables
pca_ind <- as.data.frame(env_pca_result$ind$coord)  # Individual scores
pca_var <- as.data.frame(env_pca_result$var$coord)  # Variable loadings

# Add additional columns for plotting
pca_ind$Site <- env_pca$Site
pca_ind$Treatment <- env_pca$Treatment

# Extract variable loadings (coordinates) and contributions
pca_var <- as.data.frame(env_pca_result$var$coord)  # Loadings (coordinates)
var <- get_pca_var(env_pca_result)
pca_var$contrib <- var$contrib[, 1] + var$contrib[, 2]  # Summing contributions to Dim.1 and Dim.2

# Save PCA contribution plots
# Define a scaling factor to adjust arrow sizes
scaling_factor <- 3.2  # Adjust this factor to increase/decrease arrow size --> Default Fviz biplot arrows are much larger than woth extracted values, unsure why or what sclaing factor (if any they use) 
site_colors <- c(
  "Migneint"  = "#e41a1c",
  "Bowness" = "#377eb8",
  "Crocach" = "#4daf4a",
  "Langwell" = "#984ea3",
  "Migneint" = "#ff7f00",
  "Moors_House" = "#ffff33",
  "Stean" = "#a65628"
)


pca_var$custom_labels <- c( "Moisture" ,  "pH", "Moss cover",  "O2" , "C:N", "OM chemistry" )

pca_ind$Treatment <- factor(pca_ind$Treatment, levels = c("NAT", "REST", "DAM"))

ggplot(data = pca_ind, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(fill = Site, shape = Treatment), 
             size = 3, color = "black", alpha = 0.7) +  # Ensure color is only for outline
  # Define shape legend for Treatment
  scale_shape_manual(values = c(21, 22, 23), 
                     labels = c("Natural", "Restored", "Degraded")) +
  stat_ellipse(
    aes(group = Treatment, fill = Treatment), 
    geom = "polygon", alpha = 0.1, color = NA, level = 0.9) +
  
  # Define fill legend for Site
  scale_fill_manual(values = c("Balmoral" = "#e41a1c", 
                               "Bowness" = "#377eb8", 
                               "Crocach" = "#4daf4a", 
                               "Langwell" = "#984ea3", 
                               "Migneint" = "#ff7f00", 
                               "Moors_House" = "#ffff33", 
                               "Stean" = "#a65628"), 
                    labels = c("Balmoral", "Bowness", "Crocach", 
                               "Langwell", "Migneint", "Moor House", "Stean")) +
  
  xlab("Dim1 (35.4%)") +
  ylab("Dim2 (29.4%)") +
  coord_cartesian(ylim = c(-3, 3)) + 
  
  geom_segment(data = pca_var, 
               aes(x = 0, y = 0, xend = Dim.1 * scaling_factor, yend = Dim.2 * scaling_factor, color = contrib), 
               arrow = arrow(length = unit(0.2, "cm")), size = 0.8, alpha =0.8) +
  
  geom_text(data = pca_var, 
            aes(x = Dim.1 * scaling_factor *1.15, y = Dim.2 * scaling_factor *1.15, label = custom_labels), 
            color = "black", hjust = 0.5, vjust = 0.5, fontface = "bold") +
  
  scale_color_gradient2(low = "blue", mid = "grey", high = "red", 
                        midpoint = mean(pca_var$contrib, na.rm = TRUE)) +
  
  # Ensure Site colors appear in the legend
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment"), 
    color = guide_colorbar(order = 3, title = "Contribution") 
  ) +
  
  labs(color = "Contribution") +
  theme_light() +
  theme(legend.position = "right")


ggsave("Figure_1b_PCA.png",  dpi=1000, width = 7, height = 5.5, units = "in") 

# Figure 1c
# Ensure rows match
metadata <- env_pca %>%
  select(Sample, Site, Treatment)

# Optional: Check alignment
all(rownames(env_pca_m_imputed) == metadata$Sample) # Should be TRUE

# Run PERMANOVA
permanova <- adonis2(env_pca_m_imputed ~ Site * Treatment, data = metadata, method = "euclidean", permutations = 999)

# View results
print(permanova)


#Ecosystem health index

env_pca_s<- env_pca %>%
  select(-c("Dim.2" ,"bulk_density","vascular", "X.C","X.N")) %>%
  dplyr::rename(pH = ph) %>%
  dplyr::rename(C_N = CN) %>%
  dplyr::rename(FTIR_axis_1 = Dim.1) %>%
  dplyr::rename(O2 = oxic) %>%
  dplyr::rename(Moss_cover = sphag) %>%
  dplyr::rename(Moisture = moisture) 


env_pca_s <- env_pca %>%
  select(-c("Dim.2", "bulk_density", "vascular", "X.C", "X.N")) %>%
  rename(pH = ph, C_N = CN, FTIR_axis_1 = Dim.1, O2 = oxic, Moss_cover = sphag, Moisture = moisture)

# Standardize variables (but scale them by a factor to reduce range)
scaling_factor <- 0.01  # You can adjust this value as needed

# Standardize variables
env_data_standardized <- env_pca_s %>%
  mutate(across(c("Moisture", "pH", "Moss_cover", "FTIR_axis_1", "O2", "C_N"), ~scale(.) * scaling_factor)) %>%
  mutate(O2 = -O2, FTIR_axis_1 = -FTIR_axis_1, C_N = -C_N)  # Reverse selected variables

# Define weights
weights <- c(Moisture = 17, pH = 12.4, Moss_cover = 13.5, FTIR_axis_1 = 18.1, O2 = 19.8, C_N = 19.1)

# Adjust weights where Moss_cover is missing and add adjusted_weights column
env_data_standardized <- env_data_standardized %>%
  rowwise() %>%
  mutate(
    adjusted_weights = list(
      if (is.na(Moss_cover)) {
        # Reweight variables proportionally if Moss_cover is missing
        new_weights <- weights[names(weights) != "Moss_cover"]
        new_weights / sum(new_weights) * sum(weights)  
      } else {
        weights
      }
    )
  ) %>%
  ungroup()

# Compute degradation index, ensuring adjusted_weights is correctly referenced
env_data_standardized <- env_data_standardized %>%
  rowwise() %>%
  mutate(
    # Ensure adjusted_weights is correctly used by unlisting and matching variable names
    degradation_index = sum(c_across(all_of(names(weights))) * unlist(adjusted_weights)[names(weights)], na.rm = TRUE)
  ) %>%
  ungroup()

env_data_standardized$Site <- factor(env_data_standardized$Site,
                                     levels = rev(c("Crocach", "Langwell", "Balmoral", 
                                                    "Bowness", "Moors_House", "Stean", "Migneint")))


env_data_standardized$Treatment <- factor(env_data_standardized$Treatment, 
                                          levels = c("NAT", "REST", "DAM"))


site_means <- env_data_standardized %>%
  group_by(Site) %>%
  summarise(mean_index = mean(degradation_index, na.rm = TRUE))

mean_points <- env_data_standardized %>%
  group_by(Site, Treatment) %>%
  summarise(mean_index = mean(degradation_index, na.rm = TRUE), .groups = "drop")

ggplot(env_data_standardized, aes(y = Site, x = degradation_index, fill = Site)) +
  
  # Violin plot
  geom_violin(alpha = 0.7) +  
  # Add mean values as horizontal dashes
  stat_summary(fun = mean, geom = "crossbar", width = 0.15, color = "white", size = 0.3) +
  
  # Axis Labels
  labs(x = "Ecosystem health index") +
  # Custom Site Colors
  scale_fill_manual(values = c("Balmoral" = "#e41a1c", 
                               "Bowness" = "#377eb8", 
                               "Crocach" = "#4daf4a", 
                               "Langwell" = "#984ea3", 
                               "Migneint" = "#ff7f00", 
                               "Moors_House" = "#ffff33", 
                               "Stean" = "#a65628"), 
                    labels = c("Balmoral", "Bowness", "Crocach", 
                               "Langwell", "Migneint", "Moor House", "Stean")) +
  scale_y_discrete(labels = function(y) gsub("Moors_House", "Moor House", y)) +
  
  # Add mean points by Treatment shape
  geom_point(data = mean_points, 
             aes(x = mean_index, y = Site, shape = Treatment), 
             size = 3, fill = "white", color = "black", stroke = 1) +
  
  scale_shape_manual(values = c(21, 22, 23), 
                     labels = c("NAT", "REST", "DAM")) +
  
  # Remove x-axis labels, title, minor gridlines, and legend
  theme_light() +
  theme(
    axis.title.x = element_text(size = 13, colour = "black"),  # Adjust size for y-axis title
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 13, colour = "black"),   # Adjust size for y-axis text
    axis.text.x = element_text(size = 13, colour = "black"),   # Adjust size for x-axis text
    panel.grid.minor = element_blank(),
    legend.position = "none")

ggsave("Figure_1c_EcosystemHealthIndex.png",  dpi=1000, width = 4.5, height =5, units = "in") 

# Figure 1d
### Plot oxygen concentration
o2_plot<- ggplot(env_pca_s, mapping=aes(y=O2, x=Treatment))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=O2, x=Treatment, colour=Site), alpha=0.7, width = 0.3, size=3)+
  #ylab("Oxygen Concentration (mg/L)")+
  ylab(expression("Oxygen concentration (mg L  "^{ -1}*")"))+ 
  xlab(" ")+
  theme_light(base_size = 14)+
  theme(legend.position = "None")+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black"))


ggsave("Figure_1d_oxygen.png",  dpi=1000, width = 4, height = 4, units = "in")

# Figure 1e
### Plot moisture 

moisture_plot<- ggplot(env_pca_s,aes(x=Treatment,y=Moisture))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=Moisture, x=Treatment, color=Site), alpha=0.7, width = 0.3, size=3)+
  scale_y_continuous(name="Moisture (%)",breaks=c(75,80,85,90,95),limits = c(75, 95))+
  xlab(" ")+
  theme_light(base_size = 14)+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black"))


ggsave("Figure_1e_moisture.png",  dpi=1000, width = 4, height = 4, units = "in")

# Figure 1f
### Plot FTIR PCA axis 1 across sites and treatments 

ftir_plot<- ggplot(env_pca_s, mapping=aes(y=FTIR_axis_1, x=Treatment))+
  geom_boxplot(lwd=1,outlier.shape = NA)+
  geom_jitter(aes(y=FTIR_axis_1, x=Treatment, colour=Site), alpha=0.7, width = 0.3, size=3)+
  #ylab("Oxygen Concentration (mg/L)")+
  ylab(expression("OM chemistry (FTIR PC1)"))+ 
  xlab(" ")+
  theme_light(base_size = 14)+
  theme(legend.position = "None")+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) +
  scale_y_continuous(limits = c(-10, 12)) 


ggsave("Figure_1f_FTIR.png",  dpi=1000, width = 4, height = 4, units = "in")


# Figure 1g
###### Plot C:N ratio 

env_pca_s <- env_pca_s %>%
  mutate(CN_ratio = TC / TN)  # Calculate C:N ratio

CN_ratio_plot<- ggplot(env_pca_s, mapping=aes(y=C_N, x=Treatment))+
  geom_boxplot(lwd=1,outlier.shape = NA)+
  geom_jitter(aes(y=C_N, x=Treatment, colour=Site), alpha=0.7, width = 0.3, size=3)+
  #ylab("Oxygen Concentration (mg/L)")+
  ylab(expression("C:N ratio"))+ 
  xlab(" ")+
  theme_light(base_size = 14)+
  theme(legend.position = "None")+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

ggsave("Figure_1g_CN_ratio.png",  dpi=1000, width = 4, height = 4, units = "in")

# Figure 1h
### Plot pH 

ph_plot<- ggplot(env_pca_s,aes(x=Treatment,y=pH))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=pH, x=Treatment, color=Site), alpha=0.7, width = 0.3, size=3)+
  ylab("pH")+
  xlab(" ")+
  theme_light(base_size = 14)+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  scale_y_continuous(limits = c(3.5,4.3)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black"))

ggsave("Figure_1h_pH.png",  dpi=1000, width = 4, height = 4, units = "in")

# Figure 1i
### Plot moss cover 
 
sphag_plot<- ggplot(env_pca_s,aes(x=Treatment,y=Moss_cover))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=Moss_cover, x=Treatment, color=Site), alpha=0.7, width = 0.3, size=3)+
  ylab("Moss cover (%)")+
  xlab(" ")+
  theme_light(base_size = 14)+
  scale_colour_manual(values=c( "#e41a1c",
                                "#377eb8",
                                "#4daf4a",
                                "#984ea3",
                                "#ff7f00",
                                "#ffff33",
                                "#a65628"), 
                      labels = c("Balmoral"  = "Balmoral",
                                 "Bowness" = "Bowness",
                                 "Crocach" = "Crocach",
                                 "Langwell" = "Langwell",
                                 "Migneint" = "Migneint",
                                 "Moors_House" = "Moor House",
                                 "Stean" = "Stean" )) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black"))

ggsave("Figure_1i_moss_cover.png",  dpi=1000, width = 4, height = 4, units = "in")

