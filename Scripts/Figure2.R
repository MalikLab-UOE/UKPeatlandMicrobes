#############  Figure 2 ###################
setwd("./Data/Figure2/")

# Load required libraries
library(tidyverse)
library(scales) 
library(nlme)
library(lme4)
library(car)
library(ggplot2)
library(ggcorrplot)
library(dplyr)

# Figure 2a
dna_prod<- read.csv("CUE_calc.csv")

dna_prod<- dna_prod %>%
  filter(Sample %in% c(  "MGr3I"  , "MGr3H",   "MGr3G"  , "MGr2F"  , "MGr2E"  , "MGr2D" ,  "MGr1C"  , "MGr1B" , 
                         "MGr1A" ,  "CRr3I" ,  "CRr3H" ,  "CRr3G" ,  "CRr2F"  , "CRr2E"  , "CRr2D"  , "CRr1C"  ,
                         "CRr1B" ,  "CRr1A"  , "MHr3I" ,  "MHr3H" ,  "MHr3G"  , "MHr2F"  , "MHr2E" ,  "MHr2D",  
                         "MHr1C" ,  "MHr1B" ,  "MHr1A",   "BOr3I" ,  "BOr3H"  , "BOr3G" ,  "BOr2F" ,  "BOr2E"  ,
                         "BOr2D" ,  "BOr1C"  , "BOr1B" ,  "BOr1A",   "BAr3I" ,  "BAr3H"  , "BAr3G" ,  "BAr2F",  
                         "BAr2E" ,  "BAr2D" ,  "BAr1C"  , "BAr1B" ,  "BAr1A" ,  "SEr3I"  , "SEr3H"  , "SEr3G",  
                         "SEr2F"  , "SEr2E"  , "SEr2D"  , "LASCr2F", "LASCr2E",
                         "LASCr2D", "LABRr1C" ,"LABRr1B" ,"LABRr1A" ,"LASAr3I", "LASAr3H" ,"LASAr3G" )
         
         
  )

growth_rate_treatment_plot<- ggplot(dna_prod,aes(x=Treatment, y=Growth))+
  geom_boxplot(aes(color = Treatment), linewidth = 1,lwd=1,outlier.shape= NA)+
  geom_jitter(aes(x=Treatment, y=Growth, colour=Treatment), alpha=0.7, width = 0.3, size=3)+
  ylim(log2(1),log2(2050))+
  ylab(expression("Microbial growth rate (ng g"^{-1}*"h"^{-1}*")"))+ 
  xlab(" ")+
  scale_y_continuous(trans='log2',labels = label_number(accuracy = 0.1))+
  theme_light(base_size = 14)+
  #scale_colour_manual(values = c( "Red","Blue"))+
  #scale_colour_discrete(name = "Depth (cm)", labels = c("0-10", "40-50"))+
  #theme(legend.position = "None")+
  scale_colour_manual(values = c("DAM" = "#e41a1c","REST" = "#377eb8","NAT" = "#4daf4a"),
                      labels = c("DAM" = "Degraded", "REST" = "Restored", "NAT" = "Natural")) +
  scale_x_discrete(limits=c("DAM", "REST", "NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black"))

ggsave("Figure_2a_MicrobialGrowthRate_Treatment.png",  dpi=1000, width = 4, height = 4, units = "in") 


# Figure 2b
# Compute mean Growth for each Site and Treatment
mean_growth <- dna_prod %>%
  group_by(Site, Treatment) %>%
  summarise(mean_growth = mean(Growth, na.rm = TRUE), .groups = "drop")

site_order <- c("Balmoral","Bowness","Stean","Moors_House","Langwell","Crocach","Migneint")

dna_prod$Site <- factor(dna_prod$Site, levels = site_order)
mean_growth$Site <- factor(mean_growth$Site, levels = site_order)

# Create the plot
growth_rate_plot <- ggplot(dna_prod, aes(y = Growth, x = Site, colour = Treatment)) +
  geom_point(size = 4,shape=21) +  # Individual points
  geom_point(data = mean_growth, aes(y = mean_growth, x = Site), size = 5, alpha = 0.7) +  # Mean points
  geom_line(data = mean_growth %>% filter(Treatment %in% c("REST", "DAM")),  
            aes(y = mean_growth, x = Site, group = Site), colour = "black", linetype = "dashed") +  # Connect REST & DAM means
  scale_y_continuous(trans = 'log2', labels = label_number(accuracy = 0.1)) +
  ylab(expression("Microbial biomass production (ng g"^{-1}*"h"^{-1}*")")) +
  xlab("Site") +
  theme_light(base_size = 14) +
  scale_colour_manual(values = c("NAT" = "#4daf4a", "REST" = "#377eb8", "DAM" = "#e41a1c"),
                      labels = c("NAT" = "Natural", "REST" = "Restored", "DAM" = "Degraded")) +
  theme(axis.text = element_text(size = 13, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) ifelse(x == "Moors_House", "Moor House", x))  # Change label in plot only

ggsave("Figure_2b_MicrobialGrowthRate_Site.png",  dpi=1000, width = 6, height =5, units = "in") 


# Figure 2c
##### Prepare data
dna_prod_test <- dna_prod
dna_prod_test$Treatment <- factor(dna_prod_test$Treatment, levels = c("NAT", "DAM", "REST"))

##### Fit linear mixed-effects models with different variance structures

# Model 1: Homogeneous variance
growth_rate_lme_homo <- lme(
  log(Growth) ~ Treatment,
  random = ~1 | Site,
  data = dna_prod_test
)

summary(growth_rate_lme_homo)
plot(growth_rate_lme_homo)

# Model 2: Heterogeneous variance by Site
growth_rate_lme_hetero <- lme(
  log(Growth) ~ Treatment,
  random = ~1 | Site,
  weights = varIdent(form = ~1 | Site),
  data = dna_prod_test
)

summary(growth_rate_lme_hetero)
plot(growth_rate_lme_hetero)

# Model 3: Heterogeneous variance by Treatment
growth_rate_lme_hetero_treat <- lme(
  log(Growth) ~ Treatment,
  random = ~1 | Site,
  weights = varIdent(form = ~1 | Treatment),
  data = dna_prod_test
)

# Compare models using AIC
AIC(growth_rate_lme_homo, growth_rate_lme_hetero, growth_rate_lme_hetero_treat)

# ANOVA for best model
summary(growth_rate_lme_homo)
Anova(growth_rate_lme_homo, type = "III")

##### Model with nested random effect: Treatment within Site
growth_rate_lme_nested <- lme(
  log(Growth) ~ Treatment,
  random = ~1 | Site/Treatment,
  data = dna_prod_test
)

summary(growth_rate_lme_nested)
Anova(growth_rate_lme_nested, type = "III")
plot(growth_rate_lme_nested)

##### Fit LMMs with additional random effects
# Multiple random effects: MAT, MAP, elevation
growth_rate_lme_multi <- lmer(
  log(Growth) ~ Treatment + (1 | MAT) + (1 | MAP) + (1 | elevation),
  data = dna_prod_test
)

summary(growth_rate_lme_multi)

# Single random effect: MAT
growth_rate_lme_mat <- lmer(
  log(Growth) ~ Treatment + (1 | MAT),
  data = dna_prod_test
)

summary(growth_rate_lme_mat)

# Compare models
anova(growth_rate_lme_mat, growth_rate_lme_multi)
AIC(growth_rate_lme_mat)
AIC(growth_rate_lme_multi)

##### Fit model with Treatment × Site interaction
growth_rate_lmer_interaction <- lmer(
  log(Growth) ~ Treatment * Site + (1 | Site),
  data = dna_prod_test
)

summary(growth_rate_lmer_interaction)
Anova(growth_rate_lmer_interaction, type = "III")

##### Model with degradation_index as fixed effect
growth_rate_lmer_degradation <- lmer(
  log(Growth) ~ degradation_index + (1 | Site),
  data = dna_prod
)

summary(growth_rate_lmer_degradation)
Anova(growth_rate_lmer_degradation, type = "III")

##### Calculate R² for degradation_index model
m1 <- growth_rate_lmer_degradation

# Extract variance components
var_comp <- as.data.frame(VarCorr(m1))
var_random <- var_comp$vcov[var_comp$grp == "Site"]
var_resid <- attr(VarCorr(m1), "sc")^2
var_total <- var_random + var_resid

# Variance explained by fixed effects
y_hat_fixed <- predict(m1, re.form = NA)
var_fixed <- var(y_hat_fixed)

# Marginal and conditional R²
R2_marginal <- var_fixed / var_total
R2_conditional <- (var_fixed + var_random) / var_total

cat("Marginal R² (fixed effects only):", round(R2_marginal, 3), "\n")
cat("Conditional R² (fixed + random effects):", round(R2_conditional, 3), "\n")

##### Correlation matrix among numeric variables
cor_data <- dna_prod_test %>%
  select(Growth, Moisture, pH, Moss_cover, FTIR_axis_1, O2, TC, TN, degradation_index)

cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Heatmap of correlations
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

##### Plot: Growth vs Ecosystem Health Index by Treatment
dna_prod$Treatment <- factor(dna_prod$Treatment, levels = c("DAM", "REST", "NAT"))

ggplot(dna_prod, aes(x = degradation_index, y = Growth, colour = Treatment)) +
  geom_point(size = 4, alpha = 0.7) +  # scatter points
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # regression line
  scale_y_continuous(
    trans = 'log2',
    labels = scales::label_number(accuracy = 0.1),
    limits = c(0.05, 11)
  ) +
  scale_colour_manual(
    values = c("NAT" = "#4daf4a", "REST" = "#377eb8", "DAM" = "#e41a1c"),
    labels = c("NAT" = "Natural", "REST" = "Restored", "DAM" = "Degraded")
  ) +
  ylab(expression("Microbial growth rate (ng g"^{-1}*" h"^{-1}*")")) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14) +
  theme(axis.text = element_text(size = 13, colour = "black"))

# Save the figure
ggsave("Figure_2c_MicrobialGrowthRate_EcosystemHealth.png", dpi = 1000, width = 6, height = 4.5, units = "in")
