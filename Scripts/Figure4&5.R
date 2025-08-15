############ Figure 4 & 5 ################


## load packages 
library(tidyverse)
library(vegan)
library(ggpubr)
library(ggcorrplot)
library(scales)
library(ggplot2)
library(emmeans)


#load raw data KEGG annotations: configs based gene abundances 
kegg_abund_raw<- read.csv("kegg_abund_raw_t.csv")
Sample_id<- read.csv("Sample_ID.csv")
env_data<- read.csv("env_data.csv")

# load prodigal output of total protein coding genes to normalise data 
prodigal_output_kegg_abund<- read.csv("prodigal_output_kegg.csv")

kegg_abund<- kegg_abund_raw %>%
  left_join(x=.,y=prodigal_output_kegg_abund, by=c("Sample")) 

kegg_abund_gather<- kegg_abund %>%
  gather(key = "KO", value = "reads", 
         K00003:K24042 )

kegg_abund_gather$reads<- as.numeric(kegg_abund_gather$reads)

#normalise: gene abundances as reads divided by total per sample prodigal output multiplied by mean all sample prodigal output
kegg_abund_gather<- kegg_abund_gather %>%
  mutate(reads_norm = (reads/prodigal_output)*57474057.89)

kegg_abund_gather<- kegg_abund_gather  %>%
  left_join(x=.,y=Sample_id, by=c("Sample"))

# pivot data by site (for faster execution)

kegg_abund_mig<- kegg_abund_gather %>%
  filter(Site == "Migneint") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)


kegg_abund_moor<- kegg_abund_gather %>%
  filter(Site == "Moors_House") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)

kegg_abund_cro<- kegg_abund_gather %>%
  filter(Site == "Crocach") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)

kegg_abund_bal<- kegg_abund_gather %>%
  filter(Site == "Balmoral") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)

kegg_abund_bow<- kegg_abund_gather %>%
  filter(Site == "Bowness") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)

kegg_abund_lang<- kegg_abund_gather %>%
  filter(Site == "Langwell" ) %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)

kegg_abund_stean<- kegg_abund_gather %>%
  filter(Site == "Stean") %>%
  select(Sample,KO,reads_norm) %>%
  pivot_wider(names_from = KO, values_from = reads_norm)


# dataframe to select genes of interest 
kegg_abund <- rbind(kegg_abund_bal, kegg_abund_bow,kegg_abund_cro,kegg_abund_mig,kegg_abund_lang, kegg_abund_moor, kegg_abund_stean)

kegg_abund <- kegg_abund  %>%
  left_join(x=.,y=Sample_id, by=c("Sample"))


#load raw data CAZy annotations: configs based gene abundances 
cazy_abund_raw<- read.csv("cazy_abund_raw.csv")

rownames(cazy_abund_raw) <- cazy_abund_raw[,1]

cazy_abund_raw<- cazy_abund_raw %>%
  select(-c("Sample"))

cazy_abund_raw_t<- t(cazy_abund_raw)

cazy_abund <- cbind(Sample = rownames(cazy_abund_raw_t), cazy_abund_raw_t)

cazy_abund<- as.data.frame(cazy_abund, stringsAsFactors=FALSE)

prodigal_output_caz_abund<- read.csv("prodigal_output_caz_abund.csv")

cazy_abund<- cazy_abund %>%
  left_join(x=.,y=prodigal_output_caz_abund, by=c("Sample")) 

cazy_abund_gather<- cazy_abund %>%
  gather(key = "cazyme", value = "reads", 
         AA1_1:SLH )

cazy_abund_gather$reads<- as.numeric(cazy_abund_gather$reads)

cazy_abund_gather<- cazy_abund_gather %>%
  mutate(reads_norm = (reads/prodigal_output)*57474057.89)    

cazy_abund_gather<- cazy_abund_gather  %>%
  left_join(x=.,y=Sample_id, by=c("Sample"))

### CAZy genes by substrate classes 

cazy_substrate<- read.csv("cazyme_by_substrate.csv")

# create dataframe in same format
cazy_abund_substrate<- cazy_abund_gather %>%
  left_join(x=., y=cazy_substrate, by=("cazyme"))

cazy_substrate_sum<- cazy_abund_substrate %>%
  dplyr::group_by(Sample, Treatment, Site,Substrate) %>%
  dplyr::summarise(reads_sum= sum(reads_norm, na.rm=TRUE))

cazy_substrate_sum<- cazy_substrate_sum %>%
  filter(Sample %in% c(  "MGr3I"  , "MGr3H",   "MGr3G"  , "MGr2F"  , "MGr2E"  , "MGr2D" ,  "MGr1C"  , "MGr1B" , 
                         "MGr1A" ,  "CRr3I" ,  "CRr3H" ,  "CRr3G" ,  "CRr2F"  , "CRr2E"  , "CRr2D"  , "CRr1C"  ,
                         "CRr1B" ,  "CRr1A"  , "MHr3I" ,  "MHr3H" ,  "MHr3G"  , "MHr2F"  , "MHr2E" ,  "MHr2D",  
                         "MHr1C" ,  "MHr1B" ,  "MHr1A",   "BOr3I" ,  "BOr3H"  , "BOr3G" ,  "BOr2F" ,  "BOr2E"  ,
                         "BOr2D" ,  "BOr1C"  , "BOr1B" ,  "BOr1A",   "BAr3I" ,  "BAr3H"  , "BAr3G" ,  "BAr2F",  
                         "BAr2E" ,  "BAr2D" ,  "BAr1C"  , "BAr1B" ,  "BAr1A" ,  "SEr3I"  , "SEr3H"  , "SEr3G",  
                         "SEr2F"  , "SEr2E"  , "SEr2D"  , "LASCr2F", "LASCr2E",
                         "LASCr2D", "LABRr1C" ,"LABRr1B" ,"LABRr1A" ,"LASAr3I" ,"LASAr3G","LASAr3H" ))

## Figure 4a ##
# Fermentation
kegg_abund_ferm<- kegg_abund %>%
  select( 
    "K00163", "K00627", "K00169" , "K03737", "K00174",  "K00656", "K01568" , "K14028", "K00114",
    "K00163" , "K00627" , "K00169" , "K03737" , "K00174" , "K00656", "K13788" , "K01512" , "K00925",
    "K00101",  "K00016",  "K00102" , "K03777",
    "K01595", "K00024", "K01676", "K00239",  "K00244",  "K18209" ,"K01899" ,  "K01902", "K01847", "K05606" ,"K11264",  "K03416",Sample
  )

####

# convert to as.numeric

kegg_abund_ferm<- kegg_abund_ferm %>% 
  mutate(across(1:28, as.numeric))
# sum 
kegg_abund_ferm_sum<- kegg_abund_ferm %>%
  left_join(x=.,y=Sample_id, by=c("Sample")) %>%
  mutate(ferm_sum = K00163+ K00627+ K00169 +K03737+ K00174+  K00656+ K01568 + K14028+ K00114+
           K00163 + K00627 + K00169 + K03737 + K00174 + K00656+ K13788 + K01512 + K00925+
           K00101+  K00016+  K00102 + K03777+
           K01595+ K00024+ K01676+ K00239+  K00244+  K18209 +K01899 +  K01902+ K01847+ K05606 +K11264+  K03416 ) %>%
  mutate(ethanol_ferm = K00163+ K00627+ K00169 +K03737+ K00174+  K00656+ K01568 + K14028+ K00114) %>%
  mutate(acetato = K00163 + K00627 + K00169 + K03737 + K00174 + K00656+ K13788 + K01512 + K00925) %>%
  mutate(lactate =  K00101+  K00016+  K00102 + K03777) %>%
  mutate(prop_ferm = K01595+ K00024+ K01676+ K00239+  K00244+  K18209 +K01899 +  K01902+ K01847+ K05606 +K11264+  K03416) 

kegg_abund_ferm_sum<- kegg_abund_ferm_sum %>%
  filter(Sample %in% c(  "MGr3I"  , "MGr3H",   "MGr3G"  , "MGr2F"  , "MGr2E"  , "MGr2D" ,  "MGr1C"  , "MGr1B" , 
                         "MGr1A" ,  "CRr3I" ,  "CRr3H" ,  "CRr3G" ,  "CRr2F"  , "CRr2E"  , "CRr2D"  , "CRr1C"  ,
                         "CRr1B" ,  "CRr1A"  , "MHr3I" ,  "MHr3H" ,  "MHr3G"  , "MHr2F"  , "MHr2E" ,  "MHr2D",  
                         "MHr1C" ,  "MHr1B" ,  "MHr1A",   "BOr3I" ,  "BOr3H"  , "BOr3G" ,  "BOr2F" ,  "BOr2E"  ,
                         "BOr2D" ,  "BOr1C"  , "BOr1B" ,  "BOr1A",   "BAr3I" ,  "BAr3H"  , "BAr3G" ,  "BAr2F",  
                         "BAr2E" ,  "BAr2D" ,  "BAr1C"  , "BAr1B" ,  "BAr1A" ,  "SEr3I"  , "SEr3H"  , "SEr3G",  
                         "SEr2F"  , "SEr2E"  , "SEr2D"  , "LASCr2F", "LASCr2E",
                         "LASCr2D", "LABRr1C" ,"LABRr1B" ,"LABRr1A" ,"LASAr3I" ,"LASAr3G", "LASAr3H"  ))


Figure_4a <- ggplot(kegg_abund_ferm_sum, aes(x=Treatment, y=ferm_sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=ferm_sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Fermentation", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000 , limits = c(44000, 71000)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_4a

ggsave("Figure_4a.png",  dpi=1000, width = 4, height = 4, units = "in")

#stats for figure 4a 
ferm_sum_lm <- lm(ferm_sum ~ Treatment, data = kegg_abund_ferm_sum)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(ferm_sum_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

ferm_sum_lm_interaction <- lm(ferm_sum ~ Treatment * Site, data = kegg_abund_ferm_sum)
anova(ferm_sum_lm_interaction)

merged_ferm <- merge(kegg_abund_ferm_sum, env_data, by = "Sample")

site_order <- c("Balmoral","Bowness","Stean","Moors_House","Langwell","Crocach","Migneint")
merged_ferm$Site.x <- factor(merged_ferm$Site.x, levels = site_order)
site_labels <- c(
  "Balmoral" = "Balmoral",
  "Bowness" = "Bowness",
  "Crocach" = "Crocach",
  "Langwell" = "Langwell",
  "Migneint" = "Migneint",
  "Moors_House" = "Moor House",  # <- rename here
  "Stean" = "Stean"
)

#correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_ferm %>%
  select(ferm_sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

ferm_sum_env <- lm(ferm_sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_ferm)
summary(ferm_sum_env)

ferm_sum_env <- lm(ferm_sum ~ moisture+ oxic, data = merged_ferm)
summary(ferm_sum_env)

ferm_sum_eco <- lm(ferm_sum ~ eco_index , data = merged_ferm)
summary(ferm_sum_eco)

Figure_4f <- ggplot(merged_ferm, aes(x = eco_index, y = ferm_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # Regression line with confidence interval
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Fermentation", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_4f

ggsave("Figure_4f.png",  dpi=1000, width = 6, height = 4.5, units = "in") 



## Figure S9a ##
merged_ferm <- merged_ferm %>% filter(Site.x != "Stean")

ggplot(merged_ferm, aes(x = eco_index, y = ferm_sum)) +
  geom_point(color = "black", size = 4, alpha = 0.7, 
             aes(shape = as.character(Treatment.x),
                 fill = as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, 
              linetype = "dashed", size = 0.7, alpha = 0.3) +
  scale_fill_manual(values = c(
    "#e41a1c", "#377eb8", "#4daf4a", 
    "#984ea3", "#ff7f00", "#ffff33", "#a65628")) +
  scale_shape_manual(name = "Treatment", values = c(21, 22, 23),
                     limits = c("NAT", "REST", "DAM"),
                     labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
  labs(y = "Fermentation") +
  scale_y_continuous(labels = function(x) x / 1000) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14) +
  theme(
    axis.text = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  facet_grid(. ~ Site.x, scales = "free", labeller = as_labeller(site_labels))

ggsave("Figure_S9a.png",  dpi=1000, width = 8, height = 3, units = "in") 


## Figure S10c ##
kegg_abund_ferm_sum <- merge(kegg_abund_ferm_sum,env_data,by="Sample")
lm_eco <- lm(ferm_sum ~ eco_index, data = kegg_abund_ferm_sum)
eco_summary <- summary(lm_eco)

# Extract model fit stats
r_squared <- round(eco_summary$r.squared, 3)
adj_r_squared <- round(eco_summary$adj.r.squared, 3)
f_stat <- eco_summary$fstatistic
model_p <- round(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE), 4)

# Coefficients table
eco_df <- as.data.frame(eco_summary$coefficients)
eco_df$Term <- rownames(eco_df)

# Add model stats as additional rows
stats_df <- data.frame(
  Estimate = NA,
  `Std. Error` = NA,
  `t value` = NA,
  `Pr(>|t|)` = NA,
  Term = c(
    paste0("R² = ", r_squared),
    paste0("Adjusted R² = ", adj_r_squared),
    paste0("Model p-value = ", model_p)
  )
)

# Combine the coefficient table and model statistics
names(stats_df) <- names(eco_df)
eco_results <- rbind(eco_df, stats_df)

# Coefficients table
eco_df <- as.data.frame(eco_summary$coefficients)
eco_df$Term <- rownames(eco_df)
rownames(eco_df) <- NULL

# Reorder columns
eco_df <- eco_df[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "Term")]

# Model stats as rows (ensure same column names and order)
stats_df <- data.frame(
  Estimate = NA,
  `Std. Error` = NA,
  `t value` = NA,
  `Pr(>|t|)` = NA,
  Term = c(
    paste0("R² = ", r_squared),
    paste0("Adjusted R² = ", adj_r_squared),
    paste0("Model p-value = ", model_p)
  ),
  check.names = FALSE
)

# Bind rows
eco_results <- rbind(eco_df, stats_df)

# Write to CSV
write.csv(eco_results, file = "eco_index_model_results.csv", row.names = FALSE)

# Write to CSV
write.csv(eco_results, file = "fer_eco_index_results.csv", row.names = FALSE)

#Plot ferm vs growth 

ggplot(merged_ferm, aes(y = Growth, x = ferm_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # Regression line with confidence interval
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  xlab(expression(atop("Fermentation", paste("(normalised reads ×1000)"))))+ 
  ylab(expression("Microbial growth rate (ng g"^{-1}*"h"^{-1}*")"))+ 
  scale_x_continuous(labels = function(x) x / 1000) +
  scale_y_continuous(trans='log2',labels = label_number(accuracy = 0.1), limits = c(0.05, 11))+
  theme(
    axis.text = element_text(size = 13, colour = "black"),
    legend.position = "None"
  ) +
  theme_light(base_size = 14) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) 
  
 
ggsave("Figure_S10c.png",  dpi=1000, width = 7, height = 4.5, units = "in") 

model <- lm(log(Growth) ~ ferm_sum, data = merged_ferm)
summary(model)

## Figure 5a ##
# lignin
cazy_substrate_Lignin <- cazy_substrate_sum %>%
  filter(Substrate == "Lignin") 

Figure_5a <- ggplot(cazy_substrate_Lignin, aes(x=Treatment, y=reads_sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=reads_sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Lignin", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                          "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000 , limits = c(10500, 17500)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_5a

ggsave("Figure_5a.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure 5f ##
Lignin_lm <- lm(reads_sum ~ Treatment, data = cazy_substrate_Lignin)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(Lignin_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

Lignin_lm_interaction <- lm(reads_sum ~ Treatment * Site, data = cazy_substrate_Lignin)
anova(Lignin_lm_interaction)

merged_Lignin <- merge(cazy_substrate_Lignin, env_data, by = "Sample")  %>%
  filter(Sample != "SEr2E")

# correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_Lignin %>%
  select(reads_sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Optional: Visualize correlations with a heatmap
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

Lignin_lm_env <- lm(reads_sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_Lignin)
summary(Lignin_lm_env)

Lignin_lm_env <- lm(reads_sum ~ Dim.1, data = merged_Lignin)
summary(Lignin_lm_env)

Lignin_lm_eco <- lm(reads_sum ~ eco_index, data = merged_Lignin)
summary(Lignin_lm_eco)

Figure_5f <- ggplot(merged_Lignin, aes(x = eco_index, y = reads_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # Regression line with confidence interval
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Lignin", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_5f

ggsave("Figure_5f.png",  dpi=1000, width = 6, height = 4.5, units = "in") 

## Figure 5b ##
# Cellulose

cazy_substrate_Cellulose <- cazy_substrate_sum %>%
  filter(Substrate == "Cellulose") 

Figure_5b <- ggplot(cazy_substrate_Cellulose, aes(x=Treatment, y=reads_sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=reads_sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Cellulose", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                          "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000, 
                     limits = c(32000, 56000), 
                     breaks = seq(35000, 55000, by = 5000)) +  # Adjust breaks and limits
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_5b

ggsave("Figure_5b.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure 5g ##
Cellulose_lm <- lm(reads_sum ~ Treatment, data = cazy_substrate_Cellulose)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(Cellulose_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

Cellulose_lm_interaction <- lm(reads_sum ~ Treatment * Site, data = cazy_substrate_Cellulose)
anova(Cellulose_lm_interaction)

merged_Cellulose <- merge(cazy_substrate_Cellulose, env_data, by = "Sample")  %>%
  filter(Sample != "SEr2E")

# correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_Cellulose %>%
  select(reads_sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Optional: Visualize correlations with a heatmap
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

Cellulose_lm_env <- lm(reads_sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_Cellulose)
summary(Cellulose_lm_env)

Cellulose_lm_env <- lm(reads_sum ~ Dim.1, data = merged_Cellulose)
summary(Cellulose_lm_env)

Cellulose_lm_eco <- lm(reads_sum ~ eco_index, data = merged_Cellulose)
summary(Cellulose_lm_eco)

Figure_5g <- ggplot(merged_Cellulose, aes(x = eco_index, y = reads_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # Regression line with confidence interval
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Cellulose", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_5g

ggsave("Figure_5g.png",  dpi=1000, width = 6, height = 4.5, units = "in") 

## Figure 5c ##
# oligosacharides 
cazy_substrate_Oligosaccharides <- cazy_substrate_sum %>%
  filter(Substrate == "Oligosaccharides"  ) 

Figure_5c <- ggplot(cazy_substrate_Oligosaccharides, aes(x=Treatment, y=reads_sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=reads_sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Oligosaccharides", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                          "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000, limits = c(2000, 13500)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_5c

ggsave("Figure_5c.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure 5h ##
Oligosaccharides_lm <- lm(reads_sum ~ Treatment, data = cazy_substrate_Oligosaccharides)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(Oligosaccharides_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

Oligosaccharides_lm_interaction <- lm(reads_sum ~ Treatment * Site, data = cazy_substrate_Oligosaccharides)
anova(Oligosaccharides_lm_interaction)

merged_Oligosaccharides <- merge(cazy_substrate_Oligosaccharides, env_data, by = "Sample")  %>%
  filter(Sample != "SEr2E")

#correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_Oligosaccharides %>%
  select(reads_sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

Oligosaccharides_lm_env <- lm(reads_sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_Oligosaccharides)
summary(Oligosaccharides_lm_env)

Oligosaccharides_lm_env <- lm(reads_sum ~ moisture+ sphag+ Dim.1, data = merged_Oligosaccharides)
summary(Oligosaccharides_lm_env)

Oligosaccharides_lm_eco <- lm(reads_sum ~ eco_index, data = merged_Oligosaccharides)
summary(Oligosaccharides_lm_eco)

Figure_5h <- ggplot(merged_Oligosaccharides, aes(x = eco_index, y = reads_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Oligosaccharides", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_5h

ggsave("Figure_5h.png",  dpi=1000, width = 6, height = 4.5, units = "in") 


## Figure 5d ##
# Aromatics compounds
ko_levels <- read.csv("KO_levels.csv") 

kegg_levels <- kegg_abund_gather %>%
  left_join(x=.,y=ko_levels, by=c("KO"))

kegg_levels_aromatic <- kegg_levels %>%
  filter(Category3 == "01220_Degradation_of_aromatic_compounds_[PATH:ko01220]")

#00190_Oxidative_phosphorylation_[PATH:ko00190]

kegg_levels_aromatic <- kegg_levels_aromatic  %>%
  dplyr::group_by(Sample, Site, Treatment) %>%
  dplyr::summarise(sum = sum(reads_norm)) 

Figure_5d <- ggplot(kegg_levels_aromatic, aes(x=Treatment, y=sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Aromatic compounds", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                          "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000, limits = c(18000, 52000)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_5d

ggsave("Figure_5d.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure 5i##
aromatic_lm <- lm(sum ~ Treatment, data = kegg_levels_aromatic)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(aromatic_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

aromatic_lm_interaction <- lm(sum ~ Treatment * Site, data = kegg_levels_aromatic)
anova(aromatic_lm_interaction)

merged_aromatic <- merge(kegg_levels_aromatic, env_data, by = "Sample") %>%
  filter(Sample != "SEr2E")

#correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_aromatic %>%
  select(sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Optional: Visualize correlations with a heatmap
#library(ggcorrplot)
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

aromatic_lm_env <- lm(sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_aromatic)
summary(aromatic_lm_env)

aromatic_lm_env <- lm(sum ~ Dim.1, data = merged_aromatic)
summary(aromatic_lm_env)

aromatic_lm_eco <- lm(sum ~ eco_index, data = merged_aromatic)
summary(aromatic_lm_eco)

Figure_5i <- ggplot(merged_aromatic, aes(x = eco_index, y = sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Aromatic compounds", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_5i

ggsave("Figure_5i.png",  dpi=1000, width = 6, height = 4.5, units = "in")

## Figure 5e ##
# cell wall 

cazy_substrate_Cell_Wall<- cazy_substrate_sum %>%
  filter(Substrate == "Cell Wall"  ) 
 
Figure_5e <- ggplot(cazy_substrate_Cell_Wall, aes(x=Treatment, y=reads_sum))+
  geom_boxplot(lwd=1,outlier.shape= NA)+
  geom_jitter(aes(y=reads_sum, x=Treatment, color=Site), alpha=0.7, width = 0.3,size=3)+
  ylab(expression(atop("Microbial cell wall", paste("(normalised reads ×1000)"))))+ ### [reported] = subscript 
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
  scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                           "REST" = "Restored"))+
  scale_y_continuous(labels = function(x) x / 1000, limits = c(1600, 5400)) +
  theme(legend.text=element_text(size=16))+
  theme(legend.position = "None")+
  theme(axis.text=element_text(size=13,colour="black")) 

Figure_5e

ggsave("Figure_5e.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure 5j##
cell_wall_lm <- lm(reads_sum ~ Treatment, data = cazy_substrate_Cell_Wall)

# Perform pairwise comparisons for the Treatment variable
emmeans_results <- emmeans(cell_wall_lm, ~ Treatment)

# Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

# Print the pairwise comparisons
summary(pairwise_comparisons)

cell_wall_lm_interaction <- lm(reads_sum ~ Treatment * Site, data = cazy_substrate_Cell_Wall)
anova(cell_wall_lm_interaction)

merged_cell_wall <- merge(cazy_substrate_Cell_Wall, env_data, by = "Sample") %>%
  filter(Sample != "SEr2E")

#correaltion 
# Select relevant columns (numeric factors)
cor_data <- merged_cell_wall %>%
  select(reads_sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

# Compute Pearson correlation matrix
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Optional: Visualize correlations with a heatmap
#library(ggcorrplot)
ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

cell_wall_lm_env <- lm(reads_sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_cell_wall)
summary(cell_wall_lm_env)

cell_wall_lm_env <- lm(reads_sum ~ Dim.1, data = merged_cell_wall)
summary(cell_wall_lm_env)

cell_wall_lm_eco <- lm(reads_sum ~ eco_index, data = merged_cell_wall)
summary(cell_wall_lm_eco)

Figure_5j <- ggplot(merged_cell_wall, aes(x = eco_index, y = reads_sum)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                       fill=as.character(Site.x))) +  
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  
  scale_fill_manual(values=c( "#e41a1c",
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
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  ylab(expression(atop("Microbial cell wall", paste("(normalised reads ×1000)"))))+ 
  scale_y_continuous(labels = function(x) x / 1000) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_5j

ggsave("Figure_5j.png",  dpi=1000, width = 6, height = 4.5, units = "in") 


# Figure 4b-e & 4g-j & Figure S11
# Define your list of KEGG pathways and labels
pathways <- list(
  list("01220_Degradation_of_aromatic_compounds_[PATH:ko01220]", "Degradation of Aromatic Compounds", "Aromatic compounds", "Aromatic"),
  list("00680_Methane_metabolism_[PATH:ko00680]", "Methane metabolism", "Methane metabolism", "Methane"),
  list("00720_Carbon_fixation_pathways_in_prokaryotes_[PATH:ko00720]", "Carbon fixation pathways in prokaryotes", "Carbon fixation in prokaryotes", "CarbonFixation"),
  list("00910_Nitrogen_metabolism_[PATH:ko00910]", "Nitrogen metabolism", "Nitrogen metabolism", "Nitrogen"),
  list("00920_Sulfur_metabolism_[PATH:ko00920]", "Sulfur metabolism", "Sulfur metabolism", "Sulfur")
)

# Loop through each pathway
for (pathway in pathways) {
  kegg_id <- pathway[[1]]
  title <- pathway[[2]]
  x_axis_label <- pathway[[3]]
  filename_prefix <- pathway[[4]]
  
  # Filter and summarise
  df_pathway <- kegg_levels %>%
    filter(Category3 == kegg_id) %>%
    group_by(Sample, Site, Treatment) %>%
    summarise(sum = sum(reads_norm), .groups = "drop") %>%
    filter(Sample %in% c("MGr3I", "MGr3H", "MGr3G", "MGr2F", "MGr2E", "MGr2D", "MGr1C", "MGr1B", 
                         "MGr1A", "CRr3I", "CRr3H", "CRr3G", "CRr2F", "CRr2E", "CRr2D", "CRr1C", 
                         "CRr1B", "CRr1A", "MHr3I", "MHr3H", "MHr3G", "MHr2F", "MHr2E", "MHr2D", 
                         "MHr1C", "MHr1B", "MHr1A", "BOr3I", "BOr3H", "BOr3G", "BOr2F", "BOr2E", 
                         "BOr2D", "BOr1C", "BOr1B", "BOr1A", "BAr3I", "BAr3H", "BAr3G", "BAr2F", 
                         "BAr2E", "BAr2D", "BAr1C", "BAr1B", "BAr1A", 
                         "LASCr2F", "LASCr2E", "LASCr2D", 
                         "LABRr1C", "LABRr1B", "LABRr1A", "LASAr3I", "LASAr3G", "LASAr3H"))
  
  # Merge with environmental data
  merged_pathway <- merge(df_pathway, env_data, by = "Sample")
  
  # ---- Linear model ----
  model <- lm(log(Growth) ~ sum, data = merged_pathway)
  cat("\nPathway:", title, "\n")
  print(summary(model))
  
  # Reorder Site.x
  site_order <- c("Balmoral", "Bowness", "Moors_House", "Langwell", "Crocach", "Migneint")
  merged_pathway$Site.x <- factor(merged_pathway$Site.x, levels = (site_order))

 
  Figure_5_substrate_plot1 <- ggplot(merged_pathway, aes(x=Treatment.x, y=sum))+
     geom_boxplot(lwd=1,outlier.shape= NA)+
     geom_jitter(aes(y=sum, x=Treatment.x, color=Site.x), alpha=0.7, width = 0.3,size=3)+
     ylab(bquote(atop(.(title), "(normalised reads ×1000)")))+
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
    scale_x_discrete(limits=c("DAM", "REST","NAT"),labels=c("DAM" = "Degraded", "NAT" = "Natural",
                                                            "REST" = "Restored"))+
    scale_y_continuous(labels = function(x) x / 1000 , limits = c(floor(min(merged_pathway$sum) / 1000) * 1000, ceiling(max(merged_pathway$sum) / 1000) * 1000)) +
    theme(legend.text=element_text(size=16))+
    theme(legend.position = "None")+
    theme(axis.text=element_text(size=13,colour="black"))
 
  Figure_5_substrate_plot1

  ggsave(paste0("Figure_5_",filename_prefix, "_plot1.png"),  dpi=1000, width = 4, height = 4, units = "in")

  # Figure 5f ##
  substrate_lm <- lm(sum ~ Treatment.x, data = merged_pathway)

  # Perform pairwise comparisons for the Treatment variable
  emmeans_results <- emmeans(substrate_lm, ~ Treatment.x)

  # Perform pairwise comparisons for DAM vs. REST, REST vs. NAT, DAM vs. NAT
  pairwise_comparisons <- pairs(emmeans_results, adjust = "tukey")

  # Print the pairwise comparisons
  summary(pairwise_comparisons)

  substrate_lm_interaction <- lm(sum ~ Treatment.x * Site.x, data = merged_pathway)
  anova(substrate_lm_interaction)

  # correaltion
  # Select relevant columns (numeric factors)
  cor_data <- merged_pathway %>%
    select(sum, moisture, ph, sphag, Dim.1, oxic, CN, eco_index)

  # Compute Pearson correlation matrix
  cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

  # Print correlation matrix
  print(cor_matrix)

  # Optional: Visualize correlations with a heatmap
  ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

  substrate_lm_env <- lm(sum ~ moisture+ ph+ sphag+ Dim.1+ oxic+ CN, data = merged_pathway)
  summary(substrate_lm_env)

  substrate_lm_env <- lm(sum ~ Dim.1, data = merged_pathway)
  summary(substrate_lm_env)

  substrate_lm_eco <- lm(sum ~ eco_index, data = merged_pathway)
  summary(substrate_lm_eco)

  Figure_5_substrate_plot2 <- ggplot(merged_pathway, aes(x = eco_index, y = sum)) +
    geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(Treatment.x),
                                                         fill=as.character(Site.x))) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +  # Regression line with confidence interval
    scale_fill_manual(values=c( "#e41a1c",
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
    scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                       labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
    ylab(bquote(atop(.(title), "(normalised reads ×1000)")))+
    scale_y_continuous(labels = function(x) x / 1000) +

    guides(
      fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(order = 2, title = "Treatment") ) +
    labs(x = "Ecosystem health index") +
    theme_light(base_size = 14)+
    theme(axis.text = element_text(size = 13, colour = "black"))

  Figure_5_substrate_plot2

  ggsave(paste0("Figure_5_",filename_prefix, "_plot2.png"),  dpi=1000, width = 6, height = 4, units = "in")

  # Correlation with ecosystem health index
  p2 <- ggplot(merged_pathway, aes(x = sum, y = Growth)) +
    geom_point(color = "black", size = 4, alpha = 0.7, aes(shape = as.character(Treatment.x), fill = as.character(Site.x))) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +
    scale_y_continuous(trans='log2',labels = label_number(accuracy = 0.1), limits = c(0.05, 11))+
    scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")) +
    scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                       labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    xlab(paste0(x_axis_label))  +
    scale_x_continuous(labels = function(x) x / 1000) +
    ylab(expression("Microbial growth rate (ng g"^{-1}*"h"^{-1}*")"))+ 
    guides(fill = "none", shape = "none", color = "none") +
    theme_light(base_size = 14) +
    theme(
      axis.text = element_text(size = 13, colour = "black"),
      legend.position = "none"
    ) +
    facet_grid(. ~ Site.x, scales = "free")
  
  ggsave(paste0(filename_prefix, "_growth_facet_site.png"), plot = p2, dpi = 1000, width = 8, height = 3, units = "in")
}  



#LOOP for CAZy 

library(dplyr)
library(ggplot2)
library(emmeans)
library(ggcorrplot)

# Define CAZy substrates: name used in data, y-axis label, and filename prefix
substrates <- list(
  list("Cell Wall", "Microbial cell wall", "CellWall"),
  list("Lignin", "Lignin", "Lignin"),
  list("Oligosaccharides", "Oligosaccharides", "Oligosaccharides"),
  list("Cellulose", "Cellulose", "Cellulose"),
  list("Hemicellulose", "Hemicellulose", "Hemicellulose")
)

# Loop through each CAZy substrate
for (entry in substrates) {
  substrate_name <- entry[[1]]
  y_label <- entry[[2]]
  filename_prefix <- entry[[3]]
  
  # Filter data
  df_sub <- cazy_substrate_sum %>%
    group_by(Sample, Site, Treatment) %>%
    filter(Substrate == substrate_name) %>%
    filter(Sample %in% c("MGr3I", "MGr3H", "MGr3G", "MGr2F", "MGr2E", "MGr2D", "MGr1C", "MGr1B", 
                       "MGr1A", "CRr3I", "CRr3H", "CRr3G", "CRr2F", "CRr2E", "CRr2D", "CRr1C", 
                       "CRr1B", "CRr1A", "MHr3I", "MHr3H", "MHr3G", "MHr2F", "MHr2E", "MHr2D", 
                       "MHr1C", "MHr1B", "MHr1A", "BOr3I", "BOr3H", "BOr3G", "BOr2F", "BOr2E", 
                       "BOr2D", "BOr1C", "BOr1B", "BOr1A", "BAr3I", "BAr3H", "BAr3G", "BAr2F", 
                       "BAr2E", "BAr2D", "BAr1C", "BAr1B", "BAr1A", 
                       "LASCr2F", "LASCr2E", "LASCr2D", 
                       "LABRr1C", "LABRr1B", "LABRr1A", "LASAr3I", "LASAr3G", "LASAr3H"))
  
  # Merge with environmental data
  merged_df <- merge(df_sub, env_data, by = "Sample") 

  # Plot
  p <- ggplot(df_sub, aes(x = Treatment, y = reads_sum)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    geom_jitter(aes(color = Site), alpha = 0.7, width = 0.3, size = 3) +
    ylab(paste0(y_label))  +
    xlab("") +
    theme_light(base_size = 14) +
    theme(legend.position = "none", axis.text = element_text(size = 13, colour = "black")) +
    scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
    scale_x_discrete(
      limits = c("DAM", "REST", "NAT"),
      labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")
    ) +
    scale_y_continuous(labels = \(x) x / 1000)
  
  ggsave(paste0("Figure_", filename_prefix, ".png"), plot = p, dpi = 1000, width = 4, height = 4, units = "in")
  
  # Correlation with ecosystem health index
  p2 <- ggplot(merged_df, aes(x = eco_index, y = reads_sum)) +
    geom_point(color = "black", size = 4, alpha = 0.7, aes(shape = as.character(Treatment.x), fill = as.character(Site.x))) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +
    scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
    scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                       labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    ylab(paste0(y_label))  +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(x = "Ecosystem health index") +
    guides(fill = "none", shape = "none", color = "none") +
    theme(legend.position = "none") +
    theme_light(base_size = 14) +
    theme(axis.text = element_text(size = 13, colour = "black"))
  
  ggsave(paste0(filename_prefix, "_eco.png"), plot = p2, dpi = 1000, width = 4, height = 4, units = "in")
  # Models and results
  
  # Define output file name
  output_file <- paste0("results_", filename_prefix, ".csv")
  
  # Helper function to write each section
  write_section <- function(header, df, file, append = TRUE) {
    write(paste0("### ", header, " ###"), file = file, append = append)
    write.table(df, file = file, sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE)
    write("\n\n", file = file, append = TRUE)
  }
  
  # 1. Pairwise comparisons
  # Pairwise contrasts using emmeans
  library(emmeans)
  emm <- emmeans(lm(reads_sum ~ Treatment.x * Site.x, data = merged_df), ~ Treatment.x)
  pairwise_results <- summary(pairs(emm))
  pairwise_df <- as.data.frame(pairwise_results)
  write_section("Pairwise comparisons (Tukey-adjusted)", pairwise_df, output_file, append = FALSE)
  
  # 2. ANOVA: Treatment × Site
  # ANOVA interaction analysis
  anova_model <- lm(reads_sum ~ Treatment.x * Site.x, data = merged_df)
  anova_df <- as.data.frame(anova(anova_model))
  anova_df$Effect <- rownames(anova_df)
  write_section("ANOVA: Treatment × Site interaction", anova_df, output_file)
  
  # 3. Environmental model: sum ~ moisture + ph + sphag + Dim.1 + oxic + CN
  lm_env <- lm(reads_sum ~ moisture + ph + sphag + Dim.1 + oxic + CN, data = merged_df)
  lm_env_df <- as.data.frame(summary(lm_env)$coefficients)
  lm_env_df$Term <- rownames(lm_env_df)
  write_section("Linear model: Environmental variables", lm_env_df, output_file)
  
  # 4. Ecological index model: sum ~ eco_index
  lm_eco <- lm(reads_sum ~ eco_index, data = merged_df)
  eco_summary <- summary(lm_eco)
  
  # Extract model fit stats
  r_squared <- round(eco_summary$r.squared, 3)
  adj_r_squared <- round(eco_summary$adj.r.squared, 3)
  f_stat <- eco_summary$fstatistic
  model_p <- round(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE), 4)
  
  # Coefficients table
  eco_df <- as.data.frame(eco_summary$coefficients)
  eco_df$Term <- rownames(eco_df)
  
  # Add model stats as additional rows
  stats_df <- data.frame(
    Estimate = NA,
    `Std. Error` = NA,
    `t value` = NA,
    `Pr(>|t|)` = NA,
    Term = c(
      paste0("R² = ", r_squared),
      paste0("Adjusted R² = ", adj_r_squared),
      paste0("Model p-value = ", model_p)
    )
  )
  
  # Combine and write
  eco_df <- bind_rows(eco_df, stats_df)
  write_section("Linear model: Ecosystem health index", eco_df, output_file)
  
}


# Define CAZy substrates: name used in data, y-axis label, and filename prefix
substrates <- list(
  list("Cell Wall", "Microbial cell wall", "CellWall"),
  list("Lignin", "Lignin", "Lignin"),
  list("Oligosaccharides", "Oligosaccharides", "Oligosaccharides"),
  list("Cellulose", "Cellulose", "Cellulose"),
  list("Hemicellulose", "Hemicellulose", "Hemicellulose")
)

# Loop through each CAZy substrate
for (entry in substrates) {
  substrate_name <- entry[[1]]
  y_label <- entry[[2]]
  filename_prefix <- entry[[3]]
  
  # Filter data
  df_sub <- cazy_substrate_sum %>%
    group_by(Sample, Site, Treatment) %>%
    filter(Substrate == substrate_name) %>%
    filter(Sample %in% c("MGr3I", "MGr3H", "MGr3G", "MGr2F", "MGr2E", "MGr2D", "MGr1C", "MGr1B", 
                         "MGr1A", "CRr3I", "CRr3H", "CRr3G", "CRr2F", "CRr2E", "CRr2D", "CRr1C", 
                         "CRr1B", "CRr1A", "MHr3I", "MHr3H", "MHr3G", "MHr2F", "MHr2E", "MHr2D", 
                         "MHr1C", "MHr1B", "MHr1A", "BOr3I", "BOr3H", "BOr3G", "BOr2F", "BOr2E", 
                         "BOr2D", "BOr1C", "BOr1B", "BOr1A", "BAr3I", "BAr3H", "BAr3G", "BAr2F", 
                         "BAr2E", "BAr2D", "BAr1C", "BAr1B", "BAr1A", 
                         "LASCr2F", "LASCr2E", "LASCr2D", 
                         "LABRr1C", "LABRr1B", "LABRr1A", "LASAr3I", "LASAr3G", "LASAr3H"))
  
  # Merge with environmental data
  merged_df <- merge(df_sub, env_data, by = "Sample") 
  
  # Reorder Site.x
  site_order <- c("Balmoral", "Bowness", "Moors_House", "Langwell", "Crocach", "Migneint")
  merged_df$Site.x <- factor(merged_df$Site.x, levels = (site_order))
  
  # Correlation with ecosystem health index
  p2 <- ggplot(merged_df, aes(x = eco_index, y = reads_sum)) +
    geom_point(color = "black", size = 4, alpha = 0.7, aes(shape = as.character(Treatment.x), fill = as.character(Site.x))) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7, alpha = 0.3) +
    scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")) +
    scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                       labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    ylab(paste0(y_label))  +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(x = "Ecosystem health index") +
    guides(fill = "none", shape = "none", color = "none") +
    theme_light(base_size = 14) +
    theme(
      axis.text = element_text(size = 13, colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    facet_grid(. ~ Site.x, scales = "free")
  
  ggsave(paste0(filename_prefix, "_eco.png"), plot = p2, dpi = 1000, width = 8, height = 3, units = "in")
}

