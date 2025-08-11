#############  Figure 6 ###################


# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(ggsci)
library(patchwork)
library(gridExtra)

# Read MAG metadata & HMM results
mag_metadata <- read_csv("mag_metabolism.csv")
hmm_hits <- read_csv("mags_HMMhit.csv")

### HMM ###
# Reshape HMM results to long format
hmm_long <- melt(hmm_hits, id.vars = colnames(hmm_hits)[1:10])

# Filter where all values in group are "Present"
hmm_present <- hmm_long %>%
  group_by(Function, variable) %>%
  filter(all(value == "Present")) %>%
  ungroup() %>%
  dplyr::select(Category, Function, variable) %>%
  distinct()


# Select specific metabolic categories of interest
categories_of_interest <- c("Ethanol fermentation", "Fermentation", "Sulfur cycling enzymes (detailed)",
                            "Nitrogen cycling", "Methane metabolism", "Carbon fixation", "Sulfur cycling")

hmm_present$variable <- gsub(".Hmm.presence", "", hmm_present$variable)

# Filter to selected categories
hmm_selected <- hmm_present %>%
  filter(Category %in% categories_of_interest)

# Merge with metadata
mag_taxa <- unique(mag_metadata[, 1:12])
hmm_selected <- merge(hmm_selected, mag_taxa, by.x = "variable", by.y = "Host")

# Clean Host Trend Group
hmm_selected_clean <- hmm_selected %>%
  filter(!is.na(`Host Trend Group`) & `Host Trend Group` != "NA") %>%
  mutate(`Host Trend Group` = recode(`Host Trend Group`,
                                     "Natural-abundant" = "NAT",
                                     "Restored-abundant" = "REST",
                                     "Damaged-abundant" = "DAM"),
         `Host Trend Group` = factor(`Host Trend Group`, levels = c("NAT", "REST", "DAM")),
         Host_Abundance = as.numeric(Host_Abundance))

# Count unique hosts per category & group
hmm_unique_hosts_by_category <- hmm_selected_clean %>%
  group_by(Category, `Host Trend Group`) %>%
  summarize(unique_hosts = n_distinct(variable), .groups = "drop")

# Count unique hosts per function & group
hmm_unique_hosts <- hmm_selected_clean %>%
  group_by(Function, `Host Trend Group`, Category) %>%
  summarize(unique_hosts = n_distinct(variable), .groups = "drop")

# Prepare fermentation data & add alcohol utilization placeholder
fermentation_data <- hmm_unique_hosts %>%
  filter(Category %in% c("Ethanol fermentation", "Fermentation")) %>%
  mutate(Category = "Fermentation")

# Summarize unique_hosts per Host Trend Group
fermentation_summary <- fermentation_data %>%
  group_by(`Host Trend Group`) %>%
  summarise(total_unique_hosts = sum(unique_hosts), .groups = "drop")

# Fermentation bar plot
fermentation_plot_overall <- ggplot(hmm_unique_hosts_by_category[hmm_unique_hosts_by_category$Category == "Fermentation",], 
                           aes(x = `Host Trend Group`, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "", y = "MAGs with metabolic pathway") +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("fermentation_plot_overall.png", plot = fermentation_plot_overall, dpi = 300, width = 3, height = 1.5)

# Plot fermentation data
fermentation_data$Category <- "Fermentation\n"

# Fermentation scatter plot
fermentation_plot <- ggplot(fermentation_data, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 4, color = "black", alpha = 0.7, shape=21) +  # 设置点的边缘为黑色，填充颜色根据`fill`设置
  scale_shape_manual(values = c(21, 22, 23)) + 
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(x = NULL, y = 'MAGs with metabolic pathway', fill = 'Cluster') +
  coord_flip() +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  theme_light(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
    #axis.title.x = element_blank()
  )+
  scale_y_continuous(limits = c(0, 190), breaks = seq(0, 180, by = 30))

# Save fermentation plot
# ggsave("fermentation_plot.png", plot = fermentation_plot, dpi = 300, width = 10, height = 3)

### KEGG ###

# Load KEGG data
mags_kegg <- read_csv("mags_keggmodulehit.csv")

# Reshape long format
kegg_long <- melt(mags_kegg, id.vars = colnames(mags_kegg)[1:3])
kegg_long$variable <- gsub(".Module.presence", "", kegg_long$variable)

kegg_long_catagory <- kegg_long %>%
  filter(value == "Present")

# Filter to interested categories
categories_of_interest <- c("Carbon fixation", "Nitrogen metabolism", 
                            "Methane metabolism", "Sulfur metabolism")

# Add host metadata
taxa_info <- unique(df[, c(1, 6)])  # assume `df` defined earlier
kegg_summary <- merge(kegg_long_catagory, taxa_info, 
                      by.x = "variable", by.y = "Host", all.x = TRUE)

# Clean and recode Host Trend Group
kegg_summary <- kegg_summary %>%
  filter(!is.na(`Host Trend Group`) & `Host Trend Group` != "NA") %>%
  mutate(`Host Trend Group` = recode(`Host Trend Group`,
                                     "Natural-abundant" = "NAT",
                                     "Restored-abundant" = "REST",
                                     "Damaged-abundant" = "DAM"),
         `Host Trend Group` = factor(`Host Trend Group`, levels = c("NAT", "REST", "DAM")))


# Count unique hosts per category & group
keggresult_by_category <- kegg_summary %>%
  group_by(`Host Trend Group`, `Module.Category`) %>%
  summarise(unique_hosts = n_distinct(variable), .groups = "drop")

# Count unique hosts per function & group

kegg_long_function <- kegg_long %>%
  mutate(presence = ifelse(value == "Present", 1, 0))

# Add host metadata
taxa_info <- unique(df[, c(1, 6)])  # assume `df` defined earlier
kegg_summary_function <- merge(kegg_long_function, taxa_info, 
                      by.x = "variable", by.y = "Host", all.x = TRUE)

# Clean and recode Host Trend Group
kegg_summary_function <- kegg_summary_function %>%
  filter(!is.na(`Host Trend Group`) & `Host Trend Group` != "NA") %>%
  mutate(`Host Trend Group` = recode(`Host Trend Group`,
                                     "Natural-abundant" = "NAT",
                                     "Restored-abundant" = "REST",
                                     "Damaged-abundant" = "DAM"),
         `Host Trend Group` = factor(`Host Trend Group`, levels = c("NAT", "REST", "DAM")))

# Aggregate: if any "Present" in the same group => 1
keggresult <- kegg_summary_function %>%
  group_by(`Module`, `Host Trend Group`, `Module.Category`) %>%
  summarise(unique_hosts = sum(presence), .groups = "drop")  %>%
  group_by(`Module`) %>%
  filter(sum(unique_hosts) > 0) %>%
  ungroup()


# Rename columns for consistency with downstream (if needed)
colnames(keggresult) <- c( "Function", "Host Trend Group", "Category","unique_hosts")

keggresult <- keggresult[keggresult$Category %in% categories_of_interest,]
keggresult$Category <- factor(keggresult$Category,
                              levels = c("Methane metabolism", "Nitrogen metabolism", 
                                         "Sulfur metabolism", "Carbon fixation", "Fermentation"),
                              labels = c("Methane\nmetabolism", "Nitrogen\nmetabolism", 
                                         "Sulfur\nmetabolism", "Carbon\nfixation", "Fermentation"))

# global presence data
genecount <- read.csv("genecountinbac&arc_kegg.csv")
genecount_long <- melt(genecount)

colnames(genecount_long) <- c("Category","Function","variable","genecount")

########## Methane ###############

# Subset data where Category is "Methane metabolism"
methane_result <- keggresult[keggresult$Category == "Methane\nmetabolism",]

# Methane scatter plot
methane_plot1 <- ggplot(methane_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 4, color = "black", alpha = 0.7, shape = 21) +  # set black border and fill by `fill`
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(x = NULL, y = 'MAGs with metabolic pathway', fill = 'Cluster') +
  coord_flip() +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  theme_light(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 50, by = 10))

# Get the order of Function in methane_result
methane_function_order <- unique(methane_result$Function)

# Reorder Function in genecount_long to match methane_result
methane_genecount_long_sorted <- genecount_long[genecount_long$Category == "Methane metabolism",]
methane_genecount_long_sorted$Function <- factor(methane_genecount_long_sorted$Function,
                                                 levels = methane_function_order)
# Global methane presence
methane_plot2 <- ggplot(methane_genecount_long_sorted[!is.na(methane_genecount_long_sorted$Function),], 
                        aes(x = Function, y = variable)) +
  geom_tile(aes(fill = log10(genecount + 1)), width = 1, height = 1, color = "gray") +  # heatmap
  geom_text(aes(label = round(genecount, 2)), color = "black", size = 3) +  # add text labels
  scale_fill_gradient(low = "white", high = "lightgrey") +  # gradient fill
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # remove y-axis text
    axis.text.x = element_blank(),  # remove x-axis text
    strip.background = element_blank(),
    panel.border = element_blank(),  # remove border
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_blank(),  # remove y-axis title
    panel.background = element_blank(),  # remove panel background
    legend.position = "none",  # remove legend
    panel.grid = element_blank()  # remove grid lines
  ) +   
  coord_flip()  # flip axes

# Convert heatmap to grob
methane_grob_plot2 <- ggplotGrob(methane_plot2)

# Methane bar plot
methane_plot3 <- ggplot(keggresult_by_category[keggresult_by_category$Module.Category %in% "Methane metabolism",], 
                        aes(x = `Host Trend Group`, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "", y = "MAGs with metabolic pathway") +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("methane_plot_overall.png", plot = methane_plot3, dpi = 300, width = 3, height = 1.5)


########## Nitrogen ###############

# Subset data where Category is "Nitrogen metabolism"
nitrogen_result <- keggresult[keggresult$Category == "Nitrogen\nmetabolism",]

# Nitrogen scatter plot
nitrogen_plot1 <- ggplot(nitrogen_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 4, color = "black", alpha = 0.7, shape = 21) +  # set black border and fill by `fill`
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(x = NULL, y = 'MAGs with metabolic pathway', fill = 'Cluster') +
  coord_flip() +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  theme_light(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 50, by = 10))

# Get the order of Function in nitrogen_result
nitrogen_function_order <- unique(nitrogen_result$Function)

# Reorder Function in genecount_long to match nitrogen_result
nitrogen_genecount_long_sorted <- genecount_long[genecount_long$Category == "Nitrogen metabolism",]
nitrogen_genecount_long_sorted$Function <- factor(nitrogen_genecount_long_sorted$Function,
                                                 levels = nitrogen_function_order)

# Global nitrogen presence
nitrogen_plot2 <- ggplot(nitrogen_genecount_long_sorted[!is.na(nitrogen_genecount_long_sorted$Function),], 
                        aes(x = Function, y = variable)) +
  geom_tile(aes(fill = log10(genecount + 1)), width = 1, height = 1, color = "gray") +  # heatmap
  geom_text(aes(label = round(genecount, 2)), color = "black", size = 3) +  # add text labels
  scale_fill_gradient(low = "white", high = "lightgrey") +  # gradient fill
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # remove y-axis text
    axis.text.x = element_blank(),  # remove x-axis text
    strip.background = element_blank(),
    panel.border = element_blank(),  # remove border
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_blank(),  # remove y-axis title
    panel.background = element_blank(),  # remove panel background
    legend.position = "none",  # remove legend
    panel.grid = element_blank()  # remove grid lines
  ) +   
  coord_flip()  # flip axes

# Convert heatmap to grob
nitrogen_grob_plot2 <- ggplotGrob(nitrogen_plot2)

# Nitrogen bar plot
nitrogen_plot3 <- ggplot(keggresult_by_category[keggresult_by_category$Module.Category %in% "Nitrogen metabolism",], 
                        aes(x = `Host Trend Group`, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "", y = "MAGs with metabolic pathway") +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("nitrogen_plot_overall.png", plot = nitrogen_plot3, dpi = 300, width = 3, height = 1.5)


########## Sulfur ###############

# Subset data where Category is "sulfur metabolism"
sulfur_result <- keggresult[keggresult$Category == "Sulfur\nmetabolism",]

# Sulfur scatter plot
sulfur_plot1 <- ggplot(sulfur_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 4, color = "black", alpha = 0.7, shape = 21) +  # set black border and fill by `fill`
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(x = NULL, y = 'MAGs with metabolic pathway', fill = 'Cluster') +
  coord_flip() +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  theme_light(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 13, colour = "black"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 50, by = 10))

# Get the order of Function in sulfur_result
sulfur_function_order <- unique(sulfur_result$Function)

# Reorder Function in genecount_long to match sulfur_result
sulfur_genecount_long_sorted <- genecount_long[genecount_long$Category == "Sulfur metabolism",]
sulfur_genecount_long_sorted$Function <- factor(sulfur_genecount_long_sorted$Function,
                                                  levels = sulfur_function_order)

# Global sulfur presence
sulfur_plot2 <- ggplot(sulfur_genecount_long_sorted[!is.na(sulfur_genecount_long_sorted$Function),], 
                         aes(x = Function, y = variable)) +
  geom_tile(aes(fill = log10(genecount + 1)), width = 1, height = 1, color = "gray") +  # heatmap
  geom_text(aes(label = round(genecount, 2)), color = "black", size = 3) +  # add text labels
  scale_fill_gradient(low = "white", high = "lightgrey") +  # gradient fill
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # remove y-axis text
    axis.text.x = element_blank(),  # remove x-axis text
    strip.background = element_blank(),
    panel.border = element_blank(),  # remove border
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_blank(),  # remove y-axis title
    panel.background = element_blank(),  # remove panel background
    legend.position = "none",  # remove legend
    panel.grid = element_blank()  # remove grid lines
  ) +   
  coord_flip()  # flip axes

# Convert heatmap to grob
sulfur_grob_plot2 <- ggplotGrob(sulfur_plot2)

# Sulfur bar plot
sulfur_plot3 <- ggplot(keggresult_by_category[keggresult_by_category$Module.Category %in% "Sulfur metabolism",], 
                         aes(x = `Host Trend Group`, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "", y = "MAGs with metabolic pathway") +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("sulfur_plot_overall.png", plot = sulfur_plot3, dpi = 300, width = 3, height = 1.5)

########## Carbon fixation ###############

# Subset data where Category is "carbon metabolism"
carbon_result <- keggresult[keggresult$Category == "Carbon\nfixation",]

# Carbon fixation scatter plot
carbon_plot1 <- ggplot(carbon_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 4, color = "black", alpha = 0.7, shape = 21) +  # set black border and fill by `fill`
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(x = NULL, y = 'MAGs with metabolic pathway', fill = 'Cluster') +
  coord_flip() +
  facet_grid(Category ~ ., scales = 'free_y', space = 'free') +
  theme_light(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_blank(),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 190), breaks = seq(0, 180, by = 30))

# Get the order of Function in carbon_result
carbon_function_order <- unique(carbon_result$Function)

# Reorder Function in genecount_long to match carbon_result
carbon_genecount_long_sorted <- genecount_long[genecount_long$Category == "Carbon fixation",]
carbon_genecount_long_sorted$Function <- factor(carbon_genecount_long_sorted$Function,
                                                levels = carbon_function_order)

# Global carbon fixation presence
carbon_plot2 <- ggplot(carbon_genecount_long_sorted[!is.na(carbon_genecount_long_sorted$Function),], 
                       aes(x = Function, y = variable)) +
  geom_tile(aes(fill = log10(genecount + 1)), width = 1, height = 1, color = "gray") +  # heatmap
  geom_text(aes(label = round(genecount, 2)), color = "black", size = 3) +  # add text labels
  scale_fill_gradient(low = "white", high = "lightgrey") +  # gradient fill
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # remove y-axis text
    axis.text.x = element_blank(),  # remove x-axis text
    strip.background = element_blank(),
    panel.border = element_blank(),  # remove border
    axis.title.x = element_blank(),  # remove x-axis title
    axis.title.y = element_blank(),  # remove y-axis title
    panel.background = element_blank(),  # remove panel background
    legend.position = "none",  # remove legend
    panel.grid = element_blank()  # remove grid lines
  ) +   
  coord_flip()  # flip axes

# Convert heatmap to grob
carbon_grob_plot2 <- ggplotGrob(carbon_plot2)

# Carbon fixation bar plot
carbon_plot3 <- ggplot(keggresult_by_category[keggresult_by_category$Module.Category %in% "Carbon fixation",], 
                       aes(x = `Host Trend Group`, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "", y = "MAGs with metabolic pathway") +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  coord_flip()

ggsave("carbon_plot_overall.png", plot = carbon_plot3, dpi = 300, width = 3, height = 1.5)


# Combine plots
methane_plot <- methane_plot1 + 
  annotation_custom(methane_grob_plot2, xmin = -.2, xmax = 10, ymin = 56, ymax = 63.8)

nitrogen_plot <- nitrogen_plot1 + 
  annotation_custom(nitrogen_grob_plot2, xmin = -.2, xmax = 7, ymin = 56, ymax = 63.8) 

sulfur_plot <- sulfur_plot1 + 
  annotation_custom(sulfur_grob_plot2, xmin = -.2, xmax = 3, ymin = 56, ymax = 63.8) 

carbon_plot <- carbon_plot1 + 
  annotation_custom(carbon_grob_plot2, xmin = -.2, xmax = 7, ymin = 177, ymax = 203) 


# Adjusting position
nitrogen_plot <- nitrogen_plot + 
  theme(plot.margin = margin(-2, 7, 70,-4)) 

sulfur_plot <- sulfur_plot + 
  theme(plot.margin = margin(-63, 6.8, 165, 18)) 

carbon_plot <- carbon_plot + 
  theme(plot.margin = margin(-163, 6.5, 190,-76)) 

fermentation_plot <- fermentation_plot  + 
  theme(plot.margin = margin(-183, 6.5, 200, 162)) 

# Combine all pathway plots into one figure
combined_plot <- grid.arrange(
  methane_plot, nitrogen_plot, sulfur_plot, carbon_plot, fermentation_plot,
  ncol = 1
)

ggsave("combined_plot.pdf", plot = combined_plot, dpi = 300, width = 15, height = 12)


# Optional: Plot legend only
legend_plot <- ggplot(carbon_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
  geom_point(size = 0, alpha = 0, shape = 21) +
  scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
  labs(fill = 'Treatment') +
  theme_void() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 5)))

ggsave("legend_only.png", plot = legend_plot, width = 4, height = 4)


fermentation_plot+theme(plot.margin = margin(-146, 5.5, 194, 320))