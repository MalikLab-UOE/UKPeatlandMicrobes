############   Figure 3 ################

library("tidyverse")
library("vegan")
library("ape")
library("RColorBrewer")
library("cowplot")
library("lme4")
library("lmerTest")    
library("emmeans")     
library("performance") 
library("nlme")
library("dplyr")
library("purrr")
library("broom")
library("svglite")

## Figure 3a ##
# Load metadata
metadata <- readRDS("metadata_simple.RDS") %>%
  mutate(treatment = factor(treatment, levels = c("NAT", "REST", "DAM")))

env_data<- read.csv("eco_index.csv")

# input data, hosts = prokaryotes 
host.metabolism <- readRDS("host_abundance_long_all_and_all_sites_and_metabolism.RDS")
#write.csv(host.metabolism, "host_metabolism.csv", row.names = FALSE)
host.contig.genome <- readRDS("host_contig_to_genome.RDS")

#Create dissimilarity matrix
tmeans.host <- readRDS("host_tmeans_norm_50.RDS")
tmeans.xform.host <- decostand(tmeans.host, method="hellinger")
bray_curtis_dist_host <- as.matrix(vegdist(t(tmeans.xform.host), method='bray'))

pcoa_host <- pcoa(as.dist(bray_curtis_dist_host))
axes_host <- as.data.frame(pcoa_host$vectors) # make a dataframe named axes, put pcoa values in there
axes_host$SampleID <- rownames(axes_host) # Give df extra column with the rownames in it 
axes_host <- merge(metadata[,c("site", "treatment")] %>% mutate("SampleID" = row.names(.)), axes_host, by.x = "SampleID", by.y = "SampleID")
head(axes_host)

# Store eigenvalues
eigval_host <- round(pcoa_host$values$Relative_eig * 100, digits = 2)
eigval_host <- data.frame(PC = 1:length(eigval_host), Eigval = eigval_host)

# ANOSIM for prokaryoteic MAGs
anosim_result_hosts <- anosim(bray_curtis_dist_host, metadata$site)
r_statistic_host <- anosim_result_hosts$statistic
p_value_host <- anosim_result_hosts$signif
summary(anosim_result_hosts)

### Permanova for MAGs Figure 3a 

set.seed(123)  # Ensure reproducibility

permanova <- adonis2(bray_curtis_dist_host ~ site*treatment, data = metadata, permutations = 9999)

# Print results
print(permanova)

# Plot PCoA for prokaryoteic MAGs
plot.pcoa.host <- ggplot(axes_host, aes(Axis.1, Axis.2)) +
  geom_point(aes(shape=as.character(treatment),
                 fill=as.character(site)),
             color = "black",
             size = 4,
             alpha=0.7,
             stroke=0.5) +
  xlab(paste("PCo1 (", eigval_host$Eigval[1], " %)", sep = "")) +
  ylab(paste("PCo2 (", eigval_host$Eigval[2], " %)", sep = "")) +
  scale_fill_manual(name="Site", limits=c("Migneint",  "Moor_House" ,"Crocach" ,    "Balmoral" ,   "Bowness" ,    "Langwell", "Stean"),
                    labels=c("Moor_House" = "Moor House"), values = c("#ff7f00" , "#ffff33" , "#4daf4a" ,"#e41a1c" ,  "#377eb8" , "#984ea3","#a65628"))+
  scale_shape_manual(name="Treatment", values=c(21, 22, 23), limits=c("NAT", "REST", "DAM"),
                     labels=c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored"))+
  # Ensure Site colors appear in the legend
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  theme_light() +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.position = "right")

plot.pcoa.host

ggsave("Figure_3a.png", dpi = 1000, width = 7, height = 5.5, units = "in")


## Figure 3b ##
axes_host2<- read.csv("pcoa_axes_host_all_eco_index.csv")

# Fit the global interaction model
global_model <- lm(Axis.1 ~ eco_index * site, data = axes_host2)

# Summary of the model
summary(global_model)

# R² (marginal = fixed only, conditional = fixed + random)
performance::r2(global_model)

site_slopes <- emtrends(global_model, ~ site, var = "eco_index")
summary(site_slopes)

# Example: per-site linear models and get R² and adjusted R²
site_r2 <- axes_host2 %>%
  group_by(site) %>%
  group_split() %>%
  map_df(~ {
    model <- lm( Axis.1 ~ eco_index, data = .x)
    summary_model <- summary(model)
    tibble(
      Site = unique(.x$site),
      R2 = summary_model$r.squared,
      adj_R2 = summary_model$adj.r.squared,
      p_value = coef(summary_model)[2, "Pr(>|t|)"]
    )
  })

print(site_r2)


#Stats for Figure 3b

axes_host <- merge(axes_host,env_data,by="SampleID")
m1 <- lmer(Axis.1 ~ eco_index + (1 | site), data = axes_host)

# Extract variance components
var_comp <- as.data.frame(VarCorr(m1))
var_random <- var_comp$vcov[var_comp$grp == "Site"]         # Random intercept variance
var_resid  <- attr(VarCorr(m1), "sc")^2                     # Residual variance

# Calculate total variance
var_total <- var_random + var_resid

# Marginal R²: fixed effects only
# Conditional R²: fixed + random effects
# We'll use the method from Nakagawa & Schielzeth (2013)

# Calculate fitted values based on fixed effects
y_hat_fixed <- predict(m1, re.form = NA)

# Variance of fitted values from fixed effects
var_fixed <- var(y_hat_fixed)

# R² values
R2_marginal <- var_fixed / var_total
R2_conditional <- (var_fixed + var_random) / var_total

# Show results
cat("Marginal R² (fixed effects only):", round(R2_marginal, 3), "\n")
cat("Conditional R² (fixed + random effects):", round(R2_conditional, 3), "\n")

# PCo1 vs. ecosystem health index for all samples
Figure_3b <- ggplot(axes_host2, aes(x = eco_index, y = Axis.1)) +
  geom_point(color = "black",size = 4, alpha = 0.7,aes(shape=as.character(treatment),
                                                       fill=as.character(site))) +  
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
  ylab("PCo1")+ 
  # Ensure Site colors appear in the legend
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  labs(x = "Ecosystem health index") +
  theme_light(base_size = 14)+
  theme(axis.text = element_text(size = 13, colour = "black")) 

Figure_3b 

ggsave("Figure_3b.png",  dpi=1000, width = 4, height = 4, units = "in")

## Figure S8 ##
site_order <- c("Balmoral","Bowness","Stean","Moor_House","Langwell","Crocach","Migneint")
axes_host2$site <- factor(axes_host2$site, levels = site_order)
site_labels <- c(
  "Balmoral" = "Balmoral",
  "Bowness" = "Bowness",
  "Crocach" = "Crocach",
  "Langwell" = "Langwell",
  "Migneint" = "Migneint",
  "Moor_House" = "Moor House",  # <- rename here
  "Stean" = "Stean"
)

ggplot(axes_host2, aes(x = eco_index, y = Axis.1)) +
  geom_point(aes(fill = site, shape = treatment), size = 4, color = "black", alpha = 0.7,stroke=0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", aes(group = 1)) + 
  scale_fill_manual(
    name = "Site",
    limits = c("Migneint", "Moor_House", "Crocach", "Balmoral", "Bowness", "Langwell", "Stean"),
    values = c("#ff7f00", "#ffff33", "#4daf4a", "#e41a1c", "#377eb8", "#984ea3", "#a65628")
  ) +
  scale_shape_manual(
    name = "Treatment",
    values = c("NAT" = 21, "REST" = 22, "DAM" = 23),
    labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")
  ) +
  labs(
    x = "Ecosystem health index",
    y = "PCo1"
  ) +
  theme_light() +
  theme(legend.position = "None")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank()) +
  guides(
    fill = guide_legend(order = 1, title = "Site", override.aes = list(shape = 21, color = "black")), 
    shape = guide_legend(order = 2, title = "Treatment") ) +
  facet_grid(. ~ site, scales = "free", labeller = as_labeller(site_labels))

ggsave("Figure_S8.png",  dpi=1000, width = 8, height = 3, units = "in")

## Figure 3c ##

# Load and format data
trimmed.mean.cov <- read.csv(
  "MAG_trimmed_mean_cov.tsv",
  sep = "\t") %>%
  rename_with(~str_remove(., '.filtered.Trimmed.Mean'))
rownames(trimmed.mean.cov) <- trimmed.mean.cov$Genome
trimmed.mean.cov <- trimmed.mean.cov[, !(colnames(trimmed.mean.cov) %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F"))]
trimmed.mean.cov <- trimmed.mean.cov %>%
  select(-Genome)

coverage <- read.csv("MAG_covered_fraction.tsv", sep = "\t", header = TRUE, row.names = 1)
colnames(coverage) <- sub(".filtered.Covered.Fraction", "", colnames(coverage))
coverage <- coverage %>%
  select(-contains(".filtered"))
coverage <- coverage[, !(colnames(coverage) %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F"))]

trimmed.mean.cov[coverage < 0.50] <- 0
trimmed.mean.cov <- trimmed.mean.cov[rowSums(trimmed.mean.cov>0) > 1,] # Remove singletons
trimmed.mean.cov <- trimmed.mean.cov[, colSums(trimmed.mean.cov)>0] # Remove singletons

seqdepth <- read_csv(file="seq_depth.csv", show_col_types = FALSE)
seqdepth_R1 <- subset(seqdepth, Pair=="R1")
trimmed.mean.cov.norm <- trimmed.mean.cov
for(col in colnames(trimmed.mean.cov.norm)) {
  trimmed.mean.cov.norm[[col]] <- trimmed.mean.cov.norm[[col]] / seqdepth_R1$Hundred.Millions.Reads[seqdepth_R1$Sample == col]
}

metadata <- read.csv("./metadata.tsv", sep = "\t") %>%
  filter(
    !sample %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F")
  )

# Transform data, compute Bray-Curtis, run PCoA
eigvals  <- c()
pcoa.data <- c()
anosims <- c()
for(Site in unique(metadata$site)) {
  meta.site <- metadata %>%
    filter(site == Site)
  data.site <- trimmed.mean.cov %>%
    select(any_of(meta.site$sample))
  data.hellinger <- decostand(data.site, method = "hellinger")
  
  bray <- as.matrix(vegdist(t(data.hellinger), method = "bray"))
  
  pcoa <- pcoa(as.dist(bray))
  pcoa.axes <- as.data.frame(pcoa$vectors) %>%
    mutate(sample = rownames(.)) %>%
    left_join(
      meta.site %>%
        select(sample, site, treatment),
    ) %>%
    select(sample, site, treatment, everything())
  pcoa.data[[Site]] <- pcoa.axes
  
  eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
  eigval <- data.frame(PCoA = 1:length(eigval), Eigval = eigval)
  eigval$site <- Site
  eigvals[[Site]] <- eigval
  
  anosim.res <- anosim(bray, meta.site$treatment)
  anosim.r <- anosim.res$statistic
  anosim.p <- anosim.res$signif
  anosims[[Site]] <- c("R" = anosim.r, "P" = anosim.p, "site" = Site)
}
axes.df <- bind_rows(pcoa.data)
eigvals.df <- bind_rows(eigvals)
anosims.df <- bind_rows(anosims)
anosims.df$P.adj <- p.adjust(anosims.df$P, method = "BH")
anosims.df <- anosims.df  %>%
  mutate(R = as.numeric(R), P = as.numeric(P), P.adj = as.numeric(P.adj))


# Plot
plots <- lapply(unique(axes.df$site), function(Site) {
  site.data <- axes.df %>% filter(site == Site)
  site.anosim <- anosims.df %>%
    filter(site == Site)
  site.eigval <- eigvals.df %>%
    filter(site == Site)
  
  x.pos <- Inf
  y.pos <- min(site.data$Axis.2, na.rm = TRUE) - 0.2 * diff(range(site.data$Axis.2, na.rm = TRUE))
  
  ggplot(site.data, aes(Axis.1, Axis.2)) +
    geom_point(aes(fill = as.character(treatment)),
               shape = 21,
               color = "black",
               size = 3,
               alpha = 0.7,
               stroke = 0.5) +
    geom_text(
      data = data.frame(
        x = x.pos, y = y.pos,
        label = paste0("italic(R)==", round(site.anosim$R, 2), "*','~italic(P)==", signif(site.anosim$P.adj, 2))
      ),
      aes(x = x, y = y, label = label),
      parse = TRUE, hjust = 1.1, vjust = -0.1, size = 2.75, inherit.aes = FALSE
    ) +
    labs(
      x = paste("PCo1 (", site.eigval$Eigval[1], " %)", sep = ""),
      y = paste("PCo2 (", site.eigval$Eigval[2], " %)", sep = ""),
      title = if(Site == "Moor_House"){title <- "Moor House"}else{title <- Site}
    ) +
    scale_fill_manual(name = "Treatment",
                      values = c("NAT" = "#4daf4a", "REST" = "#377eb8", "DAM" = "#e41a1c"),
                      labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_light() +
    theme(text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
})
plot.with.legend <- plots[[1]] + theme(legend.position = "right")
legend <- get_legend(plot.with.legend)
plots[[8]] <- ggplot() + theme_light() + theme(panel.grid.minor = element_blank(), panel.border = element_blank())
plots[[9]] <- legend
plot.pcoa <- plot_grid(plotlist = plots, ncol = 3, align = "hv", axis = "tblr")

plot.pcoa

ggsave(plot.pcoa,
       filename="fig3c.svg",
       device="svg",
       width=8,
       height=8,
       bg = "white",
       units="in")




