#############  Figure 6 ###################


# Load required libraries
library(ggplot2)
library(readr)
library(ggsci)
library(patchwork)
library(gridExtra)
library("tidyverse");packageVersion("tidyverse")
library("dplyr");packageVersion("dplyr")
library("reshape2");packageVersion("reshape2")
library("tidyverse");packageVersion("tidyverse")
library("DESeq2");packageVersion("DESeq2")
library(broom)

# MAG trend groups
# Load and format data
## Metadata
metadata <- read.csv("./metadata.tsv", sep = "\t") %>%
  filter(
    !sample %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F")
  )
row.names(metadata) <- metadata$sample


## Contig to genome table
contig2genome <- read.csv("./contig_to_genome.tsv", header=FALSE, sep="\t")
colnames(contig2genome) <- c("contig", "Genome")
contig2genome$Genome <- gsub(".fasta", "", contig2genome$Genome)
head(contig2genome)

## Genome counts
genome_counts <- read.csv("./genome_counts.tsv", sep = "\t", header = TRUE, row.names = 1)
colnames(genome_counts) <- sub(".filtered.Read.Count", "", colnames(genome_counts))
genome_counts <- genome_counts %>% select(-contains(".filtered"))
genome_counts <- genome_counts[, !(colnames(genome_counts) %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F"))] # remove extra Langwell restored replicates
head(genome_counts)


## Trimmed mean genome coverage
tmeans <- read.csv("./genome_trimmed_mean_cov.tsv", sep = "\t", header = TRUE, row.names = 1)
colnames(tmeans) <-  sub(".filtered.Trimmed.Mean", "", colnames(tmeans))
tmeans <- tmeans %>% select(-contains(".filtered"))
tmeans <- tmeans[, !(colnames(tmeans) %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F"))] # remove extra Langwell restored replicates
head(tmeans)

## Genome coverage
coverage <- read.csv("./genome_covered_fraction.tsv", sep = "\t", header = TRUE, row.names = 1)
colnames(coverage) <-  sub(".filtered.Covered.Fraction", "", colnames(coverage))
coverage <- coverage %>% select(-contains(".filtered"))
coverage <- coverage[, !(colnames(coverage) %in% c("LASYr2D", "LASYr2E", "LASYr2F", "LAWAr2D", "LAWAr2E", "LAWAr2F"))] # remove extra Langwell restored replicates

## Pre-filter trimmed means and values to 0 if coverage < 50%
#But not the genome counts, since DESeq2 requires raw, unfiltered count matrices
tmeans[coverage < 0.50] <- 0
tmeans <- tmeans[rowSums(tmeans>0) > 1,] # Remove singletons
tmeans <- tmeans[, colSums(tmeans)>0] # Remove singletons
head(tmeans)

## And normalize trimmed means by the sequencing depth per sample (per 100M reads)
seqdepth <- read_csv(file="./seq_depth.csv", show_col_types = FALSE)
seqdepth_R1 <- subset(seqdepth, Pair=="R1")
tmeans_norm <- tmeans
for(col in colnames(tmeans_norm)) {
  tmeans_norm[[col]] <- tmeans_norm[[col]] / seqdepth_R1$Hundred.Millions.Reads[seqdepth_R1$Sample == col]
}
head(tmeans_norm)

## MAG GTDB taxonomy
gtdb <- read_csv(file="./MAGs_gtdb_taxonomy.csv", show_col_types = FALSE)
head(gtdb)


## Metabolic functions
metabolic_modules <- read.csv("./METABOLIC_KEGG_module_presence.csv")
head(metabolic_modules)
metabolic_result <- read.csv("./METABOLIC_result.tsv", sep="\t")
head(metabolic_result)

## Reformat METABOLIC results
#Only retain metabolic categories/functions of interest, remove pathways that are absent in all MAGs, rename categories, and move columns.

modules_of_interest <- unique(subset(metabolic_modules, Module.Category %in% c("Methane metabolism", "Nitrogen metabolism", "Sulfur metabolism", "Carbon fixation"))$Module.ID)

# Filter the table for the specified categories and rename
metabolic_filtered <- metabolic_modules %>%
  dplyr::rename(Category = Module.Category) %>%
  filter(Module.ID %in% modules_of_interest)

# Identify columns representing bins
bin_columns <- grep("\\Module\\.presence$", colnames(metabolic_filtered), value = TRUE)

# Remove rows where the function is absent in all bins
metabolic_filtered <- metabolic_filtered %>%
  rowwise() %>%
  filter(any(across(all_of(bin_columns), ~ . == "Present")))

metabolic_filtered <- metabolic_filtered %>%
  select(Category, Module, Module.ID, everything())

# Rename Module to Pathway and split the Pathway into a Module and Function
metabolic_filtered <- metabolic_filtered %>%
  dplyr::rename(Pathway = Module)
metabolic_filtered[c("Module", "Function")] <- do.call(rbind, lapply(metabolic_filtered$Pathway, function(x) {
  parts <- unlist(strsplit(x, ", ", fixed = TRUE))
  if (length(parts) > 1) {
    c(parts[1], paste(parts[-1], collapse = ", "))
  } else {
    c(parts[1], NA)
  }
}))

metabolic_filtered <- metabolic_filtered %>%
  select(Module, everything()) %>%
  select(-Function) %>%
  arrange(Category, Module)

# Convert back to data frame if necessary
metabolic_filtered <- as.data.frame(metabolic_filtered)

head(metabolic_filtered)


## Additional fermentation functions
# Filter the subset of fermentation data
fermentation_subset <- subset(metabolic_result, Category %in% c("Fermentation"))

# Rename the columns of the subset to match `metabolic_filtered`
fermentation_subset <- fermentation_subset %>%
  dplyr::rename(Module = Category, Pathway = Function) %>%
  mutate(Module = case_when(Module == "Fermentation" ~ "Fermentation"),
         Category = case_when(Module == "Fermentation" ~ "Fermentation"),
         Module_orig = Pathway)

fermentation_subset <- fermentation_subset %>%
  select(Module, Category, Pathway, everything()) %>%
  select(-Module_orig) %>%
  arrange(Category, Module)

# Identify bin columns in both dataframes
bin_columns_result <- grep("\\.Hmm\\.presence$", colnames(fermentation_subset), value = TRUE)
bin_columns_filtered <- grep("\\.Module\\.presence$", colnames(metabolic_filtered), value = TRUE)

# Ensure presence/absence column names match by renaming Hmm.presence to Module.presence
colnames(fermentation_subset) <- gsub("\\.Hmm\\.presence$", ".Module.presence", colnames(fermentation_subset))

# Remove any columns with ".Hits" or ".Hit.numbers" or other unwanted columns
fermentation_subset <- fermentation_subset %>%
  select(-matches("\\.Hits$|\\.Hit\\.numbers$")) %>%
  select(-Gene.abbreviation, -Corresponding.KO)

# Verify that the columns now match
bin_columns_result_renamed <- grep("\\.Module\\.presence$", colnames(fermentation_subset), value = TRUE)

# Now, append the complex carbon degradation subset to `metabolic_filtered`
metabolic_filtered <- bind_rows(metabolic_filtered, fermentation_subset)

# Clean up and reorder columns in `metabolic_filtered` (optional)
metabolic_filtered <- metabolic_filtered %>%
  select(Module, Pathway, everything()) %>%
  arrange(Category, Module)

# Convert back to data frame if necessary
metabolic_filtered <- as.data.frame(metabolic_filtered)

# Remove rows where the function is absent in all bins
metabolic_filtered <- metabolic_filtered %>%
  rowwise() %>%
  filter(any(across(all_of(bin_columns), ~ . == "Present")))

write_csv(metabolic_filtered, file = "./genome_metabolic_functions.csv")

# Preview the result
head(metabolic_filtered)


## Make a long version of the filtered METABOLIC results
metabolic_filtered_long <- metabolic_filtered %>%
  pivot_longer(
    cols = ends_with(".Module.presence"),
    names_to = "Genome", 
    values_to = "Module_presence",
  ) %>%
  filter(Module_presence == "Present") %>%
  select(Module, Category, Pathway, Module.ID, Genome, -Module_presence)

metabolic_filtered_long <- metabolic_filtered_long[with(metabolic_filtered_long, order(Category, Module, Pathway, Genome)),]

metabolic_filtered_long$Genome <- gsub('.Module.presence', '', metabolic_filtered_long$Genome)

head(metabolic_filtered_long)

# Genome differential abundance with DESeq2
## Run models and retrieve results, site by site
### Balmoral
#### DESeq
metadata_balmoral <- subset(metadata, site=="Balmoral")
dds_balmoral <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                         select(all_of(rownames(metadata_balmoral))),
                                       colData = metadata_balmoral,
                                       design =  ~ treatment)
dds_balmoral <- DESeq(dds_balmoral, test = "LRT", reduced = ~1)
res_balmoral <- results(dds_balmoral) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_balmoral

#### Balmoral results
res_balmoral <- res_balmoral %>%
  filter(!is.na(pvalue))
head(res_balmoral)


#### Subset to return genomes with padj < 0.05 in Balmoral
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_balmoral <- res_balmoral %>%
  filter(padj < padj.cutoff)
insig_balmoral <- res_balmoral %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Balmoral
tmeans_balmoral <- tmeans_norm %>%
  select(all_of(rownames(metadata_balmoral)))
tmeans_sig_balmoral <- tmeans_balmoral[rownames(tmeans_balmoral) %in% sig_balmoral$Genome, ]

### Bowness
#### DESeq
metadata_bowness <- subset(metadata, site=="Bowness")
dds_bowness <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                        select(all_of(rownames(metadata_bowness))),
                                      colData = metadata_bowness,
                                      design =  ~ treatment)
dds_bowness <- DESeq(dds_bowness, test = "LRT", reduced = ~1)
res_bowness <- results(dds_bowness) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_bowness


#### Bowness results
res_bowness <- res_bowness %>%
  filter(!is.na(pvalue))
head(res_bowness)

#### Subset to return genomes with padj < 0.05 in Bowness
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_bowness <- res_bowness %>%
  filter(padj < padj.cutoff)
insig_bowness <- res_bowness %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Bowness
tmeans_bowness <- tmeans_norm %>%
  select(all_of(rownames(metadata_bowness)))
tmeans_sig_bowness <- tmeans_bowness[rownames(tmeans_bowness) %in% sig_bowness$Genome, ]

### Crocach
#### DESeq
metadata_crocach <- subset(metadata, site=="Crocach")
dds_crocach <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                        select(all_of(rownames(metadata_crocach))),
                                      colData = metadata_crocach,
                                      design =  ~ treatment)
dds_crocach <- DESeq(dds_crocach, test = "LRT", reduced = ~1)
res_crocach <- results(dds_crocach) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_crocach

#### Crocach results
res_crocach <- res_crocach %>%
  filter(!is.na(pvalue))
head(res_crocach)

#### Subset to return genomes with padj < 0.05 in Crocach
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_crocach <- res_crocach %>%
  filter(padj < padj.cutoff)
insig_crocach <- res_crocach %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Crocach
tmeans_crocach <- tmeans_norm %>%
  select(all_of(rownames(metadata_crocach)))
tmeans_sig_crocach <- tmeans_crocach[rownames(tmeans_crocach) %in% sig_crocach$Genome, ]

### Langwell
### DESeq
metadata_langwell <- subset(metadata, site=="Langwell")
dds_langwell <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                         select(all_of(rownames(metadata_langwell))),
                                       colData = metadata_langwell,
                                       design =  ~ treatment)
dds_langwell <- DESeq(dds_langwell, test = "LRT", reduced = ~1)
res_langwell <- results(dds_langwell) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_langwell

#### Langwell results
res_langwell <- res_langwell %>%
  filter(!is.na(pvalue))
head(res_langwell)

#### Subset to return genomes with padj < 0.05 in Langwell
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_langwell <- res_langwell %>%
  filter(padj < padj.cutoff)
insig_langwell <- res_langwell %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Langwell
tmeans_langwell <- tmeans_norm %>%
  select(all_of(rownames(metadata_langwell)))
tmeans_sig_langwell <- tmeans_langwell[rownames(tmeans_langwell) %in% sig_langwell$Genome, ]

### Migneint
#### DESeq
metadata_migneint <- subset(metadata, site=="Migneint")
dds_migneint <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                         select(all_of(rownames(metadata_migneint))),
                                       colData = metadata_migneint,
                                       design =  ~ treatment)
dds_migneint <- DESeq(dds_migneint, test = "LRT", reduced = ~1)
res_migneint <- results(dds_migneint) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_migneint

#### Migneint results
res_migneint <- res_migneint %>%
  filter(!is.na(pvalue))
head(res_migneint)

#### Subset to return genomes with padj < 0.05 in Migneint
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_migneint <- res_migneint %>%
  filter(padj < padj.cutoff)
insig_migneint <- res_migneint %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Migneint
tmeans_migneint <- tmeans_norm %>%
  select(all_of(rownames(metadata_migneint)))
tmeans_sig_migneint <- tmeans_migneint[rownames(tmeans_migneint) %in% sig_migneint$Genome, ]

### Moor House
#### DESeq
metadata_moor_house <- subset(metadata, site=="Moor_House")
dds_moor_house <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                           select(all_of(rownames(metadata_moor_house))),
                                         colData = metadata_moor_house,
                                         design =  ~ treatment)
dds_moor_house <- DESeq(dds_moor_house, test = "LRT", reduced = ~1)
res_moor_house <- results(dds_moor_house) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_moor_house


#### Moor House results
res_moor_house <- res_moor_house %>%
  filter(!is.na(pvalue))
head(res_moor_house)

#### Subset to return genomes with padj < 0.05 in Moor House
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_moor_house <- res_moor_house %>%
  filter(padj < padj.cutoff)
insig_moor_house <- res_moor_house %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Moor House
tmeans_moor_house <- tmeans_norm %>%
  select(all_of(rownames(metadata_moor_house)))
tmeans_sig_moor_house <- tmeans_moor_house[rownames(tmeans_moor_house) %in% sig_moor_house$Genome, ]

### Stean
#### DESeq
metadata_stean <- subset(metadata, site=="Stean")
dds_stean <- DESeqDataSetFromMatrix(countData = genome_counts %>%
                                      select(all_of(rownames(metadata_stean))),
                                    colData = metadata_stean,
                                    design =  ~ treatment)
dds_stean <- DESeq(dds_stean, test = "LRT", reduced = ~1)
res_stean <- results(dds_stean) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename("Genome" = "gene")
res_stean

#### Stean results
res_stean <- res_stean %>%
  filter(!is.na(p.value))
head(res_stean)

#### Subset to return genomes with padj < 0.05 in Stean
padj.cutoff <- 0.05 # Set alpha to 0.05
sig_stean <- res_stean %>%
  filter(padj < padj.cutoff)
insig_stean <- res_stean %>%
  filter(padj >= padj.cutoff)

#### Subset tmeans with just significant genomes in Stean
tmeans_stean <- tmeans_norm %>%
  select(all_of(rownames(metadata_stean)))
tmeans_sig_stean <- tmeans_stean[rownames(tmeans_stean) %in% sig_stean$Genome, ]

## Combine DESeq2 result tables and write
combined_deseq_results <- rbind(res_balmoral %>% mutate(site = "Balmoral"),
                                res_bowness %>% mutate(site = "Bowness"),
                                res_crocach %>% mutate(site = "Crocach"),
                                res_langwell %>% mutate(site = "Langwell"),
                                res_migneint %>% mutate(site = "Migneint"),
                                res_moor_house %>% mutate(site = "Moor House"),
                                res_stean %>% mutate(site = "Stean")
)
write_csv(combined_deseq_results, file = "./deseq_results_combined.csv")
head(combined_deseq_results)

## Cluster genomes by their abundance patterns across treatments
### General functions and format metadata
tidy_otu <- function(otu) {
  as.data.frame(otu) %>%
    mutate(Genome = row.names(otu)) %>%
    tidyr::gather(key = "Sample", value = "Count", -Genome)
}
rel_ab <- function(otu, total = 100) {
  t(t(otu)/colSums(otu)) * 100
}

metadata$Sample <- rownames(metadata)

### Balmoral
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.balmoral <- tmeans_sig_balmoral %>%
  filter(row.names(.) %in% sig_balmoral$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.balmoral <- genome.zs.balmoral %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.balmoral) <- genome.zs.matrix.balmoral$Treatment
genome.zs.matrix.balmoral <- genome.zs.matrix.balmoral[,-1]
genome.zs.matrix.balmoral <- as.matrix(genome.zs.matrix.balmoral)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.balmoral[is.na(genome.zs.matrix.balmoral)] <- 0
genome.zs.matrix.balmoral[is.nan(genome.zs.matrix.balmoral)] <- 0
genome.zs.matrix.balmoral[is.infinite(genome.zs.matrix.balmoral)] <- 0

genome.dist.balmoral <- dist(t(genome.zs.matrix.balmoral))

#### Perform the hierarchical clustering
genome.hc.balmoral <- hclust(as.dist(genome.dist.balmoral), method = "ward.D")
genome.ord.balmoral <- genome.hc.balmoral$labels[genome.hc.balmoral$order]
genome.ord.balmoral <- data.frame(Genome = genome.ord.balmoral, order = 1:length(genome.ord.balmoral))
genome.cut.balmoral <- cutree(genome.hc.balmoral[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.balmoral <- data.frame(Genome = names(genome.cut.balmoral),
                                       Cluster = genome.cut.balmoral) %>%
  inner_join(genome.ord.balmoral, by = "Genome") %>%
  inner_join(genome.zs.balmoral, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.balmoral <- genome.clusters.balmoral %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.balmoral$Title <- paste("Group: ", genome.clusters.balmoral$Cluster, " - ", genome.clusters.balmoral$Cluster.size, " genomes", sep = "")

head(genome.clusters.balmoral)

#### Plot Z-Score clusters
plot.clusters.balmoral <-
  ggplot(data = genome.clusters.balmoral, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Balmoral")
plot.clusters.balmoral

#### Manually rename cluster titles based on abundance trends
genome.clusters.balmoral <- genome.clusters.balmoral %>%
  mutate(Title = case_when(Cluster == 1 ~ "Natural-abundant",
                           Cluster == 2 ~ "Restored-abundant",
                           Cluster == 3 ~ "Damaged-abundant"))
genome.clusters.balmoral$site <- "Balmoral"
head(genome.clusters.balmoral)

### Bowness
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.bowness <- tmeans_sig_bowness %>%
  filter(row.names(.) %in% sig_bowness$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.bowness <- genome.zs.bowness %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.bowness) <- genome.zs.matrix.bowness$Treatment
genome.zs.matrix.bowness <- genome.zs.matrix.bowness[,-1]
genome.zs.matrix.bowness <- as.matrix(genome.zs.matrix.bowness)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.bowness[is.na(genome.zs.matrix.bowness)] <- 0
genome.zs.matrix.bowness[is.nan(genome.zs.matrix.bowness)] <- 0
genome.zs.matrix.bowness[is.infinite(genome.zs.matrix.bowness)] <- 0

genome.dist.bowness <- dist(t(genome.zs.matrix.bowness))

#### Perform the hierarchical clustering
genome.hc.bowness <- hclust(as.dist(genome.dist.bowness), method = "ward.D")
genome.ord.bowness <- genome.hc.bowness$labels[genome.hc.bowness$order]
genome.ord.bowness <- data.frame(Genome = genome.ord.bowness, order = 1:length(genome.ord.bowness))
genome.cut.bowness <- cutree(genome.hc.bowness[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.bowness <- data.frame(Genome = names(genome.cut.bowness),
                                      Cluster = genome.cut.bowness) %>%
  inner_join(genome.ord.bowness, by = "Genome") %>%
  inner_join(genome.zs.bowness, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.bowness <- genome.clusters.bowness %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.bowness$Title <- paste("Group: ", genome.clusters.bowness$Cluster, " - ", genome.clusters.bowness$Cluster.size, " genomes", sep = "")

head(genome.clusters.bowness)

#### Plot Z-Score clusters
plot.clusters.bowness <-
  ggplot(data = genome.clusters.bowness, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Bowness")
plot.clusters.bowness

#### Manually rename cluster titles based on abundance trends
genome.clusters.bowness <- genome.clusters.bowness %>%
  mutate(Title = case_when(Cluster == 1 ~ "Restored-abundant",
                           Cluster == 2 ~ "Damaged-abundant",
                           Cluster == 3 ~ "Natural-abundant"))
genome.clusters.bowness$site <- "Bowness"
head(genome.clusters.bowness)

### Crocach
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.crocach <- tmeans_sig_crocach %>%
  filter(row.names(.) %in% sig_crocach$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.crocach <- genome.zs.crocach %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.crocach) <- genome.zs.matrix.crocach$Treatment
genome.zs.matrix.crocach <- genome.zs.matrix.crocach[,-1]
genome.zs.matrix.crocach <- as.matrix(genome.zs.matrix.crocach)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.crocach[is.na(genome.zs.matrix.crocach)] <- 0
genome.zs.matrix.crocach[is.nan(genome.zs.matrix.crocach)] <- 0
genome.zs.matrix.crocach[is.infinite(genome.zs.matrix.crocach)] <- 0

genome.dist.crocach <- dist(t(genome.zs.matrix.crocach))

#### Perform the hierarchical clustering
genome.hc.crocach <- hclust(as.dist(genome.dist.crocach), method = "ward.D")
genome.ord.crocach <- genome.hc.crocach$labels[genome.hc.crocach$order]
genome.ord.crocach <- data.frame(Genome = genome.ord.crocach, order = 1:length(genome.ord.crocach))
genome.cut.crocach <- cutree(genome.hc.crocach[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.crocach <- data.frame(Genome = names(genome.cut.crocach),
                                      Cluster = genome.cut.crocach) %>%
  inner_join(genome.ord.crocach, by = "Genome") %>%
  inner_join(genome.zs.crocach, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.crocach <- genome.clusters.crocach %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.crocach$Title <- paste("Group: ", genome.clusters.crocach$Cluster, " - ", genome.clusters.crocach$Cluster.size, " genomes", sep = "")

head(genome.clusters.crocach)

#### Plot Z-Score clusters
plot.clusters.crocach <-
  ggplot(data = genome.clusters.crocach, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Crocach")
plot.clusters.crocach

#### Manually rename cluster titles based on abundance trends
genome.clusters.crocach <- genome.clusters.crocach %>%
  mutate(Title = case_when(Cluster == 1 ~ "Natural-abundant",
                           Cluster == 2 ~ "Restored-abundant",
                           Cluster == 3 ~ "Damaged-abundant"))
genome.clusters.crocach$site <- "Crocach"
head(genome.clusters.crocach)

### Langwell
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.langwell <- tmeans_sig_langwell %>%
  filter(row.names(.) %in% sig_langwell$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.langwell <- genome.zs.langwell %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.langwell) <- genome.zs.matrix.langwell$Treatment
genome.zs.matrix.langwell <- genome.zs.matrix.langwell[,-1]
genome.zs.matrix.langwell <- as.matrix(genome.zs.matrix.langwell)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.langwell[is.na(genome.zs.matrix.langwell)] <- 0
genome.zs.matrix.langwell[is.nan(genome.zs.matrix.langwell)] <- 0
genome.zs.matrix.langwell[is.infinite(genome.zs.matrix.langwell)] <- 0

genome.dist.langwell <- dist(t(genome.zs.matrix.langwell))

#### Perform the hierarchical clustering
genome.hc.langwell <- hclust(as.dist(genome.dist.langwell), method = "ward.D")
genome.ord.langwell <- genome.hc.langwell$labels[genome.hc.langwell$order]
genome.ord.langwell <- data.frame(Genome = genome.ord.langwell, order = 1:length(genome.ord.langwell))
genome.cut.langwell <- cutree(genome.hc.langwell[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.langwell <- data.frame(Genome = names(genome.cut.langwell),
                                       Cluster = genome.cut.langwell) %>%
  inner_join(genome.ord.langwell, by = "Genome") %>%
  inner_join(genome.zs.langwell, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.langwell <- genome.clusters.langwell %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.langwell$Title <- paste("Group: ", genome.clusters.langwell$Cluster, " - ", genome.clusters.langwell$Cluster.size, " genomes", sep = "")

head(genome.clusters.langwell)

#### Plot Z-Score clusters
plot.clusters.langwell <-
  ggplot(data = genome.clusters.langwell, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Langwell")
plot.clusters.langwell

#### Manually rename cluster titles based on abundance trends
genome.clusters.langwell <- genome.clusters.langwell %>%
  mutate(Title = case_when(Cluster == 1 ~ "Restored-abundant",
                           Cluster == 2 ~ "Natural-abundant",
                           Cluster == 3 ~ "Damaged-abundant"))
genome.clusters.langwell$site <- "Langwell"
head(genome.clusters.langwell)

### Migneint
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.migneint <- tmeans_sig_migneint %>%
  filter(row.names(.) %in% sig_migneint$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.migneint <- genome.zs.migneint %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.migneint) <- genome.zs.matrix.migneint$Treatment
genome.zs.matrix.migneint <- genome.zs.matrix.migneint[,-1]
genome.zs.matrix.migneint <- as.matrix(genome.zs.matrix.migneint)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.migneint[is.na(genome.zs.matrix.migneint)] <- 0
genome.zs.matrix.migneint[is.nan(genome.zs.matrix.migneint)] <- 0
genome.zs.matrix.migneint[is.infinite(genome.zs.matrix.migneint)] <- 0

genome.dist.migneint <- dist(t(genome.zs.matrix.migneint))

#### Perform the hierarchical clustering
genome.hc.migneint <- hclust(as.dist(genome.dist.migneint), method = "ward.D")
genome.ord.migneint <- genome.hc.migneint$labels[genome.hc.migneint$order]
genome.ord.migneint <- data.frame(Genome = genome.ord.migneint, order = 1:length(genome.ord.migneint))
genome.cut.migneint <- cutree(genome.hc.migneint[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.migneint <- data.frame(Genome = names(genome.cut.migneint),
                                       Cluster = genome.cut.migneint) %>%
  inner_join(genome.ord.migneint, by = "Genome") %>%
  inner_join(genome.zs.migneint, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.migneint <- genome.clusters.migneint %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.migneint$Title <- paste("Group: ", genome.clusters.migneint$Cluster, " - ", genome.clusters.migneint$Cluster.size, " genomes", sep = "")

head(genome.clusters.migneint)

#### Plot Z-Score clusters
plot.clusters.migneint <-
  ggplot(data = genome.clusters.migneint, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Migneint")
plot.clusters.migneint

#### Manually rename cluster titles based on abundance trends
genome.clusters.migneint <- genome.clusters.migneint %>%
  mutate(Title = case_when(Cluster == 1 ~ "Damaged-abundant",
                           Cluster == 2 ~ "Natural-abundant",
                           Cluster == 3 ~ "Restored-abundant"))
genome.clusters.migneint$site <- "Migneint"
head(genome.clusters.migneint)

### Moor House
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.moor_house <- tmeans_sig_moor_house %>%
  filter(row.names(.) %in% sig_moor_house$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.moor_house <- genome.zs.moor_house %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.moor_house) <- genome.zs.matrix.moor_house$Treatment
genome.zs.matrix.moor_house <- genome.zs.matrix.moor_house[,-1]
genome.zs.matrix.moor_house <- as.matrix(genome.zs.matrix.moor_house)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.moor_house[is.na(genome.zs.matrix.moor_house)] <- 0
genome.zs.matrix.moor_house[is.nan(genome.zs.matrix.moor_house)] <- 0
genome.zs.matrix.moor_house[is.infinite(genome.zs.matrix.moor_house)] <- 0

genome.dist.moor_house <- dist(t(genome.zs.matrix.moor_house))

#### Perform the hierarchical clustering
genome.hc.moor_house <- hclust(as.dist(genome.dist.moor_house), method = "ward.D")
genome.ord.moor_house <- genome.hc.moor_house$labels[genome.hc.moor_house$order]
genome.ord.moor_house <- data.frame(Genome = genome.ord.moor_house, order = 1:length(genome.ord.moor_house))
genome.cut.moor_house <- cutree(genome.hc.moor_house[c(1,2,4)],k = 3) # k = 3 for three treatments
genome.clusters.moor_house <- data.frame(Genome = names(genome.cut.moor_house),
                                         Cluster = genome.cut.moor_house) %>%
  inner_join(genome.ord.moor_house, by = "Genome") %>%
  inner_join(genome.zs.moor_house, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.moor_house <- genome.clusters.moor_house %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM", "REST")))
genome.clusters.moor_house$Title <- paste("Group: ", genome.clusters.moor_house$Cluster, " - ", genome.clusters.moor_house$Cluster.size, " genomes", sep = "")

head(genome.clusters.moor_house)

#### Plot Z-Score clusters
plot.clusters.moor_house <-
  ggplot(data = genome.clusters.moor_house, aes(x = factor(treatment, levels = c("NAT", "DAM", "REST")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged", "Restored")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Moor House")
plot.clusters.moor_house

#### Manually rename cluster titles based on abundance trends
genome.clusters.moor_house <- genome.clusters.moor_house %>%
  mutate(Title = case_when(Cluster == 1 ~ "Natural-abundant",
                           Cluster == 2 ~ "Damaged-abundant",
                           Cluster == 3 ~ "Restored-abundant"))
genome.clusters.moor_house$site <- "Moor_House"
head(genome.clusters.moor_house)

### Stean
#### Calculate the zscores of each genome across samples and then calculate the mean zvalue for each treatment
genome.zs.stean <- tmeans_sig_stean %>%
  filter(row.names(.) %in% sig_stean$Genome) %>% # Only include genomes significant from DESeq
  rel_ab() %>% 
  tidy_otu %>%
  group_by(Genome) %>%
  mutate(zValue = (Count - mean(Count))/sd(Count)) %>%
  inner_join(metadata, by = "Sample") %>%
  group_by(treatment, Genome) %>%
  summarise(MeanZS = mean(zValue)) 

#### Format zscores as a matrix
genome.zs.matrix.stean <- genome.zs.stean %>%
  spread(key = Genome, value = MeanZS) %>%
  as.data.frame()
row.names(genome.zs.matrix.stean) <- genome.zs.matrix.stean$Treatment
genome.zs.matrix.stean <- genome.zs.matrix.stean[,-1]
genome.zs.matrix.stean <- as.matrix(genome.zs.matrix.stean)

# Check for NA, NaN, or Inf values and handle them
genome.zs.matrix.stean[is.na(genome.zs.matrix.stean)] <- 0
genome.zs.matrix.stean[is.nan(genome.zs.matrix.stean)] <- 0
genome.zs.matrix.stean[is.infinite(genome.zs.matrix.stean)] <- 0

genome.dist.stean <- dist(t(genome.zs.matrix.stean))

#### Perform the hierarchical clustering
genome.hc.stean <- hclust(as.dist(genome.dist.stean), method = "ward.D")
genome.ord.stean <- genome.hc.stean$labels[genome.hc.stean$order]
genome.ord.stean <- data.frame(Genome = genome.ord.stean, order = 1:length(genome.ord.stean))
genome.cut.stean <- cutree(genome.hc.stean[c(1,2,4)],k = 2) # k = 2 for TWO treatments (stean has no Restored treatment)
genome.clusters.stean <- data.frame(Genome = names(genome.cut.stean),
                                    Cluster = genome.cut.stean) %>%
  inner_join(genome.ord.stean, by = "Genome") %>%
  inner_join(genome.zs.stean, by = "Genome") %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM")))
genome.clusters.stean <- genome.clusters.stean %>%
  group_by(Cluster) %>%
  mutate(Cluster.size = n_distinct(Genome)) %>%
  mutate(Title = factor(treatment, levels = c("NAT", "DAM")))
genome.clusters.stean$Title <- paste("Group: ", genome.clusters.stean$Cluster, " - ", genome.clusters.stean$Cluster.size, " genomes", sep = "")

head(genome.clusters.stean)

#### Plot Z-Score clusters
plot.clusters.stean <-
  ggplot(data = genome.clusters.stean, aes(x = factor(treatment, levels = c("NAT", "DAM")), y = MeanZS)) +
  geom_boxplot(aes(color = treatment), outlier.shape = NA) +
  geom_smooth(aes(group = 1), method = "lm", color="black", se = F, formula = y ~ poly(x, 2),) +
  xlab("Restoration status") +
  ylab("Abundance\n(mean Z-Score)") +
  scale_x_discrete(labels = c("Natural", "Damaged")) +
  facet_wrap(~Title, ncol = 3, labeller = labeller(label_column = as.character)) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(-1,1) +
  ggtitle("Stean")
plot.clusters.stean

#### Manually rename cluster titles based on abundance trends
genome.clusters.stean <- genome.clusters.stean %>%
  mutate(Title = case_when(Cluster == 1 ~ "Damaged-abundant",
                           Cluster == 2 ~ "Natural-abundant"))
genome.clusters.stean$site <- "Stean"
head(genome.clusters.stean)

### Combine cluster dataframes into one
genome.clusters <- rbind(genome.clusters.balmoral,
                         genome.clusters.bowness,
                         genome.clusters.crocach,
                         genome.clusters.langwell,
                         genome.clusters.migneint,
                         genome.clusters.moor_house,
                         genome.clusters.stean)
head(genome.clusters)
unique(genome.clusters$site)

### Combine cluster barplots
plot.clusters.combined <- cowplot::plot_grid(plot.clusters.balmoral,
                                             plot.clusters.bowness,
                                             plot.clusters.crocach,
                                             plot.clusters.langwell,
                                             plot.clusters.migneint,
                                             plot.clusters.moor_house,
                                             plot.clusters.stean,
                                             nrow = 4,
                                             ncol = 2,
                                             label_size = 16,
                                             label_fontfamily = "sans",
                                             label_fontface = "bold")
ggsave(plot = plot.clusters.combined,
       filename = "./trend_groups.png",
       device = "png",
       dpi = 600,
       width = 10,
       height = 10,
       units = "in",
       bg = "white")
plot.clusters.combined

### Merge and Save DESeq2 Results and Cluster Results
combined_deseq_clusters <- combined_deseq_results %>%
  full_join(genome.clusters, by = join_by("Genome", "site"))
head(combined_deseq_clusters)
write_csv(combined_deseq_clusters, file = "./deseq_and_clusters_combined.csv")

# Combined genome clusters, metabolic functions, abundance, and taxonomy
## Genome abundance with taxonomy and metabolic functions
genome_clusters <- genome.clusters %>%
  group_by(Genome) %>%
  select(Genome, site,  Title) %>%
  dplyr::rename(c(`Trend Group` = Title)) %>%
  distinct()

abundance_long_all <- tmeans_norm %>%
  rownames_to_column(var = "Genome") %>%
  gather(key = "Sample", value = "Abundance", -Genome) %>%
  left_join(metadata %>%
              select(Sample, site, treatment))  %>%
  left_join(genome_clusters, by = join_by("Genome", "site")) %>% # Cluster assignment isn't the same across all sites!
  left_join(gtdb) %>%
  filter(Abundance > 0) # 0 is non-detection not actual 0

abundance_long_all$Abundance[is.na(abundance_long_all$Abundance)] <- 0

# List the columns from "Domain" through "Genus"
tax_columns <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Replace NA values with "Unknown" and refactor to make "Unknown" the last level
abundance_long_all <- abundance_long_all %>%
  mutate(across(all_of(tax_columns), ~ ifelse(is.na(.), "Unknown", .))) %>%
  mutate(across(all_of(tax_columns), ~ factor(., levels = c(setdiff(levels(factor(.)), "Unknown"), "Unknown"))))

head(abundance_long_all)

## Merge genome abundance data with metabolic predictions
abundance_long_all_and_metabolism <- abundance_long_all %>%
  inner_join(metabolic_filtered_long, by = "Genome") %>%
  filter(Abundance > 0) # 0 is non-detection not actual 0

write_csv(abundance_long_all_and_metabolism, file = "./genome_abundance_annotated.csv")
head(abundance_long_all_and_metabolism)

# Final, simplified table
# Just genomes, their cluster/trend group in each site (if available), and functions

genomes.annotated <- metabolic_filtered_long %>%
  full_join(combined_deseq_clusters %>%
              select(Genome, site, Title) %>%
              distinct()) %>%
  dplyr::rename(`Trend group` = Title)

write_csv(genomes.annotated, file = "./genomes_annotated_simple.csv")
head(genomes.annotated)

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