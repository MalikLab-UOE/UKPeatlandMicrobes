################################################################################
##  Figure 4 & 5  —  DESeq2 (VST) 
################################################################################

library(tidyverse)
library(DESeq2)
library(compositions)   # CLR transform (Figure4_5_merged_clr output set)
library(emmeans)
library(scales)
library(ggpubr)
library(KEGGREST)

## ---------------------------------------------------------------------------
## 0. Global settings: output root, site palette, shared labels, site schemes
## ---------------------------------------------------------------------------
new_root <- "Figure4_5_merged_new"
dir.create(new_root, showWarnings = FALSE)

# Second output set for figure stitching: y-axis shows the pathway name only
# (no "(Abundance VST value)" line), and the per-Site facet is abundance (y)
# vs Ecosystem health index (x) instead of growth-vs-abundance.
nolabel_root <- "Figure4_5_merged_nolabel"
dir.create(nolabel_root, showWarnings = FALSE)

# Third output set: identical layout to Figure4_5_merged_new (two-line y label,
# growth-vs-abundance facet) but the abundance values are CLR-transformed
# (centred-log-ratio on prodigal-normalised counts) instead of VST. Pathway
# definitions are unchanged (de-oxygenated KO sets for N/S/Methane).
clr_root <- "Figure4_5_merged_clr"
dir.create(clr_root, showWarnings = FALSE)

# Fourth output set: CLR values with the nolabel y-axis style (plain pathway
# name only, abundance-vs-eco-index facet) for stitching panels together.
clr_nolabel_root <- "Figure4_5_merged_clr_nolabel"
dir.create(clr_nolabel_root, showWarnings = FALSE)

# Prodigal normalisation constant (shared by the KEGG & CAZyme CLR matrices;
# a per-sample linear scale factor that does not affect the CLR geometry).
prodigal_scale <- 57474057.89

# Named palette so each site keeps the same colour regardless of which sites are
# present in a given scheme (avoids colour shifting when points are dropped).
site_colors <- c(
  "Balmoral"    = "#e41a1c",
  "Bowness"     = "#377eb8",
  "Crocach"     = "#4daf4a",
  "Langwell"    = "#984ea3",
  "Migneint"    = "#ff7f00",
  "Moors_House" = "#ffff33",
  "Stean"       = "#a65628"
)

# Display names for sites (labels only; data/palette keys stay as-is).
# "Moors_House" is shown as "Moor House" (no underscore) in facet strips & legend.
site_labels <- c(
  "Balmoral"    = "Balmoral",
  "Bowness"     = "Bowness",
  "Crocach"     = "Crocach",
  "Langwell"    = "Langwell",
  "Migneint"    = "Migneint",
  "Moors_House" = "Moor House",
  "Stean"       = "Stean"
)

# Facet strip order (left -> right) for the per-Site facet panels.
# Uses data keys; display names come from site_labels. Any site absent from a
# scheme is dropped automatically by facet_grid(drop = TRUE).
facet_site_order <- c("Balmoral", "Bowness", "Moors_House", "Langwell", "Crocach", "Migneint")

# Shared growth-rate y-axis label
growth_ylab <- expression("Microbial growth rate (ng g"^{-1}*"h"^{-1}*")")

# Two site schemes. Stean is ALWAYS excluded (removed at data load via grep on
# the "SE" Sample_ID prefix), so it is not listed here.
schemes <- list(
  list(tag = "all_sites",     drop = character(0),                label = "Stean removed"),
  list(tag = "no_MH_Crocach", drop = c("Moors_House", "Crocach"), label = "Stean + Moors_House + Crocach removed")
)

## ---------------------------------------------------------------------------
## 1. Data loading (Stean removed here — the "SE" Sample_ID prefix is Stean)
## ---------------------------------------------------------------------------
kegg_abund_raw <- read.csv("Data/Figure4&5/kegg_abund_raw_t.csv")
Sample_id      <- read.csv("Data/Figure4&5/Sample_ID.csv")
Sample_id      <- Sample_id[-grep("SE", Sample_id$Sample_ID), ]   # drop Stean
env_data       <- read.csv("Data/Figure4&5/env_data.csv")
ko_levels      <- read.csv("Data/Figure4&5/KO_levels.csv")

cazy_abund_raw <- read.csv("Data/Figure4&5/cazy_abund_raw.csv")
cazy_substrate <- read.csv("Data/Figure4&5/cazyme_by_substrate.csv")

# Environment factors used by every figure
env_subset <- env_data %>% select(Sample, eco_index, Growth)

## ---------------------------------------------------------------------------
## 2. De-oxygenated KO lists for Nitrogen / Sulfur / Methane (live from KEGG)
## ---------------------------------------------------------------------------
fetch_module_kos <- function(module_ids, label) {
  cat(paste0("Fetching ", label, " KO list from KEGG...\n"))
  raw <- list()
  for (id in module_ids) {
    tryCatch({
      res <- keggGet(id)
      raw[[id]] <- names(res[[1]]$ORTHOLOGY)
    }, error = function(e) warning(paste("Could not fetch module:", id)))
  }
  unique(unlist(raw, use.names = FALSE))
}

nitrogen_modules <- c("M00175", "M00531", "M00530", "M00804", "M00973")
sulfur_modules   <- c("M00987", "M00176", "M00596", "M00986",
                      "M00990", "M00991", "M00992", "M00993")
methane_modules  <- c("M00567", "M00357", "M00356", "M00563", "M00358",
                      "M00608", "M00346", "M00345", "M00344", "M00378",
                      "M00935", "M00422")

nitrogen_kos <- fetch_module_kos(nitrogen_modules, "Nitrogen metabolism")
sulfur_kos   <- fetch_module_kos(sulfur_modules,   "Sulfur (oxidation removed)")
methane_kos  <- fetch_module_kos(methane_modules,  "Methane (oxidation removed)")

# Fermentation KOs (original script)
ferm_kos <- c("K00163", "K00627", "K00169", "K03737", "K00174", "K00656")

## ---------------------------------------------------------------------------
## 3. VST normalisation — KEGG
## ---------------------------------------------------------------------------
counts_mat <- kegg_abund_raw %>% column_to_rownames("Sample") %>% t() %>% round()
colData    <- Sample_id %>% filter(Sample %in% colnames(counts_mat)) %>% column_to_rownames("Sample")
dds        <- DESeqDataSetFromMatrix(countData = counts_mat[, rownames(colData)],
                                     colData = colData, design = ~ Site + Treatment)
vsd        <- varianceStabilizingTransformation(dds, blind = FALSE)
kegg_vst_df <- as.data.frame(t(assay(vsd))) %>% rownames_to_column("Sample") %>%
  left_join(Sample_id, by = "Sample")

kegg_long_vst <- kegg_vst_df %>%
  gather(key = "KO", value = "vst_val", starts_with("K")) %>%
  left_join(ko_levels, by = "KO")

## ---------------------------------------------------------------------------
## 4. VST normalisation — CAZyme
## ---------------------------------------------------------------------------
cazy_counts   <- cazy_abund_raw %>% column_to_rownames("Sample") %>% round()
colData_caz   <- Sample_id %>% filter(Sample %in% colnames(cazy_counts)) %>% column_to_rownames("Sample")
cazy_counts   <- cazy_counts[, rownames(colData_caz)]
dds_caz       <- DESeqDataSetFromMatrix(countData = cazy_counts, colData = colData_caz,
                                        design = ~ Site + Treatment)
vsd_caz       <- varianceStabilizingTransformation(dds_caz, blind = FALSE)
cazy_vst_df   <- as.data.frame(t(assay(vsd_caz))) %>% rownames_to_column("Sample") %>%
  left_join(Sample_id, by = "Sample")

## ---------------------------------------------------------------------------
## 5. Build per-pathway summed-VST data frames (joined with env factors)
## ---------------------------------------------------------------------------
sum_by_kos <- function(kos) {
  kegg_vst_df %>%
    select(Sample, Site, Treatment, any_of(kos)) %>%
    mutate(sum_vst = rowSums(select(., any_of(kos)))) %>%
    select(Sample, Site, Treatment, sum_vst)
}
sum_by_cat3 <- function(cat) {
  kegg_long_vst %>% filter(Category3 == cat) %>%
    group_by(Sample, Site, Treatment) %>%
    summarise(sum_vst = sum(vst_val), .groups = "drop")
}
sum_by_cazy <- function(sub_name) {
  sub_genes <- cazy_substrate %>% filter(Substrate == sub_name) %>% pull(cazyme)
  cazy_vst_df %>%
    select(Sample, Site, Treatment, any_of(sub_genes)) %>%
    mutate(sum_vst = rowSums(select(., any_of(sub_genes)))) %>%
    select(Sample, Site, Treatment, sum_vst)
}

# KEGG pathway definitions
kegg_defs <- list(
  list(name = "Fermentation",        tag = "Ferm",           data = sum_by_kos(ferm_kos)),
  list(name = "Aromatic compounds",  tag = "Aromatic",       data = sum_by_cat3("01220_Degradation_of_aromatic_compounds_[PATH:ko01220]")),
  list(name = "Carbon fixation",     tag = "CarbonFixation", data = sum_by_cat3("00720_Carbon_fixation_pathways_in_prokaryotes_[PATH:ko00720]")),
  list(name = "Nitrogen metabolism", tag = "Nitrogen",       data = sum_by_kos(nitrogen_kos)),
  list(name = "Sulfur metabolism",   tag = "Sulfur",         data = sum_by_kos(sulfur_kos)),
  list(name = "Methane metabolism",  tag = "Methane",        data = sum_by_kos(methane_kos))
)

# CAZyme substrate definitions
cazy_defs <- list(
  list(name = "Microbial cell wall", tag = "CellWall",        data = sum_by_cazy("Cell Wall")),
  list(name = "Lignin",              tag = "Lignin",          data = sum_by_cazy("Lignin")),
  list(name = "Oligosaccharides",    tag = "Oligosaccharides",data = sum_by_cazy("Oligosaccharides")),
  list(name = "Cellulose",           tag = "Cellulose",       data = sum_by_cazy("Cellulose")),
  list(name = "Hemicellulose",       tag = "Hemicellulose",   data = sum_by_cazy("Hemicellulose"))
)

# Attach env factors once
attach_env <- function(defs) lapply(defs, function(d) { d$merged <- d$data %>% inner_join(env_subset, by = "Sample"); d })
kegg_defs <- attach_env(kegg_defs)
cazy_defs <- attach_env(cazy_defs)

## ---------------------------------------------------------------------------
## 5b. CLR normalisation — KEGG + CAZyme (parallel to the VST pipeline above)
##     Prodigal-normalise raw counts, then clr(mat + 0.5) per sample (compositions
##     package). Pathway/substrate gene sets are IDENTICAL to the VST pipeline
##     (de-oxygenated KO lists for N/S/Methane). The summed column is named
##     `sum_vst` so all downstream plotting code is reused unchanged.
## ---------------------------------------------------------------------------
prodigal_kegg <- read.csv("Data/Figure4&5/prodigal_output_kegg.csv")
prodigal_cazy <- read.csv("Data/Figure4&5/prodigal_output_caz_abund.csv")

# --- KEGG CLR matrix (Sample x KO) ---
kegg_norm_wide <- kegg_abund_raw %>%
  left_join(prodigal_kegg, by = "Sample") %>%
  gather(key = "KO", value = "reads", starts_with("K")) %>%
  mutate(reads = as.numeric(reads),
         reads_norm = (reads / prodigal_output) * prodigal_scale) %>%
  select(Sample, KO, reads_norm) %>%
  spread(KO, reads_norm) %>%
  filter(Sample %in% Sample_id$Sample) %>%          # drop Stean
  column_to_rownames("Sample")

kegg_clr_df <- as.data.frame(clr(kegg_norm_wide + 0.5)) %>%
  rownames_to_column("Sample") %>%
  left_join(Sample_id, by = "Sample")
kegg_long_clr <- kegg_clr_df %>%
  gather(key = "KO", value = "clr_val", starts_with("K")) %>%
  left_join(ko_levels, by = "KO")

# --- CAZyme CLR matrix (Sample x cazyme) ---
cazy_raw_clr <- cazy_abund_raw
names(cazy_raw_clr)[1] <- "cazyme_id"
cazy_norm_wide <- cazy_raw_clr %>%
  gather(key = "Sample", value = "reads", -cazyme_id) %>%
  mutate(reads = as.numeric(reads)) %>%
  left_join(prodigal_cazy, by = "Sample") %>%
  mutate(reads_norm = (reads / prodigal_output) * prodigal_scale) %>%
  select(Sample, cazyme_id, reads_norm) %>%
  spread(key = cazyme_id, value = reads_norm) %>%
  filter(Sample %in% Sample_id$Sample) %>%          # drop Stean
  column_to_rownames("Sample")

cazy_clr_df <- as.data.frame(clr(cazy_norm_wide + 0.5)) %>%
  rownames_to_column("Sample") %>%
  left_join(Sample_id, by = "Sample")

# Per-pathway CLR summing (mirrors the VST sum_by_* helpers; column = sum_vst)
sum_by_kos_clr <- function(kos) {
  kegg_clr_df %>%
    select(Sample, Site, Treatment, any_of(kos)) %>%
    mutate(sum_vst = rowSums(select(., any_of(kos)))) %>%
    select(Sample, Site, Treatment, sum_vst)
}
sum_by_cat3_clr <- function(cat) {
  kegg_long_clr %>% filter(Category3 == cat) %>%
    group_by(Sample, Site, Treatment) %>%
    summarise(sum_vst = sum(clr_val), .groups = "drop")
}
sum_by_cazy_clr <- function(sub_name) {
  sub_genes <- cazy_substrate %>% filter(Substrate == sub_name) %>% pull(cazyme)
  cazy_clr_df %>%
    select(Sample, Site, Treatment, any_of(sub_genes)) %>%
    mutate(sum_vst = rowSums(select(., any_of(sub_genes)))) %>%
    select(Sample, Site, Treatment, sum_vst)
}

kegg_defs_clr <- list(
  list(name = "Fermentation",        tag = "Ferm",           data = sum_by_kos_clr(ferm_kos)),
  list(name = "Aromatic compounds",  tag = "Aromatic",       data = sum_by_cat3_clr("01220_Degradation_of_aromatic_compounds_[PATH:ko01220]")),
  list(name = "Carbon fixation",     tag = "CarbonFixation", data = sum_by_cat3_clr("00720_Carbon_fixation_pathways_in_prokaryotes_[PATH:ko00720]")),
  list(name = "Nitrogen metabolism", tag = "Nitrogen",       data = sum_by_kos_clr(nitrogen_kos)),
  list(name = "Sulfur metabolism",   tag = "Sulfur",         data = sum_by_kos_clr(sulfur_kos)),
  list(name = "Methane metabolism",  tag = "Methane",        data = sum_by_kos_clr(methane_kos))
)
cazy_defs_clr <- list(
  list(name = "Microbial cell wall", tag = "CellWall",        data = sum_by_cazy_clr("Cell Wall")),
  list(name = "Lignin",              tag = "Lignin",          data = sum_by_cazy_clr("Lignin")),
  list(name = "Oligosaccharides",    tag = "Oligosaccharides",data = sum_by_cazy_clr("Oligosaccharides")),
  list(name = "Cellulose",           tag = "Cellulose",       data = sum_by_cazy_clr("Cellulose")),
  list(name = "Hemicellulose",       tag = "Hemicellulose",   data = sum_by_cazy_clr("Hemicellulose"))
)
kegg_defs_clr <- attach_env(kegg_defs_clr)
cazy_defs_clr <- attach_env(cazy_defs_clr)

## ---------------------------------------------------------------------------
## 6. Plotting function: boxplot + LM(eco) + growth facet + stats
## ---------------------------------------------------------------------------
make_pathway_plots <- function(merged_all, plot_name, file_tag, category, scheme,
                               out_root = new_root, label_mode = "vst",
                               value_mode = "vst") {
  out_dir <- file.path(out_root, category, scheme$tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- merged_all %>% filter(!Site %in% scheme$drop)
  # Order the per-Site facet strips (drop levels not present in this scheme).
  dat$Site <- factor(dat$Site, levels = intersect(facet_site_order, unique(dat$Site)))

  # Transform wording for the value axis: "VST value" or "CLR value".
  value_word <- if (value_mode == "clr") "CLR value" else "VST value"

  # y-axis label = "<pathway> gene abundance". "vst" keeps the (<transform> value)
  # second line; "plain" (for stitching multiple panels) drops that wording.
  y_label_expr <- if (label_mode == "vst") {
    bquote(atop(.(paste(plot_name, "gene abundance")), .(paste0("(", value_word, ")"))))
  } else {
    paste(plot_name, "gene abundance")
  }

  # Facet plots use the full pathway name for Carbon fixation
  facet_name <- if (plot_name == "Carbon fixation") "Carbon fixation in prokaryotes" else plot_name
  # Abundance-axis label on the per-Site facet panels: plain pathway name, matching
  # the original Figure4_5_merged_nolabel version (no "gene abundance" here — the
  # "gene abundance" wording stays only on the Boxplot / LM_Eco panels).
  facet_abund_lab <- facet_name

  # A. Boxplot by treatment
  p1 <- ggplot(dat, aes(x = Treatment, y = sum_vst)) +
    geom_boxplot(lwd = 1, outlier.shape = NA) +
    geom_jitter(aes(color = Site), width = 0.3, alpha = 0.7, size = 3) +
    ylab(y_label_expr) + xlab(" ") + theme_light(base_size = 14) +
    theme(legend.position = "none", aspect.ratio = 1,
          axis.text = element_text(size = 13, colour = "black")) +
    scale_x_discrete(limits = c("DAM", "REST", "NAT"),
                     labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    scale_colour_manual(values = site_colors)
  ggsave(file.path(out_dir, paste0("Boxplot_", file_tag, ".png")), p1,
         dpi = 1000, width = 4, height = 4, units = "in")

  # B. LM vs ecosystem health index
  p2 <- ggplot(dat, aes(x = eco_index, y = sum_vst)) +
    geom_point(aes(fill = Site, shape = Treatment), size = 4, alpha = 0.7, color = "black") +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed", size = 0.7) +
    scale_shape_manual(values = c(21, 22, 23), limits = c("NAT", "REST", "DAM")) +
    scale_fill_manual(values = site_colors) +
    ylab(y_label_expr) + xlab("Ecosystem health index") + theme_light(base_size = 14) +
    theme(legend.position = "none", axis.text = element_text(size = 13, colour = "black"))
  ggsave(file.path(out_dir, paste0("LM_Eco_", file_tag, ".png")), p2,
         dpi = 1000, width = 4, height = 4, units = "in")

  # C. Faceted-by-Site panel.
  #    - "vst"   : growth (y, log2) vs abundance (x)  [original template]
  #    - "plain" : abundance (y) vs Ecosystem health index (x)
  if (label_mode == "vst") {
    p3 <- ggplot(dat, aes(x = sum_vst, y = Growth)) +
      geom_point(color = "black", size = 4, alpha = 0.7,
                 aes(shape = as.character(Treatment), fill = as.character(Site))) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed",
                  size = 0.7, alpha = 0.3) +
      scale_y_continuous(trans = "log2", labels = label_number(accuracy = 0.1), limits = c(0.05, 11)) +
      scale_fill_manual(values = site_colors) +
      scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                         labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
      xlab(facet_abund_lab) +
      (if (value_mode == "clr")
         scale_x_continuous(breaks = scales::breaks_pretty(n = 3))
       else
         scale_x_continuous(labels = function(x) x / 1000, breaks = scales::breaks_pretty(n = 3))) +
      ylab(growth_ylab) +
      guides(fill = "none", shape = "none", color = "none") +
      theme_light(base_size = 14) +
      theme(axis.text = element_text(size = 13, colour = "black"),
            axis.text.x = element_text(size = 10), legend.position = "none") +
      facet_grid(. ~ Site, scales = "free", labeller = labeller(Site = site_labels))
    ggsave(file.path(out_dir, paste0("Growth_facet_", file_tag, ".png")), p3,
           dpi = 1000, width = 8, height = 3, units = "in")
  } else {
    p3 <- ggplot(dat, aes(x = eco_index, y = sum_vst)) +
      geom_point(color = "black", size = 4, alpha = 0.7,
                 aes(shape = as.character(Treatment), fill = as.character(Site))) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed",
                  size = 0.7, alpha = 0.3) +
      scale_fill_manual(values = site_colors) +
      scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                         labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
      xlab("Ecosystem health index") +
      ylab(facet_abund_lab) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
      guides(fill = "none", shape = "none", color = "none") +
      theme_light(base_size = 14) +
      theme(axis.text = element_text(size = 13, colour = "black"),
            axis.text.x = element_text(size = 9), legend.position = "none") +
      facet_grid(. ~ Site, scales = "free", labeller = labeller(Site = site_labels))
    ggsave(file.path(out_dir, paste0("Facet_Eco_", file_tag, ".png")), p3,
           dpi = 1000, width = 8, height = 3, units = "in")
  }

  # D. Stats: treatment contrasts (Tukey)
  emm <- emmeans(lm(sum_vst ~ Treatment, data = dat), ~ Treatment)
  write.csv(as.data.frame(pairs(emm, adjust = "tukey")),
            file.path(out_dir, paste0("Stats_Boxplot_", file_tag, ".csv")), row.names = FALSE)

  # E. Stats: LM vs eco_index
  lm_fit <- summary(lm(sum_vst ~ eco_index, data = dat))
  write.csv(data.frame(R2 = lm_fit$r.squared, Adj_R2 = lm_fit$adj.r.squared,
                       P_val = lm_fit$coefficients[2, 4]),
            file.path(out_dir, paste0("Stats_LM_Eco_", file_tag, ".csv")), row.names = FALSE)
}

## ---------------------------------------------------------------------------
## 7. Growth-correlation function: Growth vs pathway abundance + lm statistics
## ---------------------------------------------------------------------------
make_growth_corr <- function(merged_all, plot_name, x_axis_label, file_tag, scheme,
                             out_root = new_root, value_mode = "vst") {
  out_dir <- file.path(out_root, "Growth_correlation", scheme$tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- merged_all %>% filter(!Site %in% scheme$drop)

  p <- ggplot(dat, aes(x = sum_vst, y = Growth)) +
    geom_point(color = "black", size = 4, alpha = 0.7,
               aes(shape = as.character(Treatment), fill = as.character(Site))) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed",
                size = 0.7, alpha = 0.3) +
    scale_y_continuous(trans = "log2", labels = label_number(accuracy = 0.1), limits = c(0.05, 11)) +
    scale_fill_manual(values = site_colors) +
    scale_shape_manual(name = "Treatment", values = c(21, 22, 23), limits = c("NAT", "REST", "DAM"),
                       labels = c("DAM" = "Degraded", "NAT" = "Natural", "REST" = "Restored")) +
    xlab(x_axis_label) +
    (if (value_mode == "clr")
       scale_x_continuous(breaks = scales::breaks_pretty(n = 4))
     else
       scale_x_continuous(labels = function(x) x / 1000, breaks = scales::breaks_pretty(n = 4))) +
    ylab(growth_ylab) +
    guides(fill = "none", shape = "none", color = "none") +
    theme_light(base_size = 14) +
    theme(axis.text = element_text(size = 13, colour = "black"), legend.position = "none")
  ggsave(file.path(out_dir, paste0("Growth_vs_", file_tag, ".png")), p,
         dpi = 1000, width = 5, height = 4, units = "in")

  # lm on raw Growth and on log2(Growth) (matches the log2 display axis), plus Pearson
  fit_raw  <- summary(lm(Growth ~ sum_vst, data = dat))
  fit_log2 <- summary(lm(log2(Growth) ~ sum_vst, data = dat))
  pear     <- cor.test(dat$sum_vst, dat$Growth)
  stats <- data.frame(
    Pathway   = plot_name,
    Scheme    = scheme$tag,
    n         = nrow(dat),
    R2_raw    = fit_raw$r.squared,
    Adj_R2_raw= fit_raw$adj.r.squared,
    P_raw     = fit_raw$coefficients[2, 4],
    R2_log2   = fit_log2$r.squared,
    P_log2    = fit_log2$coefficients[2, 4],
    Pearson_r = unname(pear$estimate),
    Pearson_p = pear$p.value
  )
  write.csv(stats, file.path(out_dir, paste0("Stats_Growth_", file_tag, ".csv")), row.names = FALSE)
  stats
}

## ---------------------------------------------------------------------------
## 8. Generate all figures — both site schemes, KEGG + CAZyme
## ---------------------------------------------------------------------------
# Two label modes: "vst" -> Figure4_5_merged_new (original, two-line y label,
# growth-vs-abundance facet); "plain" -> Figure4_5_merged_nolabel (pathway-name
# y label only, abundance-vs-eco-index facet, for stitching).
## Each output set names its abundance transform (value_mode: "vst" or "clr"),
## its y-label style (label_mode), and which per-pathway def lists to draw from.
output_sets <- list(
  list(root = new_root,     mode = "vst",   value_mode = "vst", kegg = kegg_defs,     cazy = cazy_defs),
  list(root = nolabel_root, mode = "plain", value_mode = "vst", kegg = kegg_defs,     cazy = cazy_defs),
  list(root = clr_root,         mode = "vst",   value_mode = "clr", kegg = kegg_defs_clr, cazy = cazy_defs_clr),
  list(root = clr_nolabel_root, mode = "plain", value_mode = "clr", kegg = kegg_defs_clr, cazy = cazy_defs_clr)
)

for (oset in output_sets) {
  cat(paste0("\n########## Output set: ", oset$root, " (label_mode = ", oset$mode,
             ", value_mode = ", oset$value_mode, ") ##########\n"))
  for (scheme in schemes) {
    cat(paste0("\n=== Scheme: ", scheme$tag, " (", scheme$label, ") ===\n"))
    for (d in oset$kegg) {
      cat(paste0("  KEGG   : ", d$name, "\n"))
      make_pathway_plots(d$merged, d$name, d$tag, "KEGG", scheme,
                         out_root = oset$root, label_mode = oset$mode,
                         value_mode = oset$value_mode)
    }
    for (d in oset$cazy) {
      cat(paste0("  CAZyme : ", d$name, "\n"))
      make_pathway_plots(d$merged, d$name, d$tag, "CAZyme", scheme,
                         out_root = oset$root, label_mode = oset$mode,
                         value_mode = oset$value_mode)
    }
  }
}

## ---------------------------------------------------------------------------
## 9. Growth-correlation section: Carbon fixation / Methane / Fermentation
## ---------------------------------------------------------------------------
get_merged <- function(tag, defs = kegg_defs) Filter(function(d) d$tag == tag, defs)[[1]]$merged

corr_targets <- list(
  list(tag = "CarbonFixation", plot_name = "Carbon fixation",     xlab = "Carbon fixation in prokaryotes"),
  list(tag = "Methane",        plot_name = "Methane metabolism",  xlab = "Methane metabolism"),
  list(tag = "Ferm",           plot_name = "Fermentation",        xlab = "Fermentation")
)

# Growth-correlation plots have no VST wording (y = growth rate, x = pathway
# name), so write an identical, self-contained copy into each output set.
for (oset in output_sets) {
  growth_corr_summary <- list()
  for (scheme in schemes) {
    for (ct in corr_targets) {
      s <- make_growth_corr(get_merged(ct$tag, oset$kegg), ct$plot_name, ct$xlab, ct$tag, scheme,
                            out_root = oset$root, value_mode = oset$value_mode)
      growth_corr_summary[[length(growth_corr_summary) + 1]] <- s
    }
  }
  growth_corr_summary <- bind_rows(growth_corr_summary)
  write.csv(growth_corr_summary,
            file.path(oset$root, "Growth_correlation", "Growth_correlation_summary_all.csv"),
            row.names = FALSE)
}

## ---------------------------------------------------------------------------
## 10. Export KEGG module table used (after removing oxidation steps)
## ---------------------------------------------------------------------------
## A single comprehensive supplementary table of everything analysed, with a
## description for each entry. Three kinds of definition are covered:
##   - "KEGG module"  : Nitrogen / Sulfur / Methane (oxidation modules removed)
##   - "KEGG pathway" : Aromatic compounds (ko01220) & Carbon fixation (ko00720),
##                      summed over all KOs annotated to the pathway
##   - "KO (gene)"    : Fermentation, defined by an explicit KO list
kegg_name <- function(id) {
  nm <- tryCatch(keggGet(id)[[1]]$NAME, error = function(e) NA)
  paste(nm, collapse = "; ")
}

# Module-based pathways (Nitrogen / Sulfur / Methane)
module_rows <- bind_rows(
  data.frame(Pathway = "Nitrogen metabolism", Definition_type = "KEGG module", ID = nitrogen_modules),
  data.frame(Pathway = "Sulfur metabolism",   Definition_type = "KEGG module", ID = sulfur_modules),
  data.frame(Pathway = "Methane metabolism",  Definition_type = "KEGG module", ID = methane_modules)
) %>%
  mutate(Description = sapply(ID, kegg_name))

# Pathway-based categories (summed over all KOs of the KEGG pathway)
pathway_rows <- data.frame(
  Pathway = c("Aromatic compounds", "Carbon fixation"),
  Definition_type = "KEGG pathway",
  ID = c("ko01220", "ko00720"),
  Description = c("Degradation of aromatic compounds",
                  "Carbon fixation pathways in prokaryotes")
)

# KO-list category (Fermentation) — explicit gene list from the script
ferm_rows <- data.frame(
  Pathway = "Fermentation",
  Definition_type = "KO (gene)",
  ID = ferm_kos,
  Description = sapply(ferm_kos, kegg_name)
)

kegg_module_summary <- bind_rows(ferm_rows, pathway_rows, module_rows)
write.csv(kegg_module_summary,
          file.path(new_root, "KEGG_modules_genes_used.csv"),
          row.names = FALSE)
# Keep the original filename as well (now with the Description column) so any
# existing references to it still resolve.
write.csv(kegg_module_summary,
          file.path(new_root, "KEGG_modules_used_after_removing_oxidation.csv"),
          row.names = FALSE)
# Same gene/module definitions drive the CLR output set — copy the table there.
write.csv(kegg_module_summary,
          file.path(clr_root, "KEGG_modules_genes_used.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 11. Standalone legends (one per site scheme) for stitched figures
##       - all_sites     : Stean removed        -> 6 sites
##       - no_MH_Crocach : Stean + MH + Crocach -> 4 sites
##     Each legend shows Site (colour) + Treatment (shape). Saved into both
##     output roots and also at the top level for convenience.
## ---------------------------------------------------------------------------
make_scheme_legend <- function(scheme) {
  # Reference data (any pathway works — Fermentation) filtered by this scheme
  ref <- get_merged("Ferm") %>% filter(!Site %in% scheme$drop)

  # Keep sites in the fixed palette order, restricted to those present
  sites_present <- names(site_colors)[names(site_colors) %in% unique(ref$Site)]
  ref$Site      <- factor(ref$Site, levels = sites_present)
  ref$Treatment <- factor(ref$Treatment, levels = c("NAT", "REST", "DAM"))

  p <- ggplot(ref, aes(x = eco_index, y = sum_vst)) +
    geom_point(aes(fill = Site, shape = Treatment), size = 4, alpha = 0.7, colour = "black") +
    scale_fill_manual(name = "Site", values = site_colors[sites_present],
                      labels = site_labels[sites_present], drop = FALSE) +
    scale_shape_manual(name = "Treatment", values = c(21, 22, 23),
                       limits = c("NAT", "REST", "DAM"),
                       labels = c("NAT" = "Natural", "REST" = "Restored", "DAM" = "Degraded")) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 4, alpha = 1)),
      shape = guide_legend(order = 2, override.aes = list(fill = "grey40", size = 4, alpha = 1))
    ) +
    theme_light(base_size = 14) +
    theme(legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 13, colour = "black"),
          legend.key   = element_blank())

  leg <- ggpubr::as_ggplot(ggpubr::get_legend(p))

  # Height scales with number of legend rows so the title is never clipped
  leg_h <- 1.2 + 0.42 * (length(sites_present) + 3)
  fname <- paste0("Legend_", scheme$tag, ".png")
  for (root in c(new_root, nolabel_root, clr_root, clr_nolabel_root)) {
    ggsave(file.path(root, fname), leg, dpi = 1000, width = 2.4, height = leg_h, units = "in", bg = "white")
  }
  ggsave(fname, leg, dpi = 1000, width = 2.4, height = leg_h, units = "in", bg = "white")
  cat(paste0("  Legend written: ", fname, " (", length(sites_present), " sites)\n"))
}

cat("\nGenerating standalone legends...\n")
for (scheme in schemes) make_scheme_legend(scheme)

cat("\nAll figures and statistics written under: ", new_root, "\n")
