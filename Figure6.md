# 1. Read MAG metadata & HMM results

    mag_metadata <- read_csv("mag_metabolism.csv")

    ## Rows: 43181 Columns: 17
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (16): Host, Sample, site, treatment, Host Trend Group, Domain, Phylum, Class, Order, Family, Genus, Mod...
    ## dbl  (1): Host_Abundance
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    hmm_hits <- read_csv("mags_HMMhit.csv")

    ## Rows: 315 Columns: 2815
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1880): Category, Function, Gene.abbreviation, Gene.name, Hmm.file, Corresponding.KO, Reaction, Substra...
    ## dbl  (935): BAr1A1B1C__bin_10.Hit.numbers, BAr1A1B1C__bin_100.Hit.numbers, BAr1A1B1C__bin_103.Hit.numbers, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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

    hmm_unique_hosts

    ## # A tibble: 128 × 4
    ##    Function                `Host Trend Group` Category             unique_hosts
    ##    <chr>                   <fct>              <chr>                       <int>
    ##  1 Acetaldehyde => Ethanol NAT                Ethanol fermentation           38
    ##  2 Acetaldehyde => Ethanol REST               Ethanol fermentation           38
    ##  3 Acetaldehyde => Ethanol DAM                Ethanol fermentation           48
    ##  4 Acetate => Acetaldehyde NAT                Ethanol fermentation           47
    ##  5 Acetate => Acetaldehyde REST               Ethanol fermentation           47
    ##  6 Acetate => Acetaldehyde DAM                Ethanol fermentation           47
    ##  7 Acetate to acetyl-CoA   NAT                Fermentation                  173
    ##  8 Acetate to acetyl-CoA   REST               Fermentation                  174
    ##  9 Acetate to acetyl-CoA   DAM                Fermentation                  172
    ## 10 Acetogenesis            NAT                Fermentation                   34
    ## # ℹ 118 more rows

# 2. Fermentation

    # Prepare fermentation data & add alcohol utilization placeholder
    fermentation_data <- hmm_unique_hosts %>%
      filter(Category %in% c("Ethanol fermentation", "Fermentation")) %>%
      mutate(Category = "Fermentation")

    # Summarize unique_hosts per Host Trend Group
    fermentation_summary <- fermentation_data %>%
      group_by(`Host Trend Group`) %>%
      summarise(total_unique_hosts = sum(unique_hosts), .groups = "drop")

    # Plot fermentation data
    fermentation_data$Category <- "Fermentation\n"

### Fermentation bar plot

    fermentation_bar_plot <- ggplot(hmm_unique_hosts_by_category[hmm_unique_hosts_by_category$Category == "Fermentation",], 
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

    fermentation_bar_plot

![](Figure6_files/figure-markdown_strict/Fermentation%20bar%20plot-1.png)

### Save the plot

    ggsave("fermentation_bar_plot.png", plot = fermentation_bar_plot, dpi = 300, width = 3, height = 1.5)

### Fermentation scatter plot

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

    fermentation_plot

![](Figure6_files/figure-markdown_strict/Fermentation%20scatter%20plot-1.png)

### Save the plot

    ggsave("fermentation_scatter_plot.png", plot = fermentation_plot, dpi = 300, width = 10, height = 3)

# 3. KEGG

    # Load KEGG data
    mags_kegg <- read_csv("mags_keggmodulehit.csv")

    ## Rows: 470 Columns: 938
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (938): Module.ID, Module, Module.Category, BAr1A1B1C__bin_10.Module.presence, BAr1A1B1C__bin_100.Module...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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

    keggresult_by_category

    ## # A tibble: 88 × 3
    ##    `Host Trend Group` Module.Category                      unique_hosts
    ##    <fct>              <chr>                                       <int>
    ##  1 NAT                ATP synthesis                                 161
    ##  2 NAT                Arginine and proline metabolism               216
    ##  3 NAT                Aromatic amino acid metabolism                183
    ##  4 NAT                Aromatics degradation                          32
    ##  5 NAT                Biosynthesis of other antibiotics               2
    ##  6 NAT                Branched-chain amino acid metabolism          185
    ##  7 NAT                Carbon fixation                               177
    ##  8 NAT                Central carbohydrate metabolism               224
    ##  9 NAT                Cofactor and vitamin metabolism               224
    ## 10 NAT                Cysteine and methionine metabolism            200
    ## # ℹ 78 more rows

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

    keggresult

    ## # A tibble: 78 × 4
    ##    Function                                             `Host Trend Group` Category               unique_hosts
    ##    <chr>                                                <fct>              <fct>                         <dbl>
    ##  1 Acetyl-CoA pathway, CO2 => acetyl-CoA                NAT                "Methane\nmetabolism"             7
    ##  2 Acetyl-CoA pathway, CO2 => acetyl-CoA                REST               "Methane\nmetabolism"             8
    ##  3 Acetyl-CoA pathway, CO2 => acetyl-CoA                DAM                "Methane\nmetabolism"             6
    ##  4 Assimilatory nitrate reduction, nitrate => ammonia   NAT                "Nitrogen\nmetabolism"           20
    ##  5 Assimilatory nitrate reduction, nitrate => ammonia   REST               "Nitrogen\nmetabolism"           17
    ##  6 Assimilatory nitrate reduction, nitrate => ammonia   DAM                "Nitrogen\nmetabolism"            7
    ##  7 Assimilatory sulfate reduction, sulfate => H2S       NAT                "Sulfur\nmetabolism"              6
    ##  8 Assimilatory sulfate reduction, sulfate => H2S       REST               "Sulfur\nmetabolism"              3
    ##  9 Assimilatory sulfate reduction, sulfate => H2S       DAM                "Sulfur\nmetabolism"             10
    ## 10 C4-dicarboxylic acid cycle, NADP - malic enzyme type NAT                "Carbon\nfixation"               53
    ## # ℹ 68 more rows

    # global presence data
    genecount <- read.csv("genecountinbac&arc_kegg.csv")
    genecount_long <- melt(genecount)

    ## Using Category, Function as id variables

    genecount_long

    ##               Category
    ## 1      Carbon fixation
    ## 2      Carbon fixation
    ## 3      Carbon fixation
    ## 4      Carbon fixation
    ## 5      Carbon fixation
    ## 6      Carbon fixation
    ## 7      Carbon fixation
    ## 8      Carbon fixation
    ## 9      Carbon fixation
    ## 10     Carbon fixation
    ## 11     Carbon fixation
    ## 12     Carbon fixation
    ## 13     Carbon fixation
    ## 14  Methane metabolism
    ## 15  Methane metabolism
    ## 16  Methane metabolism
    ## 17  Methane metabolism
    ## 18  Methane metabolism
    ## 19  Methane metabolism
    ## 20  Methane metabolism
    ## 21  Methane metabolism
    ## 22  Methane metabolism
    ## 23  Methane metabolism
    ## 24  Methane metabolism
    ## 25  Methane metabolism
    ## 26  Methane metabolism
    ## 27 Nitrogen metabolism
    ## 28 Nitrogen metabolism
    ## 29 Nitrogen metabolism
    ## 30 Nitrogen metabolism
    ## 31 Nitrogen metabolism
    ## 32 Nitrogen metabolism
    ## 33 Nitrogen metabolism
    ## 34   Sulfur metabolism
    ## 35   Sulfur metabolism
    ## 36   Sulfur metabolism
    ## 37   Sulfur metabolism
    ## 38   Sulfur metabolism
    ## 39   Sulfur metabolism
    ## 40   Sulfur metabolism
    ## 41   Sulfur metabolism
    ## 42   Sulfur metabolism
    ## 43   Sulfur metabolism
    ## 44   Sulfur metabolism
    ## 45     Carbon fixation
    ## 46     Carbon fixation
    ## 47     Carbon fixation
    ## 48     Carbon fixation
    ## 49     Carbon fixation
    ## 50     Carbon fixation
    ## 51     Carbon fixation
    ## 52     Carbon fixation
    ## 53     Carbon fixation
    ## 54     Carbon fixation
    ## 55     Carbon fixation
    ## 56     Carbon fixation
    ## 57     Carbon fixation
    ## 58  Methane metabolism
    ## 59  Methane metabolism
    ## 60  Methane metabolism
    ## 61  Methane metabolism
    ## 62  Methane metabolism
    ## 63  Methane metabolism
    ## 64  Methane metabolism
    ## 65  Methane metabolism
    ## 66  Methane metabolism
    ## 67  Methane metabolism
    ## 68  Methane metabolism
    ## 69  Methane metabolism
    ## 70  Methane metabolism
    ## 71 Nitrogen metabolism
    ## 72 Nitrogen metabolism
    ## 73 Nitrogen metabolism
    ## 74 Nitrogen metabolism
    ## 75 Nitrogen metabolism
    ## 76 Nitrogen metabolism
    ## 77 Nitrogen metabolism
    ## 78   Sulfur metabolism
    ## 79   Sulfur metabolism
    ## 80   Sulfur metabolism
    ## 81   Sulfur metabolism
    ## 82   Sulfur metabolism
    ## 83   Sulfur metabolism
    ## 84   Sulfur metabolism
    ## 85   Sulfur metabolism
    ## 86   Sulfur metabolism
    ## 87   Sulfur metabolism
    ## 88   Sulfur metabolism
    ##                                                                                                     Function
    ## 1                                                           Reductive pentose phosphate cycle (Calvin cycle)
    ## 2                                                                   CAM (Crassulacean acid metabolism), dark
    ## 3                                                                  CAM (Crassulacean acid metabolism), light
    ## 4                                                       C4-dicarboxylic acid cycle, NADP - malic enzyme type
    ## 5                                                        C4-dicarboxylic acid cycle, NAD - malic enzyme type
    ## 6                                         C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type
    ## 7                                                             Reductive citrate cycle (Arnon-Buchanan cycle)
    ## 8                                                                               3-Hydroxypropionate bi-cycle
    ## 9                                                                    Hydroxypropionate-hydroxybutylate cycle
    ## 10                                                                       Dicarboxylate-hydroxybutyrate cycle
    ## 11                                                     Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)
    ## 12                                 Phosphate acetyltransferase-acetate kinase pathway, acetyl-CoA => acetate
    ## 13                                            Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate
    ## 14                                                                            Methanogenesis, CO2 => methane
    ## 15                                                                        Methanogenesis, acetate => methane
    ## 16                                                                       Methanogenesis, methanol => methane
    ## 17                                       Methanogenesis, methylamine0dimethylamine0trimethylamine => methane
    ## 18                                                                                   Coenzyme M biosynthesis
    ## 19    2-Oxocarboxylic acid chain extension, 2-oxoglutarate => 2-oxoadipate => 2-oxopimelate => 2-oxosuberate
    ## 20                                                  Methane oxidation, methanotroph, methane => formaldehyde
    ## 21                                                                 Formaldehyde assimilation, serine pathway
    ## 22                                                 Formaldehyde assimilation, ribulose monophosphate pathway
    ## 23                                                 Formaldehyde assimilation, xylulose monophosphate pathway
    ## 24                                                                                F420 biosynthesis, archaea
    ## 25                                                                                 Methanofuran biosynthesis
    ## 26                                                                     Acetyl-CoA pathway, CO2 => acetyl-CoA
    ## 27                                                                    Nitrogen fixation, nitrogen => ammonia
    ## 28                                                        Assimilatory nitrate reduction, nitrate => ammonia
    ## 29                                                       Dissimilatory nitrate reduction, nitrate => ammonia
    ## 30                                                                      Denitrification, nitrate => nitrogen
    ## 31                                                                         Nitrification, ammonia => nitrite
    ## 32                                           Complete nitrification, comammox, ammonia => nitrite => nitrate
    ## 33                                                                    Anammox, nitrite + ammonia => nitrogen
    ## 34                                                    Assimilatory sulfate reduction, plants, sulfate => H2S
    ## 35                                                            Assimilatory sulfate reduction, sulfate => H2S
    ## 36                                                           Dissimilatory sulfate reduction, sulfate => H2S
    ## 37                                                      Sulfur oxidation, SOX system, thiosulfate => sulfate
    ## 38 Sulfur oxidation, tetrathionate intermediate (S4I) pathway, thiosulfate => sulfur + sulfate + thiosulfate
    ## 39                                                                       Sulfur reduction, sulfur => sulfide
    ## 40                                                                      Sulfide oxidation, sulfide => sulfur
    ## 41                Dimethylsulfoniopropionate (DMSP) degradation, demethylation pathway, DMSP => methanethiol
    ## 42                       Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => acetyl-CoA
    ## 43           Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => acrylate => acetyl-CoA
    ## 44                    Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => propionyl-CoA
    ## 45                                                          Reductive pentose phosphate cycle (Calvin cycle)
    ## 46                                                                  CAM (Crassulacean acid metabolism), dark
    ## 47                                                                 CAM (Crassulacean acid metabolism), light
    ## 48                                                      C4-dicarboxylic acid cycle, NADP - malic enzyme type
    ## 49                                                       C4-dicarboxylic acid cycle, NAD - malic enzyme type
    ## 50                                        C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type
    ## 51                                                            Reductive citrate cycle (Arnon-Buchanan cycle)
    ## 52                                                                              3-Hydroxypropionate bi-cycle
    ## 53                                                                   Hydroxypropionate-hydroxybutylate cycle
    ## 54                                                                       Dicarboxylate-hydroxybutyrate cycle
    ## 55                                                     Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)
    ## 56                                 Phosphate acetyltransferase-acetate kinase pathway, acetyl-CoA => acetate
    ## 57                                            Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate
    ## 58                                                                            Methanogenesis, CO2 => methane
    ## 59                                                                        Methanogenesis, acetate => methane
    ## 60                                                                       Methanogenesis, methanol => methane
    ## 61                                       Methanogenesis, methylamine0dimethylamine0trimethylamine => methane
    ## 62                                                                                   Coenzyme M biosynthesis
    ## 63    2-Oxocarboxylic acid chain extension, 2-oxoglutarate => 2-oxoadipate => 2-oxopimelate => 2-oxosuberate
    ## 64                                                  Methane oxidation, methanotroph, methane => formaldehyde
    ## 65                                                                 Formaldehyde assimilation, serine pathway
    ## 66                                                 Formaldehyde assimilation, ribulose monophosphate pathway
    ## 67                                                 Formaldehyde assimilation, xylulose monophosphate pathway
    ## 68                                                                                F420 biosynthesis, archaea
    ## 69                                                                                 Methanofuran biosynthesis
    ## 70                                                                     Acetyl-CoA pathway, CO2 => acetyl-CoA
    ## 71                                                                    Nitrogen fixation, nitrogen => ammonia
    ## 72                                                        Assimilatory nitrate reduction, nitrate => ammonia
    ## 73                                                       Dissimilatory nitrate reduction, nitrate => ammonia
    ## 74                                                                      Denitrification, nitrate => nitrogen
    ## 75                                                                         Nitrification, ammonia => nitrite
    ## 76                                           Complete nitrification, comammox, ammonia => nitrite => nitrate
    ## 77                                                                    Anammox, nitrite + ammonia => nitrogen
    ## 78                                                    Assimilatory sulfate reduction, plants, sulfate => H2S
    ## 79                                                            Assimilatory sulfate reduction, sulfate => H2S
    ## 80                                                           Dissimilatory sulfate reduction, sulfate => H2S
    ## 81                                                      Sulfur oxidation, SOX system, thiosulfate => sulfate
    ## 82 Sulfur oxidation, tetrathionate intermediate (S4I) pathway, thiosulfate => sulfur + sulfate + thiosulfate
    ## 83                                                                       Sulfur reduction, sulfur => sulfide
    ## 84                                                                      Sulfide oxidation, sulfide => sulfur
    ## 85                Dimethylsulfoniopropionate (DMSP) degradation, demethylation pathway, DMSP => methanethiol
    ## 86                       Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => acetyl-CoA
    ## 87           Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => acrylate => acetyl-CoA
    ## 88                    Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => propionyl-CoA
    ##    variable value
    ## 1       Bac   458
    ## 2       Bac  4061
    ## 3       Bac  1409
    ## 4       Bac     0
    ## 5       Bac     0
    ## 6       Bac     0
    ## 7       Bac    23
    ## 8       Bac     5
    ## 9       Bac     0
    ## 10      Bac     0
    ## 11      Bac    72
    ## 12      Bac  6153
    ## 13      Bac     0
    ## 14      Bac     0
    ## 15      Bac     0
    ## 16      Bac     0
    ## 17      Bac     0
    ## 18      Bac     1
    ## 19      Bac     0
    ## 20      Bac    53
    ## 21      Bac    63
    ## 22      Bac   664
    ## 23      Bac     0
    ## 24      Bac     1
    ## 25      Bac     0
    ## 26      Bac     1
    ## 27      Bac   814
    ## 28      Bac  1708
    ## 29      Bac  1547
    ## 30      Bac   333
    ## 31      Bac    45
    ## 32      Bac     6
    ## 33      Bac     2
    ## 34      Bac     0
    ## 35      Bac  3093
    ## 36      Bac   143
    ## 37      Bac   557
    ## 38      Bac    17
    ## 39      Bac   177
    ## 40      Bac   386
    ## 41      Bac    23
    ## 42      Bac    32
    ## 43      Bac     3
    ## 44      Bac     3
    ## 45      Arc     0
    ## 46      Arc   226
    ## 47      Arc     0
    ## 48      Arc     0
    ## 49      Arc     0
    ## 50      Arc     0
    ## 51      Arc     0
    ## 52      Arc     0
    ## 53      Arc    44
    ## 54      Arc    12
    ## 55      Arc     0
    ## 56      Arc    27
    ## 57      Arc    68
    ## 58      Arc   117
    ## 59      Arc    84
    ## 60      Arc    56
    ## 61      Arc    41
    ## 62      Arc    56
    ## 63      Arc   120
    ## 64      Arc     0
    ## 65      Arc     0
    ## 66      Arc     4
    ## 67      Arc     0
    ## 68      Arc   272
    ## 69      Arc    45
    ## 70      Arc   103
    ## 71      Arc    73
    ## 72      Arc    60
    ## 73      Arc     0
    ## 74      Arc     0
    ## 75      Arc     0
    ## 76      Arc     0
    ## 77      Arc     0
    ## 78      Arc     0
    ## 79      Arc    21
    ## 80      Arc     5
    ## 81      Arc     0
    ## 82      Arc    22
    ## 83      Arc   167
    ## 84      Arc     0
    ## 85      Arc     0
    ## 86      Arc     0
    ## 87      Arc     0
    ## 88      Arc     0

    colnames(genecount_long) <- c("Category","Function","variable","genecount")

## 1) Methane

    # Subset data where Category is "Methane metabolism"
    methane_result <- keggresult[keggresult$Category == "Methane\nmetabolism",]

    # Get the order of Function in methane_result
    methane_function_order <- unique(methane_result$Function)

    # Reorder Function in genecount_long to match methane_result
    methane_genecount_long_sorted <- genecount_long[genecount_long$Category == "Methane metabolism",]
    methane_genecount_long_sorted$Function <- factor(methane_genecount_long_sorted$Function,levels = methane_function_order)

### Methane scatter plot

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

    methane_plot1

![](Figure6_files/figure-markdown_strict/Methane%20scatter%20plot-1.png)

### Global methane presence

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

    methane_plot2

![](Figure6_files/figure-markdown_strict/Global%20methane%20presence-1.png)

### Methane bar plot

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

    methane_plot3

![](Figure6_files/figure-markdown_strict/Methane%20bar%20plot-1.png)
\### Save the plot

    ggsave("methane_scatter_plot.png", plot = methane_plot1, dpi = 300, width = 10, height = 3)
    ggsave("methane_global_presence.png", plot = methane_plot2, dpi = 300, width = 4, height = 1)
    ggsave("methane_bar_plot.png", plot = methane_plot3, dpi = 300, width = 3, height = 1.5)

## 2) Nitrogen

    # Subset data where Category is "Nitrogen metabolism"
    nitrogen_result <- keggresult[keggresult$Category == "Nitrogen\nmetabolism",]

    # Get the order of Function in nitrogen_result
    nitrogen_function_order <- unique(nitrogen_result$Function)

    # Reorder Function in genecount_long to match nitrogen_result
    nitrogen_genecount_long_sorted <- genecount_long[genecount_long$Category == "Nitrogen metabolism",]
    nitrogen_genecount_long_sorted$Function <- factor(nitrogen_genecount_long_sorted$Function,
                                                     levels = nitrogen_function_order)

### Nitrogen scatter plot

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

    nitrogen_plot1

![](Figure6_files/figure-markdown_strict/Nitrogen%20scatter%20plot-1.png)

### Global nitrogen presence

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

    nitrogen_plot2

![](Figure6_files/figure-markdown_strict/Global%20nitrogen%20presence-1.png)

### Nitrogen bar plot

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

    nitrogen_plot3

![](Figure6_files/figure-markdown_strict/Nitrogen%20bar%20plot-1.png)
\### Save the plot

    ggsave("nitrogen_scatter_plot.png", plot = nitrogen_plot1, dpi = 300, width = 10, height = 3)
    ggsave("nitrogen_global_presence.png", plot = nitrogen_plot2, dpi = 300, width = 4, height = 1)
    ggsave("nitrogen_bar_plot.png", plot = nitrogen_plot3, dpi = 300, width = 3, height = 1.5)

## 3) Sulfur

    # Subset data where Category is "sulfur metabolism"
    sulfur_result <- keggresult[keggresult$Category == "Sulfur\nmetabolism",]

    # Get the order of Function in sulfur_result
    sulfur_function_order <- unique(sulfur_result$Function)

    # Reorder Function in genecount_long to match sulfur_result
    sulfur_genecount_long_sorted <- genecount_long[genecount_long$Category == "Sulfur metabolism",]
    sulfur_genecount_long_sorted$Function <- factor(sulfur_genecount_long_sorted$Function,
                                                      levels = sulfur_function_order)

### Sulfur scatter plot

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

    sulfur_plot1

![](Figure6_files/figure-markdown_strict/Sulfur%20scatter%20plot-1.png)

### Global sulfur presence

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

    sulfur_plot2

![](Figure6_files/figure-markdown_strict/Global%20sulfur%20presence-1.png)

### Sulfur bar plot

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

### Save the plot

    ggsave("sulfur_scatter_plot.png", plot = sulfur_plot1, dpi = 300, width = 10, height = 3)
    ggsave("sulfur_global_presence.png", plot = sulfur_plot2, dpi = 300, width = 4, height = 1)
    ggsave("sulfur_bar_plot.png", plot = sulfur_plot3, dpi = 300, width = 3, height = 1.5)

## 4) Carbon fixation

    # Subset data where Category is "carbon metabolism"
    carbon_result <- keggresult[keggresult$Category == "Carbon\nfixation",]

    # Get the order of Function in carbon_result
    carbon_function_order <- unique(carbon_result$Function)

    # Reorder Function in genecount_long to match carbon_result
    carbon_genecount_long_sorted <- genecount_long[genecount_long$Category == "Carbon fixation",]
    carbon_genecount_long_sorted$Function <- factor(carbon_genecount_long_sorted$Function,
                                                    levels = carbon_function_order)

### Carbon fixation scatter plot

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

### Global carbon fixation presence

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

### Carbon fixation bar plot

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

### Save the plot

    ggsave("carbonfixation_scatter_plot.png", plot = methane_plot1, dpi = 300, width = 10, height = 3)
    ggsave("carbonfixation_global_presence.png", plot = methane_plot2, dpi = 300, width = 4, height = 1)
    ggsave("carbonfixation_bar_plot.png", plot = methane_plot3, dpi = 300, width = 3, height = 1.5)

    # 4.Combine plots

    # Convert global presence plots to grob
    methane_grob_plot2 <- ggplotGrob(methane_plot2)

    nitrogen_grob_plot2 <- ggplotGrob(nitrogen_plot2)

    sulfur_grob_plot2 <- ggplotGrob(sulfur_plot2)

    carbon_grob_plot2 <- ggplotGrob(carbon_plot2)

    # Combine plots
    methane_plot <- methane_plot1 + 
      annotation_custom(methane_grob_plot2, xmin = -.2, xmax = 10, ymin = 56, ymax = 63.8)

    nitrogen_plot <- nitrogen_plot1 + 
      annotation_custom(nitrogen_grob_plot2, xmin = -.2, xmax = 7, ymin = 56, ymax = 63.8) 

    sulfur_plot <- sulfur_plot1 + 
      annotation_custom(sulfur_grob_plot2, xmin = -.2, xmax = 3, ymin = 56, ymax = 63.8) 

    carbon_plot <- carbon_plot1 + 
      annotation_custom(carbon_grob_plot2, xmin = -.2, xmax = 7, ymin = 177, ymax = 203) 

    fermentation_plot+theme(plot.margin = margin(-146, 5.5, 194, 320))

![](Figure6_files/figure-markdown_strict/Combine%20plots-1.png)

    # Adjusting position
    nitrogen_plot <- nitrogen_plot + 
      theme(plot.margin = margin(-2, 7, 70,-4)) 

    sulfur_plot <- sulfur_plot + 
      theme(plot.margin = margin(-63, 6.8, 165, 18)) 

    carbon_plot <- carbon_plot + 
      theme(plot.margin = margin(-163, 6.5, 190,-76)) 

    fermentation_plot <- fermentation_plot  + 
      theme(plot.margin = margin(-183, 6.5, 200, 162)) 

### Combine all pathway plots into one figure

    combined_plot <- grid.arrange(
      methane_plot, nitrogen_plot, sulfur_plot, carbon_plot, fermentation_plot,
      ncol = 1
    )

![](Figure6_files/figure-markdown_strict/Combined%20plot-1.png)

    combined_plot

    ## TableGrob (5 x 1) "arrange": 5 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (2-2,1-1) arrange gtable[layout]
    ## 3 3 (3-3,1-1) arrange gtable[layout]
    ## 4 4 (4-4,1-1) arrange gtable[layout]
    ## 5 5 (5-5,1-1) arrange gtable[layout]

### Save the plot

    ggsave("combined_plot.pdf", plot = combined_plot, dpi = 300, width = 15, height = 12)

### Plot legend only

    legend_plot <- ggplot(carbon_result, aes(x = Function, y = unique_hosts, fill = `Host Trend Group`)) +
      geom_point(size = 0, alpha = 0, shape = 21) +
      scale_fill_manual(values = c('NAT' = '#4daf4a', 'REST' = '#377eb8', 'DAM' = '#e41a1c')) +
      labs(fill = 'Treatment') +
      theme_void() +
      theme(legend.position = "right") +
      guides(fill = guide_legend(override.aes = list(alpha = 1, size = 5)))

    legend_plot

![](Figure6_files/figure-markdown_strict/Legend%20plot-1.png) \### Save
the plot

    ggsave("legend_only.png", plot = legend_plot, width = 4, height = 4)
