# 1. Figure 2a

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

    dna_prod

    ##     Sample        Site Treatment      Growth  MAP MAT elevation Moisture   pH Moss_cover      O2
    ## 1    MGr1A    Migneint       NAT  0.38869360 2181 8.0       453 93.60220 3.79 66.1876879  5.3940
    ## 2    MGr1B    Migneint       NAT  0.34086121 2181 8.0       453 91.79173 3.94 54.3505188  2.6115
    ## 3    MGr1C    Migneint       NAT  0.57375152 2181 8.0       453 88.85526 3.92 32.2594406  3.8380
    ## 4    MGr2D    Migneint      REST  0.35266063 2181 8.0       453 90.50027 3.83 37.1269318  7.6860
    ## 5    MGr2E    Migneint      REST  0.51621336 2181 8.0       453 89.01993 3.91 30.4905093  3.2345
    ## 6    MGr2F    Migneint      REST  0.73237813 2181 8.0       453 88.06708 3.85 11.7036794  6.0275
    ## 7    MGr3G    Migneint       DAM  1.44582310 2181 8.0       453 88.79431 3.88  3.3525996  5.8700
    ## 8    MGr3H    Migneint       DAM  0.07219564 2181 8.0       453 89.94187 3.86 13.9967492  9.2865
    ## 9    MGr3I    Migneint       DAM  0.26602057 2181 8.0       453 88.66024 3.75  5.8353934  6.8585
    ## 10   CRr1A     Crocach       NAT  1.68775735 1258 7.1       189 87.17251 3.99 11.7858922  6.2760
    ## 11   CRr1B     Crocach       NAT  0.63721974 1258 7.1       189 89.04116 4.04 18.5248394  5.5245
    ## 12   CRr1C     Crocach       NAT  0.46046682 1258 7.1       189 90.24133 4.13 31.7718780  5.7825
    ## 13   CRr2D     Crocach      REST  0.39735323 1258 7.1       189 84.08085 4.03 20.7780596  4.9900
    ## 14   CRr2E     Crocach      REST  1.14204704 1258 7.1       189 85.72452 4.02 12.9237580  7.1840
    ## 15   CRr2F     Crocach      REST  0.69912406 1258 7.1       189 85.03915 4.10 22.3215287  7.5285
    ## 16   CRr3G     Crocach       DAM  0.83886133 1258 7.1       189 89.30760 3.83 11.0334833  5.2010
    ## 17   CRr3H     Crocach       DAM  0.44089702 1258 7.1       189 89.93341 3.89 17.8346918  8.0785
    ## 18   CRr3I     Crocach       DAM  1.50661803 1258 7.1       189 89.83699 4.16  4.2322748  7.4980
    ## 19   SEr2D       Stean      REST  0.17650038 1229 8.0       530 92.52280 3.86         NA  8.2180
    ## 20   SEr2E       Stean      REST  1.01480188 1229 8.0       530 91.69001 3.96         NA 11.8200
    ## 21   SEr2F       Stean      REST  0.42573170 1229 8.0       530 89.54129 4.02         NA  8.9430
    ## 22   SEr3G       Stean       DAM  2.14676311 1229 8.0       530 87.32310 3.82         NA 10.7600
    ## 23   SEr3H       Stean       DAM  0.61830575 1229 8.0       530 90.13841 3.75         NA  9.7260
    ## 24   SEr3I       Stean       DAM  3.84213321 1229 8.0       530 87.11373 3.74         NA  9.9270
    ## 25   MHr1A Moors_House       NAT  2.05837200 1699 8.0       517 89.77337 3.75 20.8881867 10.5000
    ## 26   MHr1B Moors_House       NAT  0.30296316 1699 8.0       517 92.64590 3.98 19.3135510  7.5605
    ## 27   MHr1C Moors_House       NAT  0.15537793 1699 8.0       517 92.69704 3.95 31.7324024  8.2995
    ## 28   MHr2D Moors_House      REST  1.20253865 1699 8.0       517 87.99768 3.89  6.1741580  9.1380
    ## 29   MHr2E Moors_House      REST  0.76496729 1699 8.0       517 87.48715 3.91  1.6966111  8.6860
    ## 30   MHr2F Moors_House      REST  0.14974297 1699 8.0       517 90.40447 3.74  1.7706887 10.3225
    ## 31   MHr3G Moors_House       DAM  1.68531614 1699 8.0       517 90.54129 3.70 15.7959643  8.2860
    ## 32   MHr3H Moors_House       DAM  0.38781041 1699 8.0       517 90.30659 3.64  1.3777282  9.7405
    ## 33   MHr3I Moors_House       DAM  0.92478687 1699 8.0       517 88.86773 3.70  2.3262473  7.1335
    ## 34 LABRr1A    Langwell       NAT  0.43283820 1223 7.0       255 89.87880 4.15 15.6898898  6.2085
    ## 35 LABRr1B    Langwell       NAT  0.09182724 1223 7.0       255 92.43188 4.02 21.8623275 10.0500
    ## 36 LABRr1C    Langwell       NAT  0.79484349 1223 7.0       255 89.73957 3.99 27.6016439  6.8445
    ## 37 LASAr3G    Langwell       DAM  0.32197303 1223 7.0       255 88.97433 3.78  3.1468596  8.5450
    ## 38 LASAr3H    Langwell       DAM  1.15336929 1223 7.0       255 88.81946 3.84 17.9265668  9.1345
    ## 39 LASAr3I    Langwell       DAM  0.66822415 1223 7.0       255 89.12712 3.82  3.5386288  9.7240
    ## 40   BOr1A     Bowness       NAT  0.76778036  953 9.6        66 89.44489 3.62 18.5421161  6.9800
    ## 41   BOr1B     Bowness       NAT  0.55018238  953 9.6        66 87.74365 3.70  6.9777328  6.6940
    ## 42   BOr1C     Bowness       NAT  2.92101820  953 9.6        66 85.52546 3.69  5.0326398  6.0255
    ## 43   BOr2D     Bowness      REST  1.23910383  953 9.6        66 88.39204 4.15  4.9809653  8.0655
    ## 44   BOr2E     Bowness      REST  0.25789753  953 9.6        66 86.48992 4.04  0.7924928  7.7955
    ## 45   BOr2F     Bowness      REST  2.17106877  953 9.6        66 87.46045 4.06 13.0070921  8.8140
    ## 46   BOr3G     Bowness       DAM  5.81848179  953 9.6        66 83.35184 3.60  2.0521152  9.8890
    ## 47   BOr3H     Bowness       DAM  4.34250819  953 9.6        66 81.90553 3.59  3.6322984  9.8460
    ## 48   BOr3I     Bowness       DAM  2.87260494  953 9.6        66 85.83180 3.62  5.2124817  8.8085
    ## 49   BAr1A    Balmoral       NAT  0.79109328 1412 5.5       695 86.23216 3.73  2.5936194  9.8740
    ## 50   BAr1B    Balmoral       NAT  5.14052550 1412 5.5       695 87.57867 3.80 14.6546829  8.8665
    ## 51   BAr1C    Balmoral       NAT  4.90633352 1412 5.5       695 90.59173 3.91  5.0644358  7.0405
    ## 52   BAr2D    Balmoral      REST 10.94489133 1412 5.5       695 87.49755 3.71  8.8829984  9.2875
    ## 53   BAr2E    Balmoral      REST  7.93344148 1412 5.5       695 87.44708 3.73 14.0147498  8.2275
    ## 54   BAr2F    Balmoral      REST  6.60151559 1412 5.5       695 85.75619 3.51 20.5633611  8.3215
    ## 55   BAr3G    Balmoral       DAM 10.47585538 1412 5.5       695 79.57904 3.60 14.7352104  9.4475
    ## 56   BAr3H    Balmoral       DAM  8.35337506 1412 5.5       695 79.88972 3.59  1.5695768  9.3530
    ## 57   BAr3I    Balmoral       DAM  7.49684880 1412 5.5       695 75.38969 3.62  2.2326978 10.2955
    ## 58 LASCr2D    Langwell      REST  2.85406531 1223 7.0       255 90.79751 4.08 14.5214756  9.9630
    ## 59 LASCr2E    Langwell      REST  0.39066372 1223 7.0       255 90.99193 4.04 18.0108674  8.5460
    ## 60 LASCr2F    Langwell      REST  0.52917065 1223 7.0       255 89.54370 4.03 11.7345203 10.1550
    ##          TC        TN      C.N FTIR_axis_1 degradation_index
    ## 1  49.51977 1.8755337 26.40303   4.7201213       1.039833529
    ## 2  50.84798 2.4383523 20.85342   3.5834275       1.396351814
    ## 3  43.88880 2.1822319 20.11189   8.4835592       0.776499950
    ## 4  45.53867 1.9665606 23.15650  -4.1697732       0.688060648
    ## 5  48.84855 1.9837331 24.62456  -0.1013926       0.919580766
    ## 6  51.89397 2.3162584 22.40422   1.6136481       0.364005021
    ## 7  50.16782 1.7664410 28.40051  10.5528150      -0.026448296
    ## 8  47.01845 1.9150014 24.55270  -2.4318758       0.202461524
    ## 9  49.97535 2.0100722 24.86246   7.6933825      -0.041149766
    ## 10 47.89912 2.0214064 23.69594  10.5843964       0.137120448
    ## 11 49.00604 1.8297397 26.78307   2.8671865       0.529210487
    ## 12 47.24883 1.7600350 26.84539  -5.5195485       0.974647050
    ## 13 53.16726 1.3586918 39.13121   5.0262560      -0.033306484
    ## 14 49.32398 1.9584445 25.18528   2.2507617       0.173641711
    ## 15 49.88054 2.1696458 22.99018  -1.2004110       0.401446952
    ## 16 51.33075 1.5050452 34.10579  -2.1481348       0.168110443
    ## 17 46.77607 1.6273959 28.74290   6.9155312       0.042392608
    ## 18 43.69639 1.7116307 25.52910  -1.9512279       0.465350015
    ## 19 46.44328 1.1211059 41.42631  -4.0531888       0.055305662
    ## 20 45.62141 1.0063529 45.33341 -35.5976450       0.462907656
    ## 21 45.01654 1.2458235 36.13396  -1.5807610       0.016103806
    ## 22 52.86704 1.1503860 45.95592   2.6460162      -0.922147193
    ## 23 49.64233 0.9110925 54.48660  -3.7233433      -0.763943486
    ## 24 51.86971 1.3637865 38.03360   1.3774953      -0.627714840
    ## 25 47.97347 1.3372927 35.87357  -1.8734713      -0.252022444
    ## 26 49.47584 1.5526851 31.86469  -5.6546686       0.551376705
    ## 27 49.35648 1.3935555 35.41766  -8.3060380       0.555763355
    ## 28 50.64264 1.7610246 28.75749   3.1013350      -0.188796232
    ## 29 50.98798 1.9943821 25.56580  -2.4901959       0.024697972
    ## 30 48.64301 1.4993716 32.44227  -4.7352177      -0.240792645
    ## 31 49.85317 1.5411164 32.34874  -1.4938349       0.008831597
    ## 32 50.68416 1.5832690 32.01235   8.3938829      -0.582883862
    ## 33 47.86358 1.7537367 27.29234  -2.7526421       0.068085083
    ## 34 49.90676 1.7855492 27.95037  -4.4841524       0.709029404
    ## 35 44.34337 1.2840947 34.53279  -7.5750774       0.318794546
    ## 36 47.69652 1.7810409 26.78014  -5.8549973       0.704354976
    ## 37 46.83960 1.3907166 33.68019   5.3588569      -0.374910130
    ## 38 50.28053 1.9692106 25.53334   6.4976671      -0.044798328
    ## 39 45.57399 1.7046403 26.73526   9.7899008      -0.384444538
    ## 40 47.23617 1.4358318 32.89812  -2.5185531       0.066541715
    ## 41 48.25644 1.5578886 30.97554  -2.7836535       0.006889205
    ## 42 46.29145 1.3583310 34.07965  -4.1653549      -0.112365103
    ## 43 49.29037 1.6016167 30.77538  -6.0779021       0.299486592
    ## 44 47.98551 1.7381197 27.60771  -3.0249243       0.111499409
    ## 45 49.52766 1.6779313 29.51710  -5.0174196       0.195822181
    ## 46 46.64932 1.7134116 27.22599  -0.6555999      -0.626651551
    ## 47 45.12646 1.7978941 25.09962   9.9198881      -0.935412031
    ## 48 49.97747 2.2437608 22.27397   6.1138685      -0.381937074
    ## 49 46.27535 1.1701019 39.54814   6.3061509      -0.872388753
    ## 50 46.42556 1.1706003 39.65961  -1.9565497      -0.320731707
    ## 51 49.78674 1.2474421 39.91106  -5.5507649       0.089770005
    ## 52 51.12538 1.1211820 45.59954   6.5308311      -0.861788205
    ## 53 52.07274 1.5402536 33.80790   8.8456267      -0.438826990
    ## 54 51.15368 1.5189359 33.67731  -3.2045526      -0.325096378
    ## 55 42.35862 1.4844972 28.53398   2.2691565      -0.755040527
    ## 56 50.73779 1.8095157 28.03943   1.2498036      -0.831613140
    ## 57 47.02564 1.4212520 33.08748  10.5733020      -1.497950267
    ## 58 41.44368 0.9723937 42.62027  -7.9420040       0.011132604
    ## 59 46.26957 1.0615925 43.58506  -6.2830715       0.105479608
    ## 60 51.38515 1.2663594 40.57707   6.0899599      -0.437896483

### Microbial Growth Rate Across Treatments plot

    growth_rate_treatment_plot <- ggplot(dna_prod,aes(x=Treatment, y=Growth))+
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

    growth_rate_treatment_plot

![](Figure2_files/figure-markdown_strict/Microbial%20Growth%20Rate%20Across%20Treatments%20plot-1.png)
\### Save the plot

    ggsave("Figure_2a_MicrobialGrowthRate_Treatment.png",  dpi=1000, width = 4, height = 4, units = "in") 

# 2. Figure 2b

    # Compute mean Growth for each Site and Treatment
    mean_growth <- dna_prod %>%
      group_by(Site, Treatment) %>%
      summarise(mean_growth = mean(Growth, na.rm = TRUE), .groups = "drop")

    site_order <- c("Balmoral","Bowness","Stean","Moors_House","Langwell","Crocach","Migneint")

    dna_prod$Site <- factor(dna_prod$Site, levels = site_order)
    mean_growth$Site <- factor(mean_growth$Site, levels = site_order)

### Microbial Growth Rate Across Sites plot

    growth_rate_sites_plot <- ggplot(dna_prod, aes(y = Growth, x = Site, colour = Treatment)) +
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

    growth_rate_sites_plot

![](Figure2_files/figure-markdown_strict/Microbial%20Growth%20Rate%20Across%20Sites%20plot-1.png)
\### Save the plot

    ggsave("Figure_2b_MicrobialGrowthRate_Site.png",  dpi=1000, width = 6, height =5, units = "in") 

# 3. Figure 2c

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

    ## Linear mixed-effects model fit by REML
    ##   Data: dna_prod_test 
    ##        AIC     BIC    logLik
    ##   177.2717 187.487 -83.63587
    ## 
    ## Random effects:
    ##  Formula: ~1 | Site
    ##         (Intercept)  Residual
    ## StdDev:   0.8619272 0.8620367
    ## 
    ## Fixed effects:  log(Growth) ~ Treatment 
    ##                    Value Std.Error DF    t-value p-value
    ## (Intercept)   -0.3621749 0.3872240 51 -0.9353110  0.3540
    ## TreatmentDAM   0.6638070 0.2814201 51  2.3587761  0.0222
    ## TreatmentREST  0.2545121 0.2814201 51  0.9043849  0.3700
    ##  Correlation: 
    ##               (Intr) TrtDAM
    ## TreatmentDAM  -0.402       
    ## TreatmentREST -0.402  0.553
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.5090888 -0.5015986  0.1675783  0.6377418  1.9065916 
    ## 
    ## Number of Observations: 60
    ## Number of Groups: 7

    plot(growth_rate_lme_homo)

![](Figure2_files/figure-markdown_strict/Figure%202c-1.png)

    # Model 2: Heterogeneous variance by Site
    growth_rate_lme_hetero <- lme(
      log(Growth) ~ Treatment,
      random = ~1 | Site,
      weights = varIdent(form = ~1 | Site),
      data = dna_prod_test
    )

    summary(growth_rate_lme_hetero)

    ## Linear mixed-effects model fit by REML
    ##   Data: dna_prod_test 
    ##        AIC      BIC    logLik
    ##   186.5721 209.0457 -82.28605
    ## 
    ## Random effects:
    ##  Formula: ~1 | Site
    ##         (Intercept) Residual
    ## StdDev:   0.8806573 0.695258
    ## 
    ## Variance function:
    ##  Structure: Different standard deviations per stratum
    ##  Formula: ~1 | Site 
    ##  Parameter estimates:
    ##    Balmoral     Bowness       Stean Moors_House    Langwell     Crocach    Migneint 
    ##   1.0000000   1.3142760   1.4316946   1.3722513   1.3088909   0.8710345   1.3210794 
    ## Fixed effects:  log(Growth) ~ Treatment 
    ##                    Value Std.Error DF    t-value p-value
    ## (Intercept)   -0.3371449 0.3874116 51 -0.8702499  0.3882
    ## TreatmentDAM   0.6051352 0.2613726 51  2.3152205  0.0247
    ## TreatmentREST  0.2765186 0.2613726 51  1.0579478  0.2951
    ##  Correlation: 
    ##               (Intr) TrtDAM
    ## TreatmentDAM  -0.368       
    ## TreatmentREST -0.368  0.535
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.3130463 -0.4653400  0.1739348  0.6123585  1.7776855 
    ## 
    ## Number of Observations: 60
    ## Number of Groups: 7

    plot(growth_rate_lme_hetero)

![](Figure2_files/figure-markdown_strict/Figure%202c-2.png)

    # Model 3: Heterogeneous variance by Treatment
    growth_rate_lme_hetero_treat <- lme(
      log(Growth) ~ Treatment,
      random = ~1 | Site,
      weights = varIdent(form = ~1 | Treatment),
      data = dna_prod_test
    )

    # Compare models using AIC
    AIC(growth_rate_lme_homo, growth_rate_lme_hetero, growth_rate_lme_hetero_treat)

    ##                              df      AIC
    ## growth_rate_lme_homo          5 177.2717
    ## growth_rate_lme_hetero       11 186.5721
    ## growth_rate_lme_hetero_treat  7 181.2430

    # ANOVA for best model
    summary(growth_rate_lme_homo)

    ## Linear mixed-effects model fit by REML
    ##   Data: dna_prod_test 
    ##        AIC     BIC    logLik
    ##   177.2717 187.487 -83.63587
    ## 
    ## Random effects:
    ##  Formula: ~1 | Site
    ##         (Intercept)  Residual
    ## StdDev:   0.8619272 0.8620367
    ## 
    ## Fixed effects:  log(Growth) ~ Treatment 
    ##                    Value Std.Error DF    t-value p-value
    ## (Intercept)   -0.3621749 0.3872240 51 -0.9353110  0.3540
    ## TreatmentDAM   0.6638070 0.2814201 51  2.3587761  0.0222
    ## TreatmentREST  0.2545121 0.2814201 51  0.9043849  0.3700
    ##  Correlation: 
    ##               (Intr) TrtDAM
    ## TreatmentDAM  -0.402       
    ## TreatmentREST -0.402  0.553
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.5090888 -0.5015986  0.1675783  0.6377418  1.9065916 
    ## 
    ## Number of Observations: 60
    ## Number of Groups: 7

    Anova(growth_rate_lme_homo, type = "III")

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: log(Growth)
    ##              Chisq Df Pr(>Chisq)  
    ## (Intercept) 0.8748  1    0.34963  
    ## Treatment   5.7949  2    0.05516 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##### Model with nested random effect: Treatment within Site
    growth_rate_lme_nested <- lme(
      log(Growth) ~ Treatment,
      random = ~1 | Site/Treatment,
      data = dna_prod_test
    )

    summary(growth_rate_lme_nested)

    ## Linear mixed-effects model fit by REML
    ##   Data: dna_prod_test 
    ##        AIC    BIC    logLik
    ##   179.2717 191.53 -83.63587
    ## 
    ## Random effects:
    ##  Formula: ~1 | Site
    ##         (Intercept)
    ## StdDev:   0.8619272
    ## 
    ##  Formula: ~1 | Treatment %in% Site
    ##          (Intercept)  Residual
    ## StdDev: 0.0001988674 0.8620367
    ## 
    ## Fixed effects:  log(Growth) ~ Treatment 
    ##                    Value Std.Error DF    t-value p-value
    ## (Intercept)   -0.3621749 0.3872240 40 -0.9353110  0.3552
    ## TreatmentDAM   0.6638070 0.2814201 11  2.3587759  0.0379
    ## TreatmentREST  0.2545121 0.2814201 11  0.9043849  0.3852
    ##  Correlation: 
    ##               (Intr) TrtDAM
    ## TreatmentDAM  -0.402       
    ## TreatmentREST -0.402  0.553
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.5090887 -0.5015986  0.1675783  0.6377417  1.9065915 
    ## 
    ## Number of Observations: 60
    ## Number of Groups: 
    ##                Site Treatment %in% Site 
    ##                   7                  20

    Anova(growth_rate_lme_nested, type = "III")

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: log(Growth)
    ##              Chisq Df Pr(>Chisq)  
    ## (Intercept) 0.8748  1    0.34963  
    ## Treatment   5.7949  2    0.05516 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plot(growth_rate_lme_nested)

![](Figure2_files/figure-markdown_strict/Figure%202c-3.png)

    ##### Fit LMMs with additional random effects
    # Multiple random effects: MAT, MAP, elevation
    growth_rate_lme_multi <- lmer(
      log(Growth) ~ Treatment + (1 | MAT) + (1 | MAP) + (1 | elevation),
      data = dna_prod_test
    )

    summary(growth_rate_lme_multi)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: log(Growth) ~ Treatment + (1 | MAT) + (1 | MAP) + (1 | elevation)
    ##    Data: dna_prod_test
    ## 
    ## REML criterion at convergence: 164.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7925 -0.4174  0.0938  0.4916  1.9293 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  elevation (Intercept) 0.0000   0.000   
    ##  MAP       (Intercept) 0.0000   0.000   
    ##  MAT       (Intercept) 0.9274   0.963   
    ##  Residual              0.7413   0.861   
    ## Number of obs: 60, groups:  elevation, 7; MAP, 7; MAT, 5
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error      df t value Pr(>|t|)  
    ## (Intercept)    -0.2032     0.4770  5.1518  -0.426   0.6874  
    ## TreatmentDAM    0.7151     0.2774 53.0095   2.578   0.0127 *
    ## TreatmentREST   0.3058     0.2774 53.0095   1.102   0.2752  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) TrtDAM
    ## TreatmntDAM -0.306       
    ## TretmntREST -0.306  0.541
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    # Single random effect: MAT
    growth_rate_lme_mat <- lmer(
      log(Growth) ~ Treatment + (1 | MAT),
      data = dna_prod_test
    )

    summary(growth_rate_lme_mat)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: log(Growth) ~ Treatment + (1 | MAT)
    ##    Data: dna_prod_test
    ## 
    ## REML criterion at convergence: 164.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7925 -0.4174  0.0938  0.4916  1.9293 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  MAT      (Intercept) 0.9274   0.963   
    ##  Residual             0.7413   0.861   
    ## Number of obs: 60, groups:  MAT, 5
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error      df t value Pr(>|t|)  
    ## (Intercept)    -0.2032     0.4770  5.1518  -0.426   0.6874  
    ## TreatmentDAM    0.7151     0.2774 53.0095   2.578   0.0127 *
    ## TreatmentREST   0.3058     0.2774 53.0095   1.102   0.2752  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) TrtDAM
    ## TreatmntDAM -0.306       
    ## TretmntREST -0.306  0.541

    # Compare models
    anova(growth_rate_lme_mat, growth_rate_lme_multi)

    ## Data: dna_prod_test
    ## Models:
    ## growth_rate_lme_mat: log(Growth) ~ Treatment + (1 | MAT)
    ## growth_rate_lme_multi: log(Growth) ~ Treatment + (1 | MAT) + (1 | MAP) + (1 | elevation)
    ##                       npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
    ## growth_rate_lme_mat      5 172.62 183.09 -81.309   162.62                    
    ## growth_rate_lme_multi    7 176.62 191.28 -81.309   162.62     0  2          1

    AIC(growth_rate_lme_mat)

    ## [1] 174.3437

    AIC(growth_rate_lme_multi)

    ## [1] 178.3437

    ##### Fit model with Treatment × Site interaction
    growth_rate_lmer_interaction <- lmer(
      log(Growth) ~ Treatment * Site + (1 | Site),
      data = dna_prod_test
    )

    summary(growth_rate_lmer_interaction)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: log(Growth) ~ Treatment * Site + (1 | Site)
    ##    Data: dna_prod_test
    ## 
    ## REML criterion at convergence: 123.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.66022 -0.55639  0.01716  0.61929  1.81024 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Site     (Intercept) 0.1927   0.4390  
    ##  Residual             0.7458   0.8636  
    ## Number of obs: 60, groups:  Site, 7
    ## 
    ## Fixed effects:
    ##                               Estimate Std. Error      df t value Pr(>|t|)   
    ## (Intercept)                     0.9978     0.6643 40.0000   1.502  0.14097   
    ## TreatmentDAM                    1.1643     0.7051 40.0000   1.651  0.10653   
    ## TreatmentREST                   1.1193     0.7051 40.0000   1.587  0.12030   
    ## SiteBowness                    -0.9277     0.9395 40.0000  -0.987  0.32935   
    ## SiteStean                      -2.9750     0.9395 40.0000  -3.167  0.00295 **
    ## SiteMoors_House                -1.7758     0.9395 40.0000  -1.890  0.06600 . 
    ## SiteLangwell                   -2.1494     0.9395 40.0000  -2.288  0.02751 * 
    ## SiteCrocach                    -1.2320     0.9395 40.0000  -1.311  0.19721   
    ## SiteMigneint                   -1.8567     0.9395 40.0000  -1.976  0.05505 . 
    ## TreatmentDAM:SiteBowness        0.1939     0.9972 40.0000   0.194  0.84682   
    ## TreatmentREST:SiteBowness      -1.3112     0.9972 40.0000  -1.315  0.19602   
    ## TreatmentDAM:SiteStean          1.3560     0.9972 40.0000   1.360  0.18151   
    ## TreatmentDAM:SiteMoors_House   -0.5541     0.9972 40.0000  -0.556  0.58155   
    ## TreatmentREST:SiteMoors_House  -1.0020     0.9972 40.0000  -1.005  0.32100   
    ## TreatmentDAM:SiteLangwell      -0.4773     0.9972 40.0000  -0.479  0.63483   
    ## TreatmentREST:SiteLangwell     -0.1436     0.9972 40.0000  -0.144  0.88625   
    ## TreatmentDAM:SiteCrocach       -1.1250     0.9972 40.0000  -1.128  0.26598   
    ## TreatmentREST:SiteCrocach      -1.2677     0.9972 40.0000  -1.271  0.21096   
    ## TreatmentDAM:SiteMigneint      -1.5000     0.9972 40.0000  -1.504  0.14038   
    ## TreatmentREST:SiteMigneint     -0.9320     0.9972 40.0000  -0.935  0.35558   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## unable to evaluate scaled gradient
    ##  Hessian is numerically singular: parameters are not uniquely determined

    Anova(growth_rate_lmer_interaction, type = "III")

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(Growth)
    ##                  Chisq Df Pr(>Chisq)  
    ## (Intercept)     2.2558  1     0.1331  
    ## Treatment       3.5003  2     0.1738  
    ## Site           14.0671  6     0.0289 *
    ## Treatment:Site 10.8602 11     0.4551  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##### Model with degradation_index as fixed effect
    growth_rate_lmer_degradation <- lmer(
      log(Growth) ~ degradation_index + (1 | Site),
      data = dna_prod
    )

    summary(growth_rate_lmer_degradation)

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
    ## Formula: log(Growth) ~ degradation_index + (1 | Site)
    ##    Data: dna_prod
    ## 
    ## REML criterion at convergence: 165.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3868 -0.5856  0.2932  0.7059  1.7374 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Site     (Intercept) 0.3838   0.6195  
    ##  Residual             0.7598   0.8717  
    ## Number of obs: 60, groups:  Site, 7
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error      df t value Pr(>|t|)   
    ## (Intercept)        -0.0451     0.2602  5.2932  -0.173  0.86883   
    ## degradation_index  -0.7268     0.2564 57.0608  -2.834  0.00634 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## degrdtn_ndx 0.014

    Anova(growth_rate_lmer_degradation, type = "III")

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(Growth)
    ##                    Chisq Df Pr(>Chisq)   
    ## (Intercept)       0.0300  1   0.862396   
    ## degradation_index 8.0339  1   0.004591 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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

    ## Marginal R² (fixed effects only): 0.145

    cat("Conditional R² (fixed + random effects):", round(R2_conditional, 3), "\n")

    ## Conditional R² (fixed + random effects): 0.481

    ##### Correlation matrix among numeric variables
    cor_data <- dna_prod_test %>%
      select(Growth, Moisture, pH, Moss_cover, FTIR_axis_1, O2, TC, TN, degradation_index)

    cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")

    # Heatmap of correlations
    ggcorrplot(cor_matrix, lab = TRUE, colors = c("blue", "white", "red"))

![](Figure2_files/figure-markdown_strict/Figure%202c-4.png)

    ##### Plot: Growth vs Ecosystem Health Index by Treatment
    dna_prod$Treatment <- factor(dna_prod$Treatment, levels = c("DAM", "REST", "NAT"))

### Microbial Growth Rate Across Ecosystem health index plot

    growth_rate_ecosystemhealth_plot <- ggplot(dna_prod, aes(x = degradation_index, y = Growth, colour = Treatment)) +
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

    growth_rate_ecosystemhealth_plot

![](Figure2_files/figure-markdown_strict/Methane%20scatter%20plot-1.png)
\### Save the plot

    ggsave("Figure_2c_MicrobialGrowthRate_EcosystemHealth.png", dpi = 1000, width = 6, height = 4.5, units = "in")
