correlation_plot_unweight_weight <- function(model_stats_wt_adj,
                                             model_stats_wt_unadj)
{
  #PURPOSE: scaled cotinine adj vs scaled non-cotinine adj, both weighted
  
  library(tidyverse)
  library(ggrepel)
  # library(rstatix)
  library(broom)
  
  #TEMPORARY
  # model_stats_wt_adj <- model_stats_smk_scaled
  # model_stats_wt_unadj <- model_stats_scale_no_smk
  
  ###########################################################################################################
  ############################################## CLEAN DATASETS #############################################
  ###########################################################################################################
  
  # pull the estimates from each dataset
  
  #weighted and cotinine
  wt_est <- model_stats_wt_adj %>%
    filter(term == "chem_log_measurement") %>% 
    mutate(immune_measure = case_when(immune_measure == "Lymphocyte (1000 cells/uL)" ~ "Lymphocytes (1000 cells/uL)",
                                      immune_measure == "Monocyte (1000 cells/uL)" ~ "Monocytes (1000 cells/uL)",
                                      TRUE ~ immune_measure)) %>%
    dplyr::select(chemical_name,
                  chem_family,
                  immune_measure,
                  celltype_codename,
                  estimate,
                  FDR) %>%
    rename(estimate_wt = estimate,
           FDR_wt = FDR)
  #weighted and no cotinine
  wt_unadj_est <- model_stats_wt_unadj %>%
    filter(term == "chem_log_measurement") %>%
    mutate(immune_measure = case_when(immune_measure == "Lymphocyte (1000 cells/uL)" ~ "Lymphocytes (1000 cells/uL)",
                                      immune_measure == "Monocyte (1000 cells/uL)" ~ "Monocytes (1000 cells/uL)",
                                      TRUE ~ immune_measure)) %>%
    dplyr::select(chemical_name,
                  chem_family,
                  immune_measure,
                  celltype_codename,
                  estimate,
                  FDR) %>%
    rename(estimate_wt_nosmk = estimate,
           FDR_wt = FDR)
  
  # View(wt_est)
  # View(wt_unadj_est)
  
  #merge by these columns
  merge_by <- c("chemical_name",
                "chem_family",
                "immune_measure",
                "celltype_codename")
  

  ###########################################################################################################
  ############################################## MERGE DATASETS #############################################
  ###########################################################################################################
  
  #merge datasets
  wt_unadj_adj <- full_join(wt_est, wt_unadj_est, by = merge_by)
  
  # View(wt_unadj_adj)
  
  #############################################################################################################
  ################################# DEFINE COLORS FOR FOREST PLOT - CHEM CLASSES ##############################
  #############################################################################################################
  
  # Define a vector of chemical family names in a particular order
  chem_family_levels <- c("Acrylamide"
                          # , "Melamine"
                          , "Brominated Flame Retardants (BFR)"
                          , "Phosphate Flame Retardants (PFR)"
                          , "Polychlorinated Biphenyls (PCB)"
                          , "Dioxins"
                          , "Furans"
                          , "Metals"
                          , "Phthalates & Plasticizers"
                          , "Personal Care & Consumer Product Compounds"
                          , "Pesticides"
                          , "Aromatic Amines"
                          # , "Phytoestrogens"
                          , "Polyaromatic Hydrocarbons (PAH)"
                          , "Volatile Organic Compounds (VOC)"
                          , "Smoking Related Compounds"
                          , "Per- and Polyfluoroalkyl Substances (PFAS)"
                          , "Aldehydes"
                          # , "Dietary Components"
                          , "Other")
  
  # Define a string vector of color hexcodes for the chemical family in corresponding order
  chem_family_colors <- c("#8B0000"      # Acrylamide
                          # , "#9b870c"    # Melamine
                          , "#EE0000"    # BFRs
                          , "#FF6B00"    # PFRs
                          , "#FF69B4"    # PCBs
                          , "#FFA500"    # Dioxins
                          , "#EEEE00"    # Furans
                          , "#228B22"    # Metals
                          , "#A4D3EE"    # Phthalates & Plasticizers
                          , "#A2CD5A"    # Personal Care
                          , "#1E90FF"    # Pesticides
                          , "#be67c9"    # Aromatic Amines
                          # , "#7D26CD"    # Phytoestrogens
                          , "#cf9b76"    # PAHs
                          , "#828282"    # VOCs
                          , "#8B4513"    # Smoking
                          , "#FFB6C1"    # PFCs
                          , "#0E1171"    # Aldehydes
                          , "#BABABA" )  # Other
  
  # Define a string vector of shape codes for the chemical family in corresponding order
  chem_family_shapes <- c(16      # Acrylamide
                          # , "#9b870c"    # Melamine
                          , 16    # BFRs
                          , 16    # PFRs
                          , 16    # PCBs
                          , 16   # Dioxins
                          , 16    # Furans
                          , 18    # Metals
                          , 16    # Phthalates & Plasticizers
                          , 16    # Personal Care
                          , 16    # Pesticides
                          , 17    # Aromatic Amines
                          # , 16    # Phytoestrogens
                          , 16    # PAHs
                          , 15    # VOCs
                          , 16    # Smoking
                          , 16    # PFCs
                          , 16    # Aldehydes
                          , 25 )  # Other
  
 
  ############ Second Plot ############ 
  
  # Ensure that the levels of the chemical family are in a defined order to ensure proper color scheme
  wt_unadj_adj$chem_family <- factor(wt_unadj_adj$chem_family
                                         , levels = chem_family_levels)
  
  #this drops the units from the chemical names
  wt_unadj_adj$chemical_name <- gsub("\\s\\(([^()]+)\\)$"
                                         , ""
                                         , wt_unadj_adj$chemical_name)
  
  #shorten this terrible name
  wt_unadj_adj$chemical_name <- gsub("N-Acetyl-S-(2-hydroxy-3-methyl-3-butenyl)-L-cysteine + N-Acetyl-S-(2-hydroxy-2-methyl-3-butenyl)-L-cysteine",
                                         "N-Acetyl-S-(2-hydroxy-2/3-methyl-3-butenyl)-L-cysteine",
                                     wt_unadj_adj$chemical_name,
                                         fixed = TRUE)
  
  #set up the order of the facets
  wt_unadj_adj$immune_measure <- factor(wt_unadj_adj$immune_measure,
                                            levels = c("Lymphocytes (1000 cells/uL)",
                                                       "Neutrophils (1000 cells/uL)",
                                                       "Monocytes (1000 cells/uL)",
                                                       "Basophils (1000 cells/uL)",
                                                       "Eosinophils (1000 cells/uL)",
                                                       "WBC (1000 cells/uL)",
                                                       "RBC (million cells/uL)",
                                                       "Mean Corpuscular Volume (fL)"))
  
  ###########################################################################################################
  ########################################### MAKE 2ND SCATTERPLOT ##########################################
  ###########################################################################################################

  wt_unadj_adj_plot <-
    ggplot(data = wt_unadj_adj,
           aes(x = estimate_wt_nosmk,
               y = estimate_wt,
               color = chem_family,
               shape = chem_family,
               label = chemical_name)) +
    geom_point(aes(col = chem_family),
               size = 3) +
    geom_abline(intercept = 0,
                slope = 1)+
    ggrepel::geom_text_repel(aes(label = chemical_name,
                                 size = 3.5),
                             box.padding = unit(0.5, "lines"),
                             point.padding = unit(0.2, "lines"),
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 35),
                             show.legend = FALSE)+
    theme_bw() +
    #axes
    xlab('Weighted, Not Cotinine Adjusted Beta Coefficients')+
    ylab('Weighted, Cotinine Adjusted Beta Coefficients')+
    theme(axis.title = element_text(size = 15))+
    theme(axis.text = element_text(size = 12))+
    #legend
    scale_color_manual(name = "Chemical Family",
                       values = chem_family_colors)+
    scale_shape_manual(name = "Chemical Family",
                       values = chem_family_shapes)+
    guides(color = guide_legend(nrow = 3))+
    theme(legend.position = "top",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    #facet
    facet_wrap(vars(immune_measure),
               ncol = 2,
               nrow = 4,
               scales = "free")+
    theme(strip.text = element_text(size = 15,
                                    face = "bold"))



  setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))

  # Save the plot as a pdf for viewing at a high resolution
  print("weighted cotinine vs no cotinine estimates.pdf")
  ggsave(filename = "new weighted cotinine vs no cotinine estimates_new.pdf"
         , plot = wt_unadj_adj_plot
         , width = 14
         , height = 9)

  # Save the plot as a png for presentation
  print("weighted cotinine vs no cotinine estimates.png")
  ggsave(filename = "new weighted cotinine vs no cotinine estimates_new.png"
         , plot = wt_unadj_adj_plot
         , units = "in"
         , width = 14
         , height = 9
         , dpi = 900)
  
  setwd(current_directory)
  
  ###########################################################################################################
  ############################################### CORRELATIONS ##############################################
  ###########################################################################################################

  #calculate correlations between weighted cotinine and weighted no cotinine for each set of immune measures
  cor_chems <- wt_unadj_adj %>%
    group_by(immune_measure) %>%
    summarise(correlations = cor(estimate_wt_nosmk, estimate_wt)) %>%
    ungroup()
  print("correlations between weighted cotinine and no cotinine")
  print(cor_chems)
  print(mean(cor_chems$correlations))
  
  
  
  
  
  # CODE FROM: https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
  #get correlations of old and new models to compare the estimates by immune measure
  data_test <- wt_unadj_adj %>%
    dplyr::select(-chem_family,
                  -celltype_codename,
                  -FDR_wt.x,
                  -FDR_wt.y)
  
  #nest the data into lists
  nested <- group_by(data_test, immune_measure) %>% nest()
  
  # nested %>% 
    # mutate(test = map(data, ~ cor.test(.x$estimate_wt, .x$estimate_wt_nosmk)))
  
  cor_df <- nested %>% 
    mutate(test = map(data, ~ cor.test(.x$estimate_wt,
                                       .x$estimate_wt_nosmk,
                                       method = "pearson")), # S3 list-col
           tidied = map(test, tidy)) %>% 
    unnest(tidied) %>%
    as.data.frame()
  
  print("correlations and p-values")
  cor_df %>%
    dplyr::select(immune_measure,
                  estimate,
                  p.value) %>%
    print()
  
}