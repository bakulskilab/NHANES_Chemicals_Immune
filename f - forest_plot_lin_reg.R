#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################  FOREST PLOTS  ###############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates forest plots of 196 chemicals and cell types
#          
# Inputs:   model_stats - tidy output of linear regression stats adjusted for demographics
#           conversion - dataframe of chemicals to use and info about them
#
# Outputs:  forest plots

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

forest_plot_lin_reg <- function(model_stats,
                                conversion)
{
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  
  setwd(current_directory)

  #TEMPORARY  
  conversion <- use_these_chems
  model_stats <- model_stats_smk_scaled

  #############################################################################################################
  ########################################## FIX UP LINEAR REGRESSIONS ########################################
  #############################################################################################################

  #keep only the measurement values
  model_stats_scaled <- model_stats %>%
    filter(term == "chem_log_measurement") %>%
    mutate(immune_measure = case_when(immune_measure == "Mean Corpuscular Volume (fL)" ~ "MCV (fL)",
                                      TRUE ~ as.character(immune_measure)))
  
  #############################################################################################################
  ######################################## SCALE THE ESTIMATE AND STDEV #######################################
  #############################################################################################################

  #calculate the confidence intervals
  z_score <- 1.96
  model_stats_CI <- model_stats_scaled %>%
    mutate(lower = estimate - (z_score*std.error),
           upper = estimate + (z_score*std.error))

  #calculate FDR just to have to number for reporting
  model_stats_CI <- model_stats_CI %>%
    group_by(celltype_codename) %>%
    mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
  
#############################################################################################################
###################################### RENAME CHEMICALS AND CELL TYPES ######################################
#############################################################################################################  

  conversion_subset <- conversion %>%
    dplyr::select(chemical_codename_use,
                  chemical_name,
                  chem_family) %>%
    rename(chemical_codename = chemical_codename_use)
  
  #merge in the chemical_names
  merge_by <- c("chemical_codename", "chem_family", "chemical_name")
  model_stats_names <- left_join(model_stats_CI, conversion_subset, by = merge_by)
  
  #make a conversion for the cell types into names
  celltype_codename <- c("LBXLYPCT", #lymphocytes
                         "LBXMOPCT", #monocytes
                         "LBXNEPCT", #neutrophils
                         "LBXEOPCT", #eosinophils
                         "LBXBAPCT", #basophils
                         "LBXWBCSI", #WBC count
                         "LBXRBCSI", #RBC count
                         "LBXMCVSI" #MCV
                        )
  cell_name <- c("Lymphocytes (%)",
                 "Monocytes (%)",
                 "Neutrophils (%)",
                 "Eosinophils (%)",
                 "Basophils (%)",
                 "WBC (1000 cells/uL)",
                 "RBC (million cells/uL)",
                 "MCV (fL)")
  cell_conversion <- as.data.frame(cbind(celltype_codename,
                                         cell_name))
  cell_conversion$celltype_codename <- as.factor(cell_conversion$celltype_codename)
  
  #merge in the cell type names
  model_stats_names$celltype_codename <- as.factor(model_stats_names$celltype_codename)
  model_stats_names_cells <- left_join(model_stats_names, cell_conversion, by = "celltype_codename")
  
  #set up the order of the facets
  model_stats_names_cells$cell_name <- factor(model_stats_names_cells$cell_name,
                                              levels = c("Lymphocytes (%)",
                                                         "Neutrophils (%)",
                                                         "Monocytes (%)",
                                                         "Basophils (%)",
                                                         "Eosinophils (%)",
                                                         "WBC (1000 cells/uL)",
                                                         "RBC (million cells/uL)",
                                                         "MCV (fL)"))
  
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
  
  # Ensure that the levels of the chemical family are in a defined order to ensure proper color scheme
  model_stats_names_cells$chem_family <- factor(model_stats_names_cells$chem_family
                                         , levels = chem_family_levels)
  
  #this drops the units from the chemical names
  model_stats_names_cells$chemical_name <- gsub("\\s\\(([^()]+)\\)$"
                                                , ""
                                                , model_stats_names_cells$chemical_name)
  
  #shorten this terrible name
  model_stats_names_cells$chemical_name <- gsub("N-Acetyl-S-(2-hydroxy-3-methyl-3-butenyl)-L-cysteine + N-Acetyl-S-(2-hydroxy-2-methyl-3-butenyl)-L-cysteine",
                                           "N-Acetyl-S-(2-hydroxy-2/3-methyl-3-butenyl)-L-cysteine",
                                           model_stats_names_cells$chemical_name,
                                           fixed = TRUE)
  
  #add a numerical column to indicate the chemical families for sorting
  model_stats_names_cells <- model_stats_names_cells %>%
    mutate(num_chem_family =
             case_when(chem_family == "Acrylamide" ~ 1
                       # , "Melamine"
                       , chem_family == "Brominated Flame Retardants (BFR)" ~ 2
                       , chem_family == "Phosphate Flame Retardants (PFR)" ~ 3
                       , chem_family == "Polychlorinated Biphenyls (PCB)" ~ 4
                       , chem_family == "Dioxins" ~ 5
                       , chem_family == "Furans" ~ 6
                       , chem_family == "Metals" ~ 7
                       , chem_family == "Phthalates & Plasticizers" ~ 8
                       , chem_family == "Personal Care & Consumer Product Compounds" ~ 9
                       , chem_family == "Pesticides" ~ 10
                       , chem_family == "Aromatic Amines" ~ 11
                       # , chem_family == "Phytoestrogens" ~ 12
                       , chem_family == "Polyaromatic Hydrocarbons (PAH)" ~ 13
                       , chem_family == "Volatile Organic Compounds (VOC)" ~ 14
                       , chem_family == "Smoking Related Compounds" ~ 15
                       , chem_family == "Per- and Polyfluoroalkyl Substances (PFAS)" ~ 16
                       , chem_family == "Aldehydes" ~ 17
                       # , "Dietary Components"
                       , chem_family == "Other" ~ 18))
  
  ###########################################################################################################
  ################################# REMOVE CHEMICALS WITH LARGE ESTIMATES ###################################
  ###########################################################################################################
  
  # Vector of chemicals with large estimates
  larger_estimate_chems <- c("LBXBSE",
                             "LBXSCU",
                             "URXHCTT",
                             "URXCOTT",
                             "URXNAL",
                             "LBXFOR",
                             "LBX4AL",
                             "LBX3AL")
  
  #make a dataframe of chemicals that have large estimates for separate plotting
  model_stats_large_est <- model_stats_names_cells %>%
    filter(chemical_codename %in% larger_estimate_chems) %>%
    droplevels(.)
  
  # Re-define a string vector of color hexcodes for the chemical family in corresponding order
  chem_family_colors_large <- c("#228B22",    # Metals
                                "#8B4513",    # Smoking
                                "#0E1171")    # Aldehydes
  # Re-define a string vector of shape codes for the chemical family in corresponding order
  chem_family_shapes_large <- c(18,    # Metals
                                16,    # Smoking
                                16)    # Aldehydes
  
  #make a dataframe of the rest of the chemicals
  model_stats_small_est <- model_stats_names_cells %>%
    filter(!chemical_codename %in% larger_estimate_chems) %>%
    droplevels(.)


  ###########################################################################################################
  ############################### MAKE FOREST PLOT BY CHEMICAL AND CELL TYPE ################################
  ###########################################################################################################
  
  setwd(paste0(current_directory, "/Forest Plot - Chemicals by Cell Types"))
  
  
  forest_plot_by_chem_fam_small <-
    ggplot(data = model_stats_small_est,
           aes(x = estimate,
               y = reorder(chemical_name, -num_chem_family),
               xmin = lower, #95% CI
               xmax = upper,
               color = chem_family,
               shape = chem_family)) +
    geom_pointrange(aes(col = chem_family),
                    position = position_dodge2(width = 1),
                    size = 0.3) +
    geom_vline(xintercept = 0,
               linetype = 2) +
    theme_bw() +
    #axes
    xlab('')+
    ylab('')+
    theme(axis.title.x = element_text(size = 10))+
    theme(axis.text.x = element_text(size = 8))+
    theme(axis.text.y = element_text(size = 3))+
    #legend
    scale_color_manual(name = "Chemical Family",
                       values = chem_family_colors)+
    scale_shape_manual(name = "Chemical Family",
                       values = chem_family_shapes)+
    guides(color = guide_legend(nrow = 3))+
    theme(legend.position = "top",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))+
    guides(colour = guide_legend(override.aes = list(size=1)))+
    #facet
    facet_wrap(vars(cell_name),
               ncol = 9,
               scales = "free_x")+
    theme(strip.text = element_text(size=10,
                                    face = "bold"))+
    theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  forest_plot_by_chem_fam_large <-
    ggplot(data = model_stats_large_est,
           aes(x = estimate,
               y = reorder(chemical_name, -num_chem_family),
               xmin = lower, #95% CI
               xmax = upper,
               color = chem_family,
               shape = chem_family)) +
    geom_pointrange(aes(col = chem_family),
                    position = position_dodge2(width = 1),
                    size = 0.4) +
    geom_vline(xintercept = 0,
               linetype = 2) +
    theme_bw() +
    #axes
    xlab('Effect Estimate')+
    ylab('')+
    theme(axis.title.x = element_text(size = 10))+
    theme(axis.text.x = element_text(size = 8))+
    theme(axis.text.y = element_text(size = 3))+
    #legend
    scale_color_manual(name = "Chemical Family",
                       values = chem_family_colors_large)+
    scale_shape_manual(name = "Chemical Family",
                       values = chem_family_shapes_large)+
    guides(color = guide_legend(nrow = 3))+
    theme(legend.position = "none",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))+
    guides(colour = guide_legend(override.aes = list(size=1)))+
    #facet
    facet_wrap(vars(cell_name),
               ncol = 9,
               scales = "free_x")+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())+
    theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))

  ###########################################################################################################
  ############################################### SAVE PLOTS ################################################
  ###########################################################################################################
  
  combined_forest <- plot_grid(forest_plot_by_chem_fam_small,
                               forest_plot_by_chem_fam_large,
                               align = "v", nrow = 2, rel_heights = c(1/2, 1/28))
  
  print("forest_plot_weighted_scaled.pdf")
  save_plot(filename = "forest_plot_weighted_scaled.pdf",
            plot = combined_forest,
            base_width = 14,
            base_height = 9)
  
  print("forest_plot_weighted_scaled.png")
  save_plot(filename = "forest_plot_weighted_scaled.png",
            plot = combined_forest,
            base_width = 14,
            base_height = 9,
            dpi = 1200)
  
  #############################################################################################################
  ################################# CALCULATE MEAN EFFECT ESTIMATE BY CHEM CLASS ##############################
  #############################################################################################################
  
  model_stats_wbc <- model_stats_names_cells %>%
    filter(celltype_codename == "LBXWBCSI")
  
  model_stats_avg <- model_stats_wbc %>%
    group_by(chem_family) %>%
    summarise(mean_est = mean(estimate)) %>%
    ungroup() %>%
    mutate(converted_est = mean_est*1000)
  
  #############################################################################################################
  
  setwd(current_directory)
}