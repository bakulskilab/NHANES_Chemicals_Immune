#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####################################  LINEAR REGRESSION VOLCANO PLOT  ########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This function creates a volcano plot based on the adjusted linear regression results
#          
# Inputs:   model_stats - tidy output of linear regression stats adjusted for demographics
#           conversion - dataframe of chemicals to use and info about them
#           long_nhanes_subset - long form dataframe of nhanes dataset
#
# Outputs:  volcano_plot_lin_reg.png/pdf - volcano plot colored by chem family, saved as png and pdf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

volcano_plot_lin_reg <- function(model_stats,
                                 conversion)
{
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(dichromat)
  
  #TEMPORARY
  # conversion <- use_these_chems
  # model_stats <- model_stats_smk_scaled
  
  conversion <- conversion %>%
    rename(chemical_codename = chemical_codename_use)
  
  #grab the chemical measurement estimates
  model_adjust <- model_stats %>%
    filter(term == "chem_log_measurement")
  
  # #convert fdr into -log10
  # model_adjust_log <- model_adjust %>%
  #   mutate(fdr_log10 = -log10(FDR))
  
  #add significant id column and column of labels for the immune measures
  model_unclean <- model_adjust %>%
    mutate(Significance = ifelse(.$FDR <0.05, "FDR < 0.05", "Not Significant")) %>%
    mutate(immune_labels = 
             case_when(celltype_codename == "LBDLYMNO" ~ "Lymphocytes (1000 cells/uL)",
                       celltype_codename == "LBDMONO" ~ "Monocytes (1000 cells/uL)",
                       celltype_codename == "LBDNENO" ~ "Neutrophils (1000 cells/uL)",
                       celltype_codename == "LBDEONO" ~ "Eosinophils (1000 cells/uL)",
                       celltype_codename == "LBDBANO" ~ "Basophils (1000 cells/uL)",
                       celltype_codename == "LBXWBCSI" ~ "White Blood Cells (1000 cells/uL)",
                       celltype_codename == "LBXRBCSI" ~ "Red Blood Cells (million cells/uL)",
                       celltype_codename == "LBXMCVSI" ~ "Mean Corpuscular Volume (fL)")
           )
  # View(model_unclean)
    
  
  #merge in the chemical names
  join_by <- c("chemical_codename", "chem_family", "chemical_name")
  model_clean <- left_join(model_unclean, conversion, by = join_by) %>%
    select(-chem_family_shortened,
           -comment_codename)
  
  
  #add labels for only significant points
  model_clean <- model_clean %>%
    mutate(chem_labels = ifelse(.$FDR <0.05, chemical_name, ""))
  # View(model_clean)
  
  #############################################################################################################
  ############################################# SCALE THE ESTIMATE ############################################
  #############################################################################################################

  #calculate quartiles
  # chems_cells_iqr_summary <- long_nhanes_subset %>%
  #   group_by(chemical_codename, celltype_codename) %>%
  #   summarise(IQR = IQR(chem_log_measurement))

  #merge the IQRs into the results dataset
  # chem_cells_vector <- c("chemical_codename", "celltype_codename")
  # model_clean_iqr <- left_join(model_clean, chems_cells_iqr_summary, by = chem_cells_vector)

  #scale the estimates and standard deviations
  # model_clean_iqr <- model_clean_iqr %>%
  #   mutate(estimate_scaled = IQR*estimate)
  
  # model_clean_iqr <- model_clean %>%
  #   rename(estimate_scaled = estimate)
  
  #############################################################################################################
  ############################################## Set Up Cell Types ############################################
  #############################################################################################################
  
  #make a conversion for the cell types into names
  celltype_codename <- c("LBDLYMNO", #lymphocytes
                         "LBDNENO",  #neutrophils
                         "LBDMONO",  #monocytes
                         "LBDBANO",  #basophils
                         "LBDEONO",  #eosinophils
                         "LBXWBCSI", #WBC count
                         "LBXRBCSI", #RBC count
                         "LBXMCVSI"  #MCV
                         )
  cell_name <- c("Lymphocytes (1000 cells/uL)",
                 "Monocytes (1000 cells/uL)",
                 "Neutrophils (1000 cells/uL)",
                 "Eosinophils (1000 cells/uL)",
                 "Basophils (1000 cells/uL)",
                 "WBC (1000 cells/uL)",
                 "RBC (million cells/uL)",
                 "Mean Corpuscular Volume (fL)"
                 )
  cell_conversion <- as.data.frame(cbind(celltype_codename,
                                         cell_name))
  cell_conversion$celltype_codename <- as.factor(cell_conversion$celltype_codename)

  #merge in the cell type names
  model_clean$celltype_codename <- as.factor(model_clean$celltype_codename)
  model_clean <- left_join(model_clean, cell_conversion, by = "celltype_codename")
  
  #set up the order of the facets
  model_clean$cell_name <- factor(model_clean$cell_name,
                                      levels = c("Lymphocytes (1000 cells/uL)",
                                                 "Monocytes (1000 cells/uL)",
                                                 "Neutrophils (1000 cells/uL)",
                                                 "Eosinophils (1000 cells/uL)",
                                                 "Basophils (1000 cells/uL)",
                                                 "WBC (1000 cells/uL)",
                                                 "RBC (million cells/uL)",
                                                 "Mean Corpuscular Volume (fL)"
                                                 )
                                  )
  
  #############################################################################################################
  ######################################## Set Up Chemical Family Colors ######################################
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
  
  # Redefine the column vector containing the chemical family as a factor with the levels 
  model_clean$chem_family <- factor(model_clean$chem_family
                                        , levels = chem_family_levels)
  
  model_clean <- as.data.frame(model_clean)
  
  #this drops the units from the chemical names
  model_clean$chem_labels <- gsub("\\s\\(([^()]+)\\)$"
                                        , ""
                                        , model_clean$chem_labels)
  # View(model_clean)
  
  #############################################################################################################
  ############################################### Plot Volcano ################################################
  ############################################################################################################# 
  
  setwd(paste0(current_directory, "/Volcano Plots"))
  
  
  gg_volcano_plot <-
    ggplot(data = model_clean,
           aes(x = estimate,
               y = -log10(FDR),
               label = chem_labels,
               color = chem_family,
               shape = chem_family))+
    geom_point(size = 2) +
    scale_color_manual(name = "Chemical Family"
                       , values = chem_family_colors)+
    scale_shape_manual(name = "Chemical Family",
                       values = chem_family_shapes)+
    # geom_text(check_overlap = TRUE,
    #           size = 3,
    #           vjust = 0,
    #           nudge_y = -0.5,
    #           show.legend = FALSE)+
    ggrepel::geom_text_repel(
      # label = chem_labels,
                             size = 3,
                             box.padding = unit(0.5, "lines"),
                             point.padding = unit(0.1, "lines"),
                             max.overlaps = 40,
                             show.legend = FALSE
                             )+
    theme_bw()+
    geom_vline(xintercept = 0,
               color = "black",
               linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05),
               color = "red")+
    # scale_y_log10()+
    ylab("-log10(FDR)")+
    xlab("Beta Coefficients")+
    theme(axis.title = element_text(size = 18))+
    theme(axis.text = element_text(size = 15))+
    theme(legend.position = "top",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))+
    theme(strip.text.x = element_text(size = 15))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    facet_wrap(~ cell_name,
               ncol = 2,
               nrow = 4,
               scales = "free")
  
  # Save the plot as a pdf for viewing at a high resolution
  print("volcano_plot_lin_reg_wt_smk.pdf")
  ggsave(filename = "volcano_plot_lin_reg_wt_smk_new.pdf"
         , plot = gg_volcano_plot
         , width = 14
         , height = 9)

  # Save the plot as a png for presentation
  print("volcano_plot_lin_reg_wt_smk.png")
  ggsave(filename = "volcano_plot_lin_reg_wt_smk_new.png"
         , plot = gg_volcano_plot
         , units = "in"
         , width = 14
         , height = 9
         , dpi = 900)
  
  
  setwd(current_directory)
}