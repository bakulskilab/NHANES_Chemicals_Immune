#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##################################  Graph The Chemical Measurements By LOD  ###################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates and saves barplot and scatterplot of LOD info per chemical
#          
# Inputs:   subset_chemicals - dataframe of chemicals that passed the LOD check
#
# Outputs:  barplot and scatterplot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

limit_of_detection_bar_plots <- function(subset_chemicals)
{
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  
  setwd(current_directory)
  
  #temporary
  # subset_chemicals <- use_these_chems

  #############################################################################################################
  #################################### Make a Below Column in subset_chemicals ################################
  #############################################################################################################
  
  #make below column and move it after above column - relocate() doesn't work
  subset_chemicals_above_below <- subset_chemicals %>%
    mutate(below = total - above)
  subset_chemicals_above_below <- subset_chemicals_above_below %>%
    dplyr::select(chemical_codename_use,
                  chemical_name,
                  chem_family,
                  above,
                  below,
                  total,
                  percent_above_LOD)
  
  #remove the units from subset_chemicals for scatterplot later
  #this drops the units from the chemical names
  subset_chemicals$chemical_name <- gsub("\\s\\(([^()]+)\\)$"
                                         , ""
                                         , subset_chemicals$chemical_name)
  
  #############################################################################################################
  ############################################## Make a Long Dataset ##########################################
  #############################################################################################################
  
  long_subset_chemicals <- gather(data = subset_chemicals_above_below, #this is the wide dataset
                                  key = over_under, #this is the new column to describe the number
                                  value = count, #these are the counts per chemical
                                  above:below #these are the columns to adjust
                                 )
  # long_subset_chemicals$chemical_name <- str_wrap(long_subset_chemicals$chemical_name, width = 75)
  
  #define colors
  mycolors <- c("darkmagenta", "grey50")
  
  #this drops the units from the chemical names
  long_subset_chemicals$chemical_name <- gsub("\\s\\(([^()]+)\\)$"
                                              , ""
                                              , long_subset_chemicals$chemical_name)
  
  #shorten this terrible name, fixed=T is the "did I stutter mode"
  long_subset_chemicals$chemical_name <- gsub("N-Acetyl-S-(2-hydroxy-3-methyl-3-butenyl)-L-cysteine + N-Acetyl-S-(2-hydroxy-2-methyl-3-butenyl)-L-cysteine",
                                              "N-Acetyl-S-(2-hydroxy-2/3-methyl-3-butenyl)-L-cysteine",
                                              long_subset_chemicals$chemical_name,
                                              fixed = TRUE)
  
  #############################################################################################################
  ######################################### Make The Barplot By Chemical ######################################
  #############################################################################################################
  
  chem_by_partic_by_LOD <-
    ggplot(data = long_subset_chemicals,
           aes(x = count,
               y = reorder(chemical_name, count)))+
    geom_col(aes(fill = over_under), width = 0.6,
             position = position_stack(reverse = TRUE))+
    labs(x = "Participants", y = "")+
    scale_fill_manual(values=mycolors,
                      name = "Limit of Detection")+
    scale_x_continuous(expand = c(0,0))+
    theme(axis.text.y = element_text(colour = "black",
                                     # angle = 90,
                                     vjust = 0.5, hjust=1,
                                     size = 4),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 10)) +
    theme(panel.background = element_rect(fill = "white", #this is the background
                                          colour = "black", #this is the border
                                          size = 0.1, linetype = "solid"))+
    theme(panel.grid.major.x = element_line(color = "grey"),
          panel.grid.minor.x = element_line(color = "grey",
                                            linetype = "dashed"))
  
  
  setwd(paste0(current_directory, "/Bar Plots - Sample Size by Toxicant"))
  # Save the plot as a pdf for viewing at a high resolution
  print("chemicals by LOD numbers_weighted.pdf")
  ggsave(filename = "chemicals by LOD numbers_weighted.pdf"
         , plot = chem_by_partic_by_LOD
         , width = 9
         , height = 9)
  
  # Save the plot as a png for presentation
  print("chemicals by LOD numbers_weighted.png")
  ggsave(filename = "chemicals by LOD numbers_weighted.png"
         , plot = chem_by_partic_by_LOD
         , units = "in"
         , width = 9
         , height = 9
         , dpi = 500)
  
  setwd(current_directory)
  
  #############################################################################################################
  ####################################### Make The Scatterplot By Chemical ####################################
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
  # Redefine the column vector containing the chemical family as a factor with the levels 
  subset_chemicals$chem_family <- factor(subset_chemicals$chem_family
                                      , levels = chem_family_levels)
  
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
                          # , "#49FA0F"    # Dietary Components
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
  
  #Scatterplot of # participants vs. % over LOD
  scatter_partic_by_chem_LOD <-
    ggplot(subset_chemicals,
           aes(x = total,
               y = percent_above_LOD,
               color = chem_family,
               shape = chem_family))+
    geom_point(size = 2)+
    scale_color_manual(name = "Chemical Family",
                       values = chem_family_colors) +
    scale_shape_manual(name = "Chemical Family",
                       values = chem_family_shapes)+
    geom_text(aes(label=ifelse(total > 30000, as.character(chemical_name),'')),
              vjust = -1,
              size = 4,
              show.legend = FALSE)+
    geom_text_repel(aes(label=ifelse(percent_above_LOD < 55 & total < 20000,
                               as.character(chemical_name),'')),
              hjust = -0.1,
              size = 4,
              show.legend = FALSE)+
    labs(x = "Participants per Chemical",
         y = "Measurements Above LOD (%)")+
    theme(axis.text.x = element_text(colour = "black",
                                     size = 15),
          axis.text.y = element_text(colour = "black",
                                     size = 15))+
    theme(axis.title = element_text(size = 20))+
    theme(legend.text = element_text(size = 15),
          legend.title = element_text(size = 20))+
    theme(panel.background = element_rect(fill = "white",
                                          colour = "black", #this is the border
                                          size = 0.1, linetype = "solid"))
  
  
  setwd(paste0(current_directory, "/Bar Plots - Sample Size by Toxicant"))
  # Save the plot as a pdf for viewing at a high resolution
  print("participants and percent over LOD_weighted.pdf")
  ggsave(filename = "participants and percent over LOD_weighted.pdf"
         , plot = scatter_partic_by_chem_LOD
         , width = 14
         , height = 9)
  
  # Save the plot as a png for presentation
  print("participants and percent over LOD_weighted.png")
  ggsave(filename = "participants and percent over LOD_weighted.png"
         , plot = scatter_partic_by_chem_LOD
         , units = "in"
         , width = 14
         , height = 9
         , dpi = 300)

  #############################################################################################################
  
  setwd(current_directory)
}