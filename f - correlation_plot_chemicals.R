#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################  MAKE PLOT OF CORRELATIONS BETWEEN THE CHEMICALS  ##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates correlation plots of demographics, chemicals, and immune measures
#          
# Inputs: nhanes_subset     - dataframe containing complete demographic and immune measures data for each
#                             participant
#         conversion        - dataframe to convert between chemical_codename and chemical family
#         
#         subset_chemicals  - dataframe that has a vector of the chemicals to use
#
# Outputs: Plot showing correlations between chemicals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

correlation_plot_chemicals <- function(subset_chemicals,
                                       nhanes_subset,
                                       conversion)
{
  library(tidyverse)
  # library(gplots) #heatmap.2 package
  # library(Hmisc)
  # library(data.table)
  library(dichromat)
  library(pheatmap)
  
  # subset_chemicals <- use_these_chems
  # nhanes_subset <- nhanes_subset_dataset
  # conversion <- use_these_chems

  #############################################################################################################
  ######################### Clean The Chemical Data To Make It Usable For Correlations ########################
  #############################################################################################################

  #select the usable chemicals
  chems <- subset_chemicals$chemical_codename_use
  
  #select those chemicals from the wide nhanes dataset,
  #log transform values
  nhanes_subset_chems <- nhanes_subset %>%
    dplyr::select(all_of(chems)) %>%
    log2()
  
  # Rename chemical familiy variable for the plot legend later
  subset_chemicals$`Chemical Family` <- subset_chemicals$chem_family
  
  #############################################################################################################
  ############################################# Set Up For Plotting ###########################################
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
  
  # Ensure that the levels of the chemical family are in a defined order to ensure proper color scheme
  subset_chemicals$`Chemical Family` <- factor(subset_chemicals$`Chemical Family`
                                               , levels = chem_family_levels)

  #make sure the dataset to color the chemicals is in the right order
  chem_fam_reorder <- subset_chemicals %>%
    arrange(`Chemical Family`)
  chem_order <- chem_fam_reorder$chemical_codename_use
  
  #make the chemicals in the right order too
  chem_cor_order <- nhanes_subset_chems %>%
    dplyr::select(all_of(chem_order))
  
  #############################################################################################################
  ####################################### Calculate Chemical Correlations #####################################
  #############################################################################################################
  
  #calculate the correlations
  chem_correlations <- cor(chem_cor_order,
                           use = "pairwise.complete.obs",
                           method = "spearman")
  
  ################################### Create Supplemental Table of Results ####################################
  
  # Create a dataframe of the column names and chem names for the correlations
  identical(colnames(chem_correlations), rownames(chem_correlations))
  conversion_update <- conversion %>%
    rename(chem_colnames = chemical_codename_use) %>%
    dplyr::select(chem_colnames, chemical_name) %>%
    mutate(chemical_name = gsub("\\s\\(([^()]+)\\)$"
                                , ""
                                , chemical_name))
  chem_colnames <- colnames(chem_correlations)
  chem_colnames_df <- data.frame(chem_colnames)
  
  # Attach the chemical names to the codenames in the order they appear in the correlation matrix
  chem_names <- left_join(chem_colnames_df, conversion_update, by = "chem_colnames")
  
  # Replace the column and row names of the chemical correlation matrix with the actual chemical names
  chem_correlations_copy <- chem_correlations
  colnames(chem_correlations_copy) <- chem_names$chemical_name
  rownames(chem_correlations_copy) <- chem_names$chemical_name
  
  # Save matrix as a csv file
  setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))
  write.csv(chem_correlations_copy, "chem_correlations_matrix.csv")
  
  
  #############################################################################################################
  
  
  # Convert the chemical codnames to names for the supplemental table
  # identical(colnames(chem_correlations), conversion$chemical_codename_use)
  
  # identical(colnames(chem_correlations), colnames(chem_cor_order))
  #TRUE - the correlations are in the same order as the chem families dataset
  
  # Compute a matrix of correlation p-values - ggcorrplot package
  # p.mat <- cor_pmat(chem_correlations)

  # info from: https://sites.ualberta.ca/~ahamann/teaching/graphics/LabHM.pdf
  
  #since the correlations are outputted as a duplicated square matrix rather than triangle, 
  #need to select only one triangle of data
  # lower_tri <- function(chem_correlations){chem_correlations[upper.tri(chem_correlations)] <- NA
  # return(chem_correlations)}
  # bottom_tri_correlations <- lower_tri(chem_correlations)
  
  # setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))
  # write.csv(chem_correlations, "chem_correlations_matrix.csv")

  #############################################################################################################
  ############################################# Set Up For Plotting ###########################################
  #############################################################################################################
  
  #set up to make the chemical family color bars
  chem_fam_names <- chem_fam_reorder %>%
    dplyr::select(`Chemical Family`)
  row.names(chem_fam_names) <- chem_fam_reorder$chemical_codename_use

  chem_fam_names <- as.data.frame(chem_fam_names)

  #identify where to put the gaps between chemical families
  num_chem_fam <- chem_fam_names %>%
    group_by(`Chemical Family`) %>%
    summarise(number = n()) %>%
    ungroup()

  #calculate where the breaks in the plot should go
  num_chem_fam$breaks <- cumsum(num_chem_fam$number)

  #assign the colors to chem families
  chem_family_colors = list(`Chemical Family` = c(Acrylamide = "#8B0000",
                                            # , "Melamine"
                                            `Brominated Flame Retardants (BFR)` = "#EE0000",
                                            `Phosphate Flame Retardants (PFR)` = "#FF6B00",
                                            `Polychlorinated Biphenyls (PCB)` = "#FF69B4",
                                            `Dioxins` = "#FFA500",
                                            `Furans` = "#EEEE00",
                                            `Metals` = "#228B22",
                                            `Phthalates & Plasticizers` = "#A4D3EE",
                                            `Personal Care & Consumer Product Compounds` = "#A2CD5A",
                                            `Pesticides` = "#1E90FF",
                                            `Aromatic Amines` = "#be67c9",
                                            # `Phytoestrogens` = "#7D26CD",
                                            `Polyaromatic Hydrocarbons (PAH)` = "#cf9b76",
                                            `Volatile Organic Compounds (VOC)` = "#828282",
                                            `Smoking Related Compounds` = "#8B4513",
                                            `Per- and Polyfluoroalkyl Substances (PFAS)` = "#FFB6C1",
                                            `Aldehydes` = "#0E1171",
                                            `Other` = "#BABABA"))

  #shorten this terrible name
  chem_fam_reorder$chemical_name <- gsub("N-Acetyl-S-(2-hydroxy-3-methyl-3-butenyl)-L-cysteine + N-Acetyl-S-(2-hydroxy-2-methyl-3-butenyl)-L-cysteine (ng/mL)",
                                         "N-Acetyl-S-(2-hydroxy-2/3-methyl-3-butenyl)-L-cysteine (ng/mL)",
                                         chem_fam_reorder$chemical_name,
                                         fixed = TRUE)

  #this drops the units from the chemical names
  chem_fam_reorder$chemical_name <- gsub("\\s\\(([^()]+)\\)$"
                                         , ""
                                         , chem_fam_reorder$chemical_name)

  #############################################################################################################
  ############################################ Make The Heatmap Plot ##########################################
  #############################################################################################################

  #set directory to correlation plots folder
  setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))

  #https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
  paletteLength <- 100
  myColor <- colorRampPalette(c("blue", "white", "yellow", "red"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(chem_correlations, na.rm = T), 0, length.out=ceiling(paletteLength/2.8) + 1),
                seq(max(chem_correlations, na.rm = T)/paletteLength,
                    max(chem_correlations, na.rm = T),
                    length.out=floor(paletteLength/1.8)))


  pdf("chemical_heatmap_correlation_smk.pdf", width = 14, height = 9)
  # pheatmap(mat = chem_correlations,
  #          cluster_rows = FALSE, cluster_cols = FALSE,
  #          annotation_col = chem_fam_names,
  #          annotation_row = chem_fam_names,
  #          labels_row = chem_fam_reorder$chemical_name, #chemical names
  #          labels_col = chem_fam_reorder$chemical_name,
  #          annotation_names_row = FALSE, #x label
  #          annotation_names_col = FALSE, #y label
  #          # angle_col = 45,
  #          fontsize_row = 2,
  #          fontsize_col = 2,
  #          gaps_col = num_chem_fam$breaks,
  #          gaps_row = num_chem_fam$breaks,
  #          annotation_colors = chem_family_colors,
  #          color=myColor,
  #          breaks=myBreaks,
  #          # cellheight=5,cellwidth=4,
  #          legend = TRUE)
  pheatmap(chem_correlations,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = chem_fam_names,
           annotation_col = chem_fam_names,
           labels_row = chem_fam_reorder$chemical_name, #chemical names
           labels_col = chem_fam_reorder$chemical_name,
           gaps_col = num_chem_fam$breaks,
           gaps_row = num_chem_fam$breaks,
           color=myColor, breaks=myBreaks,
           annotation_colors = chem_family_colors,
           annotation_names_row = FALSE, #x label
           annotation_names_col = FALSE, #y label
           fontsize_row = 2,
           fontsize_col = 2,
           legend = TRUE)
  dev.off()
  

  svg(file = "chemical_heatmap_correlation_smk.svg", width = 14, height = 9)
  pheatmap(chem_correlations,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = chem_fam_names,
           annotation_col = chem_fam_names,
           labels_row = chem_fam_reorder$chemical_name, #chemical names
           labels_col = chem_fam_reorder$chemical_name,
           gaps_col = num_chem_fam$breaks,
           gaps_row = num_chem_fam$breaks,
           color=myColor, breaks=myBreaks,
           annotation_colors = chem_family_colors,
           annotation_names_row = FALSE, #x label
           annotation_names_col = FALSE, #y label
           fontsize_row = 2,
           fontsize_col = 2,
           legend = TRUE)
  dev.off()

  #############################################################################################################
  #############################################################################################################
  #############################################################################################################
  library(beepr)
  beep()
  setwd(current_directory)
}

#Extra code:
#############################################################################################################
######################### Clean The Chemical Data To Make It Usable For Correlations ########################
#############################################################################################################

# #chem vs chem: log transform data, pearson correlation
#
# #select the chemicals and SEQN, log10 transform the chemical measurements
# long_chem_log <- long_nhanes_subset %>%
#                       select(SEQN,
#                              chemical_codename,
#                              chem_measurement) %>%
#                       mutate(chem_measurement = log10(chem_measurement))
#
# #merge in the chemical names
# long_chem_dataset_names <- left_join(long_chem_log, subset_chemicals, by = "chemical_codename")
#
# #remove the unnecessary variables
# long_chem_dataset_names <- long_chem_dataset_names %>%
#                             select(-above,
#                                    -total,
#                                    -percent_above_LOD)
#
# #remove the old long datasets
# rm(long_chem_dataset)
#
# #make a dataset of the chem classes
# chem_classes <- conversion %>%
#   select(chem_family,
#          chemical_name)
#
# #merge in the chemical classes to the long dataset
# long_chem <- left_join(long_chem_dataset_names, chem_classes, by = "chemical_name")
# # chemical_order <- unique(long_chem$chemical_name)
#
# #sort by chemical class and then drop the chemical class variable
# long_chem_sorted <- long_chem %>%
#   arrange(chem_family) %>%
#   select(-chem_family)
#
# #get the order of chemicals because spread() will remove the order
# chemical_order <- unique(long_chem_sorted$chemical_name)
#


# upper_tri <- function(chem_correlations){chem_correlations[lower.tri(chem_correlations)]<- NA
# return(chem_correlations)}

# #############GGplot method##############
#
# # up_tri
#
# melted_cor <- melt(bottom_tri, na.rm = TRUE)

# ggplot(data = melted_cor, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "red", high = "blue", mid = "white",
#                        midpoint = 0, limit = c(-1,1), space = "Lab",
#                        name="Pearson\nCorrelation") +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
#   theme(legend.justification = c(1, 0),
#         legend.position = c(0.6, 0.7),
#         legend.direction = "horizontal") +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1))+
#   scale_y_discrete(position = "right")+
#   coord_fixed()


#############################################################################################################
###################################### Create Chemical Correlation Plot #####################################
#############################################################################################################

#CORRPLOT VERSION OF PLOT

# correlation_directory <- paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals")
#
# setwd(correlation_directory)
#
# png(filename = "chemical_correlation.png", width = 2000, height = 2000)
# corrplot(chem_correlations,
#          method = "ellipse",
#          type = "upper",
#          # order = "FPC",
#          tl.col = "black",
#          tl.cex = 0.9,
#          cl.cex = 1,
#          cl.align = "r",
#          diag = FALSE,
#          na.label = ".",
#          mar = c(0.1, 0.0001, 0.0001, 0.1))
# dev.off()

#############################################################################################################
###################################### Create Chemical Correlation Plot #####################################
#############################################################################################################

#HEATMAP.2 VERSION OF PLOT

# pdf("chemical_heatmap_correlation.pdf", width = 15, height = 15)
# heatmap.2(bottom_tri_correlations,
#           dendrogram = "none",
#           col=redblue(256),
#           na.color = "black",
#           keysize = 0.5,
#           # lhei = c(1,4),
#           # density.info="none", trace="none",
#           srtCol = 65,
#           ColSideColors = chem_family_colors[as.factor(subset_chemical_family$chem_family)],
#           RowSideColors = chem_family_colors[as.factor(subset_chemical_family$chem_family)],
#           labCol = subset_chemical_family$chemical_name,
#           labRow = subset_chemical_family$chemical_name,
#           # offsetRow = -60,
#           # mtext(subset_chemical_family$chemical_name, side = 4),
#           # distfun = mydist, hclustfun = myfun,
#           trace = "none",
#           cex.main = 3,
#           margins = c(15,15),
#           Rowv = FALSE, #don't reorder the rows
#           Colv = FALSE
#           )
# legend(y = 1.1, x = 0.25, xpd = TRUE,
#        legend = unique(subset_chemical_family$chem_family),
#        col = chem_family_colors,
#        lty = 1,
#        lwd = 5,
#        cex = 0.7)
# dev.off()