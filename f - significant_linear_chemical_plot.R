#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################  Heatmap of Significant Chemicals  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose:  This function creates a heatmap of the significant chemicals per immune measure and chemical family
#          
# Inputs:   model_stats - tidy output of linear regression stats adjusted for demographics
#           conversion - dataframe of chemicals to use and info about them
#
# Outputs:  significant chemicals heatmap, quadratic significant chemicals

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

significant_linear_chemical_plot <- function(conversion,
                                             model_stats)
{

  #load libraries
  library(tidyverse)
  library(pheatmap)
  library(dichromat)
  library(usefun)
  
  #TEMPORARY
  # conversion <- use_these_chems
  # model_stats <- model_stats_smk_scaled
  
  #############################################################################################################
  ############################################ Clean Up The Datasets ##########################################
  #############################################################################################################
  
  #rename datasets for use
  conversion <- conversion %>%
    rename(chemical_codename = chemical_codename_use)
  
  #grab only the chemical measurement estimates from linear and quadratic model results
  lin_model_chem <- model_stats %>%
    filter(term == "chem_log_measurement")
  
  #keep only significant chemicals
  lin_model_sig <- lin_model_chem %>%
    filter(FDR <= 0.05)
  
  #merge in the chemical names and families
  merge_by <- c("chemical_codename", "chem_family", "chemical_name")
  lin_model_clean <- left_join(lin_model_sig, conversion, by = merge_by) %>%
    select(-term,
           -above,
           -total,
           -percent_above_LOD,
           -chem_family_shortened,
           -comment_codename)
  
  #add column of immune labels for plotting later
  lin_model_lab <- lin_model_clean %>%
    mutate(immune_labels = 
             case_when(celltype_codename == "LBXLYPCT" ~ "Lymphocytes (%)",
                       celltype_codename == "LBXMOPCT" ~ "Monocytes (%)",
                       celltype_codename == "LBXNEPCT" ~ "Neutrophils (%)",
                       celltype_codename == "LBXEOPCT" ~ "Eosinophils (%)",
                       celltype_codename == "LBXBAPCT" ~ "Basophils (%)",
                       celltype_codename == "LBXWBCSI" ~ "White Blood Cells (1000 cells/uL)",
                       celltype_codename == "LBXRBCSI" ~ "Red Blood Cells (million cells/uL)",
                       celltype_codename == "LBXMCVSI" ~ "Mean Corpuscular Volume (fL)"))
  
  
  #side quest: how many metals and smoking-related compounds are significant on at least one immune measure?
  metals <- lin_model_lab %>%
    filter(chem_family == "Metals")
  print("count of significant metals on at least one immune measure")
  print(length(unique(metals$chemical_codename)))
  #22/26
  smoking <- lin_model_lab %>%
    filter(chem_family == "Smoking Related Compounds")
  print("count of significant smoking on at least one immune measure")
  print(length(unique(smoking$chemical_codename)))
  #3/4
  phosphate <- lin_model_lab %>%
    filter(chem_family == "Phosphate Flame Retardants (PFR)")
  print("count of significant phosphate on at least one immune measure")
  print(length(unique(phosphate$chemical_codename)))
  
  #How many chemicals are significant per chemical family and what percent?
  sig_counts <- lin_model_lab %>%
    group_by(chem_family) %>%
    summarise(sig_chem_counts = n_distinct(chemical_codename)) %>%
    ungroup()
  chem_fam_count <- use_these_chems %>%
    group_by(chem_family) %>%
    summarise(chemicals_per_family = n_distinct(chemical_codename_use)) %>%
    ungroup()
  sig_chem_by_fam <- left_join(sig_counts, chem_fam_count, by = "chem_family") %>%
    mutate(perc_total = (sig_chem_counts/chemicals_per_family)*100)
  
  #How many unique chemicals are significant on any measure?
  length(unique(lin_model_lab$chemical_codename))
  #122
  
  #Are any chems associated with all immune measures?
  all_assoc <- lin_model_lab %>%
    group_by(chemical_name) %>%
    summarise(n = n()) %>%
    ungroup()
  # View(all_assoc)

  
  #What is the average effect estimate by chemical family?
  avg_effect_by_fam <- model_stats %>%
    filter(immune_measure == "WBC (1000 cells/uL)") %>%
    filter(term == "chem_log_measurement") %>%
    dplyr::select(chemical_name,
                  chem_family,
                  immune_measure,
                  estimate) %>%
    group_by(chem_family) %>%
    summarise(avg_est = mean(estimate)) %>%
    ungroup() %>%
    mutate(estimate_scaled = avg_est*1000)
  
  #############################################################################################################
  ######################################## Set Up Chemical Family Colors ######################################
  #############################################################################################################
  
  # Define a vector of chemical family names in a particular order
  chem_family_levels <- c("Total Significant Chemicals",
                          "Acrylamide"
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
  
  #############################################################################################################
  ################################## Linear Significant Per Immune and Family #################################
  ############################################################################################################# 
  
  #count the number of significant chemicals per family and immune measures
  sig_lin_chems <- lin_model_lab %>%
    group_by(immune_labels,
             chem_family) %>%
    summarise(num_chems = n()) %>%
    ungroup()
  
  #create a wide dataset and replace the NAs with zeros
  wide_lin_sig <- spread(sig_lin_chems,
                         immune_labels,
                         num_chems) %>%
    replace(is.na(.), 0)
  
  #add a blank row for BFRs because they are missing
  # bfrs <- c("Brominated Flame Retardants (BFR)", 0, 0, 0, 0, 0, 0, 0, 0)
  # wide_lin <- rbind(wide_lin_sig, bfrs)
  wide_lin_tidy <- wide_lin_sig %>% column_to_rownames("chem_family")
  
  #calculate total significant chemicals per immune measure
  sum_sig_chems <- colSums(sapply(wide_lin_tidy, as.numeric))
  
  #add a row of number of significant chemicals per immune measure
  wide_lin_untidy <- rownames_to_column(wide_lin_tidy, var = "chem_family")
  lin_sig_chems <- append("Total Significant Chemicals", sum_sig_chems)
  wide_lin_sum <- rbind(wide_lin_untidy, lin_sig_chems)
  
  #reorder the chemicals and immune columns
  wide_lin_order <- wide_lin_sum %>%
    select(chem_family,
           `Lymphocytes (%)`,
           `Neutrophils (%)`,
           `Monocytes (%)`,
           `Basophils (%)`,
           `Eosinophils (%)`,
           `White Blood Cells (1000 cells/uL)`,
           `Red Blood Cells (million cells/uL)`,
           `Mean Corpuscular Volume (fL)`) %>%
    arrange(factor(chem_family,
           levels = chem_family_levels))
  
  
  #get counts of chems per family
  chem_counts <- conversion %>%
    group_by(chem_family) %>%
    summarise(count = n()) %>%
    ungroup()
  chem_counts <- rbind(chem_counts, c("Total Significant Chemicals", nrow(conversion)))
  chem_counts$count <- as.integer(chem_counts$count)
  
  #merge the chemical counts into the significant counts dataset
  chem_total <- left_join(wide_lin_order, chem_counts, by = "chem_family")
  
  #sort the chemical families by the preset chemical family levels
  chem_total$chem_family <- factor(chem_total$chem_family
                                   , levels = chem_family_levels)
  
  #make a dataset of percent significant per chemical family
  chem_percent <- chem_total %>%
    mutate(lym_pct = (as.numeric(`Lymphocytes (%)`) / count)*100,
           neu_pct = (as.numeric(`Neutrophils (%)`) / count)*100,
           mon_pct = (as.numeric(`Monocytes (%)`) / count)*100,
           bas_pct = (as.numeric(`Basophils (%)`) / count)*100,
           eos_pct = (as.numeric(`Eosinophils (%)`) / count)*100,
           wbc_pct = (as.numeric(`White Blood Cells (1000 cells/uL)`) / count)*100,
           rbc_pct = (as.numeric(`Red Blood Cells (million cells/uL)`) / count)*100,
           mcv_pct = (as.numeric(`Mean Corpuscular Volume (fL)`) / count)*100) %>%
    column_to_rownames(var = "chem_family") %>%
    select(lym_pct,
           neu_pct,
           mon_pct,
           bas_pct,
           eos_pct,
           wbc_pct,
           rbc_pct,
           mcv_pct)
  
  immune_measures <- colnames(chem_total[2:9])
  
  #edit the rownames to add the total number of chemicals per family
  chem_family_names <- paste0(rownames(chem_percent), " ", "(", chem_total$count, ")")
  
  #get the counts of chemicals to overlay the numbers in the boxes
  chem_count_matrix <- chem_total %>%
    column_to_rownames(var = "chem_family") %>%
    select(-count)
  
  #############################################################################################################
  ############################################# Make The Linear Plot ##########################################
  #############################################################################################################
  
  setwd(paste0(current_directory, "/Volcano Plots"))
  
  png("significant_lin_chems_per_fam_smk.png", units = "in", width = 14, height = 9, res = 600)
  pheatmap(mat = chem_percent,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix,
           color = colorRampPalette(c("white", "yellow", "orange", "red3"))(20))
  dev.off()
  
  pdf("significant_lin_chems_per_fam_smk.pdf", width = 14, height = 9)
  pheatmap(mat = chem_percent,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix,
           color = colorRampPalette(c("white", "yellow", "orange", "red3"))(20))
  dev.off()
  
  #############################################################################################################
  ####################################### SPLIT THE ESTIMATES BY DIRECTION ####################################
  #############################################################################################################
  
  #identify positive and negative effect estimates and split them into datasets
  effect_direction <- lin_model_lab %>%
    mutate(direction = str_detect(estimate, "-"))
  pos_effects <- effect_direction %>%
    filter(direction == "FALSE")
  neg_effects <- effect_direction %>%
    filter(direction == "TRUE")
  
  #count the number of significant chemicals per family and immune measures
  sig_lin_chems_pos <- pos_effects %>%
    group_by(immune_labels,
             chem_family) %>%
    summarise(num_chems = n()) %>%
    ungroup()
  sig_lin_chems_neg <- neg_effects %>%
    group_by(immune_labels,
             chem_family) %>%
    summarise(num_chems = n()) %>%
    ungroup()
  
  #create a wide dataset and replace the NAs with zeros, add BFRs
  wide_lin_sig_pos <- spread(sig_lin_chems_pos,
                             immune_labels,
                             num_chems) %>%
    replace(is.na(.), 0)
  wide_lin_sig_pos <- rbind(wide_lin_sig_pos, c("Brominated Flame Retardants (BFR)", 0, 0, 0, 0, 0, 0, 0, 0))
  wide_lin_sig_pos <- rbind(wide_lin_sig_pos, c("Phosphate Flame Retardants (PFR)", 0, 0, 0, 0, 0, 0, 0, 0))
  wide_lin_sig_neg <- spread(sig_lin_chems_neg,
                             immune_labels,
                             num_chems) %>%
    replace(is.na(.), 0)
  wide_lin_sig_neg <- rbind(wide_lin_sig_neg, c("Brominated Flame Retardants (BFR)", 0, 0, 0, 0, 0, 0, 0, 0))
  # wide_lin_sig_neg <- rbind(wide_lin_sig_neg, c("Furans", 0, 0, 0, 0, 0, 0, 0, 0))
  
  #turn chem_family row into row names
  wide_lin_tidy_pos <- wide_lin_sig_pos %>% column_to_rownames("chem_family")
  wide_lin_tidy_neg <- wide_lin_sig_neg %>% column_to_rownames("chem_family")
  
  #calculate total significant chemicals per immune measure
  sum_sig_chems_pos <- colSums(sapply(wide_lin_tidy_pos, as.numeric))
  sum_sig_chems_neg <- colSums(sapply(wide_lin_tidy_neg, as.numeric))
  
  #add a row of number of significant chemicals per immune measure
  wide_lin_untidy_pos <- rownames_to_column(wide_lin_tidy_pos, var = "chem_family")
  wide_lin_untidy_neg <- rownames_to_column(wide_lin_tidy_neg, var = "chem_family")
  lin_sig_chems_pos <- append("Total Significant Chemicals", sum_sig_chems_pos)
  lin_sig_chems_neg <- append("Total Significant Chemicals", sum_sig_chems_neg)
  wide_lin_sum_pos <- rbind(wide_lin_untidy_pos, lin_sig_chems_pos)
  wide_lin_sum_neg <- rbind(wide_lin_untidy_neg, lin_sig_chems_neg)
  
  #reorder the chemicals and immune columns
  wide_lin_order_pos <- wide_lin_sum_pos %>%
    select(chem_family,
           `Lymphocytes (%)`,
           `Neutrophils (%)`,
           `Monocytes (%)`,
           `Basophils (%)`,
           `Eosinophils (%)`,
           `White Blood Cells (1000 cells/uL)`,
           `Red Blood Cells (million cells/uL)`,
           `Mean Corpuscular Volume (fL)`) %>%
    arrange(factor(chem_family,
                   levels = chem_family_levels))
  wide_lin_order_neg <- wide_lin_sum_neg %>%
    select(chem_family,
           `Lymphocytes (%)`,
           `Neutrophils (%)`,
           `Monocytes (%)`,
           `Basophils (%)`,
           `Eosinophils (%)`,
           `White Blood Cells (1000 cells/uL)`,
           `Red Blood Cells (million cells/uL)`,
           `Mean Corpuscular Volume (fL)`) %>%
    arrange(factor(chem_family,
                   levels = chem_family_levels))
  
  
  #get counts of chems per family
  chem_counts <- conversion %>%
    group_by(chem_family) %>%
    summarise(count = n()) %>%
    ungroup()
  chem_counts <- rbind(chem_counts, c("Total Significant Chemicals", nrow(conversion)))
  chem_counts$count <- as.integer(chem_counts$count)
  
  #merge the chemical counts into the significant counts dataset
  chem_total_pos <- left_join(wide_lin_order_pos, chem_counts, by = "chem_family")
  chem_total_neg <- left_join(wide_lin_order_neg, chem_counts, by = "chem_family")
  
  #sort the chemical families by the preset chemical family levels
  chem_total_pos$chem_family <- factor(chem_total_pos$chem_family
                                   , levels = chem_family_levels)
  chem_total_neg$chem_family <- factor(chem_total_neg$chem_family
                                   , levels = chem_family_levels)
  
  #make a dataset of percent significant per chemical family
  chem_percent_pos <- chem_total_pos %>%
    mutate(lym_pct = (as.numeric(`Lymphocytes (%)`) / count)*100,
           neu_pct = (as.numeric(`Neutrophils (%)`) / count)*100,
           mon_pct = (as.numeric(`Monocytes (%)`) / count)*100,
           bas_pct = (as.numeric(`Basophils (%)`) / count)*100,
           eos_pct = (as.numeric(`Eosinophils (%)`) / count)*100,
           wbc_pct = (as.numeric(`White Blood Cells (1000 cells/uL)`) / count)*100,
           rbc_pct = (as.numeric(`Red Blood Cells (million cells/uL)`) / count)*100,
           mcv_pct = (as.numeric(`Mean Corpuscular Volume (fL)`) / count)*100) %>%
    column_to_rownames(var = "chem_family") %>%
    select(lym_pct,
           neu_pct,
           mon_pct,
           bas_pct,
           eos_pct,
           wbc_pct,
           rbc_pct,
           mcv_pct)
  chem_percent_neg <- chem_total_neg %>%
    mutate(lym_pct = (as.numeric(`Lymphocytes (%)`) / count)*100,
           neu_pct = (as.numeric(`Neutrophils (%)`) / count)*100,
           mon_pct = (as.numeric(`Monocytes (%)`) / count)*100,
           bas_pct = (as.numeric(`Basophils (%)`) / count)*100,
           eos_pct = (as.numeric(`Eosinophils (%)`) / count)*100,
           wbc_pct = (as.numeric(`White Blood Cells (1000 cells/uL)`) / count)*100,
           rbc_pct = (as.numeric(`Red Blood Cells (million cells/uL)`) / count)*100,
           mcv_pct = (as.numeric(`Mean Corpuscular Volume (fL)`) / count)*100) %>%
    column_to_rownames(var = "chem_family") %>%
    select(lym_pct,
           neu_pct,
           mon_pct,
           bas_pct,
           eos_pct,
           wbc_pct,
           rbc_pct,
           mcv_pct)
  
  immune_measures <- colnames(chem_total_pos[2:9])
  
  #edit the rownames to add the total number of chemicals per family
  chem_family_names_pos <- paste0(rownames(chem_percent_pos), " ", "(", chem_total_pos$count, ")")
  chem_family_names_neg <- paste0(rownames(chem_percent_neg), " ", "(", chem_total_neg$count, ")")

  #get the counts of chemicals to overlay the numbers in the boxes
  chem_count_matrix_pos <- chem_total_pos %>%
    column_to_rownames(var = "chem_family") %>%
    select(-count)
  chem_count_matrix_neg <- chem_total_neg %>%
    column_to_rownames(var = "chem_family") %>%
    select(-count)
  
  #############################################################################################################
  ############################################# Make The Linear Plots #########################################
  #############################################################################################################
  
  setwd(paste0(current_directory, "/Volcano Plots"))
  
  png("significant_lin_chems_by_pos_direction_smk.png", units = "in", width = 14, height = 9, res = 600)
  pheatmap(mat = chem_percent_pos,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names_pos,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix_pos,
           color = colorRampPalette(c("white", "yellow", "orange", "red3"))(20))
  dev.off()
  pdf("significant_lin_chems_by_pos_direction_smk.pdf", width = 14, height = 9)
  pheatmap(mat = chem_percent_pos,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names_pos,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix_pos,
           color = colorRampPalette(c("white", "yellow", "orange", "red3"))(20))
  dev.off()
  
  
  # Set up color bar for negative estimates so that blue3 = 100%, not just highest %
  #https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
  paletteLength <- 100
  myColor <- colorRampPalette(c("white", "deepskyblue", "blue3"))(paletteLength)
  breaksList = seq(0, 100, by = 1)
  
  png("significant_lin_chems_by_neg_direction_smk.png", units = "in", width = 14, height = 9, res = 600)
  pheatmap(mat = chem_percent_neg,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names_neg,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix_neg,
           color=myColor,
           breaks=breaksList)
  dev.off()
  pdf("significant_lin_chems_by_neg_direction_smk.pdf", width = 14, height = 9)
  pheatmap(mat = chem_percent_neg,
           cluster_rows = FALSE, cluster_cols = FALSE,
           labels_col = immune_measures,
           labels_row = chem_family_names_neg,
           annotation_names_row = TRUE,
           annotation_names_col = FALSE,
           angle_col = 315,
           fontsize = 15,
           fontsize_number = 20,
           legend = TRUE,
           display_numbers = chem_count_matrix_neg,
           color=myColor,
           breaks=breaksList)
  dev.off()
  
  setwd(current_directory)
  
  
  #############################################################################################################
  ########################################### Make A Comparison Barplot #######################################
  #############################################################################################################
  
  # #create datasets of significant counts per chemical family in the positive and negative directions
  # count_neg <- chem_total_neg %>%
  #   remove_rownames() %>%
  #   column_to_rownames(var = "chem_family") %>%
  #   filter(row_number()==1) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   rename("Negative" = `Total Significant Chemicals`)
  # 
  # count_pos <- chem_total_pos %>%
  #   remove_rownames() %>%
  #   column_to_rownames(var = "chem_family") %>%
  #   filter(row_number()==1) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   rename("Positive" = `Total Significant Chemicals`)
  # 
  # #check that the rownames are the same order
  # identical(rownames(count_neg), rownames(count_pos))
  # 
  # #merge the two datasets
  # sig_count_merged <- cbind(count_pos, count_neg) %>%
  #   rownames_to_column(var = "immune_measure")
  # 
  # #make a long dataset for plotting
  # long_sig_count <- pivot_longer(sig_count_merged,
  #                                cols = c("Positive", "Negative"),
  #                                values_to = "counts",
  #                                names_to = "direction") %>%
  #   filter(!immune_measure == "count") %>%
  #   mutate(counts = as.integer(counts)) %>%
  #   mutate(immune_measure = factor(immune_measure))
  # 
  # #set the order for the barplot categories
  # # effect_order <- c("Positive", "Negative")
  # immune_levels <- rev(c("Lymphocytes (%)",
  #                    "Neutrophils (%)",
  #                    "Monocytes (%)",
  #                    "Basophils (%)",
  #                    "Eosinophils (%)",
  #                    "White Blood Cells (1000 cells/uL)",
  #                    "Red Blood Cells (million cells/uL)",
  #                    "Mean Corpuscular Volume (fL)"))
  # long_sig_count$immune_measure <- factor(long_sig_count$immune_measure, levels = immune_levels)
  # # long_sig_count$direction <- factor(long_sig_count$direction, levels = effect_order)
  # 
  # 
  # 
  # #plot the barplot
  # sig_counts_plot <- 
  # ggplot(data = long_sig_count,
  #        aes(x = immune_measure,
  #            y = counts,
  #            fill = direction))+
  #   geom_bar(width = 0.8,
  #            position = "dodge",
  #            stat = "identity",
  #            color = "black")+
  #   labs(x = "",
  #        y = "Number of Chemicals")+
  #   scale_fill_manual(values=c("#003CDC", "#FC9C00"),
  #                     name = "Effect Direction")+
  #   scale_y_continuous(breaks = seq(0, 100, by = 10))+
  #   theme(axis.text.y = element_text(color = "black",
  #                                    size = 15),
  #         axis.text.x = element_text(size = 15,
  #                                    color = "black",
  #                                    vjust = 0),
  #         axis.title = element_text(size = 20,
  #                                   color = "black")) +
  #   theme(panel.background = element_rect(fill = "white", #this is the background
  #                                         colour = "black", #this is the border
  #                                         size = 0.1, linetype = "solid"))+
  #   theme(panel.grid = element_blank())+
  #   theme(legend.title = element_text(size=20),
  #         legend.text = element_text(size=15))+
  #   coord_flip()
  # 
  # 
  # 
  # setwd(paste0(current_directory, "/Volcano Plots"))
  # 
  # #save the barplot
  # print("bar_plot_sig_counts_direction_smk.pdf")
  # ggsave(filename = "bar_plot_sig_counts_direction_smk.pdf"
  #        , plot = sig_counts_plot
  #        , width = 14
  #        , height = 9)
  # 
  # # Save the plot as a png for presentation
  # print("bar_plot_sig_counts_direction_smk.png")
  # ggsave(filename = "bar_plot_sig_counts_direction_smk.png"
  #        , plot = sig_counts_plot
  #        , units = "in"
  #        , width = 14
  #        , height = 9
  #        , dpi = 300)
  # 
  # setwd(current_directory)
  
  
  #############################################################################################################
  ##################################################### Stats #################################################
  #############################################################################################################
  
  test <- model_stats_smk %>% filter(FDR <0.05) %>% filter(term == "chem_log_measurement") %>% group_by(chem_family,chemical_codename) %>% tally() %>% ungroup()
  View(test)

}
