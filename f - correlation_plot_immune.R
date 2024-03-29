#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################  MAKE PLOT OF CORRELATIONS BETWEEN IMMUNE MEASURES  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates a correlation plot of the immune measures
#          
# Inputs: nhanes_subset     - dataframe containing complete demographic and immune measures data for each
#                             participant
#
# Outputs: Plot showing correlations between immune measures ("immune_correlation.pdf")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

correlation_plot_immune <- function(nhanes_subset)
{
  library(corrplot)
  library(dichromat)
  library(tidyverse)
  
  #identify the immune measures
  celltype_codename <- c("LBDLYMNO", #lymphocytes
                         "LBDNENO",  #neutrophils
                         "LBDMONO",  #monocytes
                         "LBDBANO",  #basophils
                         "LBDEONO",  #eosinophils
                         "LBXWBCSI", #WBC count
                         "LBXRBCSI", #RBC count
                         "LBXMCVSI"  #MCV
                         )
  
  #select the immune measures from the dataset
  nhanes_immune <- nhanes_subset_dataset %>%
    dplyr::select(all_of(celltype_codename))
  
  
  
  #calculate the correlations
  immune_cor <- cor(nhanes_immune)
  
  setwd(paste0(current_directory, "/Correlation Plots - Demog, Cells, Chemicals"))
  write.csv(immune_cor, "immune_correlations_matrix_new.csv")
  
  
  # Numbers for manuscript
  abs_immune_cor <- abs(immune_cor)
  table(abs_immune_cor)
  
  
  # Create .png image
  png("immune_correlation_smk_new.png", units = "in", res = 300, width = 10, height = 9)
  M <- (immune_cor)
  colnames(M) <- c("Lymphocytes (1000 cells/uL)",
                   "Neutrophils (1000 cells/uL)",
                   "Monocytes (1000 cells/uL)",
                   "Basophils (1000 cells/uL)",
                   "Eosinophils (1000 cells/uL)",
                   "WBC (1000 cells/uL)",
                   "RBC (million cells/uL)",
                   "Mean Corpuscular Volume (fL)")
  rownames(M) <- c("Lymphocytes (1000 cells/uL)",
                   "Neutrophils (1000 cells/uL)",
                   "Monocytes (1000 cells/uL)",
                   "Basophils (1000 cells/uL)",
                   "Eosinophils (1000 cells/uL)",
                   "WBC (1000 cells/uL)",
                   "RBC (million cells/uL)",
                   "Mean Corpuscular Volume (fL)")
  corrplot(M,
           method = "ellipse",
           type = "lower",
           col = colorRampPalette(c("blue", "white", "red"))(20),
           tl.col = "black")
  dev.off()
  
  # Create .pdf image
  M <- immune_cor
  colnames(M) <- c("Lymphocytes (1000 cells/uL)",
                   "Neutrophils (1000 cells/uL)",
                   "Monocytes (1000 cells/uL)",
                   "Basophils (1000 cells/uL)",
                   "Eosinophils (1000 cells/uL)",
                   "WBC (1000 cells/uL)",
                   "RBC (million cells/uL)",
                   "Mean Corpuscular Volume (fL)")
  rownames(M) <- c("Lymphocytes (1000 cells/uL)",
                   "Neutrophils (1000 cells/uL)",
                   "Monocytes (1000 cells/uL)",
                   "Basophils (1000 cells/uL)",
                   "Eosinophils (1000 cells/uL)",
                   "WBC (1000 cells/uL)",
                   "RBC (million cells/uL)",
                   "Mean Corpuscular Volume (fL)")
  pdf("immune_correlation_new.pdf")
  corrplot(M,
           method = "ellipse",
           type = "lower",
           col = colorRampPalette(c("blue", "white", "red"))(20),
           tl.col = "black")
  dev.off()
  
  
  
  setwd(current_directory)
}