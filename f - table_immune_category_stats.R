#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################################  CALCULATE SUMMARY STATISTICS  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function calculates summary statistics for chemicals and cell-types
#          
# Inputs:   nhanes_subset      - dataframe of complete demographics, cells, and chemicals
#
# Outputs:  summary statistics - csv file containing all summary statistics for immune measure categories

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_immune_category_stats <- function(nhanes_subset)
{
  library(tidyverse)
  library(gtsummary)
  library(flextable)
  
  #############################################################################################################
  ############################################## Create Categories ############################################
  #############################################################################################################
  
  # Categorize immune measures into High, Normal, Low
  # Categories based on reference cited in manuscript to create Table 2
  nhanes_immune_categ <- nhanes_subset %>%
    dplyr::select(SEQN,
                  LBDLYMNO,
                  LBDNENO,
                  LBDMONO,
                  LBDBANO,
                  LBDEONO,
                  LBXWBCSI,
                  LBXRBCSI,
                  LBXMCVSI) %>%
    mutate(lym_categ = case_when(LBDLYMNO < 1 ~ "Low",
                                 LBDLYMNO >= 1 & LBDLYMNO <= 4 ~ "Normal",
                                 LBDLYMNO > 4 ~ "High"),
           neu_categ = case_when(LBDNENO < 2.8 ~ "Low",
                                 LBDNENO >= 2.8 & LBDNENO <= 5.6 ~ "Normal",
                                 LBDNENO > 5.6 ~ "High"),
           mon_categ = case_when(LBDMONO < 0.1 ~ "Low",
                                 LBDMONO >= 0.1 & LBDMONO <= 0.8 ~ "Normal",
                                 LBDMONO > 0.8 ~ "High"),
           bas_categ = case_when(LBDBANO < 0.025 ~ "Low",
                                 LBDBANO >= 0.025 & LBDBANO <= 0.1 ~ "Normal",
                                 LBDBANO > 0.1 ~ "High"),
           eos_categ = case_when(LBDEONO < 0.05 ~ "Low",
                                 LBDEONO >= 0.05 & LBDEONO <= 0.4 ~ "Normal",
                                 LBDEONO > 0.4 ~ "High"),
           wbc_categ = case_when(LBXWBCSI < 4.5 ~ "Low",
                                 LBXWBCSI >= 4.5 & LBXWBCSI <= 11 ~ "Normal",
                                 LBXWBCSI > 11 ~ "High"),
           rbc_categ = case_when(LBXRBCSI < 4 ~ "Low",
                                 LBXRBCSI >= 4 & LBXRBCSI <= 6.2 ~ "Normal",
                                 LBXRBCSI > 6.2 ~ "High"),
           mcv_categ = case_when(LBXMCVSI < 82 ~ "Low",
                                 LBXMCVSI >= 82 & LBXMCVSI <= 93 ~ "Normal",
                                 LBXMCVSI > 93 ~ "High")
           ) %>%
    dplyr::select(-LBDLYMNO,
                  -LBDNENO,
                  -LBDMONO,
                  -LBDBANO,
                  -LBDEONO,
                  -LBXWBCSI,
                  -LBXRBCSI,
                  -LBXMCVSI)
  
  #############################################################################################################
  ################################################# Set Levels ################################################
  #############################################################################################################
  
  # Set the order of levels for each category
  immune_levels <- nhanes_immune_categ %>%
    mutate(lym_categ = (factor(lym_categ, levels = c("Low", "Normal", "High"))),
           neu_categ = (factor(neu_categ, levels = c("Low", "Normal", "High"))),
           mon_categ = (factor(mon_categ, levels = c("Low", "Normal", "High"))),
           bas_categ = (factor(bas_categ, levels = c("Low", "Normal", "High"))),
           eos_categ = (factor(eos_categ, levels = c("Low", "Normal", "High"))),
           wbc_categ = (factor(wbc_categ, levels = c("Low", "Normal", "High"))),
           rbc_categ = (factor(rbc_categ, levels = c("Low", "Normal", "High"))),
           mcv_categ = (factor(mcv_categ, levels = c("Low", "Normal", "High")))
           )
  str(immune_levels)
  
  #############################################################################################################
  ################################################ Create Table ###############################################
  #############################################################################################################
  
  #set directory
  setwd(paste0(current_directory, "/Tables - Table 1"))

  immune_levels %>%
    dplyr::select(-SEQN) %>%
    tbl_summary(statistic = list(all_categorical() ~ "{n} ({p}%)"),
                digits = list(all_categorical() ~ c(0, 1)),
                missing_text = "Missing (n)",
                label = c(lym_categ ~ "Lymphocytes",
                          neu_categ ~ "Neutrophils",
                          mon_categ ~ "Monocytes",
                          bas_categ ~ "Basophils",
                          eos_categ ~ "Eosinophils",
                          wbc_categ ~ "White Blood Cells",
                          rbc_categ ~ "Red Blood Cells",
                          mcv_categ ~ "Mean Corpuscular Volume")) %>%
    modify_header(label ~ "**Variable**") %>%
    bold_labels() %>%
    as_flex_table() %>%
    save_as_docx(path = "categorical_immune_stats_new.docx") #Export the table to Word using flextable package


  setwd(current_directory)
    
  #############################################################################################################
  ############################################ Count Abnormal Levels ##########################################
  #############################################################################################################
  
  # Calculate how many participants have an "abnormal" immune measurement
  long_immune_levels <- immune_levels %>%
    pivot_longer(cols = lym_categ:mcv_categ,
                 names_to = "immune_type",
                 values_to = "category_level") %>%
    mutate(category_level = as.character(category_level)) %>%
    mutate(abnormal = case_when(category_level == "Low" | category_level == "High" ~ "Abnormal",
                                TRUE ~ as.character(category_level))) %>%
    dplyr::select(-category_level) %>%
    pivot_wider(names_from = "immune_type",
                values_from = "abnormal") %>%
      nest(-SEQN) %>%
      transmute(SEQN, abnor_count = map(data, ~ sum(unlist(.x) == "Abnormal", na.rm = TRUE))) %>%
    mutate(abnor_count = as.integer(abnor_count))
  
    
  View(table(long_immune_levels$abnor_count))
   
  #############################################################################################################
}