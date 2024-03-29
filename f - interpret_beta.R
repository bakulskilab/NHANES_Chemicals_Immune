#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################  INTERPRET CHEMICAL BETAS  #########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates table of values to interpret the chemical beta coefficients
#          
# Inputs:   model_stats - tidy output of linear regression stats adjusted for demographics
#
# Outputs:  supplemental table

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

interpret_beta <- function(model_stats)
{
  library(tidyverse)
  
  setwd(current_directory)

  #TEMPORARY  
  # model_stats <- model_stats_smk

  #############################################################################################################
  ###################################### Clean Up and Interpret Columns #######################################
  #############################################################################################################

  # Define chemicals that are expressed as percentages or don't change
  chem_pct <- c("LBDLYMNO", #lymphocyte counts
                "LBDMONO",  #monocyte counts
                "LBDNENO",  #neutrophil counts
                "LBDEONO",  #eosinophil counts
                "LBDBANO"   #basophil counts
                )
  
  # Keep only the chemical betas
  model_stats_chems <- model_stats %>%
    filter(term == "chem_log_measurement") %>%
    mutate(interpret = case_when(celltype_codename %in% chem_pct ~ estimate * 1000,
                                 celltype_codename == "LBXWBCSI" ~ estimate * 1000,
                                 celltype_codename == "LBXRBCSI" ~ estimate * 1000000,
                                 celltype_codename == "LBXMCVSI" ~ estimate)) %>%
    mutate(lower_ci = case_when(celltype_codename %in% chem_pct ~ lower.CI * 1000,
                                celltype_codename == "LBXWBCSI" ~ lower.CI * 1000,
                                celltype_codename == "LBXRBCSI" ~ lower.CI * 1000000,
                                celltype_codename == "LBXMCVSI" ~ lower.CI)) %>%
    mutate(upper_ci = case_when(celltype_codename %in% chem_pct ~ upper.CI * 1000,
                                celltype_codename == "LBXWBCSI" ~ upper.CI * 1000,
                                celltype_codename == "LBXRBCSI" ~ upper.CI * 1000000,
                                celltype_codename == "LBXMCVSI" ~ upper.CI)) %>%
    mutate(std_error = case_when(celltype_codename %in% chem_pct ~ std.error * 1000,
                                celltype_codename == "LBXWBCSI" ~ std.error * 1000,
                                celltype_codename == "LBXRBCSI" ~ std.error * 1000000,
                                celltype_codename == "LBXMCVSI" ~ std.error)) %>%
    mutate(units = case_when(celltype_codename %in% chem_pct ~ "cells per uL",
                             celltype_codename == "LBXWBCSI" ~ "cells per uL",
                             celltype_codename == "LBXRBCSI" ~ "cells per uL",
                             celltype_codename == "LBXMCVSI" ~ "fL")) %>%
    mutate(chemical_name_only = gsub("\\s\\(([^()]+)\\)$", "", chemical_name)) %>%
    dplyr::select(immune_measure,
                  chemical_name_only,
                  interpret,
                  std_error,
                  lower.CI,
                  upper.CI,
                  statistic,
                  p.value,
                  FDR,
                  chem_family,
                  chemical_codename,
                  nobs)  %>%
    rename('Immune Measure' = immune_measure,
           'Chemical Name' = chemical_name_only,
           'Beta Coefficient' = interpret,
           'Standard Error' = std_error,
           'Lower 95% Confidence Interval' = lower.CI,
           'Upper 95% Confidence Interval' = upper.CI,
           'Statistic' = statistic,
           'p-value' = p.value,
           'FDR' = FDR,
           'Number of Participants' = nobs,
           'Chemical Family' = chem_family,
           'Chemical Codename' = chemical_codename)
  
  #############################################################################################################
  ######################################### Save Interpretted Dataset #########################################
  #############################################################################################################
  
  setwd(paste0(current_directory, "/Regression Results"))
  
  write.csv(model_stats_chems,
            file = "model_stats_chems_interpretted_new.csv",
            row.names = F)
  
  print("Results saved in Regression Results folder")
  
  #############################################################################################################
  ###################################### Significant Immune per Chemical ######################################
  #############################################################################################################
  
  #number significant chemicals per immune measure
  model_stats_summary <- model_stats %>%
    filter(term == "chem_log_measurement") %>%
    filter(FDR < 0.05) %>%
    group_by(immune_measure) %>%
    summarise(count = length(FDR)) %>%
    ungroup()
  View(model_stats_summary)
  
  #number significant immune per chemical
  model_stats_chems <- model_stats %>%
    filter(term == "chem_log_measurement") %>%
    filter(FDR < 0.05) %>%
    group_by(chemical_codename) %>%
    summarise(count = length(FDR)) %>%
    ungroup()
  View(model_stats_chems)
  
  #############################################################################################################
  
  setwd(current_directory)
}