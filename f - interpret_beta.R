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
  chem_pct <- c("LBXLYPCT",
                "LBXNEPCT",
                "LBXMOPCT",
                "LBXBAPCT",
                "LBXEOPCT")
  
  # Keep only the chemical betas
  model_stats_chems <- model_stats %>%
    # filter(term == "chem_log_measurement") %>%
    mutate(interpret = case_when(celltype_codename %in% chem_pct ~ estimate,
                                 celltype_codename == "LBXWBCSI" ~ estimate * 1000,
                                 celltype_codename == "LBXRBCSI" ~ estimate * 1000000,
                                 celltype_codename == "LBXMCVSI" ~ estimate)) %>%
    mutate(lower_ci = case_when(celltype_codename %in% chem_pct ~ lower.CI,
                                celltype_codename == "LBXWBCSI" ~ lower.CI * 1000,
                                celltype_codename == "LBXRBCSI" ~ lower.CI * 1000000,
                                celltype_codename == "LBXMCVSI" ~ lower.CI)) %>%
    mutate(upper_ci = case_when(celltype_codename %in% chem_pct ~ upper.CI,
                                celltype_codename == "LBXWBCSI" ~ upper.CI * 1000,
                                celltype_codename == "LBXRBCSI" ~ upper.CI * 1000000,
                                celltype_codename == "LBXMCVSI" ~ upper.CI)) %>%
    mutate(units = case_when(celltype_codename %in% chem_pct ~ "%",
                             celltype_codename == "LBXWBCSI" ~ "cells per uL",
                             celltype_codename == "LBXRBCSI" ~ "cells per uL",
                             celltype_codename == "LBXMCVSI" ~ "fL")) %>%
    mutate(covariate = case_when(term == "RIDRETH11" ~ "race_mexican_american",
                                 term == "RIDRETH12" ~ "race_other_hispanic",
                                 term == "RIDRETH14" ~ "race_non_hispanic_black",
                                 term == "RIDRETH15" ~ "race_other",
                                 term == "RIDAGEYR"  ~ "age",
                                 term == "RIAGENDR2" ~ "sex_female",
                                 term == "INDFMPIR"  ~ "poverty_income_ratio",
                                 term == "BMXWAIST"  ~ "waist_circumference",
                                 term == "BMXWAIST"  ~ "waist_circumference",
                                 term == "SDDSRVYR2" ~ "cycle_2",
                                 term == "SDDSRVYR3" ~ "cycle_3",
                                 term == "SDDSRVYR4" ~ "cycle_4",
                                 term == "SDDSRVYR5" ~ "cycle_5",
                                 term == "SDDSRVYR6" ~ "cycle_6",
                                 term == "SDDSRVYR7" ~ "cycle_7",
                                 term == "SDDSRVYR8" ~ "cycle_8",
                                 term == "SDDSRVYR9" ~ "cycle_9",
                                 term == "SDDSRVYR10" ~ "cycle_10",
                                 term == "URXUCR"    ~ "creatinine",
                                 term == "(Intercept)" ~ "intercept",
                                 term == "SMOKING" ~ "cotinine",
                                 term == "chem_log_measurement" ~ "chemical")) %>%
    dplyr::select(-estimate,
                  -lower.CI,
                  -upper.CI,
                  -term) %>%
    relocate(interpret, .after = celltype_codename) %>%
    relocate(covariate, .before = interpret) %>%
    relocate(lower_ci, .after = interpret) %>%
    relocate(upper_ci, .after = lower_ci)
  
  #############################################################################################################
  ######################################### Save Interpretted Dataset #########################################
  #############################################################################################################
  
  setwd(paste0(current_directory, "/Regression Results"))
  
  write.csv(model_stats_chems,
            file = "model_stats_chems_interpretted.csv",
            row.names = F)
  
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