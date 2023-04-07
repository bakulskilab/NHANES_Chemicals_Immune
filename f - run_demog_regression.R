#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#########################################  Adjusted Linear Regression  ########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function runs the adjusted linear regressions for the demographics and immune measures
#          
# Inputs: nhanes_subset - dataframe containing complete demographic and immune data for each
#                         participant
#         conversion - dataframe of chemical names, codenames, and families
#         demog_dataset - dataframe that includes the MEC survey weights
#
# Outputs: model_stats_demog - dataframe of linear regression outputs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

run_demog_regression <- function(nhanes_subset,
                                  conversion,
                                  demog_dataset)
{
  library(tidyverse)
  library(survey)
  library(broom)
  library(gt)
  library(flextable)
  
  #TEMPORARY
  # nhanes_subset <- nhanes_subset_dataset
  # conversion <- use_these_chems
  # demog_dataset <- demographics_clean

  #############################################################################################################
  ###################################### Select Variables, Merge Datasets #####################################
  #############################################################################################################
  
  # Select the variables to use from nhanes subset and demographics dataset
  nhanes_vars <- nhanes_subset %>%
    dplyr::select(SEQN,
                  "LBXLYPCT", #lymphocytes
                  "LBXNEPCT", #neutrophils
                  "LBXMOPCT", #monocytes
                  "LBXBAPCT", #basophils
                  "LBXEOPCT", #eosinophils
                  "LBXWBCSI", #WBC count
                  "LBXRBCSI", #RBC count
                  "LBXMCVSI",  #MCV
                  RIAGENDR,
                  RIDRETH1,
                  RIDAGEYR,
                  INDFMPIR,
                  BMXWAIST,
                  SMOKING,
                  SDDSRVYR,
                  SDMVPSU,
                  SDMVSTRA)
  
  demog_vars <- demog_dataset %>%
    dplyr::select(SEQN,
                  WTMEC4YR,
                  WTMEC2YR)
  
  # Merge the two subsets of variables and keep only participants previously selected
  vars_weights <- left_join(nhanes_vars, demog_vars, by = "SEQN")
  
  
  #############################################################################################################
  ############################################### Clean Weights ###############################################
  #############################################################################################################
  
  # Make a subset of participants who were in cycles 1 and 2
  cycles1_2 <- c(1:2)
  vars_12 <- vars_weights %>%
    filter(SDDSRVYR %in% cycles1_2)
  
  # Make WTMEC2YR NA because CDC says to use the MEC4 weights for the first two cycles
  vars_12$WTMEC2YR <- NA
  
  # Check that all participants have weights
  summary(vars_12$WTMEC4YR)
  sum(is.na(vars_12$WTMEC4YR))
  sum(!is.na(vars_12$WTMEC4YR))
  
  
  # Make a subset of participants after cycle 2
  cycle_other <- c(3:10)
  vars_other <- vars_weights %>%
    filter(SDDSRVYR %in% cycle_other)
  
  # Make WTMEC4YR NA because CDC says to use the MEC2 weights for cycles after 2
  vars_other$WTMEC4YR <- NA
  
  # Check that all participants have weights
  summary(vars_other$WTMEC2YR)
  sum(is.na(vars_other$WTMEC2YR))
  sum(!is.na(vars_other$WTMEC2YR)) #37173 + 8355 = 45528
  # Everyone has a weight and neither set of weights has extra rows
  
  
  # Combine the two subsets back together
  vars_weights_merge <- bind_rows(vars_12, vars_other)
  dim(vars_weights_merge)
  
  rm(vars_weights)
  
  # Multiply MEC4 weights by 2/10 and MEC2 by 1/10
  # Merge the weights columns
  # Make factor variables
  vars_weights_clean <- vars_weights_merge %>%
    mutate(mec4_clean = WTMEC4YR * (2/10),
           mec2_clean = WTMEC2YR * (1/10)) %>%
    mutate(weights_adjusted = coalesce(mec4_clean, mec2_clean)) %>%
    dplyr::select(-WTMEC4YR,
                  -WTMEC2YR,
                  -mec4_clean,
                  -mec2_clean,
                  -SMOKING) %>%
    mutate(RIAGENDR = relevel(factor(RIAGENDR),
                              ref = 1),
           RIDRETH1 = relevel(factor(RIDRETH1),
                              ref = 3),
           SDDSRVYR = factor(SDDSRVYR))
  sum(is.na(vars_weights_clean$weights_adjusted))
  sum(vars_weights_clean$weights_adjusted)
  #188 million
  
  #############################################################################################################
  ######################################### Create Individual Datasets ########################################
  #############################################################################################################
  
  #Chirag suggests to remove the "lonely" PSUs - strata with only one PSU
  #this is because "a single-PSU stratum makes no contribution to the variance"
  # - https://r-survey.r-forge.r-project.org/survey/html/surveyoptions.html
  options(survey.lonely.psu = "remove")
  
  # Make datasets for each immune measure
  
  # Lymphocytes dataset (don't really need to pivot longer but it will make the code the same later)
  LR_data_ly <- vars_weights_clean %>%
    dplyr::select(-"LBXNEPCT", #neutrophils
                  -"LBXMOPCT", #monocytes
                  -"LBXBAPCT", #basophils
                  -"LBXEOPCT", #eosinophils
                  -"LBXWBCSI", #WBC count
                  -"LBXRBCSI", #RBC count
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXLYPCT",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # Neutrophils dataset
  LR_data_ne <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXMOPCT",
                  -"LBXBAPCT",
                  -"LBXEOPCT",
                  -"LBXWBCSI",
                  -"LBXRBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXNEPCT",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # Monocytes dataset
  LR_data_mo <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXBAPCT",
                  -"LBXEOPCT",
                  -"LBXWBCSI",
                  -"LBXRBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXMOPCT",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # Basophils dataset
  LR_data_ba <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXMOPCT",
                  -"LBXEOPCT",
                  -"LBXWBCSI",
                  -"LBXRBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXBAPCT",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # Eosinophils dataset
  LR_data_eo <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXMOPCT",
                  -"LBXBAPCT",
                  -"LBXWBCSI",
                  -"LBXRBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXEOPCT",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # WBC dataset
  LR_data_wbc <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXMOPCT",
                  -"LBXBAPCT",
                  -"LBXEOPCT",
                  -"LBXRBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXWBCSI",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # RBC dataset
  LR_data_rbc <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXMOPCT",
                  -"LBXBAPCT",
                  -"LBXEOPCT",
                  -"LBXWBCSI",
                  -"LBXMCVSI") %>%
    pivot_longer(cols = "LBXRBCSI",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  # MCV dataset
  LR_data_mcv <- vars_weights_clean %>%
    dplyr::select(-"LBXLYPCT", 
                  -"LBXNEPCT",
                  -"LBXMOPCT",
                  -"LBXBAPCT",
                  -"LBXEOPCT",
                  -"LBXWBCSI",
                  -"LBXRBCSI") %>%
    pivot_longer(cols = "LBXMCVSI",
                 names_to = "celltype_codename",
                 values_to = "cell_measurement") %>%
    dplyr::select(-SEQN) %>%
    na.omit(.)
  
  #############################################################################################################
  ############################################# Make Survey Objects ###########################################
  #############################################################################################################
  
  # Pull everything together to get the survey adjustment
  svy_nhanes_ly <- svydesign(strata = ~SDMVSTRA
                          , id = ~SDMVPSU
                          , weights = ~weights_adjusted
                          , data = LR_data_ly
                          , nest = TRUE)
  
  svy_nhanes_ne <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_ne
                             , nest = TRUE)
  
  svy_nhanes_mo <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_mo
                             , nest = TRUE)
  
  svy_nhanes_ba <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_ba
                             , nest = TRUE)
  
  svy_nhanes_eo <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_eo
                             , nest = TRUE)
  
  svy_nhanes_wb <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_wbc
                             , nest = TRUE)
  
  svy_nhanes_rb <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_rbc
                             , nest = TRUE)
  
  svy_nhanes_mc <- svydesign(strata = ~SDMVSTRA
                             , id = ~SDMVPSU
                             , weights = ~weights_adjusted
                             , data = LR_data_mcv
                             , nest = TRUE)
  
  #############################################################################################################
  ############################################# Linear Regressions ############################################
  #############################################################################################################
  
  # Run regressions
  # Add confidence intervals + FDRs
  # Keep only age term
  # Add a row label for immune measure
  z_score <- 1.96
  
  df_regression_ly <- LR_data_ly %>%
    do(svyglm(cell_measurement ~ 
             RIDRETH1+
             RIDAGEYR+
             RIAGENDR+
             INDFMPIR+
             BMXWAIST+
             SDDSRVYR,
           design = svy_nhanes_ly,
           data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Lymphocytes")
  
  df_regression_ne <- LR_data_ne %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_ne,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Neutrophils")
  
  df_regression_mo <- LR_data_mo %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_mo,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Monocytes")
  
  df_regression_ba <- LR_data_ba %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_ba,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Basophils")
  
  df_regression_eo <- LR_data_eo %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_eo,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Eosinophils")
  
  df_regression_wb <- LR_data_wbc %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_wb,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "White Blood Cells")
  
  df_regression_rb <- LR_data_rbc %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_rb,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Red Blood Cells")
  
  df_regression_mc <- LR_data_mcv %>%
    do(svyglm(cell_measurement ~ 
                RIDRETH1+
                RIDAGEYR+
                RIAGENDR+
                INDFMPIR+
                BMXWAIST+
                SDDSRVYR,
              design = svy_nhanes_mc,
              data = .) %>%
         tidy(.)) %>%
    mutate(lower.CI = estimate - (z_score*std.error),
           upper.CI = estimate + (z_score*std.error)) %>%
    # mutate(FDR = p.adjust(p.value, method = "fdr")) %>%
    filter(term == "RIDAGEYR") %>%
    mutate(immune_measure = "Mean Corpuscular Volume")

  #############################################################################################################
  ###################################### ADD CONFIDENCE INTERVALS AND FDR #####################################
  #############################################################################################################

  # Put the results together
  model_stats_demog <- bind_rows(list(df_regression_ly,
                                   df_regression_ne,
                                   df_regression_mo,
                                   df_regression_ba,
                                   df_regression_eo,
                                   df_regression_wb,
                                   df_regression_rb,
                                   df_regression_mc))
  # Multiple comparisons of age to different immune measures so adding the FDR
  model_stats_demog$FDR <- p.adjust(model_stats_demog$p.value, method='fdr')
  
  # Interpret the results for a 10-year increase in age
   chem_pct <- c("Lymphocytes",
                 "Neutrophils",
                 "Monocytes",
                 "Basophils",
                 "Eosinophils")
  
  model_interpret <- model_stats_demog %>%
    mutate(interpret = case_when(immune_measure %in% chem_pct ~ estimate*10,
                                 immune_measure == "White Blood Cells" ~ estimate * 1000*10,
                                 immune_measure == "Red Blood Cells" ~ estimate * 1000000*10,
                                 immune_measure == "Mean Corpuscular Volume" ~ estimate*10)) %>%
    mutate(lower_ci = case_when(immune_measure %in% chem_pct ~ lower.CI*10,
                                immune_measure == "White Blood Cells" ~ lower.CI * 1000*10,
                                immune_measure == "Red Blood Cells" ~ lower.CI * 1000000*10,
                                immune_measure == "Mean Corpuscular Volume" ~ lower.CI*10)) %>%
    mutate(upper_ci = case_when(immune_measure %in% chem_pct ~ upper.CI*10,
                                immune_measure == "White Blood Cells" ~ upper.CI * 1000*10,
                                immune_measure == "Red Blood Cells" ~ upper.CI * 1000000*10,
                                immune_measure == "Mean Corpuscular Volume" ~ upper.CI*10)) %>%
    mutate(units = case_when(immune_measure %in% chem_pct ~ "%",
                             immune_measure == "White Blood Cells" ~ "cells per uL",
                             immune_measure == "Red Blood Cells" ~ "cells per uL",
                             immune_measure == "Mean Corpuscular Volume" ~ "fL"))
  
  #############################################################################################################
  ############################################# Make The Table Nice ###########################################
  #############################################################################################################
  
  # Clean up results table
  model_clean <- model_interpret %>%
    dplyr::select(-statistic,
                  -p.value,
                  -lower.CI,
                  -upper.CI,
                  -std.error,
                  -estimate) %>%
    relocate(immune_measure, .before = "term") %>%
    relocate(interpret, .after = "term") %>%
    relocate(lower_ci, .after = "interpret") %>%
    relocate(upper_ci, .after = "lower_ci") %>%
    mutate(lower_ci = sprintf("%0.2f", lower_ci),
           upper_ci = sprintf("%0.2f", upper_ci)) %>%
    mutate(ci = paste0(lower_ci, " - ", upper_ci)) %>%
    relocate(ci, .after = "interpret") %>%
    dplyr::select(-lower_ci,
                  -upper_ci)
  
  
  # Set wd
  setwd(paste0(current_directory, "/Tables - Table 1"))
  
  # Make a gt table for saving as html
  model_clean %>%
    dplyr::select(-term) %>%
    gt() %>%
    cols_label(interpret = "Beta Interpretation",
               immune_measure = "Immune Measure",
               ci = "95% Confidence Interval",
               units = "Units") %>%
    fmt_number(columns = c("interpret"),
               decimals = 2) %>%
    fmt_scientific(columns = "FDR",
                   decimals = 1) %>%
    cols_align(align = "center",
               columns = c("interpret",
                           "ci",
                           "FDR")) %>%
    cols_align(align = "right",
               columns = "units") %>%
    tab_options(column_labels.font.weight = "bold") %>%
    gtsave("interpreted_age_regressions.html")
  
  
  View(model_clean)

  # Reset the working directory
  setwd(current_directory)
}