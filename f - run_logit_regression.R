#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################################  Adjusted Logistic Regression  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function runs the adjusted logistic regressions for WBCs, cotinine and blood cadmium
#          
# Inputs: long_nhanes_subset - long dataframe containing complete demographic and immune data for each
#                              participant
#         conversion - dataframe of chemical names, codenames, and families
#
# Outputs: model_stats_logit - dataframe of logistic regression outputs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

run_logit_regression <- function(long_nhanes_subset,
                                 conversion,
                                 weights_dataset)
{
  library(tidyverse)
  library(broom)
  library(survey)
  library(gtsummary)
  library(gt)
  
  #TEMPORARY
  # long_nhanes_subset <- long_nhanes_subset_dataset
  # conversion <- use_these_chems
  # weights_dataset <- survey_weights
  
  #Outline
  # Make and clean a long dataset of covariates, survey variables + cotinine/cadmium/copper + WT_LBXCOT/WT_LBXBCD/WT_LBXSCU
  #   Recode WBC values as Normal or High, drop low
  #   Adjust the weights
  # Split the dataset into two for each of the chemicals - don't need to scale because no figure
  # Make survey object for each
  # Run regression (Could use tbl_regression to display the results of each? - run reg and imput model object)
  # (Merge in chem names and generally clean up if not using tbl_regression)
  
  #############################################################################################################
  ############################################# Clean Long Dataset ############################################
  #############################################################################################################
  
  long_nhanes_clean <- long_nhanes_subset %>%
    filter(chemical_codename == "LBXCOT" |
             chemical_codename == "LBXBCD" |
             chemical_codename == "LBXSCU") %>%
    filter(celltype_codename == "LBXWBCSI") %>%
    mutate(wbc_categ = case_when(cell_measurement < 4.5 ~ "Low",
                                 cell_measurement >= 4.5 & cell_measurement <= 11 ~ "Normal",
                                 cell_measurement > 11 ~ "High")) %>%
    filter(!wbc_categ == "Low") %>%
    mutate(wbc_categ = relevel(factor(wbc_categ, levels = c("Normal", "High")),
                               ref = "Normal")) %>%
    mutate(RIDRETH1 = relevel(factor(RIDRETH1,
                             levels = c(1, 2, 3, 4, 5),
                             labels=c("Mexican American",
                                      "Other Hispanic",
                                      "Non-Hispanic White",
                                      "Non-Hispanic Black",
                                      "Other Race")),
                             ref = 3),
           RIAGENDR = relevel(factor(RIAGENDR,
                                     levels = c(1, 2),
                                     labels = c("Male", "Female")),
                              ref = 1),
           SDDSRVYR = factor(SDDSRVYR))
  
  
  str(long_nhanes_clean) #dropped WBC low
  
  #############################################################################################################
  ############################################ Merge in the Weights ###########################################
  #############################################################################################################
  
  # Select the columns of weights for cotinine and cadmium
  chem_weights <- weights_dataset %>%
    dplyr::select(SEQN,
                  WT_LBXCOT,
                  WT_LBXBCD,
                  WT_LBXSCU)
  
  # Merge weights into covariates dataset
  nhanes_weights <- left_join(long_nhanes_clean, chem_weights, by = "SEQN")
  
  #############################################################################################################
  ############################################# Adjust the Weights ############################################
  #############################################################################################################
  
  # Split dataset into cycles 1-2
  # Multiply weights from cycles 1-2 by (2/10)
  cycles1_2 <- c(1, 2)
  nhanes_weights_1_2 <- nhanes_weights %>%
    filter(SDDSRVYR %in% cycles1_2) %>%
    mutate(WT_LBXCOT_adj = WT_LBXCOT*(2/10)) %>%
    mutate(WT_LBXBCD_adj = WT_LBXBCD*(2/10)) %>%
    mutate(WT_LBXSCU_adj = WT_LBXSCU) %>% #this is just a place holder column 
    dplyr::select(-WT_LBXCOT,
                  -WT_LBXBCD,
                  -WT_LBXSCU)
  
  # and cycles 3-10
  # Multiply weights from cycles 3-10 by (1/10)
  cycles3_10 <- 3:10
  nhanes_weights_3_10 <- nhanes_weights %>%
    filter(SDDSRVYR %in% cycles3_10) %>%
    mutate(WT_LBXCOT_adj = WT_LBXCOT*(1/10)) %>%
    mutate(WT_LBXBCD_adj = WT_LBXBCD*(1/10)) %>%
    mutate(WT_LBXSCU_adj = WT_LBXSCU*(1/3)) %>% #there are only 3 cycles of copper
    dplyr::select(-WT_LBXCOT,
                  -WT_LBXBCD,
                  -WT_LBXSCU)
  
  # Put the datasets back together
  nhanes_weights_adj <- bind_rows(nhanes_weights_1_2, nhanes_weights_3_10)
  
  #############################################################################################################
  ########################################### Split Chemical Datasets #########################################
  #############################################################################################################
  
  #cotinine
  cot_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "LBXCOT")
  
  #cadmium
  cad_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "LBXBCD")
  
  #copper
  cu_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "LBXSCU")
  
  str(cad_dataset)
  
  #############################################################################################################
  ############################################# Make Survey Objects ###########################################
  #############################################################################################################
  
  # Pull everything together to get the survey adjustment
  svy_cot <- svydesign(strata = ~SDMVSTRA
                       , id = ~SDMVPSU
                       , weights = ~WT_LBXCOT_adj
                       , data = cot_dataset
                       , nest = TRUE)
  
  svy_cad <- svydesign(strata = ~SDMVSTRA
                       , id = ~SDMVPSU
                       , weights = ~WT_LBXBCD_adj
                       , data = cad_dataset
                       , nest = TRUE)
  
  svy_cu <- svydesign(strata = ~SDMVSTRA
                       , id = ~SDMVPSU
                       , weights = ~WT_LBXSCU_adj
                       , data = cu_dataset
                       , nest = TRUE)
  

  #############################################################################################################
  ############################################## Run Regressions ##############################################
  #############################################################################################################
  
  model_cot <- svyglm(wbc_categ ~
                        chem_log_measurement+
                        RIDAGEYR+
                        INDFMPIR+
                        BMXWAIST+
                        RIAGENDR+
                        RIDRETH1+
                        SDDSRVYR,
                  na.action = na.omit,
                  design = svy_cot,
                  data = cot_dataset,
                  family = "binomial")
  
  model_cad <- svyglm(wbc_categ ~
                        chem_log_measurement+
                        RIDAGEYR+
                        INDFMPIR+
                        BMXWAIST+
                        SMOKING+
                        RIAGENDR+
                        RIDRETH1+
                        SDDSRVYR,
                      na.action = na.omit,
                      design = svy_cad,
                      data = cad_dataset,
                      family = "binomial")
  
  model_cu <- svyglm(wbc_categ ~
                        chem_log_measurement+
                        RIDAGEYR+
                        INDFMPIR+
                        BMXWAIST+
                        SMOKING+
                        RIAGENDR+
                        RIDRETH1+
                        SDDSRVYR,
                      na.action = na.omit,
                      design = svy_cu,
                      data = cu_dataset,
                      family = "binomial")
  
  #############################################################################################################
  ############################################# Calculate p-values ############################################
  #############################################################################################################
  
  cad_fdr <- tidy(model_cad) %>%
    filter(term == "chem_log_measurement")
  rownames(cad_fdr) <- "cadmium"
  
  cot_fdr <- tidy(model_cot) %>%
    filter(term == "chem_log_measurement")
  rownames(cot_fdr) <- "cotinine"
  
  cu_fdr <- tidy(model_cu) %>%
    filter(term == "chem_log_measurement")
  rownames(cu_fdr) <- "copper"
  
  fdr_chems <- bind_rows(cad_fdr, cot_fdr, cu_fdr) %>%
    rownames_to_column(var = "chemical") %>%
    dplyr::select(chemical,
                  estimate,
                  std.error,
                  p.value)
  print(fdr_chems)
  
  #############################################################################################################
  ############################################### Create Tables ###############################################
  #############################################################################################################
  
  # cotinine <- model_cot %>%
  #   tbl_regression(exponentiate = TRUE,
  #                  label = c(chem_log_measurement ~ "log2(chemical)",
  #                            RIDAGEYR ~ "Age (years)",
  #                            INDFMPIR ~ "Poverty-Income Ratio",
  #                            BMXWAIST ~ "Waist Circumference (cm)",
  #                            RIAGENDR ~ "Sex",
  #                            RIDRETH1 ~ "Race/Ethnicity",
  #                            SDDSRVYR ~ "Survey Cycle")) %>%
  #   bold_labels() %>%
  #   modify_header(label ~ "**Variable**") %>%
  #   bold_p(t = 0.05)
  # 
  # cadmium <- model_cad %>%
  #   tbl_regression(exponentiate = TRUE,
  #                  label = c(chem_log_measurement ~ "log2(chemical)",
  #                            RIDAGEYR ~ "Age (years)",
  #                            INDFMPIR ~ "Poverty-Income Ratio",
  #                            BMXWAIST ~ "Waist Circumference (cm)",
  #                            SMOKING  ~ "Cotinine (ng/mL)",
  #                            RIAGENDR ~ "Sex",
  #                            RIDRETH1 ~ "Race/Ethnicity",
  #                            SDDSRVYR ~ "Survey Cycle")) %>%
  #   bold_labels() %>%
  #   modify_header(label ~ "**Variable**") %>%
  #   bold_p(t = 0.05)
  # 
  # copper <- model_cu %>%
  #   tbl_regression(exponentiate = TRUE,
  #                  label = c(chem_log_measurement ~ "log2(chemical)",
  #                            RIDAGEYR ~ "Age (years)",
  #                            INDFMPIR ~ "Poverty-Income Ratio",
  #                            BMXWAIST ~ "Waist Circumference (cm)",
  #                            SMOKING  ~ "Cotinine (ng/mL)",
  #                            RIAGENDR ~ "Sex",
  #                            RIDRETH1 ~ "Race/Ethnicity",
  #                            SDDSRVYR ~ "Survey Cycle")) %>%
  #   bold_labels() %>%
  #   modify_header(label ~ "**Variable**") %>%
  #   bold_p(t = 0.05)
  
  #############################################################################################################
  ################################################ Save Results ###############################################
  #############################################################################################################
  
  # Set directory
  # setwd(paste0(current_directory, "/Regression Results"))
  # 
  # # Combine the cotinine and cadmium results into one table
  # tbl_merge(tbls = list(cadmium, cotinine, copper),
  #           tab_spanner = c("**Blood Cadmium**", "**Cotinine**", "**Copper**")) %>%
  #   as_gt() %>%
  #   gtsave(filename = "table_logistic_reg.html")
  # 
  # # Reset directory
  # setwd(current_directory)
}