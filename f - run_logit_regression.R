#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################################  Adjusted Logistic Regression  #######################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function runs the adjusted logistic regressions for WBCs and top 3 chems
#          
# Inputs: long_nhanes_subset - long dataframe containing complete demographic and immune data for each
#                              participant
#         conversion - dataframe of chemical names, codenames, and families
#
# Outputs: model_stats_logit - dataframe of logistic regression outputs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

run_logit_regression <- function(nhanes_subset,
                                 long_nhanes_subset,
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
  # nhanes_subset <- nhanes_subset_dataset
  
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
  
  nhanes_seqn <- nhanes_subset %>% pull(SEQN)
  
  long_nhanes_clean <- long_nhanes_subset %>%
    filter(chemical_codename == "LBXCOT" |
             chemical_codename == "URXP04" |
             chemical_codename == "URXCEM") %>%
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
  
  # Select the columns of weights for cotinine, cadmium, and copper
  chem_weights <- weights_dataset %>%
    dplyr::select(SEQN,
                  SDDSRVYR,
                  WT_LBXCOT,
                  WT_URXP04,
                  WT_URXCEM) %>%
    filter(!SDDSRVYR == -1) %>%
    filter(SEQN %in% nhanes_seqn) %>%
    select(-SDDSRVYR)
  
  
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
    mutate(WT_URXP04_adj = WT_URXP04) %>%
    mutate(WT_URXCEM_adj = WT_URXCEM) %>% #this is just a place holder column 
    dplyr::select(-WT_LBXCOT,
                  -WT_URXP04,
                  -WT_URXCEM)
  
  # and cycles 3-10
  # Multiply weights from cycles 3-10 by (1/10)
  cycles3_10 <- 3:10
  nhanes_weights_3_10 <- nhanes_weights %>%
    filter(SDDSRVYR %in% cycles3_10) %>%
    mutate(WT_LBXCOT_adj = WT_LBXCOT*(1/10)) %>%
    mutate(WT_URXP04_adj = WT_URXP04*(1/8)) %>%
    mutate(WT_URXCEM_adj = WT_URXCEM*(1/4)) %>% #there are only 4 cycles of data
    dplyr::select(-WT_LBXCOT,
                  -WT_URXP04,
                  -WT_URXCEM)
  
  # Put the datasets back together
  nhanes_weights_adj <- bind_rows(nhanes_weights_1_2, nhanes_weights_3_10)
  
  #############################################################################################################
  ########################################### Split Chemical Datasets #########################################
  #############################################################################################################
  
  #cotinine
  cot_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "LBXCOT") %>%
    drop_na(WT_LBXCOT_adj)
  
  #2-fluorene
  cad_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "URXP04") %>%
    drop_na(WT_URXP04_adj)
  
  #long name
  cu_dataset <- nhanes_weights_adj %>%
    filter(chemical_codename == "URXCEM") %>%
    drop_na(WT_URXCEM_adj)
  
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
                       , weights = ~WT_URXP04_adj
                       , data = cad_dataset
                       , nest = TRUE)
  
  svy_cu <- svydesign(strata = ~SDMVSTRA
                       , id = ~SDMVPSU
                       , weights = ~WT_URXCEM_adj
                       , data = cu_dataset
                       , nest = TRUE)
  

  #############################################################################################################
  ############################################## Run Regressions ##############################################
  #############################################################################################################
  
  # Unweighted
  # nhanes_subset <- nhanes_subset_dataset %>%
  #   mutate(wbc_categ = case_when(LBXWBCSI < 4.5 ~ "Low",
  #                                LBXWBCSI >= 4.5 & LBXWBCSI <= 11 ~ "Normal",
  #                                LBXWBCSI > 11 ~ "High")) %>%
  #   filter(!wbc_categ == "Low") %>%
  #   mutate(wbc_categ = as.factor(wbc_categ))
  # glm(formula = wbc_categ ~
  #       log2(LBXCOT)+
  #       RIDAGEYR+
  #       INDFMPIR+
  #       BMXWAIST+
  #       RIAGENDR+
  #       RIDRETH1+
  #       SDDSRVYR,
  #     data = nhanes_subset,
  #     family = "binomial") %>% tidy()
  
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
                        URXUCR+
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
                       URXUCR+
                        SDDSRVYR,
                      na.action = na.omit,
                      design = svy_cu,
                      data = cu_dataset,
                      family = "binomial")
  
  #############################################################################################################
  ############################################# Calculate p-values ############################################
  #############################################################################################################
  
  cot_fdr <- tidy(model_cot) %>%
    filter(term == "chem_log_measurement")
  rownames(cot_fdr) <- "cotinine"
  
  cad_fdr <- tidy(model_cad) %>%
    filter(term == "chem_log_measurement")
  rownames(cad_fdr) <- "2-fluorene"
  
  cu_fdr <- tidy(model_cu) %>%
    filter(term == "chem_log_measurement")
  rownames(cu_fdr) <- "N-acetyl-S-(2-Carboxyethyl)-L-cysteine"
  
  fdr_chems <- bind_rows(cot_fdr, cad_fdr, cu_fdr) %>%
    rownames_to_column(var = "chemical") %>%
    dplyr::select(chemical,
                  estimate,
                  std.error,
                  p.value) %>%
    mutate(odds_ratio = exp(estimate))
  print(fdr_chems)
  
  # exp(confint(model_cot))
  # exp(confint(model_cad))
  # exp(confint(model_cu))
  
  print("the warning messages are fine")
  
}