run_if_else_glm_logistic_wt <- function(chemical_immune_chunk,
                                        weights_dataset,
                                        df_chem_master, 
                                        nhanes_subset)
{
  library(tidyverse)
  library(broom)
  library(survey)
  
  #############################################################################################################
  ################################################ Get Basic Info #############################################
  #############################################################################################################
  
  
  chem_codename <- chemical_immune_chunk$chemical_codename %>%
    unique(.)
  # print(chem_codename)
  
  #update the name to include WT to match the weights codenames for the chemicals
  if(chem_codename == "LBX138LA")
  {
    chem_weight_codename <- "WT_LBX138158LA"
    
  } else if(chem_codename == "LBX196LA") {
    
    chem_weight_codename <- "WT_LBX196203LA"
    
  } else {
    
    chem_weight_codename <- paste0("WT_", chem_codename)
    
  }
  
  immune_codename <- chemical_immune_chunk$celltype_codename %>%
    unique(.)
  print(paste(immune_codename, chem_codename))
  # print(chem_codename)
  
  #Calculate the number of cycles per chemical
  # sort(unique(chemical_immune_chunk$SDDSRVYR))
  cycle_length <- length(unique(chemical_immune_chunk$SDDSRVYR))
  print(paste("Number of cycles:", cycle_length))
  print(unique(chemical_immune_chunk$SDDSRVYR))
  
  # View(chemical_immune_chunk)
  #remove NAs from entire dataset, 
  #this will remove participants who don't have urinary creatinine measured but have measurement
  #for a chemical measured in blood
  chemical_immune_chunk <- chemical_immune_chunk 
    # na.omit()
  # print(dim(chemical_immune_chunk))

  # grabs the unique data for the cycle year (SDDSRVYR) for non-NA chemical concentrations
  cycle_unique <- unique(chemical_immune_chunk$SDDSRVYR)

  # Check that cycle length before and after removing NA are the same
  if(cycle_length != length(cycle_unique))
  {
    print("Cycle length before and after removing NA are NOT the same")
    print(paste(immune_codename, chem_codename))
    print(cycle_length)
    print(length(cycle_unique))
  }
  
  # print(str_detect(chem_codename, "^LB"))
  
  #mute the summarize grouping message to reduce clutter in the console output
  options(dplyr.summarise.inform = FALSE)
  
  #############################################################################################################
  #################################### Adjust the Weights Based on Cycle Info #################################
  #############################################################################################################

  #select SEQN and chemical/immune combo
  subset_weights_dataset <- weights_dataset %>%
    filter(SDDSRVYR != -1 ) %>%
    dplyr::select(SEQN,
                  all_of(chem_weight_codename))
  # print("dimension of subset_weights_dataset")
  # print(dim(subset_weights_dataset))
  # print(colnames(subset_weights_dataset))
  
  #merge in the weights dataset into the chemical_immune_chunk
  chunk_and_weight <- left_join(chemical_immune_chunk
                                , subset_weights_dataset
                                , by = "SEQN") %>%
    rename_at(vars(starts_with("WT_")), ~ "unadj_weight") %>%
    drop_na(unadj_weight) %>%
    filter(!unadj_weight == 0)
  # View(chunk_and_weight)
  # print(colnames(chunk_and_weight))
  # print("dimemsion of chunk_and_weight")
  # print(dim(chunk_and_weight))
  
  #select SEQN and the survey variables - they're already in chunk_and_weight
  # survey_df <- nhanes_subset %>%
  #   dplyr::select(SEQN,
  #                 SDMVPSU,
  #                 SDMVSTRA)
  # View(survey_df)
  
  #merge in the survey variables
  # chunk_weight_survey <- left_join(chunk_and_weight, survey_df, by = "SEQN")
  chunk_weight_survey <- chunk_and_weight
  # View(chunk_weight_survey)


  # #check if any chemicals were measured in only cycle 1 and not 2
  test <- long_nhanes_subset_dataset %>%
    group_by(chemical_codename) %>%
    summarise(cycle_count = unique(SDDSRVYR)) %>%
    ungroup()
  # View(test)
  test_tab <- as.data.frame.matrix(table(test$chemical_codename, test$cycle_count))
  colnames(test_tab) <- paste0("cycle", 1:10)
  # View(test_tab)
  
  df_problematic_pesticides <- find_problematic_pesticides(df_chem_master)
  # View(df_problematic_pesticides)

  identifiers <- test_tab %>%
    mutate(category_id = case_when(cycle1 == 1 & cycle2 == 0 ~ "1 or 2",
                                   cycle2 == 1 & cycle1 == 0 ~ "1 or 2",
                                   cycle1 == 1 & cycle2 == 1 ~ "1 and 2"
                                   )) %>%
    mutate(category_id = ifelse(is.na(category_id), "other", category_id))
  
  index_problematic_pesticides <- which(rownames(identifiers) %in% df_problematic_pesticides$chemical_codename_use)
  # print(rownames(identifiers)[index_problematic_pesticides])
  
  identifiers[index_problematic_pesticides,"category_id"] <- "1 or 2"
  # View(identifiers)
  # #so there are some chemicals that are only measured in cycle 1 and not 2 and also 2 and not 1

  #merge in the category_ids into chunk_weight_survey for categorization
  #make the rownames of identifiers into a column called chemical_codenames
  identifiers_chem <- rownames_to_column(identifiers, var = "chemical_codename")
  # View(identifiers_chem)
  #only keep the identifiers and chemical_codename
  id_chem <- identifiers_chem %>%
    dplyr::select(chemical_codename,
                  category_id)

  chunk_weight_id <- left_join(chunk_weight_survey, id_chem, by = "chemical_codename")
  # View(chunk_weight_id)

  # Corrects the weights according to rules for combining data across multiple NHANES cycles
  chunk_and_ad_weight <- chunk_weight_id %>%
    mutate(adjusted_weights = case_when(SDDSRVYR %in% c(1,2) & category_id == "1 and 2" ~ ((2/cycle_length)*(unadj_weight))
                                        , SDDSRVYR %in% c(1,2) & category_id == "1 or 2" ~ ((1/cycle_length)*(unadj_weight))
                                        , SDDSRVYR %in% c(3:10) ~ (1/cycle_length)*(unadj_weight)
                                        )) 
  # View(chunk_and_ad_weight)
  # print(str(chunk_and_ad_weight$RIDRETH1))
  # print("Dimension of chunk_and_ad_weight")
  # print(dim(chunk_and_ad_weight %>% na.omit(.)))
  
  # Check for the calculation of the adjusted survey weights
  checking_adjusted_weights <- chunk_and_ad_weight %>%
    select(SDDSRVYR
           , category_id
           , unadj_weight
           , adjusted_weights) %>%
    mutate(multiplier_checked = adjusted_weights/unadj_weight) %>%
    mutate(multiplier_original = case_when(SDDSRVYR %in% c(1,2) & category_id == "1 and 2" ~ ((2/cycle_length))
                                           , SDDSRVYR %in% c(1,2) & category_id == "1 or 2" ~ ((1/cycle_length))
                                           , SDDSRVYR %in% c(3:10) ~ (1/cycle_length))) %>%
    select(SDDSRVYR
           , multiplier_checked
           , multiplier_original) %>%
    unique(.) %>%
    # The difference between the original multiplier and the checked multiplier should be 0 or
    # near 0 due to roundoff issues
    mutate(checking_multiplier = (multiplier_original - multiplier_checked) %>%
             round(., digits = 0))
  
  diff_original_checked_multiplier <- checking_adjusted_weights %>%
    pull(checking_multiplier) %>%
    unique(.)
  
  # The difference is a single value and is non-zero. 
  if(length(diff_original_checked_multiplier) == 1)
  {
    
    if(diff_original_checked_multiplier != 0)
    {
      print("Incorrect definition multipler for calculating adjusted weights")
      View(checking_adjusted_weights)
      print(diff_original_checked_multiplier)
    }
    
  # The difference has multiple values. 
  } else if(length(diff_original_checked_multiplier) > 1) {
    print("Incorrect definition multipler for calculating adjusted weights")
    View(checking_adjusted_weights)
    print(diff_original_checked_multiplier)
  }
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ACCOUNT FOR NHANES SAMPLING DESIGN  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #Chirag suggests to remove the "lonely" PSUs - strata with only one PSU
  #this is because "a single-PSU stratum makes no contribution to the variance"
  # - https://r-survey.r-forge.r-project.org/survey/html/surveyoptions.html
  options(survey.lonely.psu = "remove")

  # Define a data frame containing the outcome variable and necessary covariates
  LR_data <- data.frame(SEQN = factor(chunk_and_ad_weight$SEQN)
                        , SDDSRVYR = factor(chunk_and_ad_weight$SDDSRVYR)
                        , RIDAGEYR = I(chunk_and_ad_weight$RIDAGEYR)
                        , RIAGENDR = relevel(factor(chunk_and_ad_weight$RIAGENDR), ref = "1")
                        , RIDRETH1 = relevel(factor(chunk_and_ad_weight$RIDRETH1), ref = "3")
                        , INDFMPIR = I(chunk_and_ad_weight$INDFMPIR)
                        , BMXWAIST = I(chunk_and_ad_weight$BMXWAIST)
                        , URXUCR = I(chunk_and_ad_weight$URXUCR)
                        , SMOKING = chunk_and_ad_weight$SMOKING
                        , cell_measurement = I(chunk_and_ad_weight$cell_measurement)
                        , chem_log_measurement = I(chunk_and_ad_weight$chem_log_measurement)
                        , SDMVPSU = chunk_and_ad_weight$SDMVPSU #stage 1 sampling unit in NHANES
                        , SDMVSTRA = chunk_and_ad_weight$SDMVSTRA #stratum to which the PSU belongs
                        , adjusted_weights = chunk_and_ad_weight$adjusted_weights)
  # print("Dimension of LR_data")
  # print(dim(LR_data %>% na.omit(.)))
  

  # print(LR_data$RIDRETH1 %>% levels(.))

  #pull everything together to get the survey adjustment
  NHANES.svy <- svydesign(strata = ~SDMVSTRA
                          , id = ~SDMVPSU
                          , weights = ~adjusted_weights
                          , data = LR_data
                          , nest = TRUE)

  # print("nhanes.svy design variable created")

  #############################################################################################################
  ########################################### Sort And Run Regressions ########################################
  #############################################################################################################

  #survey cycle >1 and blood
  if (cycle_length > 1 & str_detect(chem_codename, "^LB")) {
    # min_cycle <- min(LR_data$SDDSRVYR, na.rm = TRUE)
    # print(1)
    # print(paste0("reference cycle: ", min_cycle))

    LR_data$SDDSRVYR <- factor(LR_data$SDDSRVYR)

    # print(levels(LR_data$SDDSRVYR))

    if(chem_codename == "LBXCOT")
    {
      model <- svyglm(cell_measurement ~
                        chem_log_measurement+
                        RIDRETH1+
                        RIDAGEYR+
                        RIAGENDR+
                        INDFMPIR+
                        BMXWAIST+
                        SDDSRVYR,
                      na.action = na.omit,
                      design = NHANES.svy,
                      data = LR_data,
                      family = "binomial")
  
      
      
    } else {
      model <- svyglm(cell_measurement ~
                        chem_log_measurement+
                        RIDRETH1+
                        RIDAGEYR+
                        RIAGENDR+
                        INDFMPIR+
                        BMXWAIST+
                        SDDSRVYR+
                        SMOKING,
                      na.action = na.omit,
                      design = NHANES.svy,
                      data = LR_data,
                      family = "binomial")
    }
  }
    # print(tidy(model))
    print("chemical in more than one survey cycle")
  #   #
  #   #   #survey cycle ==1 and blood
  # } else
  # {
  #   if (cycle_length == 1 & str_detect(chem_codename, "^LB"))
  #   {
  #     # print(2)
  #     print(paste0("survey cycle: ", unique(LR_data$SDDSRVYR)))
  # 
  #     model <- svyglm(cell_measurement ~
  #                       chem_log_measurement+
  #                       RIDRETH1+
  #                       RIDAGEYR+
  #                       RIAGENDR+
  #                       INDFMPIR+
  #                       BMXWAIST+
  #                       SMOKING,
  #                  na.action = na.omit,
  #                  design = NHANES.svy,
  #                  data = LR_data)
  #     # print(tidy(model))
  #     print("chemical only in one survey cycle")
  #     #
  #     #     #survey cycle >1 and urinary
  #   } else {
  #     if (cycle_length > 1 & str_detect(chem_codename, "^LB", negate = TRUE))  {
  #       # min_cycle <- min(LR_data$SDDSRVYR, na.rm = TRUE)
  #       # print(3)
  #       # print(paste0("reference cycle: ", min_cycle))
  # 
  #       LR_data$SDDSRVYR <- factor(LR_data$SDDSRVYR)
  # 
  #       # print(levels(LR_data$SDDSRVYR))
  # 
  #       model <- svyglm(cell_measurement ~
  #                         chem_log_measurement+
  #                         RIDRETH1+
  #                         RIDAGEYR+
  #                         RIAGENDR+
  #                         INDFMPIR+
  #                         BMXWAIST+
  #                         URXUCR+
  #                         SDDSRVYR+
  #                         SMOKING,
  #                    na.action = na.omit,
  #                    design = NHANES.svy,
  #                    data = LR_data)
  #       # print(tidy(model))
  #       print("chemical in more than one survey cycle")
  #       #
  #       #
  #     } else { #survey cycle ==1 and urinary
  #       # print(4)
  #       # print(paste0("survey cycle: ", unique(LR_data$SDDSRVYR)))
  # 
  #       model <- svyglm(cell_measurement ~
  #                         chem_log_measurement+
  #                         RIDRETH1+
  #                         RIDAGEYR+
  #                         RIAGENDR+
  #                         INDFMPIR+
  #                         BMXWAIST+
  #                         URXUCR+
  #                         SMOKING,
  #                    na.action = na.omit,
  #                    design = NHANES.svy,
  #                    data = LR_data)
  #       # print(tidy(model))
  #       print("chemical only in one survey cycle")
  #     }
  #   }
  # }

  #############################################################################################################
  #################################### Calculate sample size of each model ####################################
  #############################################################################################################
  
  nobs <- glance(model) %>%
    pull(nobs) %>%
    unlist(.)
  # print(nobs)
  
  
  #############################################################################################################
  ######################################### Turn Regressions Into Table #######################################
  #############################################################################################################

  #compile regression models using broom
  df_tidy <- tidy(model)
  df_tidy <- as.data.frame(df_tidy) %>%
    mutate(nobs = nobs)

  return(df_tidy)
}