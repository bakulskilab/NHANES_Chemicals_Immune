run_if_else_glm_weighted <- function(chemical_immune_chunk,
                                     weights_dataset,
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
  chem_weight_codename <- paste0("WT_", chem_codename)
  
  immune_codename <- chemical_immune_chunk$celltype_codename %>%
    unique(.)
  print(paste(immune_codename, chem_codename))
  print(chem_codename)
  
  #Calculate the number of cycles per chemical
  # sort(unique(chemical_immune_chunk$SDDSRVYR))
  cycle_length <- length(unique(chemical_immune_chunk$SDDSRVYR))
  # print(paste("Number of cycles:", cycle_length))
  print(unique(chemical_immune_chunk$SDDSRVYR))
  
  #remove NAs from entire dataset
  chemical_immune_chunk <- chemical_immune_chunk %>%
    na.omit()

  # grabs the unique data for the cycle year (SDDSRVYR) for non-NA chemical concentrations
  cycle_unique <- unique(chemical_immune_chunk$SDDSRVYR)

  # print(str_detect(chem_codename, "^LB"))
  
  #mute the summarize grouping message to reduce clutter in the console output
  options(dplyr.summarise.inform = FALSE)
  
  #############################################################################################################
  #################################### Adjust the Weights Based on Cycle Info #################################
  #############################################################################################################

  #select SEQN and chemical/immune combo
  subset_weights_dataset <- weights_dataset %>%
    dplyr::select(SEQN,
                  all_of(chem_weight_codename))
  
  #merge in the weights dataset into the chemical_immune_chunk
  chunk_and_weight <- left_join(chemical_immune_chunk, subset_weights_dataset, by = "SEQN") %>%
    rename_at(vars(starts_with("WT_")), ~ "unadj_weight") %>%
    drop_na(unadj_weight) %>%
    filter(!unadj_weight == 0)
  # View(chunk_and_weight)
  
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
  identifiers <- test_tab %>%
    mutate(category_id = case_when(cycle1 == 1 & cycle2 == 0 ~ "1 or 2",
                                   cycle2 == 1 & cycle1 == 0 ~ "1 or 2",
                                   cycle1 == 1 & cycle2 == 1 ~ "1 and 2"
                                   )) %>%
    mutate(category_id = ifelse(is.na(category_id), "other", category_id))
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
  # print(str(chunk_and_ad_weight$RIDRETH1))
  # View(chunk_and_ad_weight)

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
                      data = LR_data)
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
                      data = LR_data)
    }
    # print(tidy(model))
    print("chemical in more than one survey cycle")
    #
    #   #survey cycle ==1 and blood
  } else
  {
    if (cycle_length == 1 & str_detect(chem_codename, "^LB"))
    {
      # print(2)
      print(paste0("survey cycle: ", unique(LR_data$SDDSRVYR)))

      model <- svyglm(cell_measurement ~
                        chem_log_measurement+
                        RIDRETH1+
                        RIDAGEYR+
                        RIAGENDR+
                        INDFMPIR+
                        BMXWAIST+
                        SMOKING,
                   na.action = na.omit,
                   design = NHANES.svy,
                   data = LR_data)
      # print(tidy(model))
      print("chemical only in one survey cycle")
      #
      #     #survey cycle >1 and urinary
    } else {
      if (cycle_length > 1 & str_detect(chem_codename, "^LB", negate = TRUE))  {
        # min_cycle <- min(LR_data$SDDSRVYR, na.rm = TRUE)
        # print(3)
        # print(paste0("reference cycle: ", min_cycle))

        LR_data$SDDSRVYR <- factor(LR_data$SDDSRVYR)

        # print(levels(LR_data$SDDSRVYR))

        model <- svyglm(cell_measurement ~
                          chem_log_measurement+
                          RIDRETH1+
                          RIDAGEYR+
                          RIAGENDR+
                          INDFMPIR+
                          BMXWAIST+
                          URXUCR+
                          SDDSRVYR+
                          SMOKING,
                     na.action = na.omit,
                     design = NHANES.svy,
                     data = LR_data)
        # print(tidy(model))
        print("chemical in more than one survey cycle")
        #
        #
      } else { #survey cycle ==1 and urinary
        # print(4)
        # print(paste0("survey cycle: ", unique(LR_data$SDDSRVYR)))

        model <- svyglm(cell_measurement ~
                          chem_log_measurement+
                          RIDRETH1+
                          RIDAGEYR+
                          RIAGENDR+
                          INDFMPIR+
                          BMXWAIST+
                          URXUCR+
                          SMOKING,
                     na.action = na.omit,
                     design = NHANES.svy,
                     data = LR_data)
        # print(tidy(model))
        print("chemical only in one survey cycle")
      }
    }
  }

  #############################################################################################################
  ######################################### Turn Regressions Into Table #######################################
  #############################################################################################################

  #compile regression models using broom
  df_tidy <- tidy(model)
  df_tidy <- as.data.frame(df_tidy)

  return(df_tidy)
}