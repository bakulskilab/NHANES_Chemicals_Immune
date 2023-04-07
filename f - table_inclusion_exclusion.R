#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################  MAKE TABLES OF INCLUDED AND EXCLUDED PARTICIPANT INFORMATION  #######################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function makes supp table 2 using demographic and immune measure data by inclusion status
#          
# Inputs: nhanes_subset - dataframe containing complete demographic and immune measures data for each
#                         participant
# 
#         nhanes_full_dataset - dataframe containing all participants, demographic, chemical, and other
#                               variables
#
# Outputs: Table of means, medians, and p-values by inclusion status

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_inclusion_exclusion <- function(nhanes_subset,
                                      nhanes_full_dataset,
                                      demog_dataset)
{
  library(sjlabelled)
  library(tidyverse)
  library(survey)
  library(gtsummary)
  library(flextable)
  
  #set the number of digits after the decimal point to be 2
  options(digits=4)
  
  #TEMPORARY
  # nhanes_subset <- nhanes_subset_dataset
  # nhanes_full_dataset <- nhanes_merged_dataset
  # demog_dataset <- demographics_clean
  
  #############################################################################################################
  ###################################### Select Variables, Merge Datasets #####################################
  #############################################################################################################
  
  # Select the variables to use from nhanes subset and demographics dataset
  # Add a smoking variable
  nhanes_vars <- nhanes_full_dataset %>%
    mutate(SMOKING = LBXCOT) %>%
    dplyr::select(SEQN,
                  RIDRETH1,
                  RIDAGEYR,
                  RIAGENDR,
                  INDFMPIR,
                  BMXWAIST,
                  SDDSRVYR,
                  SMOKING,
                  URXUCR,
                  "LBXLYPCT",
                  "LBXMOPCT",
                  "LBXNEPCT",
                  "LBXEOPCT",
                  "LBXBAPCT",
                  "LBXWBCSI",
                  "LBXRBCSI",
                  "LBXMCVSI",
                  SDMVPSU,
                  SDMVSTRA)
  demog_vars <- demog_dataset %>%
    dplyr::select(SEQN,
                  WTMEC4YR,
                  WTMEC2YR)
  
  # Merge the two subsets of variables and keep only participants previously selected
  vars_weights <- left_join(nhanes_vars, demog_vars, by = "SEQN")
  

  #############################################################################################################
  ########################################## Identify The Variables ###########################################
  #############################################################################################################
  
  demog <- c("SEQN",     #ID
             "RIDAGEYR", #age
             "RIDRETH1", #race
             "RIAGENDR", #gender
             "INDFMPIR", #poverty-income ratio
             "BMXWAIST", #waist circumference
             "URXUCR",   #creatinine
             "SMOKING",  #smoking
             "SDDSRVYR"  #survey cycle
            )
  
  #cells
  cells <- c("LBXLYPCT", #lymphocytes percent
             "LBXMOPCT", #monocytes
             "LBXNEPCT", #neutrophils percent
             "LBXEOPCT", #eosinophils
             "LBXBAPCT", #basophils
             "LBXWBCSI", #WBC count
             "LBXRBCSI", #RBC count
             "LBXMCVSI"  #MCV
            )
  
  # weights
  weights <- c("SDMVPSU",
               "SDMVSTRA",
               "WTMEC4YR",
               "WTMEC2YR")
  
  #############################################################################################################
  ################################# Create The Excluded And Included Datasets #################################
  #############################################################################################################
  
  # Get the survey variables
  psu <- nhanes_full_dataset %>%
    dplyr::select(SEQN,
                  SDMVPSU,
                  SDMVSTRA)
  
  # Get the included participants
  included <- nhanes_subset %>%
    dplyr::select(all_of(demog),
                  all_of(cells)) %>%
    left_join(., psu, by = "SEQN") %>%
    left_join(., demog_vars, by = "SEQN") %>%
    mutate(status = "Included") %>%
    dplyr::select(all_of(demog),
                  all_of(cells),
                  all_of(weights),
                  status)
  dim(included)
  #45528    22
  # str(included$RIDRETH1)
  
  # Get the excluded participants
  excluded_unclean <- anti_join(vars_weights, nhanes_subset, by = "SEQN")
  dim(excluded_unclean)
  #dim: 55788  21
  excluded <- excluded_unclean %>%
    # mutate(RIDRETH1 = as.integer(RIDRETH1)) %>%
    # mutate(RIAGENDR = as.integer(RIAGENDR)) %>%
    mutate(status = "Excluded") %>%
    dplyr::select(all_of(demog),
                  all_of(cells),
                  all_of(weights),
                  status)
  dim(excluded)
  #55788    22
  # str(excluded$RIDRETH1)
  
  # Check order
  identical(colnames(included), colnames(excluded))
  
  # Put the included and excluded participants back together now that they're labeled
  full_nhanes_demog <- bind_rows(included, excluded)
  dim(full_nhanes_demog)
  #101316     22
  
  
  #############################################################################################################
  ############################################### Clean Weights ###############################################
  #############################################################################################################
  
  # Make a subset of participants who were in cycles 1 and 2
  cycles1_2 <- c(1:2)
  vars_12 <- full_nhanes_demog %>%
    filter(SDDSRVYR %in% cycles1_2)
  
  # Make WTMEC2YR NA because CDC says to use the MEC4 weights for the first two cycles
  vars_12$WTMEC2YR <- NA
  
  # Check that all participants have weights
  summary(vars_12$WTMEC4YR)
  sum(is.na(vars_12$WTMEC4YR))
  sum(!is.na(vars_12$WTMEC4YR))
  
  
  # Make a subset of participants after cycle 2
  cycle_other <- c(3:10)
  vars_other <- full_nhanes_demog %>%
    filter(SDDSRVYR %in% cycle_other)
  
  # Make WTMEC4YR NA because CDC says to use the MEC2 weights for cycles after 2
  vars_other$WTMEC4YR <- NA
  
  # Check that all participants have weights
  summary(vars_other$WTMEC2YR)
  sum(is.na(vars_other$WTMEC2YR))
  sum(!is.na(vars_other$WTMEC2YR))
  # Everyone has a weight and neither set of weights has extra rows
  
  
  # Combine the two subsets back together
  vars_weights_merge <- bind_rows(vars_12, vars_other)
  dim(vars_weights_merge)
  
  # Multiply MEC4 weights by 2/10 and MEC2 by 1/10
  # Merge the weights columns
  vars_weights_clean <- vars_weights_merge %>%
    mutate(mec4_clean = WTMEC4YR * (2/10),
           mec2_clean = WTMEC2YR * (1/10)) %>%
    mutate(weights_adjusted = coalesce(mec4_clean, mec2_clean)) %>%
    mutate(SDDSRVYR = factor(SDDSRVYR)) %>%
    dplyr::select(-WTMEC4YR,
                  -WTMEC2YR,
                  -mec4_clean,
                  -mec2_clean)
  sum(is.na(vars_weights_clean$weights_adjusted))
  sum(vars_weights_clean$weights_adjusted)
  #299 million - some people are weighted 0
  
  #############################################################################################################
  ######################################### Set The Labels And Levels #########################################
  #############################################################################################################
  
  # Set the levels for gender and race/eth
  vars_weights_clean <- vars_weights_clean %>%
    mutate(RIAGENDR = factor(RIAGENDR,
                             levels=c(1, 2),
                             labels=c("Male", "Female")),
           RIDRETH1 = factor(RIDRETH1,
                             levels = c(1, 2, 3, 4, 5),
                             labels=c("Mexican American",
                                      "Other Hispanic",
                                      "Non-Hispanic White",
                                      "Non-Hispanic Black",
                                      "Other Race")),
           status = factor(status,
                           levels = c("Included", "Excluded")))
  
  str(vars_weights_clean)
  
  # Set the labels for all the variables
  vars_weights_clean$RIDAGEYR <- set_label(vars_weights_clean$RIDAGEYR, "Age (years)")
  vars_weights_clean$RIAGENDR <- set_label(vars_weights_clean$RIAGENDR, "Sex")
  vars_weights_clean$RIDRETH1 <- set_label(vars_weights_clean$RIDRETH1, "Race/Ethnicity")
  vars_weights_clean$INDFMPIR <- set_label(vars_weights_clean$INDFMPIR, "Family PIR")
  vars_weights_clean$BMXWAIST <- set_label(vars_weights_clean$BMXWAIST, "Waist Circumference (cm)")
  vars_weights_clean$URXUCR   <- set_label(vars_weights_clean$URXUCR,   "Urinary Creatinine")
  vars_weights_clean$SMOKING  <- set_label(vars_weights_clean$SMOKING,  "Cotinine (ng/mL)")
  vars_weights_clean$LBXWBCSI <- set_label(vars_weights_clean$LBXWBCSI, "WBC count (1000 cells/uL)")
  vars_weights_clean$LBXRBCSI <- set_label(vars_weights_clean$LBXRBCSI, "RBC count (million cells/uL)")
  vars_weights_clean$LBXMCVSI <- set_label(vars_weights_clean$LBXMCVSI, "Mean corpuscular volume (fL)")
  vars_weights_clean$LBXLYPCT <- set_label(vars_weights_clean$LBXLYPCT, "Lymphocyte (%)")
  vars_weights_clean$LBXMOPCT <- set_label(vars_weights_clean$LBXMOPCT, "Monocyte (%)") 
  vars_weights_clean$LBXNEPCT <- set_label(vars_weights_clean$LBXNEPCT, "Segmented neutrophil (%)")
  vars_weights_clean$LBXEOPCT <- set_label(vars_weights_clean$LBXEOPCT, "Eosinophil (%)")
  vars_weights_clean$LBXBAPCT <- set_label(vars_weights_clean$LBXBAPCT, "Basophil (%)")
  vars_weights_clean$SDDSRVYR <- set_label(vars_weights_clean$SDDSRVYR, "Survey Cycle")
  
  str(vars_weights_clean)
  
  #############################################################################################################
  ######################################## Create Survey Design Object ########################################
  #############################################################################################################
  
  #Chirag suggests to remove the "lonely" PSUs - strata with only one PSU
  #this is because "a single-PSU stratum makes no contribution to the variance"
  # - https://r-survey.r-forge.r-project.org/survey/html/surveyoptions.html
  options(survey.lonely.psu = "remove")
  
  # Pull everything together to get the survey adjustment
  NHANES.svy <- svydesign(strata = ~SDMVSTRA
                          , id = ~SDMVPSU
                          , weights = ~weights_adjusted
                          , data = vars_weights_clean
                          , nest = TRUE)

  #############################################################################################################
  ########################################## Calculate The Statistics #########################################
  #############################################################################################################

  #set directory
  setwd(paste0(current_directory, "/Tables - Table 1"))
  
  NHANES.svy %>%
    tbl_svysummary(by = "status", #stratify by included/excluded
                   include = c(-SEQN,
                               -SDDSRVYR),
                   statistic = list(all_continuous() ~ "{mean} ({sd})",
                                    all_categorical() ~ "{n} ({p}%)"),
                   digits = list(all_categorical() ~ c(0, 1),
                                 all_continuous() ~ 1),
                   missing_text = "Missing (n)",
                   sort = all_categorical() ~ "frequency"
                   ) %>%
    add_p() %>%
    add_overall() %>%
    modify_header(label ~ "**Variable**") %>%
    bold_labels() %>%
    as_flex_table() %>%
    save_as_docx(path = "included-excluded_demog_stats.docx") #Export the table to Word using flextable package

  #############################################################################################################
  
  setwd(current_directory)
}