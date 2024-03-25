#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####################################  CALCULATE WEIGHTED SUMMARY STATISTICS  ##################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function calculates summary statistics for participants and survey weights them
#          
# Inputs:   nhanes_subset      - dataframe of complete demographics, cells, and chemicals
#
# Outputs:  weighted table 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_1_weighted <- function(nhanes_subset,
                             demog_dataset)
{
  library(tidyverse)
  library(survey)
  library(sjlabelled)
  library(gtsummary)
  library(gt)
  library(flextable)
  
  #Temporary
  # nhanes_subset <- nhanes_subset_dataset
  # demog_dataset <- demographics_clean
  
  #############################################################################################################
  ###################################### Select Variables, Merge Datasets #####################################
  #############################################################################################################
  
  # Select the variables to use from nhanes subset and demographics dataset
  nhanes_vars <- nhanes_subset %>%
    dplyr::select(SEQN,
                  RIDRETH1,
                  RIDAGEYR,
                  RIAGENDR,
                  INDFMPIR,
                  BMXWAIST,
                  SDDSRVYR,
                  SMOKING,
                  URXUCR,
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
  sum(!is.na(vars_12$WTMEC4YR)) #8355
  
  
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
  
  # Multiply MEC4 weights by 2/10 and MEC2 by 1/10
  # Merge the weights columns
  vars_weights_clean <- vars_weights_merge %>%
    mutate(mec4_clean = WTMEC4YR * (2/10),
           mec2_clean = WTMEC2YR * (1/10)) %>%
    mutate(weights_adjusted = coalesce(mec4_clean, mec2_clean)) %>%
    dplyr::select(-WTMEC4YR,
                  -WTMEC2YR,
                  -mec4_clean,
                  -mec2_clean)
  sum(is.na(vars_weights_clean$weights_adjusted))
  sum(vars_weights_clean$weights_adjusted)
  #188 million
  
  #############################################################################################################
  ######################################## Create Survey Design Object ########################################
  #############################################################################################################
  
  #Chirag suggests to remove the "lonely" PSUs - strata with only one PSU
  #this is because "a single-PSU stratum makes no contribution to the variance"
  # - https://r-survey.r-forge.r-project.org/survey/html/surveyoptions.html
  options(survey.lonely.psu = "remove")
  
  # Clean up the dataset, set the variable order
  LR_data <- vars_weights_clean %>%
    mutate(RIDRETH1 = factor(RIDRETH1,
                             levels = c(1, 2, 3, 4, 5),
                             labels=c("Mexican American",
                                      "Other Hispanic",
                                      "Non-Hispanic White",
                                      "Non-Hispanic Black",
                                      "Other Race")),
           RIAGENDR = factor(RIAGENDR,
                             levels = c(1, 2),
                             labels = c("Male", "Female")),
           SDDSRVYR = factor(SDDSRVYR)) %>%
    dplyr::select(RIAGENDR,
                  RIDRETH1,
                  RIDAGEYR,
                  INDFMPIR,
                  BMXWAIST,
                  URXUCR,
                  SMOKING,
                  SDMVPSU,
                  SDMVSTRA,
                  weights_adjusted)
  
  # Set the labels for all the variables
  LR_data$RIDAGEYR <- set_label(LR_data$RIDAGEYR, "Age (years)")
  LR_data$RIAGENDR <- set_label(LR_data$RIAGENDR, "Sex")
  LR_data$RIDRETH1 <- set_label(LR_data$RIDRETH1, "Race/Ethnicity")
  LR_data$INDFMPIR <- set_label(LR_data$INDFMPIR, "Family PIR")
  LR_data$BMXWAIST <- set_label(LR_data$BMXWAIST, "Waist Circumference (cm)")
  LR_data$URXUCR   <- set_label(LR_data$URXUCR,   "Urinary Creatinine")
  LR_data$SMOKING  <- set_label(LR_data$SMOKING,  "Cotinine (ng/mL)")
  
  str(LR_data)
    
  
  # Pull everything together to get the survey adjustment
  NHANES.svy <- svydesign(strata = ~SDMVSTRA
                          , id = ~SDMVPSU
                          , weights = ~weights_adjusted
                          , data = LR_data
                          , nest = TRUE)
  
  #############################################################################################################
  ########################################## Create Unweighted Tables #########################################
  #############################################################################################################

  # COUNT (%)
  unwt_count <- LR_data %>%
    tbl_summary(include = c(RIAGENDR,
                            RIDRETH1),
                statistic = all_categorical() ~ "{n} ({p}%)",
                digits = list(all_categorical() ~ c(1, 1)),
                missing = "no",
                sort = all_categorical() ~ "frequency"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Count (%)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)

  # MEAN (SD)
  unwt_mean <- LR_data %>%
    tbl_summary(include = c(RIDAGEYR,
                            INDFMPIR,
                            BMXWAIST,
                            URXUCR,
                            SMOKING),
                statistic = list(all_continuous() ~ "{mean} ({sd})"),
                digits = list(all_continuous() ~ c(1, 1)),
                missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Mean (sd)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  # MEDIAN (IQR)
  unwt_median <- LR_data %>%
    tbl_summary(include = c(RIDAGEYR,
                            INDFMPIR,
                            BMXWAIST,
                            URXUCR,
                            SMOKING),
                statistic = list(all_continuous() ~ "{median} ({IQR})"),
                digits = list(all_continuous() ~ c(1, 1)),
                missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Median (IQR)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  # RANGE
  unwt_range <- LR_data %>%
    tbl_summary(include = c(RIDAGEYR,
                            INDFMPIR,
                            BMXWAIST,
                            URXUCR,
                            SMOKING),
                statistic = list(all_continuous() ~ "{min} - {max}"),
                digits = list(all_continuous() ~ c(1, 1)),
                missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Range**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  
  #############################################################################################################
  ########################################### Create Weighted Tables ##########################################
  #############################################################################################################
  
  # COUNT (%)
  wt_count <- NHANES.svy %>%
    tbl_svysummary(include = c(RIAGENDR,
                               RIDRETH1),
                   statistic = all_categorical() ~ "{n} ({p}%)",
                   digits = list(all_categorical() ~ c(1, 1)),
                   missing = "no",
                   sort = all_categorical() ~ "frequency"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Count (%)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  # MEAN (SD)
  wt_mean <- NHANES.svy %>%
    tbl_svysummary(include = c(RIDAGEYR,
                               INDFMPIR,
                               BMXWAIST,
                               URXUCR,
                               SMOKING),
                   statistic = list(all_continuous() ~ "{mean} ({sd})"),
                   digits = list(all_continuous() ~ c(1, 1)),
                   missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Mean (sd)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  # MEDIAN (IQR)
  wt_median <- NHANES.svy %>%
    tbl_svysummary(include = c(RIDAGEYR,
                               INDFMPIR,
                               BMXWAIST,
                               URXUCR,
                               SMOKING),
                   statistic = list(all_continuous() ~ "{median} ({p75} - {p25})"),
                   digits = list(all_continuous() ~ c(1, 1)),
                   missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Median (IQR)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  #Calculate IQR in the most convoluted way possible, but only because I got to the end and realized I can't convert a tibble back into a gtsummary object
  #and I'm not remaking code to calculate the survey adjusted IQR
  wt_iqr_calculations <- NHANES.svy %>%
    tbl_svysummary(include = c(RIDAGEYR,
                               INDFMPIR,
                               BMXWAIST,
                               URXUCR,
                               SMOKING),
                   statistic = list(all_continuous() ~ "{median} ({p75} - {p25})"),
                   digits = list(all_continuous() ~ c(1, 1)),
                   missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Median (IQR)**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA) %>%
    as_tibble() %>%
    mutate(iqr_full = str_extract(string = .$`**Median (IQR)**`,
                                  pattern = "(?<=\\().*(?=\\))"),
           p_75 = as.numeric(gsub(" .*$", "", iqr_full)),
           p_25 = as.numeric(gsub(".*- ", "", iqr_full)),
           var_median = as.numeric(gsub(" .*$", "", `**Median (IQR)**`))) %>%
    mutate(calc_iqr = (p_75 - p_25)) %>%
    pull(calc_iqr) %>%
    as.numeric()
  
  print("paste these weighted IQR numbers into the table manually")
  print(wt_iqr_calculations)
  
  
  
  # RANGE
  wt_range <- NHANES.svy %>%
    tbl_svysummary(include = c(RIDAGEYR,
                               INDFMPIR,
                               BMXWAIST,
                               URXUCR,
                               SMOKING),
                   statistic = list(all_continuous() ~ "{min} - {max}"),
                   digits = list(all_continuous() ~ c(1, 1)),
                   missing = "no"
    ) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_header(stat_0 ~ "**Range**") %>%
    bold_labels() %>%
    modify_footnote(update = everything() ~ NA)
  
  
  #############################################################################################################
  ################################################ Merge Tables ###############################################
  #############################################################################################################
  
  #set directory
  setwd(paste0(current_directory, "/Tables - Table 1"))
  
  # Merge unweighted tables
  unwt_stats <- tbl_merge(tbls = list(unwt_count, unwt_mean, unwt_median, unwt_range)) %>%
    modify_spanning_header(update = everything() ~ NA)
  
  # Merge weighted tables
  wt_stats <- tbl_merge(tbls = list(wt_count, wt_mean, wt_median, wt_range)) %>%
    modify_spanning_header(update = everything() ~ NA)
  
  # Merge unweighted and weighted - save as html
  # tbl_merge(tbls = list(unwt_stats, wt_stats),
  #           tab_spanner = c("**Unweighted**", "**Weighted**")) %>%
  #   as_gt() %>%
  #   gtsave(filename = "table_1_unweighted-weighted_new.html")
  # Import table into Word and make any formatting edits (Insert, Text:Object:Text from file)
  
  # Merge unweighted and weighted - save as word doc
  # tbl_merge(tbls = list(unwt_stats, wt_stats),
  #           tab_spanner = c("**Unweighted**", "**Weighted**")) %>%
  #   as_flex_table() %>%
  #   save_as_docx(path = "table_1_unweighted-weighted_new.docx")
  
  
  
  # Save Weighted and Unweighted tables separately
  unwt_stats %>%
    as_flex_table() %>%
    save_as_docx(path = "table_1_unweighted_new.docx")
  wt_stats %>%
    as_flex_table() %>%
    save_as_docx(path = "table_1_weighted_new.docx")
  
  #############################################################################################################
  
  setwd(current_directory)
}