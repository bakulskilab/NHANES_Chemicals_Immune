#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###################################  MAKE TABLES OF DEMOGRAPHIC INFORMATION  ##################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function makes Table 1 using demographic and immune measure data
#          
# Inputs: nhanes_subset - dataframe containing complete demographic and immune measures data for each
#                         participant
# 
#         nhanes_full_dataset - dataframe containing all participants, demographic, chemical, and other
#                               variables
#
# Outputs: Table of means and medians + p-value by inclusion status

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_demographics <- function(nhanes_subset)
{
  library(tidyverse)
  library(sjlabelled)
  
  #set the number of digits after the decimal point to be 2
  options(digits=4)
  
  #TEMPORARY
  # nhanes_subset <- nhanes_subset_dataset
  
  #############################################################################################################
  #################################### Calculate The Continuous Demographics ##################################
  #############################################################################################################

  #select the 9 immune measures + age
  continuous <- c("RIDAGEYR", #age
                  "INDFMPIR", #PIR
                  "BMXWAIST", #waist
                  "SMOKING",  #smoking
                  "LBXLYPCT", #lymphocytes
                  "LBXNEPCT", #neutrophils
                  "LBXMOPCT", #monocytes
                  "LBXBAPCT", #basophils
                  "LBXEOPCT", #eosinophils
                  "LBXWBCSI", #WBC count
                  "LBXRBCSI", #RBC count
                  "LBXMCVSI"  #MCV
                 )
  
  #find the number of elements in vector: continuous
  num_contin_vars <- length(continuous)


  #calculate continuous measures
  demog_stats_continuous <- nhanes_subset %>%
    dplyr::select(all_of(continuous)) %>%
    summarise_all(list(mean,
                       sd,
                       median,
                       IQR,
                       min,
                       max)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "variable") %>%
    mutate(stat_type = rep(c("mean",
                             "sd",
                             "median",
                             "IQR",
                             "min",
                             "max"),
                           each = num_contin_vars)) %>%
    mutate(variable = gsub("_.*","", variable)) %>%
    arrange(variable) %>%
    pivot_wider(names_from = "stat_type",
                values_from = "V1") %>%
    mutate(mean_round = format(round(mean, digits = 1), nsmall = 1),
           sd_round = format(round(sd, digits = 1), nsmall = 1),
           median_round = format(round(median, digits = 1), nsmall = 1),
           iqr_round = format(round(IQR, digits = 1), nsmall = 1),
           min_round = format(round(min, digits = 1), nsmall = 1),
           max_round = format(round(max, digits = 1)), nsmall = 1) %>%
    dplyr::select(-mean, -sd, -median, -IQR, -min, -max, -nsmall) %>%
    mutate("mean (sd)" = paste0(mean_round, " (", sd_round, ")"),
           "median (IQR)" = paste0(median_round, " (", iqr_round, ")"),
           "range" = paste0(min_round, " - ", max_round)) %>%
    dplyr::select(variable,
                  "mean (sd)",
                  "median (IQR)",
                  "range")
  
  #calculate urinary creatinine separately because of the NAs
  creatinine <- nhanes_subset %>%
    drop_na(URXUCR) %>%
    summarise(mean = mean(URXUCR),
              sd = sd(URXUCR),
              median = median(URXUCR),
              IQR = IQR(URXUCR),
              min = min(URXUCR),
              max = max(URXUCR)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "stat_type") %>%
    mutate(variable = rep("URXUCR")) %>%
    # relocate(variable, .before = stat_type) %>%
    pivot_wider(names_from = "stat_type",
                values_from = "V1") %>%
    mutate(mean_round = format(round(mean, digits = 1), nsmall = 1),
           sd_round = format(round(sd, digits = 1), nsmall = 1),
           median_round = format(round(median, digits = 1), nsmall = 1),
           iqr_round = format(round(IQR, digits = 1), nsmall = 1),
           min_round = format(round(min, digits = 1), nsmall = 1),
           max_round = format(round(max, digits = 1)), nsmall = 1) %>%
    dplyr::select(-mean, -sd, -median, -IQR, -min, -max, -nsmall) %>%
    mutate("mean (sd)" = paste0(mean_round, " (", sd_round, ")"),
           "median (IQR)" = paste0(median_round, " (", iqr_round, ")"),
           "range" = paste0(min_round, " - ", max_round)) %>%
    dplyr::select(variable,
                  "mean (sd)",
                  "median (IQR)",
                  "range")
  
  #merge creatinine into continuous variables dataset
  contin_creat <- rbind(demog_stats_continuous,
                        creatinine)

  #split up the continuous variables into tables 1 and 2
  table_1_contin <- contin_creat %>%
    filter(variable == "RIDAGEYR" |
           variable == "INDFMPIR" |
           variable == "BMXWAIST" |
           variable == "URXUCR" |
           variable == "SMOKING")
  
  #Save table 1 continuous
  setwd(paste0(current_directory, "/Tables - Table 1"))
  write.csv(table_1_contin, "table_1_contin.csv")
  
  
  #############################################################################################################
  ################################### Calculate The Categorical Demographics ##################################
  #############################################################################################################
  
  #calculate categorical measures
  demog_stats_gender <- as.data.frame(table(nhanes_subset$RIAGENDR)) %>%
    mutate(Proportion = Freq / sum(Freq))
  
  demog_stats_race <- as.data.frame(table(nhanes_subset$RIDRETH1)) %>%
    mutate(Proportion = Freq / sum(Freq))
  
  #put them together for easier viewing
  demog_stats_categorical <- bind_rows(demog_stats_gender, demog_stats_race)
  
  # View(demog_stats_categorical)
  
  #add labels
  table_1_categ <- demog_stats_categorical %>%
    mutate(variable = NA)
  table_1_categ[1:2,4] <- "RIAGENDR"
  table_1_categ[3:7,4] <- "RIDRETH1"
  
  #combine columns
  table_1_categ <- table_1_categ %>%
    mutate(perc_categ = round(Proportion*100, digits = 1)) %>%
    mutate(count_perc = paste0(Freq, " (", perc_categ, ")"))
  View(table_1_categ)
  
  #Save
  write.csv(table_1_categ, "table_1_categ.csv")
  
  #############################################################################################################
  ########################################## Table 2 - Immune Measures ########################################
  #############################################################################################################
  
  #table 2
  table_2 <- contin_creat %>%
    filter(variable == "LBXLYPCT" |
             variable == "LBXNEPCT" |
             variable == "LBXMOPCT" |
             variable == "LBXBAPCT" |
             variable == "LBXEOPCT" |
             variable == "LBXWBCSI" |
             variable == "LBXRBCSI" |
             variable == "LBXMCVSI") %>%
    mutate(cell_order = case_when(variable == "LBXLYPCT" ~ 1,
                                  variable == "LBXNEPCT" ~ 2,
                                  variable == "LBXMOPCT" ~ 3,
                                  variable == "LBXBAPCT" ~ 4,
                                  variable == "LBXEOPCT" ~ 5,
                                  variable == "LBXWBCSI" ~ 6,
                                  variable == "LBXRBCSI" ~ 7,
                                  variable == "LBXMCVSI" ~ 8)) %>%
    arrange(cell_order) %>%
    dplyr::select(-cell_order)
  
  # Save Table 2
  write.csv(table_2, "table_2.csv")

  
  #############################################################################################################
  print("save csv files of stats in Tables folder")
  setwd(current_directory)
}