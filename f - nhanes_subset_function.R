#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########################################  CREATE WORKING DATASET  ###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Purpose: This function creates and saves the nhanes_subset dataset that I'll be using for analysis
#          
# Inputs:   nhanes_full_dataset   - dataframe containing all variables including cell types, chemicals,
#                                   and covariates (nhanes_merged_dataset)
#           subset_chemicals      - vector of chemicals that passed the LOD check
#
# Outputs:  nhanes_subset_dataset - dataframe of subsetted covariates, cell types, and chemicals (nhanes_subset)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

nhanes_subset_function <- function(nhanes_full_dataset,
                                   subset_chemicals)
{
  library(tidyverse)
  library(magrittr)
  
  # setwd(current_directory)

  #TEMPORARY
  # nhanes_full_dataset <- nhanes_merged_dataset
  # subset_chemicals <- use_these_chems
  
  #############################################################################################################
  ################################## Identify Variables to Include For Analysis ###############################
  #############################################################################################################
  
  #demographics
  demog <- c("SEQN",     #ID
             "RIDAGEYR", #age
             "RIDRETH1", #race
             "RIAGENDR", #gender
             "INDFMPIR", #poverty-income ratio
             "BMXWAIST", #waist circumference
             "SDDSRVYR", #cycle
             "URXUCR",   #creatinine
             "LBXCOT"   #cotinine
             ) 
  
  
  #immune measures
  immune <- c("LBDLYMNO", #lymphocyte counts
              "LBDMONO",  #monocyte counts
              "LBDNENO",  #neutrophil counts
              "LBDEONO",  #eosinophil counts
              "LBDBANO",  #basophil counts
              "LBXWBCSI", #WBC count
              "LBXRBCSI", #RBC count
              "LBXMCVSI"  #MCV
              )
  
  #chemicals
  chems <- subset_chemicals %>%
    pull(chemical_codename_use)
  
  #############################################################################################################
  ################################# Subset Initial Dataset To Have Complete Data ##############################
  #############################################################################################################

  #subset nhanes_full_dataset down to only these variables and drop participants with missing data
  #filter by missing cotinine as well because it will become the smoking covariate
  nhanes_subset_unclean <- nhanes_full_dataset %>%
                              dplyr::select(all_of(demog),
                                            all_of(immune),
                                            all_of(chems)) %T>%
                              {print("total dataset")} %T>%
                              {print(dim(.))} %>%
                              #101316    214
    
                              #drop participants under 18 (children)
                              filter(RIDAGEYR >= 18)  %T>%
                              {print("drop <18 years")} %T>%
                              {print(dim(.))} %>%
                              #59204   214
                              
                              #dropping participants with missing age, race, gender
                              drop_na(RIDAGEYR,
                                      RIDRETH1,
                                      RIAGENDR) %T>%
                              {print("drop age, race, gender")} %T>%
                              {print(dim(.))} %>%
                              #59204   214
                              drop_na(INDFMPIR) %T>%
                              {print("drop poverty-income ratio missings")} %T>%
                              {print(dim(.))} %>%
                              #53446   214
                              drop_na(BMXWAIST) %T>%
                              {print("drop waist circumference missings")} %T>%
                              {print(dim(.))} %>%
                              #48484   214
                              drop_na(LBXCOT) %T>%
                              {print("drop cotinine missings")} %T>%
                              {print(dim(.))} %>%
                              #46080   214
          
                              #dropping participants with missing cell counts
                              drop_na(LBXWBCSI) %T>%
                              {print("drop WBC missings")} %T>%
                              {print(dim(.))} %>%
                              #46004   214
                              drop_na(LBXRBCSI) %T>%
                              {print("drop RBC missings")} %T>%
                              {print(dim(.))} %>%
                              #46004   214
                              drop_na(LBDLYMNO, 
                                      LBDMONO, 
                                      LBDNENO, 
                                      LBDEONO, 
                                      LBDBANO,
                                      LBXMCVSI) %>%
                              filter(!SEQN == "102389")  %T>%
                              {print("drop cell type missings")} %T>%
                              {print(dim(.))}
                              #45870   214
                              
  print("dataset pre measurement cleaning")
  print(dim(nhanes_subset_unclean))
  #45870   214
  #this creates a dataset of variables and all complete participants (good)
  
  #############################################################################################################
  ################################### Drop Urinary Cadmium Measurements of 0 ##################################
  #############################################################################################################
  
  # Step 3. Drop participants with urinary cadmium measurements of 0
  # Measurements are incorrect in NHANES due to urinary cadmium measurements being crossed with molybdenum
  
  #how many measurements
  sum(!is.na(nhanes_subset_unclean$URXUCD)) #13592 total measurements
  nhanes_subset_unclean %>%
    dplyr::select(URXUCD) %>%
    filter(URXUCD == 0) %>%
    nrow() #7 measurements are 0
  
  # Convert measurements of 0 to NA
  nhanes_subset_cleaning <- nhanes_subset_unclean %>%
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD))
  # 45870   214 (198 chems + 16 other variables)
  
  #check new measurement counts
  sum(!is.na(nhanes_subset_cleaning$URXUCD)) #13585 total measurements

  
  #############################################################################################################
  ############################## Drop Participants With No Chemical Measurements ##############################
  #############################################################################################################
  
  # Drop participants with no chemical measurements
  
  print("______________")
  
  # Make a long dataset and drop participants
  long_cleaning <- nhanes_subset_cleaning %>%
    pivot_longer(cols = all_of(chems),
                 names_to = "chemical_codename",
                 values_to = "measurements") %>%
    drop_na(measurements) #no participants dropped
  
  #############################################################################################################
  ############################## Drop Participants With No Creatinine Measurements ############################
  #############################################################################################################
  
  # Step 6. Drop participants with urinary chemical measurements but no creatinine measurements

  nhanes_creatinine <- long_cleaning %>%
    mutate(urinary = case_when(str_detect(chemical_codename, "^LB", negate = TRUE) ~ 99999999)) %>%
    mutate(creat_na = ifelse(urinary == 99999999 & is.na(URXUCR), NA, measurements)) %>%
    drop_na(creat_na) %>%
    dplyr::select(-urinary,
                  -creat_na) %>%
    pivot_wider(names_from = chemical_codename,
                values_from = measurements)

  dim(nhanes_creatinine) #45528   214
  rm(long_cleaning)
  
  #also drop URXNAL measurement from SEQN: 65564 because their creatinine measurement is log2(0) = -inf
  #Drop the creatinine=0 too because this participant has no urinary measures in this subset
  #Only 1 person has creatinine=0
  nhanes_creatinine <- nhanes_creatinine %>%
    mutate(URXNAL = case_when(SEQN == 65564 & URXNAL == 0.0553 ~ 5000000000,
                              TRUE ~ as.numeric(URXNAL))) %>%
    mutate(URXNAL = na_if(URXNAL, 5000000000)) %>%
    mutate(URXUCR = ifelse(URXUCR == 0, NA, URXUCR))

  # print(colnames(nhanes_creatinine))
  #############################################################################################################
  ################################## Create A Smoking Variable From Cotinine ##################################
  #############################################################################################################

  #create a smoking variable for a covariate
  nhanes_subset_smk <- nhanes_creatinine %>%
    mutate(SMOKING = LBXCOT)

  #############################################################################################################
  ################################# Add Survey Variables Into Cleaned Dataset #################################
  #############################################################################################################

  #pull survey variables from nhanes_full_dataset
  survey_vars <- nhanes_full_dataset %>%
    dplyr::select(SEQN,
                  SDMVPSU,  #30 primary sampling units per cycle - mostly single counties, selected from strata
                  SDMVSTRA) #census region/metropolitan area/pop demographics
  # print(colnames(survey_vars))

  #merge into nhanes_subset_smk
  nhanes_subset_survey <- left_join(nhanes_subset_smk, survey_vars, by = "SEQN") #%>%
    # relocate(c("SDMVPSU",
    #            "SDMVSTRA"),
    #          .after = SDDSRVYR)

  # print(colnames(nhanes_subset_survey))

  #############################################################################################################
  ######################################### Fix The Variable Classes ##########################################
  #############################################################################################################

  # str(nhanes_subset_survey)

  nhanes_subset_relevel <- nhanes_subset_survey %>%
    mutate(RIAGENDR = relevel(factor(RIAGENDR), ref = 1),
           RIDRETH1 = relevel(factor(RIDRETH1), ref = 3))
  # print(str(nhanes_subset_relevel))

  print("final dataset dimensions")
  print(dim(nhanes_subset_relevel))
  #45528   217

  #############################################################################################################

  print("return(nhanes_subset_dataset)")
  return(nhanes_subset_survey)
}