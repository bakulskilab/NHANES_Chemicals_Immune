identify_chemicals_degrees_freedom <- function(nhanes_full_dataset,
                                               nhanes_comments,
                                               chemical_dataset,
                                               chem_master,
                                               weights_dataset)
{
  library(tidyverse)
  library(sjlabelled)
  library(magrittr) # %T>% package
  library(usefun)
  
  total_chems <- chemical_dataset %>%
    dplyr::select(-SEQN,
                  -SDDSRVYR) %>%
    colnames(.)
  
  chems <- chem_master %>%
    filter(!chem_family %in% c("Dietary Components",
                               "Phytoestrogens")) %>%
    filter(!chemical_codename_use %in% c("LBXEST",
                                         "SSCYA",
                                         "SSMEL",
                                         "URXUCR",
                                         "URXANBT",
                                         "URXANTT",
                                         "URXCOXT",
                                         "URXHPBT",
                                         "URXNICT",
                                         "URXNNCT",
                                         "URXNOXT")) %>%
    distinct(chemical_codename_use) %>%
    pull(chemical_codename_use)
  
  #make a vector of the 8 immune measures
  immune <- c("LBXLYPCT", #lymphocytes
              "LBXMOPCT", #monocytes
              "LBXNEPCT", #neutrophils
              "LBXEOPCT", #eosinophils
              "LBXBAPCT", #basophils
              "LBXWBCSI", #WBC count
              "LBXRBCSI", #RBC count
              "LBXMCVSI")  #MCV
  
  #make a vector of the 5 demographics + SEQN + creatinine
  demog <- c("SEQN",     #ID
             "RIDAGEYR", #age
             "RIDRETH1", #race
             "RIAGENDR", #gender
             "INDFMPIR", #poverty-income ratio
             "BMXWAIST", #waist circumference
             "URXUCR")   #creatinine
  
  nhanes_subset_unclean <- nhanes_full_dataset %>%
    dplyr::select(all_of(demog),
                  all_of(immune),
                  all_of(chems),
                  "SDMVSTRA",
                  "SDMVPSU",
                  "SDDSRVYR") %T>%
    {print("total dataset")} %T>%
    {print(dim(.))} %>% #101316    501 (7 demog + 8 immune + 486 chems)
    
    #drop participants under 18 (children)
    filter(RIDAGEYR >= 18)  %T>%
    {print("drop <18 years")} %T>%
    {print(dim(.))} %>% #59204
    
    #dropping participants with missing demographics
    drop_na(RIDAGEYR,
            RIDRETH1,
            RIAGENDR) %T>%
    {print("drop age, race, gender")} %T>%
    {print(dim(.))} %>% #59204
    drop_na(INDFMPIR) %T>%
    {print("drop poverty-income ratio missings")} %T>%
    {print(dim(.))} %>% #53446
    drop_na(BMXWAIST) %T>%
    {print("drop waist circumference missings")} %T>%
    {print(dim(.))} %>% #48484
    drop_na(LBXCOT) %T>%
    {print("drop cotinine missings")} %T>%
    {print(dim(.))} %>% #46080
    
    #dropping participants with missing cell counts
    drop_na(LBXWBCSI) %T>%
    {print("drop WBC missings")} %T>%
    {print(dim(.))} %>% #46004
    drop_na(LBXRBCSI) %T>%
    {print("drop RBC missings")} %T>%
    {print(dim(.))} %>% #46004
    drop_na(LBXLYPCT,
            LBXMOPCT,
            LBXNEPCT,
            LBXEOPCT,
            LBXBAPCT,
            LBXMCVSI) %>%
    filter(!SEQN == "102389")  %T>% #WBC count of 400 which is > upper LOD
    {print("drop cell type missings")} %T>%
    {print(dim(.))} #45870
  
  sum(!is.na(nhanes_subset_unclean$URXUCD)) #13592 total measurements
  nhanes_subset_unclean %>%
    dplyr::select(URXUCD) %>%
    filter(URXUCD == 0) %>%
    nrow() #7
  
  #figure out who these participants are for dropping their measurements later in the comments dataset
  ucd_partic <- nhanes_subset_unclean %>%
    filter(URXUCD == 0) %>%
    pull(SEQN)
  
  # Convert measurements of 0 to NA
  nhanes_subset_cleaning <- nhanes_subset_unclean %>%
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD))
  # 45870   501 (486 chems + 15 other variables)
  
  #check new measurement counts
  sum(!is.na(nhanes_subset_cleaning$URXUCD)) #13585 total measurements
  # rm(nhanes_subset_unclean)
  
  stats_degrees_freedom <- chems %>% 
    map(.
        , calculate_degrees_of_freedom
        , nhanes_subset_cleaning) %>%
    bind_rows(.)
  
  return(stats_degrees_freedom)
}