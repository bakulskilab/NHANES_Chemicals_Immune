identify_chemicals_all_df <- function(nhanes_full_dataset,
                                        nhanes_comments,
                                        chemical_dataset,
                                        chem_master,
                                        weights_dataset)
{
  library(tidyverse)
  library(sjlabelled)
  library(magrittr) # %T>% package
  library(usefun)
  
  weights_dataset <- weights_dataset %>%
    filter(SDDSRVYR != -1)
  
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
  
  #make a vector of the 13 immune measures
  immune <- c("LBDLYMNO", #lymphocyte counts
              "LBDMONO",  #monocyte counts
              "LBDNENO",  #neutrophil counts
              "LBDEONO",  #eosinophil counts
              "LBDBANO",  #basophil counts
              "LBXWBCSI", #WBC count
              "LBXRBCSI", #RBC count
              "LBXMCVSI"  #MCV
              )  
  
  #make a vector of the 5 demographics + SEQN + creatinine + cotinine + PSU + strata
  demog <- c("SEQN",     #ID
             "RIDAGEYR", #age
             "RIDRETH1", #race
             "RIAGENDR", #gender
             "INDFMPIR", #poverty-income ratio
             "BMXWAIST", #waist circumference
             "URXUCR",   #creatinine
             "LBXCOT",   #cotinine
             "SDMVPSU",  #PSU
             "SDMVSTRA") #strata
  
  nhanes_subset_unclean <- nhanes_full_dataset %>%
    dplyr::select(all_of(demog),
                  all_of(immune),
                  all_of(chems)) %T>%
    {print("total dataset")} %T>%
    {print(dim(.))} %>% #101316    501 (7 demog + 13 immune + 486 chems)
    
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
    drop_na(LBDLYMNO, 
            LBDMONO, 
            LBDNENO, 
            LBDEONO, 
            LBDBANO,
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
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD)) %T>% 
    {print("convert urinary cadmium measurements of 0 to NA")} %T>%
    {print(dim(.))}
  # 45870   501 (486 chems + 15 other variables)
  
  long_cleaning <- nhanes_subset_cleaning %>%
    pivot_longer(cols = all_of(chems),
                 names_to = "chemical_codename",
                 values_to = "measurements") %>%
    drop_na(measurements) #no participants dropped
  
  long_creatinine_seqn <- long_cleaning %>%
    mutate(urinary = case_when(str_detect(chemical_codename, "^LB", negate = TRUE) ~ 99999999)) %>%
    mutate(creat_na = ifelse(urinary == 99999999 & is.na(URXUCR), NA, measurements)) %>%
    drop_na(creat_na) %>%
    dplyr::select(-urinary,
                  -creat_na) %>%
    distinct(SEQN) %>%
    pull(SEQN)
  
  rm(long_cleaning)
  print("number of participants with urinary chemical measurements but no creatinine measurements")
  print(length(long_creatinine_seqn)) #45528 - dropped 342 participants
  
  #keep this subset of participants
  nhanes_creatinine <- nhanes_subset_cleaning %>%
    filter(SEQN %in% long_creatinine_seqn)
  
  nhanes_creatinine <- nhanes_creatinine %>%
    mutate(URXNAL = case_when(SEQN == 65564 & URXNAL == 0.0553 ~ 5000000000,
                              TRUE ~ as.numeric(URXNAL))) %>%
    mutate(URXNAL = na_if(URXNAL, 5000000000)) %>%
    mutate(URXUCR = ifelse(URXUCR == 0, NA, URXUCR))
  
  #check new measurement counts
  # sum(!is.na(nhanes_subset_unclean$URXUCD)) #13585 total measurements

  subset_comments <- nhanes_subset_unclean %>%
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD)) %>%
    dplyr::select("SEQN"
                  , all_of(demog)
                  , "URXUCD") %>%
    left_join(.
              , nhanes_comments
              , by = "SEQN") %>%
    filter(SEQN %in% long_creatinine_seqn) %T>%
    {print("dimension of comments dataset with included participants")} %T>%
    {print(dim(.))}
  # 45528   436

  comments_codenames <- colnames(nhanes_comments)
  # Remove SEQN, SDDSRVYR, and estradiol to prevent calculations of their detection frequency
  # estradiol is not an exogenous chemical
  excess_codenames <- which(comments_codenames %in% c("SEQN"
                                                      , "SDDSRVYR"
                                                      , "LBDESTLC"
                                                      ))
  comments_codenames <- comments_codenames[-excess_codenames]

  subset_weights <- nhanes_subset_unclean %>%
    mutate(URXUCD = ifelse(URXUCD == 0, NA, URXUCD)) %>%
    dplyr::select("SEQN"
                  , "URXUCR") %>%
    left_join(.
              , weights_dataset
              , by = "SEQN") %>%
    filter(SEQN %in% long_creatinine_seqn) %T>%
    {print("dimension of weights dataset with included participants")} %T>%
    {print(dim(.))}
  # 45528   542

  stats_weight <- comments_codenames %>% #[1:2] %>%
    #c(comments_codenames[1:6], "URD4FPLC") %>%
    # c("URDUCDLC"
    #   # , "LBDBCDLC"
    #   # , "LBD196LC"
    #   # , "LBD138LC"
    #   # , "LBDBPBLC"
    #   # , "SSMONPL"
    #   ) %>%
    map(.
        , calculate_detection_frequency_degrees_freedom
        , subset_comments
        , subset_weights
        , chem_master
        , demog) %>%
    bind_rows(.)
  # View(stats_weight)

  stats_weight <- get_counts_for_cu_zn(df_inclusion_criteria_stats = stats_weight
                                       , nhanes_subset = nhanes_subset_unclean
                                       , df_weights = subset_weights
                                       , demographics = demog)

  weird_smk_chems <- c("URXANBT",
                       "URXANTT",
                       "URXCOXT",
                       "URXHPBT",
                       "URXNICT",
                       "URXNNCT",
                       "URXNOXT")

  stats_weight <- stats_weight %>%
    filter(!(chemical_codename_use %in% weird_smk_chems)) %>%
    filter(chem_family != "Dietary Components") %>%
    filter(chem_family != "Phytoestrogens") %>%
    mutate(include = ifelse(above_percentage_unweighted >= 50 &
                              above_percentage_weighted >= 50 &
                              degrees_of_freedom >= 8
                            , "yes"
                            , "no")) %>%
    mutate(chemical_codename_use = ifelse(chemical_codename_use == "LBX138158LA"
                                          , "LBX138LA"
                                          , chemical_codename_use)) %>%
    mutate(chemical_codename_use = ifelse(chemical_codename_use == "LBX196203LA"
                                          , "LBX196LA"
                                          , chemical_codename_use)) %>%
    filter(include == "yes")
  # View(stats_weight)

  return(stats_weight)
}